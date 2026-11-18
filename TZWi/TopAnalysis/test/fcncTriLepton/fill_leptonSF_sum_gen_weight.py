#!/usr/bin/env python3
"""
Compute sum_gen_weight = sum(Runs.genEventSumw) for every MC sample in
leptonSF_samples_2024.yaml and write a filled copy of the YAML file.

The sum is read from the ORIGINAL NanoAOD files listed under the dataset
catalogue. It must not be computed from selected lepton-SF NTuples.
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

try:
    import ROOT
except ImportError as exc:
    raise RuntimeError(
        "PyROOT is unavailable. Run this inside the CMSSW environment after "
        "'cmsenv' or 'eval `scram runtime -sh`'."
    ) from exc

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)


# Keys are the sample 'name' values in leptonSF_samples_2024.yaml.
# Values are the dataset-group directory names in the NanoAOD catalogue.
DEFAULT_DATASET_GROUPS: Dict[str, str] = {
    "WZto3LNu": "MC2024.WZto3LNu",
    "DYGto2LG_4to50": "MC2024.DYGto2LG-1Jets_Bin-MLL-4to50",
    "DYGto2LG_50": "MC2024.DYGto2LG-1Jets_Bin-MLL-50",
    "WGtoLNuG": "MC2024.WGtoLNuG-1Jets",
    "ZZto4L": "MC2024.ZZto4L",
    "TTLL_4to50": "MC2024.TTLL_Bin-MLL-4to50",
    "TTLL_50": "MC2024.TTLL_Bin-MLL-50",
    "TZQB": "MC2024.TZQB-Zto2L-4FS_Bin-MLL-30",
    "TTto2L2Nu": "MC2024.TTto2L2Nu",
    "WWZ": "MC2024.WWZ",
    "WZZ": "MC2024.WZZ",
    "ZZZ": "MC2024.ZZZ",
    "TTH_HtoNon2B": "MC2024.TTH-HtoNon2B_Par-M-125",
    "TWZ_2Q_LNu_2L": "MC2024.TWZ-Tto2Q-WtoLNu-Zto2L-DR1",
    "TWZ_LNu_2Q_2L": "MC2024.TWZ-TtoLNu-Wto2Q-Zto2L-DR1",
    "TWZ_LNu_LNu_2L": "MC2024.TWZ-TtoLNu-WtoLNu-Zto2L-DR1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate sum(Runs.genEventSumw) from original NanoAOD files and "
            "fill sum_gen_weight in a lepton-SF YAML configuration."
        )
    )
    parser.add_argument(
        "--config",
        default="leptonSF_samples_2024.yaml",
        help="Input lepton-SF YAML configuration.",
    )
    parser.add_argument(
        "--dataset-base",
        default="",
        help=(
            "Base directory containing MC.RunIISummer24.central/MC2024.*/...txt. "
            "If omitted, common locations under $CMSSW_BASE are tried."
        ),
    )
    parser.add_argument(
        "--output",
        default="",
        help=(
            "Output YAML path. Default: <config stem>_filled.yaml. "
            "The input file is not overwritten."
        ),
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Overwrite --config after creating <config>.bak.",
    )
    parser.add_argument(
        "--lumi-pb",
        type=float,
        default=None,
        help="Optionally set lumi_pb in the output YAML at the same time.",
    )
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        help="Compute only this sample name. May be specified more than once.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Number of ROOT files opened concurrently. Default 1 is safest; "
            "use 2-4 only when the storage endpoint permits it."
        ),
    )
    parser.add_argument(
        "--xrootd-prefix",
        default="root://cms-xrdr.private.lo:2094//xrd/",
        help=(
            "Prefix used when a file-list entry starts with /xrootd/ but is "
            "not locally mounted."
        ),
    )
    parser.add_argument(
        "--allow-bad-files",
        action="store_true",
        help=(
            "Continue when an original NanoAOD file cannot be read. This is "
            "not recommended for final normalization."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve file lists and print file counts without opening ROOT files.",
    )
    return parser.parse_args()


def resolve_dataset_base(value: str) -> Path:
    candidates: List[Path] = []
    if value:
        candidates.append(Path(os.path.expandvars(os.path.expanduser(value))))

    cmssw_base = os.environ.get("CMSSW_BASE", "")
    if cmssw_base:
        top = Path(cmssw_base) / "src" / "TZWi" / "TopAnalysis"
        candidates.extend(
            [
                top / "data" / "datasets" / "NanoAOD" / "2024",
                top / "test" / "fcncTriLepton" / "NanoAOD" / "2024",
            ]
        )

    candidates.extend(
        [
            Path("NanoAOD/2024"),
            Path("data/datasets/NanoAOD/2024"),
            Path("../../data/datasets/NanoAOD/2024"),
        ]
    )

    seen = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_dir():
            return resolved

    formatted = "\n  - ".join(str(path) for path in candidates)
    raise RuntimeError(
        "Could not locate the NanoAOD dataset catalogue. Tried:\n  - " + formatted
    )


def collect_mc_samples(config: dict) -> List[dict]:
    samples: List[dict] = []
    signal = config.get("signal") or {}
    if signal:
        samples.append(signal)
    samples.extend(list(config.get("backgrounds") or []))

    for sample in samples:
        name = str(sample.get("name", "")).strip()
        if not name:
            raise RuntimeError("Every signal/background sample needs a non-empty name.")
    return samples


def find_file_lists(dataset_base: Path, dataset_group: str) -> List[Path]:
    # Normal layout:
    #   <base>/MC.RunIISummer24.central/<dataset_group>/*.txt
    patterns = [
        str(dataset_base / "*" / dataset_group / "*.txt"),
        str(dataset_base / dataset_group / "*.txt"),
    ]
    matches: List[str] = []
    for pattern in patterns:
        matches.extend(glob.glob(pattern))

    if not matches:
        matches.extend(
            glob.glob(str(dataset_base / "**" / dataset_group / "*.txt"), recursive=True)
        )

    return [Path(path) for path in sorted(dict.fromkeys(matches))]


def read_file_list(path: Path) -> List[str]:
    files: List[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            files.append(line)
    return files


def normalize_root_url(raw: str, xrootd_prefix: str) -> str:
    value = os.path.expandvars(os.path.expanduser(raw.strip()))
    if value.startswith(("root://", "file://")):
        return value

    local = Path(value)
    if local.exists():
        return str(local)

    if value.startswith("/xrootd/"):
        suffix = value[len("/xrootd/") :]
        return xrootd_prefix.rstrip("/") + "/" + suffix

    # Some catalogues already expose /new_xrootd as a local mount. If it is not
    # visible here, retain the path so ROOT gives a clear open error.
    return value


def read_runs_sumw(url: str) -> Tuple[str, float, float, int, Optional[str]]:
    """Return (url, sumw, sumw2, run_entries, error)."""
    root_file = ROOT.TFile.Open(url, "READ")
    if not root_file or root_file.IsZombie():
        if root_file:
            root_file.Close()
        return url, 0.0, 0.0, 0, "TFile.Open failed or returned a zombie file"

    try:
        runs = root_file.Get("Runs")
        if not runs:
            return url, 0.0, 0.0, 0, "missing Runs tree"
        if not runs.GetBranch("genEventSumw"):
            return url, 0.0, 0.0, 0, "Runs tree is missing genEventSumw"

        has_sumw2 = bool(runs.GetBranch("genEventSumw2"))
        sumw = 0.0
        sumw2 = 0.0
        n_entries = int(runs.GetEntries())
        for entry in runs:
            sumw += float(entry.genEventSumw)
            if has_sumw2:
                sumw2 += float(entry.genEventSumw2)

        if not math.isfinite(sumw):
            return url, 0.0, 0.0, n_entries, "non-finite genEventSumw sum"
        return url, sumw, sumw2, n_entries, None
    except Exception as exc:  # PyROOT errors vary by release.
        return url, 0.0, 0.0, 0, "%s: %s" % (type(exc).__name__, exc)
    finally:
        root_file.Close()


def calculate_sample(
    sample_name: str,
    dataset_group: str,
    file_lists: Sequence[Path],
    workers: int,
    xrootd_prefix: str,
    allow_bad_files: bool,
    dry_run: bool,
) -> dict:
    raw_files: List[str] = []
    for file_list in file_lists:
        raw_files.extend(read_file_list(file_list))

    # Avoid double counting if the same ROOT URL appears in more than one list.
    unique_raw = list(dict.fromkeys(raw_files))
    urls = [normalize_root_url(path, xrootd_prefix) for path in unique_raw]

    print(
        "[sumw] %-24s group=%-52s lists=%2d files=%5d"
        % (sample_name, dataset_group, len(file_lists), len(urls))
    )
    for path in file_lists:
        print("       file list:", path)

    if not urls:
        raise RuntimeError("No original NanoAOD files found for %s" % sample_name)

    if dry_run:
        return {
            "name": sample_name,
            "dataset_group": dataset_group,
            "file_lists": [str(path) for path in file_lists],
            "n_files": len(urls),
            "n_good_files": 0,
            "n_bad_files": 0,
            "n_runs_entries": 0,
            "sum_gen_weight": None,
            "sum_gen_weight2": None,
            "bad_files": [],
        }

    total_sumw = 0.0
    total_sumw2 = 0.0
    total_runs = 0
    bad: List[Tuple[str, str]] = []
    completed = 0

    def consume(result: Tuple[str, float, float, int, Optional[str]]) -> None:
        nonlocal total_sumw, total_sumw2, total_runs, completed
        url, sumw, sumw2, n_runs, error = result
        completed += 1
        if error:
            bad.append((url, error))
        else:
            total_sumw += sumw
            total_sumw2 += sumw2
            total_runs += n_runs
        if completed % 100 == 0 or completed == len(urls):
            print("       progress: %d/%d files" % (completed, len(urls)))

    if workers <= 1:
        for url in urls:
            consume(read_runs_sumw(url))
    else:
        # PyROOT file access from threads needs ROOT thread-safety enabled.
        if hasattr(ROOT, "EnableThreadSafety"):
            ROOT.EnableThreadSafety()
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(read_runs_sumw, url) for url in urls]
            for future in as_completed(futures):
                consume(future.result())

    if bad:
        print("[sumw] ERROR: %s has %d unreadable files:" % (sample_name, len(bad)))
        for url, error in bad[:20]:
            print("       BAD:", url)
            print("            ", error)
        if len(bad) > 20:
            print("       ... and %d more" % (len(bad) - 20))
        if not allow_bad_files:
            raise RuntimeError(
                "%s has unreadable files. Refusing to produce an incomplete "
                "normalization. Fix the files or use --allow-bad-files only "
                "for diagnostics." % sample_name
            )

    if not math.isfinite(total_sumw) or total_sumw == 0.0:
        raise RuntimeError(
            "Invalid total sum_gen_weight for %s: %r" % (sample_name, total_sumw)
        )

    result = {
        "name": sample_name,
        "dataset_group": dataset_group,
        "file_lists": [str(path) for path in file_lists],
        "n_files": len(urls),
        "n_good_files": len(urls) - len(bad),
        "n_bad_files": len(bad),
        "n_runs_entries": total_runs,
        "sum_gen_weight": total_sumw,
        "sum_gen_weight2": total_sumw2,
        "bad_files": bad,
    }
    print(
        "[sumw] RESULT %-24s sum_gen_weight = %.17g"
        % (sample_name, total_sumw)
    )
    return result


def format_float(value: float) -> str:
    # 17 significant digits preserve a C double when parsed again.
    return "%.17g" % float(value)


def patch_yaml_text(
    original_text: str,
    sumw_by_name: Dict[str, float],
    lumi_pb: Optional[float],
) -> str:
    lines = original_text.splitlines(keepends=True)
    current_name: Optional[str] = None
    replaced = set()
    lumi_replaced = False

    name_re = re.compile(r"^(\s*)name\s*:\s*([^#\n]+?)\s*(?:#.*)?$")
    sumw_re = re.compile(r"^(\s*)sum_gen_weight\s*:\s*[^#\n]*(\s*(?:#.*)?)$")
    lumi_re = re.compile(r"^(\s*)lumi_pb\s*:\s*[^#\n]*(\s*(?:#.*)?)$")

    output: List[str] = []
    for raw_line in lines:
        line = raw_line.rstrip("\r\n")
        newline = raw_line[len(line) :]

        name_match = name_re.match(line)
        if name_match:
            raw_name = name_match.group(2).strip().strip('"\'')
            current_name = raw_name

        if lumi_pb is not None:
            lumi_match = lumi_re.match(line)
            if lumi_match:
                line = "%slumi_pb: %s%s" % (
                    lumi_match.group(1),
                    format_float(lumi_pb),
                    lumi_match.group(2),
                )
                lumi_replaced = True

        sumw_match = sumw_re.match(line)
        if sumw_match and current_name in sumw_by_name:
            line = "%ssum_gen_weight: %s%s" % (
                sumw_match.group(1),
                format_float(sumw_by_name[current_name]),
                sumw_match.group(2),
            )
            replaced.add(current_name)

        output.append(line + newline)

    missing = sorted(set(sumw_by_name) - replaced)
    if missing:
        raise RuntimeError(
            "Could not locate sum_gen_weight lines for: %s" % ", ".join(missing)
        )
    if lumi_pb is not None and not lumi_replaced:
        raise RuntimeError("Could not locate lumi_pb in the YAML text.")

    return "".join(output)


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise RuntimeError("--workers must be at least 1")
    if args.lumi_pb is not None and args.lumi_pb <= 0.0:
        raise RuntimeError("--lumi-pb must be positive")
    if args.in_place and args.output:
        raise RuntimeError("Use either --in-place or --output, not both.")

    config_path = Path(args.config).resolve()
    if not config_path.is_file():
        raise RuntimeError("Configuration file not found: %s" % config_path)

    dataset_base = resolve_dataset_base(args.dataset_base)
    print("[sumw] config       =", config_path)
    print("[sumw] dataset base =", dataset_base)

    original_text = config_path.read_text()
    config = yaml.safe_load(original_text) or {}
    samples = collect_mc_samples(config)

    requested = set(args.sample)
    known_names = {str(sample["name"]) for sample in samples}
    unknown = sorted(requested - known_names)
    if unknown:
        raise RuntimeError("Unknown --sample name(s): %s" % ", ".join(unknown))

    results = []
    for sample in samples:
        name = str(sample["name"])
        if requested and name not in requested:
            continue

        # A future YAML can override the built-in mapping by adding:
        #   dataset_group: MC2024.Whatever
        dataset_group = str(
            sample.get("dataset_group") or DEFAULT_DATASET_GROUPS.get(name, "")
        )
        if not dataset_group:
            raise RuntimeError(
                "No dataset-group mapping for %s. Add 'dataset_group' to the "
                "sample block or extend DEFAULT_DATASET_GROUPS." % name
            )

        file_lists = find_file_lists(dataset_base, dataset_group)
        if not file_lists:
            raise RuntimeError(
                "No .txt file lists found for %s (%s) under %s"
                % (name, dataset_group, dataset_base)
            )

        results.append(
            calculate_sample(
                sample_name=name,
                dataset_group=dataset_group,
                file_lists=file_lists,
                workers=args.workers,
                xrootd_prefix=args.xrootd_prefix,
                allow_bad_files=args.allow_bad_files,
                dry_run=args.dry_run,
            )
        )

    print("\n# Copyable YAML values")
    for result in results:
        if result["sum_gen_weight"] is None:
            print("%-24s : <dry-run>" % result["name"])
        else:
            print(
                "%-24s : %s"
                % (result["name"], format_float(result["sum_gen_weight"]))
            )

    if args.dry_run:
        print("\n[sumw] Dry run complete; no YAML file was written.")
        return 0

    sumw_by_name = {
        result["name"]: float(result["sum_gen_weight"]) for result in results
    }
    patched = patch_yaml_text(original_text, sumw_by_name, args.lumi_pb)

    if args.in_place:
        backup = config_path.with_name(config_path.name + ".bak")
        shutil.copy2(config_path, backup)
        output_path = config_path
        print("[sumw] backup YAML =", backup)
    else:
        if args.output:
            output_path = Path(args.output).resolve()
        else:
            output_path = config_path.with_name(
                config_path.stem + "_filled" + config_path.suffix
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(patched)
    print("[sumw] wrote YAML   =", output_path)

    if args.allow_bad_files and any(result["n_bad_files"] for result in results):
        print(
            "[sumw] WARNING: output contains incomplete normalizations because "
            "--allow-bad-files was used.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print("[sumw] ERROR:", exc, file=sys.stderr)
        sys.exit(1)
