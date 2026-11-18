#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import math
import os
from array import array

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.TH1.AddDirectory(False)

CHANNELS_DEFAULT = ["ElElEl", "ElMuMu", "MuElEl", "MuMuMu"]
CL_1SIGMA = 0.682689492137


def parse_binning(text):
    vals = [float(x) for x in text.replace(" ", "").split(",") if x != ""]
    if len(vals) < 2:
        raise RuntimeError("Binning must contain at least two edges: %s" % text)
    if any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise RuntimeError("Binning edges must be strictly increasing: %s" % text)
    return vals


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute 2024 exact-3-lepton trigger scale factors from trigger-SF NTuples."
    )
    parser.add_argument("--data-glob", default="ntuple_2024_trigSF/reco/{channel}/*Run2024*/*.root",
                        help="Glob pattern for data ROOT files. Use {channel} as placeholder.")
    parser.add_argument("--mc-glob", default="ntuple_2024_trigSF/reco/{channel}/*WZto3LNu*/*.root",
                        help="Glob pattern for MC ROOT files. Use {channel} as placeholder.")
    parser.add_argument("--channels", nargs="+", default=CHANNELS_DEFAULT)
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--analysis-branch", default="HLT_analysis")
    parser.add_argument("--reference-branch", default="HLT_reference")
    parser.add_argument("--cutstep-branch", default="CutStep")
    parser.add_argument("--min-cutstep", type=int, default=1)
    parser.add_argument("--mc-weight-branch", default="",
                        help="Optional MC weight branch. If empty or missing, MC events are unweighted.")
    parser.add_argument("--min-nbjet", type=int, default=None,
                        help="Optional extra selection: nBjet >= min-nbjet.")
    parser.add_argument("--max-nbjet", type=int, default=None,
                        help="Optional extra selection: nBjet <= max-nbjet.")
    parser.add_argument("--deduplicate-data", action="store_true",
                        help="Remove duplicate data events using (run, luminosityBlock, event).")
    parser.add_argument("--make-2d", action="store_true",
                        help="Also build 2D efficiency and SF maps.")
    parser.add_argument("--pt-bins", default="20,25,30,40,50,70,100,150,200",
                        help="Comma-separated pT bin edges for 2D maps.")
    parser.add_argument("--eta-bins", default="0,0.8,1.479,2.0,2.5",
                        help="Comma-separated |eta| bin edges for 2D maps.")
    parser.add_argument("--min-2d-den-data", type=float, default=1.0,
                        help="Minimum data denominator to quote a 2D SF bin.")
    parser.add_argument("--min-2d-den-mc", type=float, default=1.0,
                        help="Minimum MC denominator to quote a 2D SF bin.")
    parser.add_argument("--output-root", default="triggerSF_2024.root")
    parser.add_argument("--output-csv", default="triggerSF_2024.csv")
    parser.add_argument("--output-json", default="triggerSF_2024.json")
    parser.add_argument("--bad-file-list", default="",
                        help="Optional text file to write input ROOT files that are zombie or missing the requested tree.")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def expand_files(pattern, channel):
    return sorted(glob.glob(pattern.format(channel=channel)))


def validate_files(files, tree_name):
    good = []
    bad = []
    for path in files:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            bad.append({"file": path, "reason": "zombie_or_unopenable"})
            if f:
                f.Close()
            continue
        t = f.Get(tree_name)
        if not t:
            bad.append({"file": path, "reason": "missing_tree_%s" % tree_name})
            f.Close()
            continue
        good.append(path)
        f.Close()
    return good, bad


def make_chain(files, tree_name):
    chain = ROOT.TChain(tree_name)
    for path in files:
        chain.Add(path)
    return chain


def require_branch(chain, branch):
    if not chain.GetBranch(branch):
        raise RuntimeError("Required branch '%s' is missing in chain with %d entries." % (branch, chain.GetEntries()))


def has_branch(chain, branch):
    return bool(chain.GetBranch(branch))


def get_bool(event, branch):
    return bool(getattr(event, branch))


def get_float(event, branch, default=1.0):
    try:
        return float(getattr(event, branch))
    except Exception:
        return float(default)


def get_int(event, branch, default=0):
    try:
        return int(getattr(event, branch))
    except Exception:
        return int(default)


def cp_interval(num_events, den_events, cl=CL_1SIGMA):
    """Return central efficiency and asymmetric Clopper-Pearson errors.

    This is used for unweighted event counts.  It remains non-zero for the
    one-sided case num == den, which is essential for trigger efficiencies that
    are observed to be 100% in finite data samples.
    """
    den_events = int(den_events)
    num_events = int(num_events)
    if den_events <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    eff = float(num_events) / float(den_events)
    low = float(ROOT.TEfficiency.ClopperPearson(den_events, num_events, cl, False))
    high = float(ROOT.TEfficiency.ClopperPearson(den_events, num_events, cl, True))
    low = max(0.0, min(1.0, low))
    high = max(0.0, min(1.0, high))
    return eff, low, high, eff - low, high - eff


def weighted_eff_error(den, num, sumw2_pass, sumw2_fail):
    """Approximate symmetric error for weighted subset efficiencies."""
    if den == 0.0:
        return 0.0
    var_eff = (((den - num) ** 2) * sumw2_pass + (num ** 2) * sumw2_fail) / (den ** 4)
    return math.sqrt(max(0.0, var_eff))


def efficiency_from_acc(acc, use_cp=True):
    if acc["den"] <= 0.0:
        return {
            "eff": 0.0, "eff_low": 0.0, "eff_high": 0.0,
            "eff_err_down": 0.0, "eff_err_up": 0.0, "eff_err": 0.0,
        }
    if use_cp and not acc.get("weighted", False):
        eff, low, high, err_down, err_up = cp_interval(acc["n_num_events"], acc["n_ref_events"])
        return {
            "eff": eff, "eff_low": low, "eff_high": high,
            "eff_err_down": err_down, "eff_err_up": err_up,
            "eff_err": max(err_down, err_up),
        }
    eff = acc["num"] / acc["den"]
    err = weighted_eff_error(acc["den"], acc["num"], acc["sumw2_pass"], acc["sumw2_fail"])
    return {
        "eff": eff,
        "eff_low": max(0.0, eff - err),
        "eff_high": min(1.0, eff + err),
        "eff_err_down": err,
        "eff_err_up": err,
        "eff_err": err,
    }


def passes_common_selection(event, args, use_btag):
    if int(getattr(event, args.cutstep_branch)) < args.min_cutstep:
        return False
    if use_btag:
        nb = int(getattr(event, "nBjet"))
        if args.min_nbjet is not None and nb < args.min_nbjet:
            return False
        if args.max_nbjet is not None and nb > args.max_nbjet:
            return False
    return True


def event_key(event):
    return (get_int(event, "run", -1), get_int(event, "luminosityBlock", -1), get_int(event, "event", -1))


def accumulate(files, args, is_mc=False, seen_data_events=None):
    good_files, bad_files = validate_files(files, args.tree)
    chain = make_chain(good_files, args.tree)
    n_entries = int(chain.GetEntries())

    base = {
        "n_files": len(files),
        "n_good_files": len(good_files),
        "n_bad_files": len(bad_files),
        "bad_files": bad_files[:20],
        "n_entries": n_entries,
        "den": 0.0,
        "num": 0.0,
        "sumw2_pass": 0.0,
        "sumw2_fail": 0.0,
        "n_ref_events": 0,
        "n_num_events": 0,
        "n_dup_events": 0,
        "weighted": False,
        "used_weight_branch": "",
    }
    if n_entries == 0:
        base.update(efficiency_from_acc(base, use_cp=True))
        return base

    require_branch(chain, args.analysis_branch)
    require_branch(chain, args.reference_branch)
    require_branch(chain, args.cutstep_branch)

    use_weight = bool(is_mc and args.mc_weight_branch and has_branch(chain, args.mc_weight_branch))
    use_btag = bool(args.min_nbjet is not None or args.max_nbjet is not None)
    if use_btag:
        require_branch(chain, "nBjet")

    use_dedup = bool((not is_mc) and args.deduplicate_data)
    if use_dedup:
        for br in ["run", "luminosityBlock", "event"]:
            require_branch(chain, br)
        if seen_data_events is None:
            seen_data_events = set()

    for iev, event in enumerate(chain):
        if args.verbose and iev and iev % 200000 == 0:
            print("  processed %d / %d" % (iev, n_entries))

        if use_dedup:
            key = event_key(event)
            if key in seen_data_events:
                base["n_dup_events"] += 1
                continue
            seen_data_events.add(key)

        if not passes_common_selection(event, args, use_btag):
            continue
        if not get_bool(event, args.reference_branch):
            continue

        w = get_float(event, args.mc_weight_branch, 1.0) if use_weight else 1.0
        passed = get_bool(event, args.analysis_branch)

        base["den"] += w
        base["n_ref_events"] += 1
        if passed:
            base["num"] += w
            base["sumw2_pass"] += w * w
            base["n_num_events"] += 1
        else:
            base["sumw2_fail"] += w * w

    base["weighted"] = bool(use_weight)
    base["used_weight_branch"] = args.mc_weight_branch if use_weight else ""
    base.update(efficiency_from_acc(base, use_cp=(not use_weight)))
    return base


def ratio_and_error(data, mc):
    ed = data["eff"]
    em = mc["eff"]
    if ed <= 0.0 or em <= 0.0:
        return {
            "sf": 0.0,
            "sf_err": 0.0,
            "sf_err_down": 0.0,
            "sf_err_up": 0.0,
            "sf_low": 0.0,
            "sf_high": 0.0,
        }
    sf = ed / em
    # Up variation: data up, MC down.  Down variation: data down, MC up.
    rel_up2 = 0.0
    rel_down2 = 0.0
    if ed > 0.0:
        rel_up2 += (data["eff_err_up"] / ed) ** 2
        rel_down2 += (data["eff_err_down"] / ed) ** 2
    if em > 0.0:
        rel_up2 += (mc["eff_err_down"] / em) ** 2
        rel_down2 += (mc["eff_err_up"] / em) ** 2
    sf_err_up = sf * math.sqrt(max(0.0, rel_up2))
    sf_err_down = sf * math.sqrt(max(0.0, rel_down2))
    return {
        "sf": sf,
        "sf_err": max(sf_err_down, sf_err_up),
        "sf_err_down": sf_err_down,
        "sf_err_up": sf_err_up,
        "sf_low": max(0.0, sf - sf_err_down),
        "sf_high": sf + sf_err_up,
    }


def make_map_specs(channel):
    same_flavor = "El" if channel in ("ElElEl", "MuElEl") else "Mu"
    dom_flavor = "El" if channel in ("ElElEl", "MuElEl") else "Mu"
    return [
        {
            "name": "leadpt_subleadpt",
            "title": "leading lepton pT vs subleading lepton pT",
            "x_branch": "trigSF2D_leadLep_pt",
            "y_branch": "trigSF2D_subleadLep_pt",
            "x_kind": "pt",
            "y_kind": "pt",
        },
        {
            "name": "leadpt_absleadeta",
            "title": "leading lepton pT vs leading lepton |eta|",
            "x_branch": "trigSF2D_leadLep_pt",
            "y_branch": "trigSF2D_leadLep_absEta",
            "x_kind": "pt",
            "y_kind": "eta",
        },
        {
            "name": "sameflavor_lead_sublead_pt",
            "title": "same-flavor leading pT vs same-flavor subleading pT",
            "x_branch": "trigSF2D_lead%s_pt" % same_flavor,
            "y_branch": "trigSF2D_sublead%s_pt" % same_flavor,
            "x_kind": "pt",
            "y_kind": "pt",
        },
        {
            "name": "dominantflavor_leadpt_subleadpt",
            "title": "dominant-flavor leading pT vs inclusive subleading lepton pT",
            "x_branch": "trigSF2D_lead%s_pt" % dom_flavor,
            "y_branch": "trigSF2D_subleadLep_pt",
            "x_kind": "pt",
            "y_kind": "pt",
        },
    ]


def new_th2(name, title, x_edges, y_edges):
    h = ROOT.TH2D(
        name,
        title,
        len(x_edges) - 1,
        array("d", x_edges),
        len(y_edges) - 1,
        array("d", y_edges),
    )
    h.SetDirectory(0)
    return h


def fill_2d_histograms(files, args, channel, is_mc, map_specs, pt_bins, eta_bins, seen_data_events=None):
    good_files, bad_files = validate_files(files, args.tree)
    chain = make_chain(good_files, args.tree)
    n_entries = int(chain.GetEntries())

    hists = {}
    for spec in map_specs:
        x_edges = pt_bins if spec["x_kind"] == "pt" else eta_bins
        y_edges = pt_bins if spec["y_kind"] == "pt" else eta_bins
        prefix = "mc" if is_mc else "data"
        h_den = new_th2(
            "h2_%s_den_%s_%s" % (prefix, channel, spec["name"]),
            "%s %s denominator;%s;%s" % (channel, spec["title"], spec["x_branch"], spec["y_branch"]),
            x_edges,
            y_edges,
        )
        h_num = new_th2(
            "h2_%s_num_%s_%s" % (prefix, channel, spec["name"]),
            "%s %s numerator;%s;%s" % (channel, spec["title"], spec["x_branch"], spec["y_branch"]),
            x_edges,
            y_edges,
        )
        hists[spec["name"]] = {"den": h_den, "num": h_num, "bad_files": bad_files[:20]}

    if n_entries == 0:
        return hists

    require_branch(chain, args.analysis_branch)
    require_branch(chain, args.reference_branch)
    require_branch(chain, args.cutstep_branch)
    use_weight = bool(is_mc and args.mc_weight_branch and has_branch(chain, args.mc_weight_branch))
    use_btag = bool(args.min_nbjet is not None or args.max_nbjet is not None)
    if use_btag:
        require_branch(chain, "nBjet")

    use_dedup = bool((not is_mc) and args.deduplicate_data)
    if use_dedup:
        for br in ["run", "luminosityBlock", "event"]:
            require_branch(chain, br)
        if seen_data_events is None:
            seen_data_events = set()

    # Only keep map specs whose branches exist.  Missing derived 2D branches are
    # not fatal for the inclusive 1D SF.
    valid_specs = []
    for spec in map_specs:
        if has_branch(chain, spec["x_branch"]) and has_branch(chain, spec["y_branch"]):
            valid_specs.append(spec)
        else:
            print("[WARN] %s: missing 2D branches for map %s: %s, %s" % (
                channel, spec["name"], spec["x_branch"], spec["y_branch"]
            ))

    for event in chain:
        if use_dedup:
            key = event_key(event)
            if key in seen_data_events:
                continue
            seen_data_events.add(key)

        if not passes_common_selection(event, args, use_btag):
            continue
        if not get_bool(event, args.reference_branch):
            continue

        passed = get_bool(event, args.analysis_branch)
        w = get_float(event, args.mc_weight_branch, 1.0) if use_weight else 1.0

        for spec in valid_specs:
            x = get_float(event, spec["x_branch"], -999.0)
            y = get_float(event, spec["y_branch"], -999.0)
            if x < 0.0 or y < 0.0:
                continue
            hists[spec["name"]]["den"].Fill(x, y, w)
            if passed:
                hists[spec["name"]]["num"].Fill(x, y, w)

    return hists


def build_eff_sf_2d(channel, spec, h_den_data, h_num_data, h_den_mc, h_num_mc, args):
    h_eff_data = h_den_data.Clone("h2_eff_data_%s_%s" % (channel, spec["name"]))
    h_eff_mc = h_den_mc.Clone("h2_eff_mc_%s_%s" % (channel, spec["name"]))
    h_sf = h_den_data.Clone("h2_triggerSF_%s_%s" % (channel, spec["name"]))
    for h in [h_eff_data, h_eff_mc, h_sf]:
        h.SetDirectory(0)
        h.Reset("ICES")

    bins_json = []
    nx = h_den_data.GetNbinsX()
    ny = h_den_data.GetNbinsY()
    for ix in range(1, nx + 1):
        for iy in range(1, ny + 1):
            dd = float(h_den_data.GetBinContent(ix, iy))
            nd = float(h_num_data.GetBinContent(ix, iy))
            dm = float(h_den_mc.GetBinContent(ix, iy))
            nm = float(h_num_mc.GetBinContent(ix, iy))

            eff_data = nd / dd if dd > 0.0 else 0.0
            eff_mc = nm / dm if dm > 0.0 else 0.0
            sf = eff_data / eff_mc if eff_data > 0.0 and eff_mc > 0.0 else 0.0

            # For 2D ROOT display, use symmetric binomial/weighted-like errors.
            # Full asymmetric CP intervals are saved in JSON for unweighted data.
            if dd > 0.0:
                _, data_low, data_high, data_down, data_up = cp_interval(int(round(nd)), int(round(dd)))
            else:
                data_low = data_high = data_down = data_up = 0.0
            if dm > 0.0:
                _, mc_low, mc_high, mc_down, mc_up = cp_interval(int(round(nm)), int(round(dm)))
            else:
                mc_low = mc_high = mc_down = mc_up = 0.0

            sf_info = ratio_and_error(
                {"eff": eff_data, "eff_err_down": data_down, "eff_err_up": data_up},
                {"eff": eff_mc, "eff_err_down": mc_down, "eff_err_up": mc_up},
            )

            if dd < args.min_2d_den_data or dm < args.min_2d_den_mc:
                sf_to_store = 0.0
                sf_err_to_store = 0.0
            else:
                sf_to_store = sf
                sf_err_to_store = sf_info["sf_err"]

            h_eff_data.SetBinContent(ix, iy, eff_data)
            h_eff_data.SetBinError(ix, iy, max(data_down, data_up))
            h_eff_mc.SetBinContent(ix, iy, eff_mc)
            h_eff_mc.SetBinError(ix, iy, max(mc_down, mc_up))
            h_sf.SetBinContent(ix, iy, sf_to_store)
            h_sf.SetBinError(ix, iy, sf_err_to_store)

            bins_json.append({
                "ix": ix,
                "iy": iy,
                "x_low": h_den_data.GetXaxis().GetBinLowEdge(ix),
                "x_high": h_den_data.GetXaxis().GetBinUpEdge(ix),
                "y_low": h_den_data.GetYaxis().GetBinLowEdge(iy),
                "y_high": h_den_data.GetYaxis().GetBinUpEdge(iy),
                "den_data": dd,
                "num_data": nd,
                "eff_data": eff_data,
                "eff_data_low": data_low,
                "eff_data_high": data_high,
                "den_mc": dm,
                "num_mc": nm,
                "eff_mc": eff_mc,
                "eff_mc_low": mc_low,
                "eff_mc_high": mc_high,
                "sf": sf_to_store,
                "sf_err": sf_err_to_store,
                "sf_err_down": sf_info["sf_err_down"] if sf_to_store > 0.0 else 0.0,
                "sf_err_up": sf_info["sf_err_up"] if sf_to_store > 0.0 else 0.0,
                "valid": bool(dd >= args.min_2d_den_data and dm >= args.min_2d_den_mc and sf_to_store > 0.0),
            })

    return h_eff_data, h_eff_mc, h_sf, bins_json


def write_hist(root_file, directory, hist, write_name):
    root_file.cd()
    d = root_file.GetDirectory(directory)
    if not d:
        d = root_file.mkdir(directory)
    d.cd()
    hist.Write(write_name, ROOT.TObject.kOverwrite)
    root_file.cd()


def make_json_serializable_maps2d(maps2d):
    """Return a JSON-safe copy of the 2D payload.

    The in-memory maps2d object also contains ROOT TH2 objects used for writing
    the output ROOT file.  ROOT objects are not JSON-serializable, so only plain
    Python metadata and per-bin numeric payloads are copied here.
    """
    serializable = {}
    for channel, channel_maps in maps2d.items():
        serializable[channel] = {}
        for map_name, payload in channel_maps.items():
            serializable[channel][map_name] = {
                "title": payload.get("title", ""),
                "x_branch": payload.get("x_branch", ""),
                "y_branch": payload.get("y_branch", ""),
                "x_bins": list(payload.get("x_bins", [])),
                "y_bins": list(payload.get("y_bins", [])),
                "bins": list(payload.get("bins", [])),
            }
    return serializable


def collect_bad_files(results, maps2d):
    """Collect bad input files seen during 1D and 2D processing."""
    out = []
    seen = set()

    def add(item, context):
        if not item:
            return
        path = item.get("file", "")
        reason = item.get("reason", "unknown")
        key = (path, reason, context)
        if path and key not in seen:
            seen.add(key)
            out.append({"file": path, "reason": reason, "context": context})

    for row in results:
        channel = row.get("channel", "unknown")
        for item in row.get("bad_data_files", []):
            add(item, "%s:data:1d" % channel)
        for item in row.get("bad_mc_files", []):
            add(item, "%s:mc:1d" % channel)

    for channel, channel_maps in maps2d.items():
        for map_name, payload in channel_maps.items():
            for item in payload.get("bad_data_files", []):
                add(item, "%s:%s:data:2d" % (channel, map_name))
            for item in payload.get("bad_mc_files", []):
                add(item, "%s:%s:mc:2d" % (channel, map_name))
    return out


def write_bad_file_list(results, maps2d, path):
    if not path:
        return
    bad = collect_bad_files(results, maps2d)
    with open(path, "w") as handle:
        for item in bad:
            handle.write("%s\t%s\t%s\n" % (item["context"], item["reason"], item["file"]))
    print("Wrote %s" % path)

def write_outputs(results, maps2d, args):
    fieldnames = [
        "channel",
        "sf", "sf_err", "sf_err_down", "sf_err_up", "sf_low", "sf_high",
        "eff_data", "eff_data_err", "eff_data_err_down", "eff_data_err_up", "eff_data_low", "eff_data_high",
        "den_data", "num_data", "n_ref_data", "n_num_data",
        "eff_mc", "eff_mc_err", "eff_mc_err_down", "eff_mc_err_up", "eff_mc_low", "eff_mc_high",
        "den_mc", "num_mc", "n_ref_mc", "n_num_mc",
        "n_data_files", "n_good_data_files", "n_bad_data_files", "n_dup_data",
        "n_mc_files", "n_good_mc_files", "n_bad_mc_files",
    ]
    with open(args.output_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    # Optional flat 2D-bin CSV.  Keep the main CSV compact and stable.
    if args.make_2d:
        stem, ext = os.path.splitext(args.output_csv)
        out_2d_csv = stem + "_2d_bins" + (ext if ext else ".csv")
        fields_2d = [
            "channel", "map", "x_branch", "y_branch", "ix", "iy", "x_low", "x_high", "y_low", "y_high",
            "den_data", "num_data", "eff_data", "eff_data_low", "eff_data_high",
            "den_mc", "num_mc", "eff_mc", "eff_mc_low", "eff_mc_high",
            "sf", "sf_err", "sf_err_down", "sf_err_up", "valid",
        ]
        with open(out_2d_csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields_2d)
            writer.writeheader()
            for channel, channel_maps in maps2d.items():
                for map_name, payload in channel_maps.items():
                    for b in payload["bins"]:
                        row = {"channel": channel, "map": map_name,
                               "x_branch": payload["x_branch"], "y_branch": payload["y_branch"]}
                        row.update(b)
                        writer.writerow({key: row.get(key, "") for key in fields_2d})
        print("Wrote %s" % out_2d_csv)

    # JSON.  Keep a top-level object so 2D maps are not silently lost.
    # Do not put ROOT TH2 objects into the JSON payload.
    json_payload = {
        "summary": results,
        "maps2d": make_json_serializable_maps2d(maps2d) if args.make_2d else {},
    }
    with open(args.output_json, "w") as handle:
        json.dump(json_payload, handle, indent=2, sort_keys=True)

    write_bad_file_list(results, maps2d, args.bad_file_list)

    root_file = ROOT.TFile(args.output_root, "RECREATE")
    nbins = len(results)
    h_sf = ROOT.TH1D("h_triggerSF_2024", "2024 trigger SF;channel;SF", nbins, 0, nbins)
    h_eff_data = ROOT.TH1D("h_triggerEffData_2024", "2024 data trigger efficiency;channel;efficiency", nbins, 0, nbins)
    h_eff_mc = ROOT.TH1D("h_triggerEffMC_2024", "2024 MC trigger efficiency;channel;efficiency", nbins, 0, nbins)
    h_den_data = ROOT.TH1D("h_triggerDenData_2024", "2024 data denominator;channel;events", nbins, 0, nbins)
    h_num_data = ROOT.TH1D("h_triggerNumData_2024", "2024 data numerator;channel;events", nbins, 0, nbins)
    h_den_mc = ROOT.TH1D("h_triggerDenMC_2024", "2024 MC denominator;channel;weighted events", nbins, 0, nbins)
    h_num_mc = ROOT.TH1D("h_triggerNumMC_2024", "2024 MC numerator;channel;weighted events", nbins, 0, nbins)
    for h in [h_sf, h_eff_data, h_eff_mc, h_den_data, h_num_data, h_den_mc, h_num_mc]:
        h.SetDirectory(0)

    for ibin, row in enumerate(results, start=1):
        channel = row["channel"]
        for hist in [h_sf, h_eff_data, h_eff_mc, h_den_data, h_num_data, h_den_mc, h_num_mc]:
            hist.GetXaxis().SetBinLabel(ibin, channel)

        h_sf.SetBinContent(ibin, row["sf"])
        h_sf.SetBinError(ibin, row["sf_err"])
        h_eff_data.SetBinContent(ibin, row["eff_data"])
        h_eff_data.SetBinError(ibin, row["eff_data_err"])
        h_eff_mc.SetBinContent(ibin, row["eff_mc"])
        h_eff_mc.SetBinError(ibin, row["eff_mc_err"])
        h_den_data.SetBinContent(ibin, row["den_data"])
        h_num_data.SetBinContent(ibin, row["num_data"])
        h_den_mc.SetBinContent(ibin, row["den_mc"])
        h_num_mc.SetBinContent(ibin, row["num_mc"])

    root_file.cd()
    for hist in [h_sf, h_eff_data, h_eff_mc, h_den_data, h_num_data, h_den_mc, h_num_mc]:
        hist.Write(hist.GetName(), ROOT.TObject.kOverwrite)

    if args.make_2d:
        for channel, channel_maps in maps2d.items():
            for map_name, payload in channel_maps.items():
                directory = "%s/%s" % (channel, map_name)
                for key, write_name in [
                    ("h_den_data", "h_den_data"),
                    ("h_num_data", "h_num_data"),
                    ("h_eff_data", "h_eff_data"),
                    ("h_den_mc", "h_den_mc"),
                    ("h_num_mc", "h_num_mc"),
                    ("h_eff_mc", "h_eff_mc"),
                    ("h_triggerSF", "h_triggerSF"),
                ]:
                    write_hist(root_file, directory, payload[key], write_name)

    root_file.Close()


def build_2d_for_channel(channel, data_files, mc_files, args):
    pt_bins = parse_binning(args.pt_bins)
    eta_bins = parse_binning(args.eta_bins)
    specs = make_map_specs(channel)

    data_seen_2d = set() if args.deduplicate_data else None
    data_hists = fill_2d_histograms(data_files, args, channel, False, specs, pt_bins, eta_bins, data_seen_2d)
    mc_hists = fill_2d_histograms(mc_files, args, channel, True, specs, pt_bins, eta_bins, None)

    out = {}
    for spec in specs:
        name = spec["name"]
        hd = data_hists[name]["den"]
        hn = data_hists[name]["num"]
        hdm = mc_hists[name]["den"]
        hnm = mc_hists[name]["num"]
        h_eff_data, h_eff_mc, h_sf, bins_json = build_eff_sf_2d(channel, spec, hd, hn, hdm, hnm, args)
        out[name] = {
            "title": spec["title"],
            "x_branch": spec["x_branch"],
            "y_branch": spec["y_branch"],
            "x_bins": pt_bins if spec["x_kind"] == "pt" else eta_bins,
            "y_bins": pt_bins if spec["y_kind"] == "pt" else eta_bins,
            "bad_data_files": data_hists[name].get("bad_files", []),
            "bad_mc_files": mc_hists[name].get("bad_files", []),
            "h_den_data": hd,
            "h_num_data": hn,
            "h_eff_data": h_eff_data,
            "h_den_mc": hdm,
            "h_num_mc": hnm,
            "h_eff_mc": h_eff_mc,
            "h_triggerSF": h_sf,
            "bins": bins_json,
        }
    return out


def main():
    args = parse_args()
    results = []
    maps2d = {}
    for channel in args.channels:
        data_files_all = expand_files(args.data_glob, channel)
        mc_files_all = expand_files(args.mc_glob, channel)

        print("[%s] data files: %d" % (channel, len(data_files_all)))
        print("[%s] MC files  : %d" % (channel, len(mc_files_all)))

        if not data_files_all:
            print("[WARN] No data files for %s; result will be zero." % channel)
        if not mc_files_all:
            print("[WARN] No MC files for %s; result will be zero." % channel)

        seen_data_channel = set() if args.deduplicate_data else None
        data_acc = accumulate(data_files_all, args, is_mc=False, seen_data_events=seen_data_channel) if data_files_all else {
            "n_files": 0, "n_good_files": 0, "n_bad_files": 0, "bad_files": [], "n_entries": 0,
            "n_ref_events": 0, "n_num_events": 0, "n_dup_events": 0,
            "den": 0.0, "num": 0.0, "sumw2_pass": 0.0, "sumw2_fail": 0.0, "weighted": False,
            "eff": 0.0, "eff_low": 0.0, "eff_high": 0.0,
            "eff_err": 0.0, "eff_err_down": 0.0, "eff_err_up": 0.0,
        }
        mc_acc = accumulate(mc_files_all, args, is_mc=True, seen_data_events=None) if mc_files_all else {
            "n_files": 0, "n_good_files": 0, "n_bad_files": 0, "bad_files": [], "n_entries": 0,
            "n_ref_events": 0, "n_num_events": 0, "n_dup_events": 0,
            "den": 0.0, "num": 0.0, "sumw2_pass": 0.0, "sumw2_fail": 0.0, "weighted": False,
            "eff": 0.0, "eff_low": 0.0, "eff_high": 0.0,
            "eff_err": 0.0, "eff_err_down": 0.0, "eff_err_up": 0.0,
        }

        sf_info = ratio_and_error(data_acc, mc_acc)
        row = {
            "channel": channel,
            **sf_info,
            "eff_data": data_acc["eff"],
            "eff_data_err": data_acc["eff_err"],
            "eff_data_err_down": data_acc["eff_err_down"],
            "eff_data_err_up": data_acc["eff_err_up"],
            "eff_data_low": data_acc["eff_low"],
            "eff_data_high": data_acc["eff_high"],
            "den_data": data_acc["den"],
            "num_data": data_acc["num"],
            "n_ref_data": data_acc.get("n_ref_events", 0),
            "n_num_data": data_acc.get("n_num_events", 0),
            "eff_mc": mc_acc["eff"],
            "eff_mc_err": mc_acc["eff_err"],
            "eff_mc_err_down": mc_acc["eff_err_down"],
            "eff_mc_err_up": mc_acc["eff_err_up"],
            "eff_mc_low": mc_acc["eff_low"],
            "eff_mc_high": mc_acc["eff_high"],
            "den_mc": mc_acc["den"],
            "num_mc": mc_acc["num"],
            "n_ref_mc": mc_acc.get("n_ref_events", 0),
            "n_num_mc": mc_acc.get("n_num_events", 0),
            "n_data_files": data_acc["n_files"],
            "n_good_data_files": data_acc["n_good_files"],
            "n_bad_data_files": data_acc["n_bad_files"],
            "bad_data_files": data_acc.get("bad_files", []),
            "n_dup_data": data_acc.get("n_dup_events", 0),
            "n_mc_files": mc_acc["n_files"],
            "n_good_mc_files": mc_acc["n_good_files"],
            "n_bad_mc_files": mc_acc["n_bad_files"],
            "bad_mc_files": mc_acc.get("bad_files", []),
        }
        results.append(row)
        print(
            "[%s] eff_data=%.6f -%.6f/+%.6f, eff_mc=%.6f -%.6f/+%.6f, "
            "SF=%.6f -%.6f/+%.6f, skipped data duplicates=%d, bad data files=%d, bad MC files=%d" % (
                channel,
                row["eff_data"], row["eff_data_err_down"], row["eff_data_err_up"],
                row["eff_mc"], row["eff_mc_err_down"], row["eff_mc_err_up"],
                row["sf"], row["sf_err_down"], row["sf_err_up"],
                row["n_dup_data"], row["n_bad_data_files"], row["n_bad_mc_files"],
            )
        )
        if row["bad_data_files"]:
            print("[%s] first bad data file: %s (%s)" % (
                channel, row["bad_data_files"][0]["file"], row["bad_data_files"][0]["reason"]
            ))

        if args.make_2d:
            maps2d[channel] = build_2d_for_channel(channel, data_files_all, mc_files_all, args)

    write_outputs(results, maps2d, args)
    print("Wrote %s" % args.output_csv)
    print("Wrote %s" % args.output_json)
    print("Wrote %s" % args.output_root)


if __name__ == "__main__":
    main()
