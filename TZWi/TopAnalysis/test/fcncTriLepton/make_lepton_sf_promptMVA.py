#!/usr/bin/env python3
"""
Build Run-3 2024 promptMVA lepton-selection SF maps from dedicated SF NTuples.

Main corrections relative to the original counting script:
  1. Reapply the analysis denominator/numerator definitions explicitly.
  2. Apply configurable WZ-enrichment cuts using SF_met_pt and SF_w_mt.
  3. Extract the Data efficiency after subtracting normalized background MC
     separately in the pass and fail samples.
  4. Use prompt WZ probes only for the MC efficiency. Optionally subtract the
     nonprompt-probe component of WZ from Data as another background.
  5. Propagate pass/fail statistical errors; do not use a binomial formula on
     background-subtracted yields.

The output map convention is x = probe pT and y = signed eta (SC eta for e).
"""

import argparse
import array
import csv
import glob
import math
import os
import re
import sys

import ROOT
import yaml

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.TH1.AddDirectory(False)

#DEFAULT_PT_BINS = "10,15,20,25,30,35,40,50,70,100,200"
DEFAULT_PT_BINS = "10,20,30,40,50,100,200"
DEFAULT_ETA_BINS = {
    ("mu", "abs"): "0,0.9,1.2,2.1,2.4",
    ("el", "abs"): "0,0.8,1.479,2.0,2.5",
    ("mu", "signed"): "-2.4,-2.1,-1.2,-0.9,0,0.9,1.2,2.1,2.4",
    ("el", "signed"): "-2.5,-2.0,-1.479,-0.8,0,0.8,1.479,2.0,2.5",
}


def parse_bins(text):
    values = [float(x) for x in str(text).split(",") if str(x).strip()]
    if len(values) < 2:
        raise ValueError("At least two bin edges are required: %s" % text)
    if any(values[i] >= values[i + 1] for i in range(len(values) - 1)):
        raise ValueError("Bin edges must be strictly increasing: %s" % text)
    return array.array("d", values)


def as_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def expand_patterns(patterns):
    matched = []
    for raw in as_list(patterns):
        pattern = os.path.expanduser(os.path.expandvars(str(raw)))
        paths = sorted(glob.glob(pattern))
        if not paths and os.path.isfile(pattern):
            paths = [pattern]
        matched.extend(paths)
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(matched))


def build_chain(patterns, tree_name="Events", label="sample"):
    files = expand_patterns(patterns)
    if not files:
        raise RuntimeError("No ROOT files matched for %s: %s" % (label, patterns))
    chain = ROOT.TChain(tree_name)
    for path in files:
        chain.Add(path)
    return chain, files


def make_hist(name, title, pt_bins, eta_bins, eta_title):
    hist = ROOT.TH2D(
        name,
        title,
        len(pt_bins) - 1,
        pt_bins,
        len(eta_bins) - 1,
        eta_bins,
    )
    hist.Sumw2()
    hist.GetXaxis().SetTitle("p_{T} [GeV]")
    hist.GetYaxis().SetTitle(eta_title)
    return hist


def clone_reset(hist, name, title=None):
    out = hist.Clone(name)
    out.Reset("ICES")
    if title is not None:
        out.SetTitle(title)
    out.SetDirectory(0)
    return out


def clamp_to_axis(value, axis):
    if value < axis.GetXmin():
        return axis.GetXmin() + 1.0e-6
    if value >= axis.GetXmax():
        return axis.GetXmax() - 1.0e-6
    return value


def has_branch(chain, name):
    return bool(chain.GetBranch(name))


def require_branches(chain, names, label):
    missing = [name for name in names if not has_branch(chain, name)]
    if missing:
        raise RuntimeError(
            "%s is missing required branches: %s" % (label, ", ".join(missing))
        )


def event_key(ev):
    try:
        return (int(ev.SF_run), int(ev.SF_luminosityBlock), int(ev.SF_event))
    except Exception:
        return None


def choose_eta(ev, flavor, eta_mode, eta_var):
    if eta_var == "auto":
        value = float(ev.SF_probe_scEta if flavor == "el" else ev.SF_probe_eta)
    elif eta_var == "eta":
        value = float(ev.SF_probe_eta)
    elif eta_var == "scEta":
        value = float(ev.SF_probe_scEta)
    else:
        raise RuntimeError("Unsupported eta_var: %s" % eta_var)
    return abs(value) if eta_mode == "abs" else value


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")


def sample_scale(sample, lumi_pb, require_absolute=True):
    """Return the multiplicative factor applied on top of the event weight branch."""
    if "scale" in sample and sample["scale"] is not None:
        return float(sample["scale"])

    xsec = sample.get("xsec_pb")
    sumw = sample.get("sum_gen_weight")
    if xsec is None or sumw is None:
        if require_absolute:
            raise RuntimeError(
                "Sample %s needs either 'scale' or both 'xsec_pb' and "
                "'sum_gen_weight' for absolute normalization."
                % sample.get("name", "unnamed")
            )
        return 1.0

    sumw = float(sumw)
    if not math.isfinite(sumw) or abs(sumw) <= 0.0:
        raise RuntimeError("Invalid sum_gen_weight for %s: %r" % (sample.get("name"), sumw))
    return float(lumi_pb) * float(xsec) / sumw


def get_selection(config, args):
    cfg = dict(config.get("selection", {}) or {})

    def pick(name, cli_value, default):
        return cli_value if cli_value is not None else cfg.get(name, default)

    return {
        "min_probe_pt": float(pick("min_probe_pt", args.min_probe_pt, 10.0)),
        "min_met": float(pick("min_met", args.min_met, 30.0)),
        "min_w_mt": float(pick("min_w_mt", args.min_w_mt, 30.0)),
        "max_w_mt": pick("max_w_mt", args.max_w_mt, None),
        "min_electron_cutbased": int(pick("min_electron_cutbased", args.min_electron_cutbased, 3)),
        "max_muon_pf_rel_iso04": float(pick("max_muon_pf_rel_iso04", args.max_muon_pf_rel_iso04, 0.40)),
        "electron_prompt_mva": float(pick("electron_prompt_mva", args.electron_prompt_mva, 0.90)),
        "muon_prompt_mva": float(pick("muon_prompt_mva", args.muon_prompt_mva, 0.70)),
        "min_electron_tightcharge": int(pick("min_electron_tightcharge", args.min_electron_tightcharge, 2)),
    }


def passes_common_selection(ev, flavor, selection, strict_object_check=True):
    target_flavor = 11 if flavor == "el" else 13
    if int(ev.SF_probe_flavor) != target_flavor:
        return False
    if not bool(ev.SF_probe_passRelax):
        return False
    if float(ev.SF_probe_pt) < selection["min_probe_pt"]:
        return False
    if float(ev.SF_met_pt) < selection["min_met"]:
        return False
    if float(ev.SF_w_mt) < selection["min_w_mt"]:
        return False
    if selection["max_w_mt"] is not None and float(ev.SF_w_mt) >= float(selection["max_w_mt"]):
        return False

    if flavor == "el":
        if int(ev.SF_probe_cutBased) < selection["min_electron_cutbased"]:
            return False
        if not bool(ev.SF_probe_convVeto):
            return False
    elif strict_object_check:
        # This branch is present in the fixed NTuple producer. Requiring it here
        # prevents silently using the older denominator definition.
        if float(ev.SF_probe_pfRelIso04_all) >= selection["max_muon_pf_rel_iso04"]:
            return False
        if not bool(ev.SF_probe_mediumId):
            return False
    return True


def passes_strict(ev, flavor, selection):
    if flavor == "el":
        return (
            int(ev.SF_probe_tightCharge) >= selection["min_electron_tightcharge"]
            and float(ev.SF_probe_promptMVA) > selection["electron_prompt_mva"]
        )
    return float(ev.SF_probe_promptMVA) > selection["muon_prompt_mva"]


def fill_pass_fail(
    chain,
    flavor,
    pt_bins,
    eta_bins,
    eta_mode,
    eta_var,
    selection,
    label,
    weight_branch="",
    scale=1.0,
    prompt_requirement=None,
    deduplicate=False,
    max_events=-1,
    strict_object_check=True,
):
    eta_title = "|#eta|" if eta_mode == "abs" else "#eta"
    sname = safe_name(label)
    h_pass = make_hist("h_pass_%s" % sname, "pass %s" % label, pt_bins, eta_bins, eta_title)
    h_fail = make_hist("h_fail_%s" % sname, "fail %s" % label, pt_bins, eta_bins, eta_title)
    h_pass_count = make_hist("h_pass_count_%s" % sname, "unweighted pass %s" % label, pt_bins, eta_bins, eta_title)
    h_fail_count = make_hist("h_fail_count_%s" % sname, "unweighted fail %s" % label, pt_bins, eta_bins, eta_title)

    required = [
        "SF_probe_flavor", "SF_probe_passRelax", "SF_probe_pt",
        "SF_probe_eta", "SF_probe_scEta", "SF_probe_promptMVA",
        "SF_met_pt", "SF_w_mt",
    ]
    if flavor == "el":
        required += ["SF_probe_cutBased", "SF_probe_convVeto", "SF_probe_tightCharge"]
    else:
        required += ["SF_probe_mediumId"]
        if strict_object_check:
            required += ["SF_probe_pfRelIso04_all"]
    if prompt_requirement is not None:
        required += ["SF_probe_isPromptMC"]
    if weight_branch:
        required += [weight_branch]
    require_branches(chain, required, label)

    seen = set()
    info = {
        "label": label,
        "n_total": 0,
        "n_selected": 0,
        "n_pass": 0,
        "n_fail": 0,
        "n_duplicate_skipped": 0,
        "n_prompt_requirement_skipped": 0,
        "weight_branch": weight_branch or "1.0",
        "scale": float(scale),
    }

    for ev in chain:
        info["n_total"] += 1
        if max_events > 0 and info["n_total"] > max_events:
            break

        if not passes_common_selection(ev, flavor, selection, strict_object_check):
            continue

        if prompt_requirement is not None:
            is_prompt = bool(ev.SF_probe_isPromptMC)
            if is_prompt != bool(prompt_requirement):
                info["n_prompt_requirement_skipped"] += 1
                continue

        if deduplicate:
            key = event_key(ev)
            if key is not None:
                key = key + (11 if flavor == "el" else 13,)
                if key in seen:
                    info["n_duplicate_skipped"] += 1
                    continue
                seen.add(key)

        event_weight = float(scale)
        if weight_branch:
            event_weight *= float(getattr(ev, weight_branch))
        if not math.isfinite(event_weight):
            continue

        pt = clamp_to_axis(float(ev.SF_probe_pt), h_pass.GetXaxis())
        eta = clamp_to_axis(choose_eta(ev, flavor, eta_mode, eta_var), h_pass.GetYaxis())
        passed = passes_strict(ev, flavor, selection)

        info["n_selected"] += 1
        if passed:
            h_pass.Fill(pt, eta, event_weight)
            h_pass_count.Fill(pt, eta, 1.0)
            info["n_pass"] += 1
        else:
            h_fail.Fill(pt, eta, event_weight)
            h_fail_count.Fill(pt, eta, 1.0)
            info["n_fail"] += 1

    return h_pass, h_fail, h_pass_count, h_fail_count, info


def make_den(pass_hist, fail_hist, name, title):
    den = pass_hist.Clone(name)
    den.SetDirectory(0)
    den.SetTitle(title)
    den.Add(fail_hist)
    return den


def efficiency_from_pass_fail(pass_hist, fail_hist, name, title):
    eff = clone_reset(pass_hist, name, title)
    valid = clone_reset(pass_hist, name + "_valid", title + " valid mask")

    for ix in range(1, pass_hist.GetNbinsX() + 1):
        for iy in range(1, pass_hist.GetNbinsY() + 1):
            p = pass_hist.GetBinContent(ix, iy)
            f = fail_hist.GetBinContent(ix, iy)
            vp = pass_hist.GetBinError(ix, iy) ** 2
            vf = fail_hist.GetBinError(ix, iy) ** 2
            d = p + f

            if p < 0.0 or f < 0.0 or d <= 0.0:
                eff.SetBinContent(ix, iy, 0.0)
                eff.SetBinError(ix, iy, 0.0)
                valid.SetBinContent(ix, iy, 0.0)
                continue

            value = p / d
            # p and f are disjoint pass/fail yields, so they can be treated as
            # independent in this approximation, including after subtraction.
            variance = (f * f * vp + p * p * vf) / (d ** 4)
            eff.SetBinContent(ix, iy, value)
            eff.SetBinError(ix, iy, math.sqrt(max(0.0, variance)))
            valid.SetBinContent(ix, iy, 1.0)

    return eff, valid


def build_sf(
    data_pass,
    data_fail,
    mc_pass,
    mc_fail,
    empty_bin_sf,
):
    h_eff_data, h_valid_data = efficiency_from_pass_fail(
        data_pass, data_fail, "h2_eff_data", "background-subtracted Data efficiency"
    )
    h_eff_mc, h_valid_mc = efficiency_from_pass_fail(
        mc_pass, mc_fail, "h2_eff_mc", "prompt WZ MC efficiency"
    )

    h_sf = clone_reset(h_eff_data, "h2_sf", "promptMVA lepton selection scale factor")
    h_up = clone_reset(h_eff_data, "h2_sf_statUp", "scale factor statistical up")
    h_down = clone_reset(h_eff_data, "h2_sf_statDown", "scale factor statistical down")
    h_valid = clone_reset(h_eff_data, "h2_sf_valid", "scale factor valid mask")

    rows = []
    for ix in range(1, h_sf.GetNbinsX() + 1):
        for iy in range(1, h_sf.GetNbinsY() + 1):
            dp = data_pass.GetBinContent(ix, iy)
            df = data_fail.GetBinContent(ix, iy)
            mp = mc_pass.GetBinContent(ix, iy)
            mf = mc_fail.GetBinContent(ix, iy)
            ed = h_eff_data.GetBinContent(ix, iy)
            em = h_eff_mc.GetBinContent(ix, iy)
            sed = h_eff_data.GetBinError(ix, iy)
            sem = h_eff_mc.GetBinError(ix, iy)

            data_ok = h_valid_data.GetBinContent(ix, iy) > 0.5
            mc_ok = h_valid_mc.GetBinContent(ix, iy) > 0.5 and em > 0.0
            if not data_ok:
                sf = float(empty_bin_sf)
                sf_err = 0.0
                status = "invalid_data_after_subtraction"
            elif not mc_ok:
                sf = float(empty_bin_sf)
                sf_err = 0.0
                status = "invalid_or_empty_prompt_mc"
            else:
                sf = ed / em
                variance = (sed / em) ** 2 + ((ed * sem) / (em * em)) ** 2
                sf_err = math.sqrt(max(0.0, variance))
                status = "ok"

            h_sf.SetBinContent(ix, iy, sf)
            h_sf.SetBinError(ix, iy, sf_err)
            h_up.SetBinContent(ix, iy, sf + sf_err)
            h_down.SetBinContent(ix, iy, max(0.0, sf - sf_err))
            h_valid.SetBinContent(ix, iy, 1.0 if status == "ok" else 0.0)

            rows.append({
                "pt_low": h_sf.GetXaxis().GetBinLowEdge(ix),
                "pt_high": h_sf.GetXaxis().GetBinUpEdge(ix),
                "eta_low": h_sf.GetYaxis().GetBinLowEdge(iy),
                "eta_high": h_sf.GetYaxis().GetBinUpEdge(iy),
                "data_pass_sub": dp,
                "data_fail_sub": df,
                "data_den_sub": dp + df,
                "mc_pass_prompt": mp,
                "mc_fail_prompt": mf,
                "mc_den_prompt": mp + mf,
                "eff_data": ed,
                "eff_mc": em,
                "err_eff_data": sed,
                "err_eff_mc": sem,
                "sf": sf,
                "sf_err": sf_err,
                "sf_stat_up": sf + sf_err,
                "sf_stat_down": max(0.0, sf - sf_err),
                "status": status,
            })

    return h_eff_data, h_eff_mc, h_sf, h_up, h_down, h_valid, rows


def write_csv(path, rows):
    if not path:
        return
    fields = [
        "pt_low", "pt_high", "eta_low", "eta_high",
        "data_pass_raw", "data_fail_raw", "bkg_pass", "bkg_fail",
        "data_pass_sub", "data_fail_sub", "data_den_sub",
        "mc_pass_prompt", "mc_fail_prompt", "mc_den_prompt",
        "eff_data", "eff_mc", "err_eff_data", "err_eff_mc",
        "sf", "sf_err", "sf_stat_up", "sf_stat_down", "status",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Build background-subtracted Run-3 2024 promptMVA lepton SF maps."
    )
    parser.add_argument("--config", required=True, help="YAML sample/normalization configuration.")
    parser.add_argument("--flavor", choices=["el", "mu"], required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--csv", default="")
    parser.add_argument("--tree", default="Events")
    parser.add_argument("--pt-bins", default=DEFAULT_PT_BINS)
    parser.add_argument("--eta-mode", choices=["signed", "abs"], default="signed")
    parser.add_argument("--eta-var", choices=["auto", "eta", "scEta"], default="auto")
    parser.add_argument("--eta-bins", default="")
    parser.add_argument("--empty-bin-sf", type=float, default=1.0)
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("--allow-no-backgrounds", action="store_true")
    parser.add_argument("--allow-old-muon-ntuple", action="store_true",
                        help="Diagnostic only: do not require the stored muon PF isolation branch.")

    # Optional command-line overrides of config selection values.
    parser.add_argument("--min-probe-pt", type=float, default=None)
    parser.add_argument("--min-met", type=float, default=None)
    parser.add_argument("--min-w-mt", type=float, default=None)
    parser.add_argument("--max-w-mt", type=float, default=None)
    parser.add_argument("--min-electron-cutbased", type=int, default=None)
    parser.add_argument("--max-muon-pf-rel-iso04", type=float, default=None)
    parser.add_argument("--electron-prompt-mva", type=float, default=None)
    parser.add_argument("--muon-prompt-mva", type=float, default=None)
    parser.add_argument("--min-electron-tightcharge", type=int, default=None)
    args = parser.parse_args()

    with open(args.config) as handle:
        config = yaml.safe_load(handle) or {}

    lumi_pb = float(config.get("lumi_pb", 0.0))
    if lumi_pb <= 0.0:
        raise RuntimeError("config must provide a positive lumi_pb")

    data_cfg = config.get("data") or {}
    signal_cfg = config.get("signal") or {}
    backgrounds_cfg = list(config.get("backgrounds") or [])
    if not backgrounds_cfg and not args.allow_no_backgrounds:
        raise RuntimeError(
            "No backgrounds are configured. Raw Data / prompt-only WZ is the biased method "
            "that produced the small SFs. Configure backgrounds or use --allow-no-backgrounds "
            "only for a diagnostic run."
        )

    selection = get_selection(config, args)
    pt_bins = parse_bins(args.pt_bins)
    eta_bins_text = args.eta_bins or DEFAULT_ETA_BINS[(args.flavor, args.eta_mode)]
    eta_bins = parse_bins(eta_bins_text)

    print("[make_lepton_sf_fixed] flavor     =", args.flavor)
    print("[make_lepton_sf_fixed] lumi_pb    =", lumi_pb)
    print("[make_lepton_sf_fixed] selection =", selection)
    print("[make_lepton_sf_fixed] pt bins    =", list(pt_bins))
    print("[make_lepton_sf_fixed] eta bins   =", list(eta_bins))

    data_chain, data_files = build_chain(data_cfg.get("inputs"), args.tree, "data")
    signal_chain, signal_files = build_chain(signal_cfg.get("inputs"), args.tree, "signal")
    signal_weight = str(signal_cfg.get("weight_branch", "SF_genWeight"))
    signal_norm = sample_scale(signal_cfg, lumi_pb, require_absolute=True)

    print("[make_lepton_sf_fixed] data files   =", len(data_files))
    print("[make_lepton_sf_fixed] signal files =", len(signal_files))
    print("[make_lepton_sf_fixed] signal norm  =", signal_norm)

    strict_object_check = not args.allow_old_muon_ntuple
    data_pass_raw, data_fail_raw, data_pass_count, data_fail_count, info_data = fill_pass_fail(
        data_chain,
        args.flavor,
        pt_bins,
        eta_bins,
        args.eta_mode,
        args.eta_var,
        selection,
        "data_raw",
        weight_branch=str(data_cfg.get("weight_branch", "")),
        scale=float(data_cfg.get("scale", 1.0)),
        prompt_requirement=None,
        deduplicate=bool(data_cfg.get("deduplicate", True)),
        max_events=args.max_events,
        strict_object_check=strict_object_check,
    )

    signal_pass, signal_fail, signal_pass_count, signal_fail_count, info_signal = fill_pass_fail(
        signal_chain,
        args.flavor,
        pt_bins,
        eta_bins,
        args.eta_mode,
        args.eta_var,
        selection,
        "signal_prompt",
        weight_branch=signal_weight,
        scale=signal_norm,
        prompt_requirement=True,
        deduplicate=False,
        max_events=args.max_events,
        strict_object_check=strict_object_check,
    )

    bkg_pass_total = clone_reset(data_pass_raw, "h_pass_background_total", "total background pass")
    bkg_fail_total = clone_reset(data_fail_raw, "h_fail_background_total", "total background fail")
    background_objects = []
    infos_background = []

    if bool(signal_cfg.get("subtract_nonprompt_component", True)):
        sp, sf, spc, sfc, info = fill_pass_fail(
            signal_chain,
            args.flavor,
            pt_bins,
            eta_bins,
            args.eta_mode,
            args.eta_var,
            selection,
            "signal_nonprompt_component",
            weight_branch=signal_weight,
            scale=signal_norm,
            prompt_requirement=False,
            deduplicate=False,
            max_events=args.max_events,
            strict_object_check=strict_object_check,
        )
        bkg_pass_total.Add(sp)
        bkg_fail_total.Add(sf)
        background_objects.extend([sp, sf, spc, sfc])
        infos_background.append(info)

    for sample in backgrounds_cfg:
        name = str(sample.get("name", "background"))
        chain, files = build_chain(sample.get("inputs"), args.tree, name)
        weight_branch = str(sample.get("weight_branch", "SF_genWeight"))
        norm = sample_scale(sample, lumi_pb, require_absolute=True)
        hp, hf, hpc, hfc, info = fill_pass_fail(
            chain,
            args.flavor,
            pt_bins,
            eta_bins,
            args.eta_mode,
            args.eta_var,
            selection,
            "background_%s" % name,
            weight_branch=weight_branch,
            scale=norm,
            prompt_requirement=None,
            deduplicate=False,
            max_events=args.max_events,
            strict_object_check=strict_object_check,
        )
        info["n_files"] = len(files)
        bkg_pass_total.Add(hp)
        bkg_fail_total.Add(hf)
        background_objects.extend([hp, hf, hpc, hfc])
        infos_background.append(info)
        print("[make_lepton_sf_fixed] background %-35s files=%5d norm=%g" % (name, len(files), norm))

    data_pass_sub = data_pass_raw.Clone("h_pass_data_subtracted")
    data_fail_sub = data_fail_raw.Clone("h_fail_data_subtracted")
    data_pass_sub.SetDirectory(0)
    data_fail_sub.SetDirectory(0)
    data_pass_sub.Add(bkg_pass_total, -1.0)
    data_fail_sub.Add(bkg_fail_total, -1.0)

    # Compatibility names used by the previous output format.
    h_num_data = data_pass_sub.Clone("h_num_data")
    h_den_data = make_den(data_pass_sub, data_fail_sub, "h_den_data", "background-subtracted Data denominator")
    h_num_mc = signal_pass.Clone("h_num_mc")
    h_den_mc = make_den(signal_pass, signal_fail, "h_den_mc", "prompt WZ MC denominator")

    h_eff_data, h_eff_mc, h_sf, h_up, h_down, h_valid, rows = build_sf(
        data_pass_sub,
        data_fail_sub,
        signal_pass,
        signal_fail,
        args.empty_bin_sf,
    )

    for row_ix, row in enumerate(rows):
        ix = row_ix // h_sf.GetNbinsY() + 1
        iy = row_ix % h_sf.GetNbinsY() + 1
        row["data_pass_raw"] = data_pass_raw.GetBinContent(ix, iy)
        row["data_fail_raw"] = data_fail_raw.GetBinContent(ix, iy)
        row["bkg_pass"] = bkg_pass_total.GetBinContent(ix, iy)
        row["bkg_fail"] = bkg_fail_total.GetBinContent(ix, iy)

    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    fout = ROOT.TFile(args.out, "RECREATE")
    objects = [
        data_pass_raw, data_fail_raw, data_pass_count, data_fail_count,
        bkg_pass_total, bkg_fail_total,
        data_pass_sub, data_fail_sub,
        signal_pass, signal_fail, signal_pass_count, signal_fail_count,
        h_num_data, h_den_data, h_num_mc, h_den_mc,
        h_eff_data, h_eff_mc, h_sf, h_up, h_down, h_valid,
    ] + background_objects
    for obj in objects:
        obj.Write()

    metadata = ROOT.TNamed(
        "metadata",
        "flavor=%s; eta_mode=%s; eta_var=%s; config=%s; selection=%s"
        % (args.flavor, args.eta_mode, args.eta_var, os.path.abspath(args.config), selection),
    )
    metadata.Write()
    fout.Close()

    csv_path = args.csv or os.path.splitext(args.out)[0] + ".csv"
    write_csv(csv_path, rows)

    print("[make_lepton_sf_fixed] data summary   =", info_data)
    print("[make_lepton_sf_fixed] signal summary =", info_signal)
    for info in infos_background:
        print("[make_lepton_sf_fixed] bkg summary    =", info)
    print("[make_lepton_sf_fixed] wrote ROOT     =", args.out)
    print("[make_lepton_sf_fixed] wrote CSV      =", csv_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("[make_lepton_sf_fixed] ERROR:", exc, file=sys.stderr)
        raise
