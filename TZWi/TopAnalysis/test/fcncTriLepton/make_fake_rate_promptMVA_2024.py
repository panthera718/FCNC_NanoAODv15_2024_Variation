#!/usr/bin/env python3
import argparse
import array
import csv
import glob
import math
import os

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

DEFAULT_PT_BINS = "10,15,20,25,30,40,60,100,200"
DEFAULT_ETA_BINS_MU = "0,0.9,1.2,2.1,2.4"
DEFAULT_ETA_BINS_EL = "0,0.8,1.479,2.0,2.5"
CL = 0.682689492137


def parse_bins(text):
    values = [float(x) for x in text.split(",") if x.strip()]
    if len(values) < 2:
        raise ValueError("At least two bin edges are required")
    if any(values[i] >= values[i + 1] for i in range(len(values) - 1)):
        raise ValueError("Bin edges must be strictly increasing: %s" % text)
    return array.array("d", values)


def build_chain(inputs, tree_name):
    chain = ROOT.TChain(tree_name)
    files = []
    for pattern in inputs:
        matches = sorted(glob.glob(pattern))
        files.extend(matches)
    files = list(dict.fromkeys(files))
    if not files:
        raise RuntimeError("No input ROOT files matched: %s" % inputs)
    for path in files:
        chain.Add(path)
    return chain, files


def make_hist(name, title, pt_bins, eta_bins):
    hist = ROOT.TH2D(name, title, len(pt_bins) - 1, pt_bins, len(eta_bins) - 1, eta_bins)
    hist.Sumw2()
    hist.GetXaxis().SetTitle("p_{T} [GeV]")
    hist.GetYaxis().SetTitle("|#eta|")
    return hist


def clamp_to_axis(value, axis):
    if value < axis.GetXmin():
        return axis.GetXmin() + 1.0e-6
    if value >= axis.GetXmax():
        return axis.GetXmax() - 1.0e-6
    return value


def cp_interval(num, den):
    if den <= 0:
        return 0.0, 0.0, 0.0
    num_i = int(round(num))
    den_i = int(round(den))
    num_i = max(0, min(num_i, den_i))
    fr = float(num_i) / float(den_i)
    low = ROOT.TEfficiency.ClopperPearson(den_i, num_i, CL, False)
    high = ROOT.TEfficiency.ClopperPearson(den_i, num_i, CL, True)
    return fr, max(0.0, fr - low), max(0.0, high - fr)


def set_rate_bin(h_fr, h_up, h_down, h_weight, ix, iy, fr, err_down, err_up, max_fr):
    fr = max(0.0, min(max_fr, fr))
    fr_up = min(max_fr, fr + err_up)
    fr_down = max(0.0, fr - err_down)
    h_fr.SetBinContent(ix, iy, fr)
    h_fr.SetBinError(ix, iy, max(err_up, err_down))
    h_up.SetBinContent(ix, iy, fr_up)
    h_down.SetBinContent(ix, iy, fr_down)

    if fr < 1.0:
        w = fr / (1.0 - fr)
        w_up = fr_up / (1.0 - fr_up) if fr_up < 1.0 else 0.0
        w_down = fr_down / (1.0 - fr_down) if fr_down < 1.0 else 0.0
    else:
        w = 0.0
        w_up = 0.0
        w_down = 0.0
    h_weight.SetBinContent(ix, iy, w)
    h_weight.SetBinError(ix, iy, max(abs(w_up - w), abs(w - w_down)))
    return fr, fr_down, fr_up, w, w_down, w_up


def main():
    parser = argparse.ArgumentParser(description="Build 2024 promptMVA fake-rate maps from FR NTuples.")
    parser.add_argument("--input", nargs="+", required=True, help="Input ROOT file glob(s).")
    parser.add_argument("--flavor", choices=["mu", "el"], required=True, help="Probe flavor.")
    parser.add_argument("--out", required=True, help="Output ROOT file.")
    parser.add_argument("--csv", default=None, help="Output CSV table. Default: replace .root with .csv")
    parser.add_argument("--tree", default="Events", help="Input tree name. Default: Events")
    parser.add_argument("--pt-bins", default=DEFAULT_PT_BINS, help="Comma-separated pT bin edges.")
    parser.add_argument("--eta-bins", default=None, help="Comma-separated |eta| bin edges.")
    parser.add_argument("--weight-branch", default="eventWeight", help="Weight branch. Use empty string for unit weights.")
    parser.add_argument("--max-fr", type=float, default=0.95, help="Maximum fake-rate value stored in output maps.")
    parser.add_argument("--progress", type=int, default=500000, help="Print progress every N events. 0 disables it.")
    args = parser.parse_args()

    pt_bins = parse_bins(args.pt_bins)
    eta_bins = parse_bins(args.eta_bins or (DEFAULT_ETA_BINS_MU if args.flavor == "mu" else DEFAULT_ETA_BINS_EL))
    target_flavor = 13 if args.flavor == "mu" else 11

    chain, files = build_chain(args.input, args.tree)
    nentries = chain.GetEntries()

    h_den = make_hist("h_den_%s" % args.flavor, "denominator %s" % args.flavor, pt_bins, eta_bins)
    h_num = make_hist("h_num_%s" % args.flavor, "numerator %s" % args.flavor, pt_bins, eta_bins)
    h_fr = make_hist("h_fr_%s" % args.flavor, "fake rate %s" % args.flavor, pt_bins, eta_bins)
    h_fr_up = make_hist("h_fr_%s_statUp" % args.flavor, "fake rate statUp %s" % args.flavor, pt_bins, eta_bins)
    h_fr_down = make_hist("h_fr_%s_statDown" % args.flavor, "fake rate statDown %s" % args.flavor, pt_bins, eta_bins)
    h_weight = make_hist("h_weight_%s" % args.flavor, "f/(1-f) %s" % args.flavor, pt_bins, eta_bins)

    n_den_events = 0
    n_num_events = 0
    for i, ev in enumerate(chain):
        if args.progress and i > 0 and i % args.progress == 0:
            print("[INFO] processed %d / %d entries" % (i, nentries))
        if not hasattr(ev, "FRRegion"):
            raise RuntimeError("Input tree does not contain FRRegion. Check that FR NTuples were produced with promptMVAFakeRateNtupleProducer_Run3.py")
        if int(ev.FRProbe_flavor) != target_flavor:
            continue
        if not bool(ev.FRProbe_isLoose):
            continue
        if int(ev.FRRegion) not in (1, 2):
            continue

        weight = 1.0
        if args.weight_branch and hasattr(ev, args.weight_branch):
            weight = float(getattr(ev, args.weight_branch))
        pt = clamp_to_axis(float(ev.FRProbe_pt), h_den.GetXaxis())
        abseta = clamp_to_axis(abs(float(ev.FRProbe_absEta)), h_den.GetYaxis())

        h_den.Fill(pt, abseta, weight)
        n_den_events += 1
        if int(ev.FRRegion) == 2 and bool(ev.FRProbe_isTight):
            h_num.Fill(pt, abseta, weight)
            n_num_events += 1

    csv_path = args.csv or os.path.splitext(args.out)[0] + ".csv"
    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    csv_dir = os.path.dirname(os.path.abspath(csv_path))
    if csv_dir:
        os.makedirs(csv_dir, exist_ok=True)

    rows = []
    for ix in range(1, h_den.GetNbinsX() + 1):
        for iy in range(1, h_den.GetNbinsY() + 1):
            den = h_den.GetBinContent(ix, iy)
            num = h_num.GetBinContent(ix, iy)
            fr, err_down, err_up = cp_interval(num, den)
            fr, fr_down, fr_up, weight, weight_down, weight_up = set_rate_bin(
                h_fr, h_fr_up, h_fr_down, h_weight, ix, iy, fr, err_down, err_up, args.max_fr
            )
            rows.append({
                "flavor": args.flavor,
                "pt_low": h_den.GetXaxis().GetBinLowEdge(ix),
                "pt_high": h_den.GetXaxis().GetBinUpEdge(ix),
                "abseta_low": h_den.GetYaxis().GetBinLowEdge(iy),
                "abseta_high": h_den.GetYaxis().GetBinUpEdge(iy),
                "denominator": den,
                "numerator": num,
                "fake_rate": fr,
                "fake_rate_err_down": err_down,
                "fake_rate_err_up": err_up,
                "fake_rate_statDown": fr_down,
                "fake_rate_statUp": fr_up,
                "fake_weight": weight,
                "fake_weight_statDown": weight_down,
                "fake_weight_statUp": weight_up,
            })

    fout = ROOT.TFile(args.out, "RECREATE")
    for obj in (h_den, h_num, h_fr, h_fr_up, h_fr_down, h_weight):
        obj.Write()
    meta = ROOT.TNamed(
        "metadata",
        "flavor=%s;nFiles=%d;nEntries=%d;nDenEvents=%d;nNumEvents=%d;ptBins=%s;etaBins=%s" % (
            args.flavor,
            len(files),
            nentries,
            n_den_events,
            n_num_events,
            ",".join(str(x) for x in pt_bins),
            ",".join(str(x) for x in eta_bins),
        ),
    )
    meta.Write()
    fout.Close()

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)

    print("[INFO] input files: %d" % len(files))
    print("[INFO] tree entries: %d" % nentries)
    print("[INFO] denominator events: %d" % n_den_events)
    print("[INFO] numerator events: %d" % n_num_events)
    print("[INFO] wrote ROOT: %s" % args.out)
    print("[INFO] wrote CSV : %s" % csv_path)


if __name__ == "__main__":
    main()
