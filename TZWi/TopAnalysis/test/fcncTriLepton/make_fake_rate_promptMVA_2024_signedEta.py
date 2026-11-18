#!/usr/bin/env python3
import argparse
import array
import csv
import glob
import os
from collections import Counter

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

DEFAULT_PT_BINS = "10,15,20,25,30,40,60,100,200"

# Nominal signed-eta binning.
DEFAULT_SIGNED_ETA_BINS_MU = "-2.4,-2.1,-1.2,-0.9,0,0.9,1.2,2.1,2.4"
DEFAULT_SIGNED_ETA_BINS_EL = "-2.5,-2.0,-1.479,-0.8,0,0.8,1.479,2.0,2.5"

# Backward-compatible absolute-eta binning.
DEFAULT_ABSETA_BINS_MU = "0,0.9,1.2,2.1,2.4"
DEFAULT_ABSETA_BINS_EL = "0,0.8,1.479,2.0,2.5"

CL = 0.682689492137

BTAG_SELECTIONS = {
    "inclusive": "all nBJet categories",
    "0b": "FR_nBJet == 0",
    "ge1b": "FR_nBJet >= 1",
    "1b": "FR_nBJet == 1",
    "ge2b": "FR_nBJet >= 2",
}


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
        files.extend(sorted(glob.glob(pattern)))
    files = list(dict.fromkeys(files))
    if not files:
        raise RuntimeError("No input ROOT files matched: %s" % inputs)
    for path in files:
        chain.Add(path)
    return chain, files


def branch_exists(tree, name):
    return bool(tree.GetBranch(name))


def make_hist(name, title, pt_bins, eta_bins, eta_mode):
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
    hist.GetYaxis().SetTitle("#eta" if eta_mode == "signed" else "|#eta|")
    return hist


def make_nbjet_hist(name, title):
    # Last bin includes nBJet >= 5 through explicit clamping in the event loop.
    hist = ROOT.TH1D(name, title, 6, -0.5, 5.5)
    hist.Sumw2()
    hist.GetXaxis().SetTitle("nBJet")
    hist.GetYaxis().SetTitle("Events")
    for ibin, label in enumerate(("0", "1", "2", "3", "4", ">=5"), start=1):
        hist.GetXaxis().SetBinLabel(ibin, label)
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

    # Nominal use is unweighted data. Convert the stored counts to integers
    # before evaluating the exact binomial Clopper-Pearson interval.
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
        fake_weight = fr / (1.0 - fr)
        weight_up = fr_up / (1.0 - fr_up) if fr_up < 1.0 else 0.0
        weight_down = fr_down / (1.0 - fr_down) if fr_down < 1.0 else 0.0
    else:
        fake_weight = 0.0
        weight_up = 0.0
        weight_down = 0.0

    h_weight.SetBinContent(ix, iy, fake_weight)
    h_weight.SetBinError(
        ix,
        iy,
        max(abs(weight_up - fake_weight), abs(fake_weight - weight_down)),
    )
    return fr, fr_down, fr_up, fake_weight, weight_down, weight_up


def default_eta_bins(flavor, eta_mode):
    if eta_mode == "signed":
        return DEFAULT_SIGNED_ETA_BINS_MU if flavor == "mu" else DEFAULT_SIGNED_ETA_BINS_EL
    return DEFAULT_ABSETA_BINS_MU if flavor == "mu" else DEFAULT_ABSETA_BINS_EL


def passes_btag_selection(n_bjet, category):
    if category == "inclusive":
        return True
    if category == "0b":
        return n_bjet == 0
    if category == "ge1b":
        return n_bjet >= 1
    if category == "1b":
        return n_bjet == 1
    if category == "ge2b":
        return n_bjet >= 2
    raise ValueError("Unsupported b-tag category: %s" % category)


def check_stored_category_consistency(ev, n_bjet, has_cat_0ge1, has_cat_01ge2):
    if has_cat_0ge1:
        stored = int(ev.FR_btagCategory0GE1)
        expected = 0 if n_bjet == 0 else 1
        if stored != expected:
            raise RuntimeError(
                "Inconsistent FR_btagCategory0GE1: stored=%d expected=%d for FR_nBJet=%d"
                % (stored, expected, n_bjet)
            )

    if has_cat_01ge2:
        stored = int(ev.FR_btagCategory01GE2)
        expected = min(n_bjet, 2)
        if stored != expected:
            raise RuntimeError(
                "Inconsistent FR_btagCategory01GE2: stored=%d expected=%d for FR_nBJet=%d"
                % (stored, expected, n_bjet)
            )


def main():
    parser = argparse.ArgumentParser(
        description="Build 2024 promptMVA fake-rate maps from b-tag-aware FR NTuples."
    )
    parser.add_argument("--input", nargs="+", required=True, help="Input ROOT file glob(s).")
    parser.add_argument("--flavor", choices=["mu", "el"], required=True, help="Probe flavor.")
    parser.add_argument("--out", required=True, help="Output ROOT file.")
    parser.add_argument(
        "--csv",
        default=None,
        help="Output CSV table. Default: replace .root with .csv",
    )
    parser.add_argument("--tree", default="Events", help="Input tree name. Default: Events")
    parser.add_argument(
        "--pt-bins",
        default=DEFAULT_PT_BINS,
        help="Comma-separated pT bin edges.",
    )
    parser.add_argument(
        "--eta-mode",
        choices=["signed", "abs"],
        default="signed",
        help="Use signed eta or absolute eta. Default: signed",
    )
    parser.add_argument(
        "--eta-bins",
        default=None,
        help="Comma-separated eta bin edges. Defaults depend on flavor and eta mode.",
    )
    parser.add_argument(
        "--btag-category",
        choices=list(BTAG_SELECTIONS.keys()),
        default="inclusive",
        help=(
            "Event-level b-tag category: inclusive, 0b, ge1b, 1b, or ge2b. "
            "Default: inclusive"
        ),
    )
    parser.add_argument(
        "--btag-branch",
        default="FR_nBJet",
        help="Event-level b-jet multiplicity branch. Default: FR_nBJet",
    )
    parser.add_argument(
        "--weight-branch",
        default="eventWeight",
        help="Weight branch. Use an empty string for unit weights.",
    )
    parser.add_argument(
        "--max-fr",
        type=float,
        default=0.95,
        help="Maximum fake-rate value stored in output maps.",
    )
    parser.add_argument(
        "--progress",
        type=int,
        default=500000,
        help="Print progress every N events. 0 disables it.",
    )
    args = parser.parse_args()

    pt_bins = parse_bins(args.pt_bins)
    eta_bins = parse_bins(args.eta_bins or default_eta_bins(args.flavor, args.eta_mode))
    target_flavor = 13 if args.flavor == "mu" else 11
    btag_selection_text = BTAG_SELECTIONS[args.btag_category]

    chain, files = build_chain(args.input, args.tree)
    nentries = chain.GetEntries()
    if nentries <= 0:
        raise RuntimeError("Input chain has no entries. Check the input glob and tree name.")

    required = [
        "FRRegion",
        "FRProbe_flavor",
        "FRProbe_isLoose",
        "FRProbe_isTight",
        "FRProbe_pt",
    ]
    required.append("FRProbe_eta" if args.eta_mode == "signed" else "FRProbe_absEta")

    has_btag_branch = branch_exists(chain, args.btag_branch)
    if args.btag_category != "inclusive":
        required.append(args.btag_branch)

    missing = [br for br in required if not branch_exists(chain, br)]
    if missing:
        raise RuntimeError(
            "Input FR NTuple is missing required branch(es): %s" % ", ".join(missing)
        )

    has_weight_branch = bool(args.weight_branch) and branch_exists(chain, args.weight_branch)
    has_cat_0ge1 = branch_exists(chain, "FR_btagCategory0GE1")
    has_cat_01ge2 = branch_exists(chain, "FR_btagCategory01GE2")
    has_btag_wp_value = branch_exists(chain, "FR_btagWPValue")
    has_btag_wp_code = branch_exists(chain, "FR_btagWPCode")

    category_title = "%s, %s" % (args.flavor, btag_selection_text)
    h_den = make_hist(
        "h_den_%s" % args.flavor,
        "denominator %s" % category_title,
        pt_bins,
        eta_bins,
        args.eta_mode,
    )
    h_num = make_hist(
        "h_num_%s" % args.flavor,
        "numerator %s" % category_title,
        pt_bins,
        eta_bins,
        args.eta_mode,
    )
    h_fr = make_hist(
        "h_fr_%s" % args.flavor,
        "fake rate %s" % category_title,
        pt_bins,
        eta_bins,
        args.eta_mode,
    )
    h_fr_up = make_hist(
        "h_fr_%s_statUp" % args.flavor,
        "fake rate statUp %s" % category_title,
        pt_bins,
        eta_bins,
        args.eta_mode,
    )
    h_fr_down = make_hist(
        "h_fr_%s_statDown" % args.flavor,
        "fake rate statDown %s" % category_title,
        pt_bins,
        eta_bins,
        args.eta_mode,
    )
    h_weight = make_hist(
        "h_weight_%s" % args.flavor,
        "f/(1-f) %s" % category_title,
        pt_bins,
        eta_bins,
        args.eta_mode,
    )

    h_nbjet_pre = None
    h_nbjet_selected = None
    if has_btag_branch:
        h_nbjet_pre = make_nbjet_hist(
            "h_nBJet_preselection_%s" % args.flavor,
            "nBJet before category selection, %s" % args.flavor,
        )
        h_nbjet_selected = make_nbjet_hist(
            "h_nBJet_selected_%s" % args.flavor,
            "nBJet after %s, %s" % (args.btag_category, args.flavor),
        )

    n_pre_btag_candidates = 0
    n_btag_rejected = 0
    n_den_events = 0
    n_num_events = 0
    n_pos_eta = 0
    n_neg_eta = 0
    nbjet_counter = Counter()
    selected_nbjet_counter = Counter()
    btag_wp_values = set()
    btag_wp_codes = set()

    for i, ev in enumerate(chain):
        if args.progress and i > 0 and i % args.progress == 0:
            print("[INFO] processed %d / %d entries" % (i, nentries))

        if int(ev.FRProbe_flavor) != target_flavor:
            continue
        if not bool(ev.FRProbe_isLoose):
            continue
        if int(ev.FRRegion) not in (1, 2):
            continue

        n_pre_btag_candidates += 1

        n_bjet = None
        if has_btag_branch:
            n_bjet = int(getattr(ev, args.btag_branch))
            if n_bjet < 0:
                raise RuntimeError(
                    "Negative %s=%d encountered" % (args.btag_branch, n_bjet)
                )

            check_stored_category_consistency(
                ev,
                n_bjet,
                has_cat_0ge1,
                has_cat_01ge2,
            )

            nbjet_counter[n_bjet] += 1
            if h_nbjet_pre:
                h_nbjet_pre.Fill(min(n_bjet, 5))

            if has_btag_wp_value:
                btag_wp_values.add(round(float(ev.FR_btagWPValue), 10))
            if has_btag_wp_code:
                btag_wp_codes.add(int(ev.FR_btagWPCode))

        if args.btag_category != "inclusive":
            if n_bjet is None:
                raise RuntimeError(
                    "Cannot apply %s without branch %s"
                    % (args.btag_category, args.btag_branch)
                )
            if not passes_btag_selection(n_bjet, args.btag_category):
                n_btag_rejected += 1
                continue

        if n_bjet is not None:
            selected_nbjet_counter[n_bjet] += 1
            if h_nbjet_selected:
                h_nbjet_selected.Fill(min(n_bjet, 5))

        weight = float(getattr(ev, args.weight_branch)) if has_weight_branch else 1.0

        pt = clamp_to_axis(float(ev.FRProbe_pt), h_den.GetXaxis())
        if args.eta_mode == "signed":
            eta_value = float(ev.FRProbe_eta)
            if eta_value >= 0.0:
                n_pos_eta += 1
            else:
                n_neg_eta += 1
        else:
            eta_value = abs(float(ev.FRProbe_absEta))
        eta = clamp_to_axis(eta_value, h_den.GetYaxis())

        h_den.Fill(pt, eta, weight)
        n_den_events += 1

        if int(ev.FRRegion) == 2 and bool(ev.FRProbe_isTight):
            h_num.Fill(pt, eta, weight)
            n_num_events += 1

    csv_path = args.csv or os.path.splitext(args.out)[0] + ".csv"
    for path in (args.out, csv_path):
        directory = os.path.dirname(os.path.abspath(path))
        if directory:
            os.makedirs(directory, exist_ok=True)

    rows = []
    for ix in range(1, h_den.GetNbinsX() + 1):
        for iy in range(1, h_den.GetNbinsY() + 1):
            den = h_den.GetBinContent(ix, iy)
            num = h_num.GetBinContent(ix, iy)
            fr, err_down, err_up = cp_interval(num, den)
            fr, fr_down, fr_up, fake_weight, weight_down, weight_up = set_rate_bin(
                h_fr,
                h_fr_up,
                h_fr_down,
                h_weight,
                ix,
                iy,
                fr,
                err_down,
                err_up,
                args.max_fr,
            )
            rows.append({
                "flavor": args.flavor,
                "eta_mode": args.eta_mode,
                "btag_category": args.btag_category,
                "btag_selection": btag_selection_text,
                "btag_branch": args.btag_branch,
                "pt_low": h_den.GetXaxis().GetBinLowEdge(ix),
                "pt_high": h_den.GetXaxis().GetBinUpEdge(ix),
                "eta_low": h_den.GetYaxis().GetBinLowEdge(iy),
                "eta_high": h_den.GetYaxis().GetBinUpEdge(iy),
                "denominator": den,
                "numerator": num,
                "fake_rate": fr,
                "fake_rate_err_down": err_down,
                "fake_rate_err_up": err_up,
                "fake_rate_statDown": fr_down,
                "fake_rate_statUp": fr_up,
                "fake_weight": fake_weight,
                "fake_weight_statDown": weight_down,
                "fake_weight_statUp": weight_up,
            })

    fout = ROOT.TFile(args.out, "RECREATE")
    output_objects = [h_den, h_num, h_fr, h_fr_up, h_fr_down, h_weight]
    if h_nbjet_pre:
        output_objects.append(h_nbjet_pre)
    if h_nbjet_selected:
        output_objects.append(h_nbjet_selected)
    for obj in output_objects:
        obj.Write()

    metadata = (
        "flavor=%s;etaMode=%s;btagCategory=%s;btagSelection=%s;"
        "btagBranch=%s;nFiles=%d;nEntries=%d;nPreBtagCandidates=%d;"
        "nBtagRejected=%d;nDenEvents=%d;nNumEvents=%d;"
        "nPosEtaDenEvents=%d;nNegEtaDenEvents=%d;"
        "btagWPValues=%s;btagWPCodes=%s;ptBins=%s;etaBins=%s"
        % (
            args.flavor,
            args.eta_mode,
            args.btag_category,
            btag_selection_text,
            args.btag_branch,
            len(files),
            nentries,
            n_pre_btag_candidates,
            n_btag_rejected,
            n_den_events,
            n_num_events,
            n_pos_eta,
            n_neg_eta,
            ",".join(str(x) for x in sorted(btag_wp_values)) or "notStored",
            ",".join(str(x) for x in sorted(btag_wp_codes)) or "notStored",
            ",".join(str(x) for x in pt_bins),
            ",".join(str(x) for x in eta_bins),
        )
    )
    ROOT.TNamed("metadata", metadata).Write()
    fout.Close()

    with open(csv_path, "w", newline="") as f:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("[INFO] input files: %d" % len(files))
    print("[INFO] tree entries: %d" % nentries)
    print("[INFO] eta mode: %s" % args.eta_mode)
    print("[INFO] b-tag category: %s (%s)" % (args.btag_category, btag_selection_text))
    print("[INFO] pre-btag denominator candidates: %d" % n_pre_btag_candidates)
    print("[INFO] rejected by b-tag category: %d" % n_btag_rejected)
    print("[INFO] denominator events: %d" % n_den_events)
    print("[INFO] numerator events: %d" % n_num_events)
    if args.eta_mode == "signed":
        print(
            "[INFO] signed-eta denominator events: eta>=0: %d, eta<0: %d"
            % (n_pos_eta, n_neg_eta)
        )
    if has_btag_branch:
        print("[INFO] nBJet distribution before category selection: %s" % dict(sorted(nbjet_counter.items())))
        print("[INFO] nBJet distribution after category selection : %s" % dict(sorted(selected_nbjet_counter.items())))
    if btag_wp_values:
        print("[INFO] FR_btagWPValue values: %s" % sorted(btag_wp_values))
    if btag_wp_codes:
        print("[INFO] FR_btagWPCode values: %s" % sorted(btag_wp_codes))
    print("[INFO] wrote ROOT: %s" % args.out)
    print("[INFO] wrote CSV : %s" % csv_path)


if __name__ == "__main__":
    main()
