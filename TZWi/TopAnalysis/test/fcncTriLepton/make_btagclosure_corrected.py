#!/usr/bin/env python32
# -*- coding: utf-8 -*-

import os
import sys
import argparse

import ROOT
ROOT.gROOT.SetBatch(True)


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def is_th1(obj):
    return obj and obj.InheritsFrom("TH1")


def ensure_sumw2(h):
    """Create Sumw2 only if it does not exist (prevents TH1D::Sumw2 warnings)."""
    if not h:
        return
    if h.GetSumw2N() == 0:
        h.Sumw2()


def main():
    parser = argparse.ArgumentParser(
        description="Make corrected ROOT files: (with-btagSF) * (ratio) -> corrected-to-no-btagSF, for BasicSelection/hnGoodJet (excluding Data, excluding hnBjet)."
    )
    parser.add_argument("--with", dest="dir_with", required=True,
                        help="Input directory for with-btagSF ROOT files (e.g. hist_BtagClosure_with-btagSF_2018)")
    parser.add_argument("--ratio", dest="dir_ratio", required=True,
                        help="Input directory for ratio ROOT files (e.g. hist_BtagClosure_ratio_2018)")
    parser.add_argument("-o", "--out", dest="dir_out",
                        default="hist_BtagClosure_corrected-to-no-btagSF_2018",
                        help="Output directory (default: hist_BtagClosure_corrected-to-no-btagSF_2018)")
    args = parser.parse_args()

    dir_with = args.dir_with
    dir_ratio = args.dir_ratio
    dir_out = args.dir_out

    if not os.path.isdir(dir_with):
        print(("[ERROR] --with directory does not exist: %s" % dir_with))
        return 1
    if not os.path.isdir(dir_ratio):
        print(("[ERROR] --ratio directory does not exist: %s" % dir_ratio))
        return 1

    ensure_dir(dir_out)

    files = sorted([f for f in os.listdir(dir_with) if f.endswith(".root")])
    if not files:
        print(("[ERROR] No .root files found in: %s" % dir_with))
        return 1

    for fname in files:
        path_with = os.path.join(dir_with, fname)
        path_ratio = os.path.join(dir_ratio, fname)
        out_path = os.path.join(dir_out, fname)

        if not os.path.isfile(path_ratio):
            print(("[WARNING] Missing corresponding file in --ratio dir, skip: %s" % path_ratio))
            continue

        f_with = ROOT.TFile.Open(path_with, "READ")
        f_ratio = ROOT.TFile.Open(path_ratio, "READ")
        if (not f_with) or f_with.IsZombie():
            print(("[ERROR] Failed to open (with): %s" % path_with))
            continue
        if (not f_ratio) or f_ratio.IsZombie():
            print(("[ERROR] Failed to open (ratio): %s" % path_ratio))
            f_with.Close()
            continue

        d_with = f_with.Get("BasicSelection/hnGoodJet")
        d_ratio = f_ratio.Get("BasicSelection/hnGoodJet")
        if not d_with:
            print(("[ERROR] Directory not found in (with): BasicSelection/hnGoodJet in %s" % path_with))
            f_with.Close()
            f_ratio.Close()
            continue
        if not d_ratio:
            print(("[ERROR] Directory not found in (ratio): BasicSelection/hnGoodJet in %s" % path_ratio))
            f_with.Close()
            f_ratio.Close()
            continue

        # Create output file and required directories
        f_out = ROOT.TFile.Open(out_path, "RECREATE")
        if (not f_out) or f_out.IsZombie():
            print(("[ERROR] Failed to create output: %s" % out_path))
            f_with.Close()
            f_ratio.Close()
            continue

        d_basic = f_out.mkdir("BasicSelection")
        d_basic.cd()

        # hnBjet is excluded as requested (do NOT create it)
        d_hnGoodJet_out = d_basic.mkdir("hnGoodJet")

        # Loop over keys in with-btagSF hnGoodJet
        n_done = 0
        n_skip = 0
        n_missing = 0
        n_notth1 = 0

        for key in d_with.GetListOfKeys():
            name = key.GetName()

            # Exclude Data explicitly
            if name == "Data":
                n_skip += 1
                continue

            h_with = d_with.Get(name)
            h_rat = d_ratio.Get(name)

            if (not is_th1(h_with)) or (not is_th1(h_rat)):
                # If either missing or not TH1, skip
                n_notth1 += 1
                continue

            # Ensure Sumw2 exists for correct error propagation
            ensure_sumw2(h_with)
            ensure_sumw2(h_rat)

            # corrected = with * ratio
            h_corr = h_with.Clone(name)
            h_corr.SetDirectory(0)
            ensure_sumw2(h_corr)

            # Multiply bin-by-bin (errors propagated if Sumw2 exists)
            h_corr.Multiply(h_rat)

            # Write to output directory
            d_hnGoodJet_out.cd()
            h_corr.Write(name, ROOT.TObject.kOverwrite)

            n_done += 1

        f_out.Close()
        f_with.Close()
        f_ratio.Close()

        print(("[INFO] Wrote: %s (corrected hists: %d, skipped(Data): %d, non-TH1/missing: %d)" %
              (out_path, n_done, n_skip, n_notth1)))

    print(("[INFO] Done. Output directory: %s" % dir_out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
