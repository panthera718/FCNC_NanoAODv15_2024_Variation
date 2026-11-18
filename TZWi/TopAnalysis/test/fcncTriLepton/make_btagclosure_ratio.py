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
        description="Make ratio ROOT files: no-btagSF / with-btagSF for BasicSelection/hnGoodJet (excluding Data)."
    )
    parser.add_argument("--no", dest="dir_no", required=True,
                        help="Input directory for no-btagSF ROOT files (e.g. hist_BtagClosure_no-btagSF_2018)")
    parser.add_argument("--with", dest="dir_with", required=True,
                        help="Input directory for with-btagSF ROOT files (e.g. hist_BtagClosure_with-btagSF_2018)")
    parser.add_argument("-o", "--out", dest="dir_out", default="hist_BtagClosure_ratio_2018",
                        help="Output directory (default: hist_BtagClosure_ratio_2018)")
    parser.add_argument("--skip-hnBjet", action="store_true",
                        help="If set, do NOT create BasicSelection/hnBjet in output (default: create it empty).")
    args = parser.parse_args()

    dir_no = args.dir_no
    dir_with = args.dir_with
    dir_out = args.dir_out

    if not os.path.isdir(dir_no):
        print(("[ERROR] --no directory does not exist: %s" % dir_no))
        return 1
    if not os.path.isdir(dir_with):
        print(("[ERROR] --with directory does not exist: %s" % dir_with))
        return 1

    ensure_dir(dir_out)

    files = sorted([f for f in os.listdir(dir_no) if f.endswith(".root")])
    if not files:
        print(("[ERROR] No .root files found in: %s" % dir_no))
        return 1

    for fname in files:
        path_no = os.path.join(dir_no, fname)
        path_with = os.path.join(dir_with, fname)
        out_path = os.path.join(dir_out, fname)

        if not os.path.isfile(path_with):
            print(("[WARNING] Missing corresponding file in --with dir, skip: %s" % path_with))
            continue

        f_no = ROOT.TFile.Open(path_no, "READ")
        f_with = ROOT.TFile.Open(path_with, "READ")
        if (not f_no) or f_no.IsZombie():
            print(("[ERROR] Failed to open (no): %s" % path_no))
            continue
        if (not f_with) or f_with.IsZombie():
            print(("[ERROR] Failed to open (with): %s" % path_with))
            f_no.Close()
            continue

        d_no = f_no.Get("BasicSelection/hnGoodJet")
        d_with = f_with.Get("BasicSelection/hnGoodJet")
        if not d_no:
            print(("[ERROR] Directory not found in (no): BasicSelection/hnGoodJet in %s" % path_no))
            f_no.Close()
            f_with.Close()
            continue
        if not d_with:
            print(("[ERROR] Directory not found in (with): BasicSelection/hnGoodJet in %s" % path_with))
            f_no.Close()
            f_with.Close()
            continue

        # Create output file and required directories
        f_out = ROOT.TFile.Open(out_path, "RECREATE")
        if (not f_out) or f_out.IsZombie():
            print(("[ERROR] Failed to create output: %s" % out_path))
            f_no.Close()
            f_with.Close()
            continue

        d_basic = f_out.mkdir("BasicSelection")
        d_basic.cd()

        # Keep structure compatible by default, but allow skipping hnBjet directory
        if not args.skip_hnBjet:
            d_basic.mkdir("hnBjet")

        d_hnGoodJet_out = d_basic.mkdir("hnGoodJet")

        # Loop over keys in no-btagSF hnGoodJet
        n_done = 0
        n_skip = 0
        n_missing = 0

        for key in d_no.GetListOfKeys():
            name = key.GetName()

            # Exclude Data explicitly
            if name == "Data":
                n_skip += 1
                continue

            obj_no = d_no.Get(name)
            obj_with = d_with.Get(name)

            if (not is_th1(obj_no)) or (not is_th1(obj_with)):
                n_skip += 1
                continue

            # Safety: if denominator histogram missing, skip
            if not obj_with:
                n_missing += 1
                continue

            # Ensure input Sumw2 exists for correct error propagation
            ensure_sumw2(obj_no)
            ensure_sumw2(obj_with)

            # Create ratio histogram as clone of numerator
            h_ratio = obj_no.Clone(name)
            h_ratio.SetDirectory(0)

            # Ensure ratio Sumw2 exists without warnings
            ensure_sumw2(h_ratio)

            # Compute ratio = no / with
            # (If denom bin is 0, ROOT sets result bin to 0; this is standard behavior.)
            h_ratio.Divide(obj_with)

            # Write to output directory
            d_hnGoodJet_out.cd()
            h_ratio.Write(name, ROOT.TObject.kOverwrite)

            n_done += 1

        f_out.Close()
        f_no.Close()
        f_with.Close()

        print(("[INFO] Wrote: %s (ratio hists: %d, skipped: %d, missing denom: %d)" %
              (out_path, n_done, n_skip, n_missing)))

    print(("[INFO] Done. Output directory: %s" % dir_out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
