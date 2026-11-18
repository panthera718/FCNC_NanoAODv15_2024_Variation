#!/usr/bin/env python32
# -*- coding: utf-8 -*-

import os
import sys
import re
import argparse

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)


def ensure_dir(path):
    if path and (not os.path.isdir(path)):
        os.makedirs(path)


def is_th1(obj):
    return obj and obj.InheritsFrom("TH1")


def sanitize_filename(s):
    # Keep only [A-Za-z0-9._-], replace others with "_"
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s)


def ensure_sumw2(h):
    if not h:
        return
    if h.GetSumw2N() == 0:
        h.Sumw2()


def draw_single_hist(h, out_pdf, out_png, label_mode="", label_hist="",
                     ymin=None, ymax=None, draw_line_at_1=True):
    c = ROOT.TCanvas("c", "c", 800, 650)
    c.SetTicks(1, 1)

    # Basic style
    h.SetLineColor(ROOT.kBlue + 1)
    h.SetMarkerColor(ROOT.kBlue + 1)
    h.SetMarkerStyle(20)

    # Y-range control
    if ymin is not None:
        h.SetMinimum(float(ymin))
    if ymax is not None:
        h.SetMaximum(float(ymax))

    h.Draw("E1")

    # Optional horizontal line at 1
    if draw_line_at_1:
        x1 = h.GetXaxis().GetXmin()
        x2 = h.GetXaxis().GetXmax()
        line = ROOT.TLine(x1, 1.0, x2, 1.0)
        line.SetLineStyle(2)
        line.Draw("SAME")

    # Labels
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.035)
    if label_mode:
        latex.DrawLatex(0.12, 0.92, "Mode: %s" % label_mode)
    if label_hist:
        latex.DrawLatex(0.12, 0.87, "Hist: %s" % label_hist)

    ensure_dir(os.path.dirname(out_pdf))
    c.SaveAs(out_pdf)
    c.SaveAs(out_png)

    c.Close()
    del c


def main():
    parser = argparse.ArgumentParser(
        description="Draw histograms from hist_BtagClosure_ratio_2018 (BasicSelection/hnGoodJet, excluding Data). Save as pdf+png."
    )
    parser.add_argument("-i", "--in", dest="dir_in", required=True,
                        help="Input directory (e.g. hist_BtagClosure_ratio_2018)")
    parser.add_argument("-o", "--out", dest="dir_out", default="plots_BtagClosure_ratio_2018",
                        help="Output directory for plots (default: plots_BtagClosure_ratio_2018)")
    parser.add_argument("--ymin", type=float, default=None, help="Set y-axis minimum (optional)")
    parser.add_argument("--ymax", type=float, default=None, help="Set y-axis maximum (optional)")
    parser.add_argument("--no-line1", action="store_true",
                        help="If set, do not draw horizontal line at y=1")
    args = parser.parse_args()

    dir_in = args.dir_in
    dir_out = args.dir_out

    if not os.path.isdir(dir_in):
        print(("[ERROR] Input directory does not exist: %s" % dir_in))
        return 1

    ensure_dir(dir_out)

    root_files = sorted([f for f in os.listdir(dir_in) if f.endswith(".root")])
    if not root_files:
        print(("[ERROR] No .root files found in: %s" % dir_in))
        return 1

    for fname in root_files:
        path = os.path.join(dir_in, fname)

        f = ROOT.TFile.Open(path, "READ")
        if (not f) or f.IsZombie():
            print(("[ERROR] Failed to open: %s" % path))
            continue

        d = f.Get("BasicSelection/hnGoodJet")
        if not d:
            print(("[ERROR] Missing BasicSelection/hnGoodJet in: %s" % path))
            f.Close()
            continue

        mode = os.path.splitext(fname)[0]
        mode_out = os.path.join(dir_out, mode)
        ensure_dir(mode_out)

        n_plots = 0
        n_skip = 0

        for key in d.GetListOfKeys():
            hname = key.GetName()

            # Exclude Data
            if hname == "Data":
                n_skip += 1
                continue

            h = d.Get(hname)
            if not is_th1(h):
                n_skip += 1
                continue

            # Clone to detach from file
            h_cl = h.Clone("ratio_%s" % hname)
            h_cl.SetDirectory(0)
            ensure_sumw2(h_cl)

            safe = sanitize_filename(hname)
            out_pdf = os.path.join(mode_out, "%s.pdf" % safe)
            out_png = os.path.join(mode_out, "%s.png" % safe)

            draw_single_hist(
                h_cl,
                out_pdf, out_png,
                label_mode=mode,
                label_hist=hname,
                ymin=args.ymin,
                ymax=args.ymax,
                draw_line_at_1=(not args.no_line1)
            )

            n_plots += 1

        f.Close()
        print(("[INFO] %s: made %d plots, skipped %d (Data/non-TH1)" % (fname, n_plots, n_skip)))

    print(("[INFO] Done. Output plots: %s" % dir_out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
