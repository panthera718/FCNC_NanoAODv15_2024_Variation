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
    if not os.path.isdir(path):
        os.makedirs(path)


def is_th1(obj):
    return obj and obj.InheritsFrom("TH1")


def ensure_sumw2(h):
    """Create Sumw2 only if missing (avoid TH1::Sumw2 warnings)."""
    if not h:
        return
    if h.GetSumw2N() == 0:
        h.Sumw2()


def sanitize_filename(s):
    # Keep only [A-Za-z0-9._-], replace others with "_"
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s)


def make_ratio_hist(h_num, h_den, name="ratio"):
    """
    Build ratio histogram = num/den with safe handling for den=0,
    and standard error propagation.
    """
    h = h_num.Clone(name)
    h.SetDirectory(0)
    h.Reset("ICES")  # reset contents/errors

    nb = h_num.GetNbinsX()
    for i in range(1, nb + 1):
        a  = h_num.GetBinContent(i)
        ea = h_num.GetBinError(i)
        b  = h_den.GetBinContent(i)
        eb = h_den.GetBinError(i)

        if b == 0.0:
            h.SetBinContent(i, 0.0)
            h.SetBinError(i, 0.0)
            continue

        r = a / b

        # error propagation:
        # r = a/b
        # dr^2 = (ea/b)^2 + (a*eb/b^2)^2
        # If a==0 -> second term is 0
        er2 = (ea / b) * (ea / b)
        if a != 0.0:
            er2 += (a * eb / (b * b)) * (a * eb / (b * b))

        h.SetBinContent(i, r)
        h.SetBinError(i, er2 ** 0.5)

    return h


def draw_compare_and_ratio(h_no, h_corr, out_pdf, out_png,
                           label_mode="", label_hist="",
                           rmin=0.5, rmax=1.5):
    """
    Top pad: overlay (no vs corrected)
    Bottom pad: ratio (corrected/no)
    Save to pdf & png.
    """
    # Styling (simple, readable)
    h_no.SetLineColor(ROOT.kBlack)
    h_no.SetMarkerColor(ROOT.kBlack)
    h_no.SetMarkerStyle(20)

    h_corr.SetLineColor(ROOT.kRed + 1)
    h_corr.SetMarkerColor(ROOT.kRed + 1)
    h_corr.SetMarkerStyle(24)

    # Canvas with two pads
    c = ROOT.TCanvas("c", "c", 800, 800)

    pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.30)

    pad1.SetBottomMargin(0.02)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.35)

    pad1.Draw()
    pad2.Draw()

    # ---------- Top pad ----------
    pad1.cd()

    maxy = max(h_no.GetMaximum(), h_corr.GetMaximum())
    # give headroom
    h_no.SetMaximum(maxy * 1.35 if maxy > 0 else 1.0)

    # Hide x labels on top pad
    h_no.GetXaxis().SetLabelSize(0)
    h_no.GetXaxis().SetTitleSize(0)

    h_no.Draw("E1")
    h_corr.Draw("E1 SAME")

    leg = ROOT.TLegend(0.60, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

#    leg.AddEntry(h_no, "no-btagSF", "lep")
#    leg.AddEntry(h_corr, "corrected-to-no", "lep")

#    leg.AddEntry(h_no, "with-btagSF", "lep")
#    leg.AddEntry(h_corr, "no-btagSF", "lep")

    leg.AddEntry(h_no, "no-btagSF", "lep")
    leg.AddEntry(h_corr, "corrected to no-btagSF", "lep")

    leg.Draw()

    # Optional text labels
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.035)
    if label_mode:
        latex.DrawLatex(0.12, 0.92, "Mode: %s" % label_mode)
    if label_hist:
        latex.DrawLatex(0.12, 0.87, "Hist: %s" % label_hist)

    # ---------- Bottom pad (ratio) ----------
    pad2.cd()

    ensure_sumw2(h_no)
    ensure_sumw2(h_corr)
    h_ratio = make_ratio_hist(h_corr, h_no, "ratio_corr_over_no")
    ensure_sumw2(h_ratio)

    h_ratio.SetLineColor(ROOT.kBlue + 1)
    h_ratio.SetMarkerColor(ROOT.kBlue + 1)
    h_ratio.SetMarkerStyle(20)

#    h_ratio.GetYaxis().SetTitle("corr / no")
#    h_ratio.GetYaxis().SetTitle("no-btagSF / with-btagSF")
    h_ratio.GetYaxis().SetTitle("corrected / no-btagSF")

    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitleSize(0.075)
    h_ratio.GetYaxis().SetTitleOffset(0.50)
    h_ratio.GetYaxis().SetLabelSize(0.09)

    h_ratio.GetXaxis().SetTitle(h_no.GetXaxis().GetTitle())
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetTitleOffset(1.10)
    h_ratio.GetXaxis().SetLabelSize(0.10)

    h_ratio.SetMinimum(rmin)
    h_ratio.SetMaximum(rmax)

    h_ratio.Draw("E1")

    # horizontal line at 1
    x1 = h_ratio.GetXaxis().GetXmin()
    x2 = h_ratio.GetXaxis().GetXmax()
    line = ROOT.TLine(x1, 1.0, x2, 1.0)
    line.SetLineStyle(2)
    line.Draw("SAME")

    # Save
    ensure_dir(os.path.dirname(out_pdf))
    c.SaveAs(out_pdf)
    c.SaveAs(out_png)

    # Clean up (avoid name collisions in ROOT global lists)
    c.Close()
    del c


def main():
    parser = argparse.ArgumentParser(
        description="Overlay no-btagSF vs corrected-to-no-btagSF and draw ratio (corrected/no) for BasicSelection/hnGoodJet (excluding Data). Save as pdf+png."
    )
    parser.add_argument("--no", dest="dir_no", required=True,
                        help="Directory: hist_BtagClosure_no-btagSF_2018")
    parser.add_argument("--corr", dest="dir_corr", required=True,
                        help="Directory: hist_BtagClosure_corrected-to-no-btagSF_2018")
    parser.add_argument("-o", "--out", dest="out_dir", default="plots_BtagClosure_compare_2018",
                        help="Output plot directory (default: plots_BtagClosure_compare_2018)")
    parser.add_argument("--rmin", type=float, default=0.5, help="Ratio plot min (default: 0.5)")
    parser.add_argument("--rmax", type=float, default=1.5, help="Ratio plot max (default: 1.5)")
    args = parser.parse_args()

    dir_no = args.dir_no
    dir_corr = args.dir_corr
    out_dir = args.out_dir

    if not os.path.isdir(dir_no):
        print(("[ERROR] --no directory does not exist: %s" % dir_no))
        return 1
    if not os.path.isdir(dir_corr):
        print(("[ERROR] --corr directory does not exist: %s" % dir_corr))
        return 1

    ensure_dir(out_dir)

    root_files = sorted([f for f in os.listdir(dir_no) if f.endswith(".root")])
    if not root_files:
        print(("[ERROR] No .root files in: %s" % dir_no))
        return 1

    for fname in root_files:
        path_no = os.path.join(dir_no, fname)
        path_corr = os.path.join(dir_corr, fname)

        if not os.path.isfile(path_corr):
            print(("[WARNING] Missing corrected file, skip: %s" % path_corr))
            continue

        f_no = ROOT.TFile.Open(path_no, "READ")
        f_corr = ROOT.TFile.Open(path_corr, "READ")
        if (not f_no) or f_no.IsZombie():
            print(("[ERROR] Failed to open no file: %s" % path_no))
            continue
        if (not f_corr) or f_corr.IsZombie():
            print(("[ERROR] Failed to open corr file: %s" % path_corr))
            f_no.Close()
            continue

        d_no = f_no.Get("BasicSelection/hnGoodJet")
        d_corr = f_corr.Get("BasicSelection/hnGoodJet")
        if not d_no:
            print(("[ERROR] Missing BasicSelection/hnGoodJet in: %s" % path_no))
            f_no.Close()
            f_corr.Close()
            continue
        if not d_corr:
            print(("[ERROR] Missing BasicSelection/hnGoodJet in: %s" % path_corr))
            f_no.Close()
            f_corr.Close()
            continue

        mode = os.path.splitext(fname)[0]  # ElElEl, ElMuMu, ...
        mode_dir = os.path.join(out_dir, mode)
        ensure_dir(mode_dir)

        n_plots = 0
        n_skip = 0

        # Loop by keys in no directory (reference list)
        for key in d_no.GetListOfKeys():
            hname = key.GetName()

            # Exclude Data as requested
            if hname == "Data":
                n_skip += 1
                continue

            h_no = d_no.Get(hname)
            h_corr = d_corr.Get(hname)

            if (not is_th1(h_no)) or (not is_th1(h_corr)):
                n_skip += 1
                continue

            # Detach from files to avoid lifetime issues after file close
            h_no_cl = h_no.Clone("no_%s" % hname)
            h_corr_cl = h_corr.Clone("corr_%s" % hname)
            h_no_cl.SetDirectory(0)
            h_corr_cl.SetDirectory(0)

            ensure_sumw2(h_no_cl)
            ensure_sumw2(h_corr_cl)

            safe = sanitize_filename(hname)
            out_pdf = os.path.join(mode_dir, "%s.pdf" % safe)
            out_png = os.path.join(mode_dir, "%s.png" % safe)

            draw_compare_and_ratio(
                h_no_cl, h_corr_cl,
                out_pdf, out_png,
                label_mode=mode,
                label_hist=hname,
                rmin=args.rmin,
                rmax=args.rmax
            )

            n_plots += 1

        f_no.Close()
        f_corr.Close()

        print(("[INFO] %s: made %d plots, skipped %d (Data/non-TH1/missing)" % (fname, n_plots, n_skip)))

    print(("[INFO] Done. Output plots: %s" % out_dir))
    return 0


if __name__ == "__main__":
    sys.exit(main())
