#!/usr/bin/env python3

import yaml
from glob import *
from ROOT import *
import os

#modes = [os.path.basename(x) for x in glob("raw_hist/*")]
#modes = ["ElElEl", "MuElEl", "ElMuMu", "MuMuMu"]
modes = [
    "ElElEl",
    "MuMuMu",
    "TaTaTa",
    "ElMuMu",
    "ElTaTa",
    "MuElEl",
    "MuTaTa",
    "TaElEl",
    "TaMuMu"
]
odName = "hist_2018"

## Load configurations
info = {}
xsecSetFile = "config/crosssection.yaml"
histSetFile = "config/histogramming.yaml"
info.update(yaml.safe_load(open(histSetFile)))
info.update(yaml.safe_load(open(xsecSetFile)))
#info.update(yaml.safe_load(open("config/systematics.yaml")))
info.update(yaml.safe_load(open("config/grouping.yaml")))
for f in glob("NanoAOD/2018/*.yaml"):
    if 'dataset' not in info: info['dataset'] = {}
    info['dataset'].update(yaml.safe_load(open(f))['dataset'])

if not os.path.exists(odName): os.makedirs(odName)

def makedirs(d, dName):
    if '/' not in dName: return d.mkdir(dName)

    d1, d2 = dName.split('/',1)
    dNext = d.GetDirectory(d1)
    if dNext == None: dNext = d.mkdir(d1)
    return makedirs(dNext, d2)

for mode in modes:
    fout = TFile("%s/%s.root" % (odName, mode), "recreate")
    ## Plan how to merge histograms
    houts = {}

    ## Collect files and their statistics
    fins = {}
    for f in glob("raw_hist_2018/%s/*/*.root" % mode):
        fin = TFile(f)
        fins[f.split('/',1)[-1]] = fin

    ## Loop over steps x histograms
    print("@@ Collecting source histograms...", mode)
    for stepInfo in info['steps']:
        print(stepInfo['name'], "/", end=' ')
        for hName in stepInfo['hists']:
            if hName not in info['hists']: continue
            hinPath = "%s/h%s" % (stepInfo['name'], hName)
            print(hName, end=' ')

            ## open raw_hist root files
            for proc, procInfo in info['processes'].items():
                title = procInfo['title']
                longTitle = procInfo['longTitle'] if 'longTitle' in procInfo else title

                houtPath = "%s/%s" % (hinPath, title)
                if houtPath not in houts: houts[houtPath] = []

                for datasetGroup in procInfo['datasets']:
                    if datasetGroup not in info['dataset']: continue
                    physProcName = datasetGroup.split('.')[1]

                    xsec = 1.0
                    if physProcName in info['crosssection']:
                        xsec = info['crosssection'][physProcName]
                    elif title != "Data":
                        print("Could not find", physProcName, "in the cross section list. setting it to be 1.0")

                    hout = None
                    nEvents = 0.
                    for ds in list(info['dataset'][datasetGroup].keys()):
                        fName = "%s/%s/%s.root" % (mode, proc, ds[1:].replace('/', '.'))
                        if fName not in fins: continue
                        fin = fins[fName]
                        hin = fin.Get(hinPath)
                        if hin == None: continue

                        nEvents += fin.Get("hCutFlow").GetBinContent(1)

                        if hout == None:
                            hout = hin.Clone()
                            hout.SetTitle(longTitle)
                            hout.Reset()
                        hout.Add(hin)
                    if hout == None: continue
                    if nEvents != 0.: hout.Scale(59800.0*xsec/nEvents)

                    houts[houtPath].append(hout)
        print("")

    mode_np = "NPL" + mode
    houts_np = {}

    ## Collect files and their statistics
    fins_np = {}
    for f in glob("raw_hist_2018/%s/*/*.root" % mode_np):
        fin = TFile(f)
        fins_np[f.split('/',1)[-1]] = fin

    ## Loop over steps x histograms
    print("@@ Collecting source histograms...", mode_np)
    for stepInfo in info['steps']:
        print(stepInfo['name'], "/", end=' ')
        for hName in stepInfo['hists']:
            if hName not in info['hists']: continue
            hinPath = "%s/h%s" % (stepInfo['name'], hName)
            print(hName, end=' ')

            houtPath = "%s/Nonprompt" % (hinPath)
            if houtPath not in houts_np: houts_np[houtPath] = []

            ## open raw_hist root files
            for proc, procInfo in info['processes'].items():
                title = procInfo['title']
                longTitle = procInfo['longTitle'] if 'longTitle' in procInfo else title

#                houtPath = "%s/%s" % (hinPath, title)
#                if houtPath not in houts_np: houts_np[houtPath] = []

                for datasetGroup in procInfo['datasets']:
                    if datasetGroup not in info['dataset']: continue
                    physProcName = datasetGroup.split('.')[1]

                    xsec = 1.0
                    if physProcName in info['crosssection']:
                        xsec_np = info['crosssection'][physProcName]
                        xsec = -1.0 * xsec_np
                    elif title != "Data":
                        print("Could not find", physProcName, "in the cross section list. setting it to be 1.0")

                    hout = None
                    nEvents = 0.
                    for ds in list(info['dataset'][datasetGroup].keys()):
                        fName = "%s/%s/%s.root" % (mode_np, proc, ds[1:].replace('/', '.'))
                        if fName not in fins_np: continue
                        fin = fins_np[fName]
                        hin = fin.Get(hinPath)
                        if hin == None: continue

                        nEvents += fin.Get("hCutFlow").GetBinContent(1)

                        if hout == None:
                            hout = hin.Clone()
                            longTitle = "Nonprompt"
                            hout.SetTitle(longTitle)
                            hout.Reset()
                        hout.Add(hin)
                    if hout == None: continue
#                    if nEvents != 0.: hout.Scale(19500.0*xsec/nEvents)
#                    else:
#                        scale_weight = 1.0
#                        hout.Scale(scale_weight)
                    if ( "Zut" in proc ) or ( "Zct" in proc ): hout.Scale(0)
                    elif ( "Data" in proc ): hout.Scale(1.0)
                    else: hout.Scale(59800.0*xsec/nEvents)

                    houts_np[houtPath].append(hout)
        print("")

    print("Writing output...")
    for houtPath, hists in sorted(iter(houts.items()), key=lambda x: x[0]):
        if len(hists) == 0: continue
        dout = fout.GetDirectory(os.path.dirname(houtPath))
        if dout == None: dout = makedirs(fout, os.path.dirname(houtPath))
        dout.cd()

        hout = None
        for h in hists:
            if hout == None: hout = h
            else: hout.Add(h)
        if hout == None: continue
        hout.SetName(os.path.basename(houtPath))
        hout.Write()

    for houtPath, hists in sorted(iter(houts_np.items()), key=lambda x: x[0]):
        if len(hists) == 0: continue
        dout = fout.GetDirectory(os.path.dirname(houtPath))
        if dout == None: dout = makedirs(fout, os.path.dirname(houtPath))
        dout.cd()

        hout = None
        for h in hists:
            if hout == None: hout = h
            else: hout.Add(h)
        if hout == None: continue
        hout.SetName(os.path.basename(houtPath))
        hout.Write()

    print("Closing...")
    fout.Close()
    for fin in list(fins.values()): fin.Close()
    print("done.")

#os.system("hadd -f %s/All.root %s" % (odName, " ".join(["%s/%s.root" % (odName, mode) for mode in modes])))

print("Finishing...")

