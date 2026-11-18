#!/usr/bin/env python3
import os
"""
modes = [
    "ElElEl",
    "MuMuMu",
    "TaTaTa",
    "ElMuMu",
    "ElTaTa",
    "MuElEl",
    "MuTaTa",
    "TaElEl",
    "TaMuMu",
    "NPLElElEl",
    "NPLMuMuMu",
    "NPLTaTaTa",
    "NPLElMuMu",
    "NPLElTaTa",
    "NPLMuElEl",
    "NPLMuTaTa",
    "NPLTaElEl",
    "NPLTaMuMu",
]

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
"""
modes = [
    "ElElEl",
    "MuMuMu",
    "ElMuMu",
    "MuElEl",
    "NPLElElEl",
    "NPLMuMuMu",
    "NPLElMuMu",
    "NPLMuElEl",
]

baseDir = "NanoAOD/2024"

#nFilePerJob = int(os.environ["NFILE"]) if "NFILE" in os.environ else 1
nFilePerJobSD = int(os.environ["NFILE"]) if "NFILE" in os.environ else 5
nFilePerJobMC = int(os.environ["NFILE"]) if "NFILE" in os.environ else 20
nFilePerJobDATA = int(os.environ["NFILE"]) if "NFILE" in os.environ else 20
if not os.path.exists("submit_2024"): os.mkdir("submit_2024")

import yaml
from glob import glob
procInfo = yaml.safe_load(open("config/grouping.yaml").read())["processes"]
datasetInfo = {}
for f in glob("%s/*.yaml" % baseDir):
    for datasetGroup, dataset in yaml.safe_load(open(f).read())['dataset'].items():
        datasetInfo[datasetGroup] = list(dataset.keys())

fLists = []
toSubmit = {}
for proc in procInfo:
    for datasetGroup in procInfo[proc]['datasets']:
        if datasetGroup not in datasetInfo: continue

        for dataset in datasetInfo[datasetGroup]:
            fLists = glob("%s/*/%s/%s.txt" % (baseDir, datasetGroup, dataset.replace('/','.')[1:]))
            for fList in fLists:
                if fList not in toSubmit: toSubmit[fList] = []
                toSubmit[fList].extend(procInfo[proc]['modes'] if 'modes' in procInfo[proc] else modes)
for k, v in list(toSubmit.items()):
    toSubmit[k] = list(dict.fromkeys(v))

from math import ceil
for fList, modes in toSubmit.items():
    for mode in modes:
        nFiles = len([x for x in open(fList).readlines() if len(x) != 0 and x[0] != '#'])
#        nJobs = ceil(1.*nFiles/nFilePerJob)
        nJobsSD = ceil(1.*nFiles/nFilePerJobSD)
        nJobsMC = ceil(1.*nFiles/nFilePerJobMC)
        nJobsDATA = ceil(1.*nFiles/nFilePerJobDATA)

        jobName = "%s.%s" % (mode, os.path.basename(fList)[:-4])
#        cmd = "cd submit_2024; create-batch3_mine bash ../01_prod_ntuple.sh %s ../%s %d --jobName %s -T --nJobs %d" % (mode, fList, nFilePerJob, jobName, nJobs)
        if ("FCNC" in jobName) or ("TZQB" in jobName) or ("TTto" in jobName) or ("WWto" in jobName) or ("WZto" in jobName) or ("ZZto" in jobName) or ("TWZ" in jobName):
            cmd = "cd submit_2024; create-batch3_mine bash ../01_prod_ntuple.sh %s ../%s %d --jobName %s -T --nJobs %d" % (mode, fList, nFilePerJobSD, jobName, nJobsSD)
        elif "NANOAODSIM" in jobName:
            cmd = "cd submit_2024; create-batch3_mine bash ../01_prod_ntuple.sh %s ../%s %d --jobName %s -T --nJobs %d" % (mode, fList, nFilePerJobMC, jobName, nJobsMC)
        else:
            cmd = "cd submit_2024; create-batch3_mine bash ../01_prod_ntuple.sh %s ../%s %d --jobName %s -T --nJobs %d" % (mode, fList, nFilePerJobDATA, jobName, nJobsDATA)
        os.system(cmd)

