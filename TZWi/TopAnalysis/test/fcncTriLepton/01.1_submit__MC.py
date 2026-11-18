#!/usr/bin/env python3
import os
from glob import glob
from math import ceil

import yaml

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

baseDir = os.environ.get("BASE_DIR", "NanoAOD/2024__MC")
submitDir = os.environ.get("SUBMIT_DIR", "submit_2024__MC")
prodScript = os.environ.get("PROD_SCRIPT", "01_prod_ntuple.sh")
dryRun = os.environ.get("DRY_RUN", "0") == "1"

# Optional explicit override.
#
# If this environment variable is omitted, 01_prod_ntuple.sh uses its own
# DEFAULT_FAKERATE_MAP_MODE setting.  If it is set, it is passed as the fourth
# explicit argument to NPL jobs; create-batch3_mine appends JOBNUMBER last.
#
# Accepted aliases:
#   split     : split, btag, btag_split, 1, true, yes
#   inclusive : inclusive, inc, 0, false, no
def normalize_fr_map_mode(value):
    value = str(value or "").strip().lower()
    if not value:
        return ""
    if value in ("split", "btag", "btag_split", "btag-split", "1", "true", "yes"):
        return "split"
    if value in ("inclusive", "inc", "0", "false", "no"):
        return "inclusive"
    raise RuntimeError(
        "Unsupported FAKERATE_MAP_MODE='%s'. Use split or inclusive." % value
    )


fakeRateMapMode = normalize_fr_map_mode(
    os.environ.get("FAKERATE_MAP_MODE", "")
)

nFilePerJobSD = int(os.environ.get("NFILE", "5"))
nFilePerJobMC = int(os.environ.get("NFILE", "10"))
nFilePerJobDATA = int(os.environ.get("NFILE", "10"))

os.makedirs(submitDir, exist_ok=True)

procInfo = yaml.safe_load(
    open("config/grouping.yaml").read()
)["processes"]

datasetInfo = {}
for yamlFile in glob("%s/*.yaml" % baseDir):
    content = yaml.safe_load(open(yamlFile).read()) or {}
    for datasetGroup, dataset in content.get("dataset", {}).items():
        datasetInfo[datasetGroup] = list(dataset.keys())

toSubmit = {}
for proc in procInfo:
    for datasetGroup in procInfo[proc]["datasets"]:
        if datasetGroup not in datasetInfo:
            continue

        for dataset in datasetInfo[datasetGroup]:
            fileLists = glob(
                "%s/*/%s/%s.txt"
                % (
                    baseDir,
                    datasetGroup,
                    dataset.replace("/", ".")[1:],
                )
            )
            for fileList in fileLists:
                if fileList not in toSubmit:
                    toSubmit[fileList] = []
                toSubmit[fileList].extend(
                    procInfo[proc]["modes"]
                    if "modes" in procInfo[proc]
                    else modes
                )

for fileList, configuredModes in list(toSubmit.items()):
    toSubmit[fileList] = list(dict.fromkeys(configuredModes))

print("[01.1_submit.py] FAKERATE_MAP_MODE=%s" % (
    fakeRateMapMode if fakeRateMapMode else
    "<use DEFAULT_FAKERATE_MAP_MODE from 01_prod_ntuple.sh>"
))
print("[01.1_submit.py] SUBMIT_DIR=%s" % submitDir)
print("[01.1_submit.py] DRY_RUN=%d" % int(dryRun))

for fileList, configuredModes in toSubmit.items():
    for mode in configuredModes:
        nFiles = len([
            line
            for line in open(fileList).readlines()
            if line.strip() and not line.lstrip().startswith("#")
        ])

        nJobsSD = int(ceil(1.0 * nFiles / nFilePerJobSD))
        nJobsMC = int(ceil(1.0 * nFiles / nFilePerJobMC))
        nJobsDATA = int(ceil(1.0 * nFiles / nFilePerJobDATA))

        jobName = "%s.%s" % (
            mode,
            os.path.basename(fileList)[:-4],
        )

        if (
            ("FCNC" in jobName)
            or ("TZQB" in jobName)
            or ("TTto" in jobName)
            or ("WWto" in jobName)
            or ("WZto" in jobName)
            or ("ZZto" in jobName)
            or ("TWZ" in jobName)
        ):
            nFilePerJob = nFilePerJobSD
            nJobs = nJobsSD
        elif "NANOAODSIM" in jobName:
            nFilePerJob = nFilePerJobMC
            nJobs = nJobsMC
        else:
            nFilePerJob = nFilePerJobDATA
            nJobs = nJobsDATA

        # Only NPL jobs use FakeRateProducer.  When an explicit mode is supplied,
        # pass it before the process number that create-batch3_mine appends.
        frModeArgument = ""
        if mode.startswith("NPL") and fakeRateMapMode:
            frModeArgument = " %s" % fakeRateMapMode
            jobName += ".FRMap_%s" % fakeRateMapMode

        cmd = (
            "cd {submit}; "
            "create-batch3_mine bash ../{script} "
            "{mode} ../{filelist} {nfile}{frmode} "
            "--jobName {jobname} -T --nJobs {njobs}"
        ).format(
            submit=submitDir,
            script=prodScript,
            mode=mode,
            filelist=fileList,
            nfile=nFilePerJob,
            frmode=frModeArgument,
            jobname=jobName,
            njobs=nJobs,
        )

        print(cmd)
        if not dryRun:
            ret = os.system(cmd)
            if ret != 0:
                raise RuntimeError(
                    "create-batch3_mine failed with status %s: %s"
                    % (ret, cmd)
                )
