#!/usr/bin/env python3
import os
from glob import glob
from math import ceil
import yaml

# Dedicated submission script for Run-3 2024 NanoAODv15 promptMVA fake-rate NTuples.
#
# Nominal fake-rate measurement strategy:
#   Run2024.JetMET -> FREl
#   Run2024.JetMET -> FRMu
#
# Rationale:
#   The fake-rate measurement should be selected with probe-lepton-independent
#   jet triggers. EGamma/Muon/MuonEG lepton streams are not used by default for
#   the central fake-rate map because their lepton triggers bias the loose->tight
#   fake-rate measurement.

BASE_DIR = os.environ.get("BASE_DIR", "NanoAOD/2024")
SUBMIT_DIR = os.environ.get("SUBMIT_DIR", "submit_2024_FR")
PROD_SCRIPT = os.environ.get("PROD_SCRIPT", "01_prod_ntuple_FR.sh")

# Use a moderate default. Increase only after a small production test succeeds.
NFILE_DATA = int(os.environ.get("NFILE_DATA", os.environ.get("NFILE", "5")))
NFILE_MC = int(os.environ.get("NFILE_MC", os.environ.get("NFILE", "5")))
INCLUDE_FR_MC = os.environ.get("INCLUDE_FR_MC", "0") == "1"
INCLUDE_LEPTON_STREAMS = os.environ.get("INCLUDE_LEPTON_STREAMS", "0") == "1"
DRY_RUN = os.environ.get("DRY_RUN", "0") == "1"

# Nominal data mapping: both electron and muon fake rates are measured from JetMET.
DATA_FR_MODES = {
    "Run2024.JetMET": ["FREl", "FRMu"],
}

# Optional fallback / cross-check only, not nominal.
if INCLUDE_LEPTON_STREAMS:
    DATA_FR_MODES.update({
        "Run2024.EGamma": ["FREl"],
        "Run2024.Muon": ["FRMu"],
    })

# Optional MC samples are useful only for validation, closure, or prompt-contamination studies.
# They are not used by default to build the central data fake-rate map.
OPTIONAL_MC_DATASETS = {
    "MC2024.TTto2L2Nu",
    "MC2024.DYGto2LG-1Jets_Bin-MLL-4to50",
    "MC2024.DYGto2LG-1Jets_Bin-MLL-50",
    "MC2024.WZto3LNu",
    "MC2024.ZZto4L",
    "MC2024.TTLL_Bin-MLL-4to50",
    "MC2024.TTLL_Bin-MLL-50",
}

os.makedirs(SUBMIT_DIR, exist_ok=True)

# Read dataset groups from NanoAOD/2024/*.yaml using the same convention as the
# analysis submitter. The YAML values are dataset names such as
# /JetMET0/Run2024C-MINIv6NANOv15-v1/NANOAOD.
dataset_info = {}
for yaml_file in glob(os.path.join(BASE_DIR, "*.yaml")):
    with open(yaml_file) as f:
        content = yaml.safe_load(f) or {}
    for dataset_group, dataset_dict in content.get("dataset", {}).items():
        dataset_info[dataset_group] = list(dataset_dict.keys())

jobs = []
for dataset_group, datasets in sorted(dataset_info.items()):
    if dataset_group in DATA_FR_MODES:
        modes = DATA_FR_MODES[dataset_group]
        nfile = NFILE_DATA
    elif INCLUDE_FR_MC and dataset_group in OPTIONAL_MC_DATASETS:
        modes = ["FREl", "FRMu"]
        nfile = NFILE_MC
    else:
        continue

    for dataset in datasets:
        file_name = dataset.replace("/", ".")[1:] + ".txt"
        pattern = os.path.join(BASE_DIR, "*", dataset_group, file_name)
        for file_list in sorted(glob(pattern)):
            with open(file_list) as f:
                nfiles = sum(1 for line in f if line.strip() and not line.lstrip().startswith("#"))
            if nfiles == 0:
                print(f"[submit_FR] skip empty file list: {file_list}")
                continue
            njobs = int(ceil(float(nfiles) / float(nfile)))
            for mode in modes:
                job_name = f"{mode}.{os.path.basename(file_list)[:-4]}"
                cmd = (
                    f"cd {SUBMIT_DIR}; "
                    f"create-batch3_mine bash ../{PROD_SCRIPT} {mode} ../{file_list} {nfile} "
                    f"--jobName {job_name} -T --nJobs {njobs}"
                )
                jobs.append(cmd)

print(f"[submit_FR] BASE_DIR={BASE_DIR}")
print(f"[submit_FR] SUBMIT_DIR={SUBMIT_DIR}")
print(f"[submit_FR] PROD_SCRIPT={PROD_SCRIPT}")
print(f"[submit_FR] NFILE_DATA={NFILE_DATA}")
print(f"[submit_FR] INCLUDE_LEPTON_STREAMS={int(INCLUDE_LEPTON_STREAMS)}")
print(f"[submit_FR] INCLUDE_FR_MC={int(INCLUDE_FR_MC)}")
print(f"[submit_FR] active DATA_FR_MODES={DATA_FR_MODES}")
print(f"[submit_FR] found {len(jobs)} batch submissions")
for cmd in jobs:
    print(cmd)
    if not DRY_RUN:
        os.system(cmd)
