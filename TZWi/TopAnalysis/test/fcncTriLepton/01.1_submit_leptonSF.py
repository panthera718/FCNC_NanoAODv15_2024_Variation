#!/usr/bin/env python3
import os
from glob import glob
from math import ceil

import yaml

"""
Dedicated submission script for Run-3 2024 NanoAODv15 promptMVA lepton-SF NTuples.

Run from TZWi/TopAnalysis/test/fcncTriLepton:
    python3 01.1_submit_leptonSF.py

Environment controls:
    BASE_DIR             default: NanoAOD/2024
    SUBMIT_DIR           default: submit_2024_leptonSF
    PROD_SCRIPT          default: 01_prod_ntuple_leptonSF.sh
    NFILE_DATA           default: 10
    NFILE_MC             default: 10
    LEPSF_MODES          default: LepSFEl,LepSFMu
    LEPSF_DATA_GROUPS    default: Run2024.EGamma,Run2024.Muon,Run2024.MuonEG
    LEPSF_MC_GROUPS      default: MC2024.WZto3LNu
    INCLUDE_LEPSF_BKG_MC default: 1
    LEPSF_BKG_GROUPS     default: MC2024.DYGto2LG-1Jets_Bin-MLL-4to50,MC2024.DYGto2LG-1Jets_Bin-MLL-50,MC2024.WGtoLNuG-1Jets,MC2024.ZZto4L,MC2024.TTLL_Bin-MLL-50,MC2024.TTLL_Bin-MLL-4to50,MC2024.TZQB-Zto2L-4FS_Bin-MLL-30,MC2024.TTto2L2Nu,MC2024.WWZ,MC2024.WZZ,MC2024.ZZZ,MC2024.TTH-HtoNon2B_Par-M-125,MC2024.TWZ-Tto2Q-WtoLNu-Zto2L-DR1,MC2024.TWZ-TtoLNu-Wto2Q-Zto2L-DR1,MC2024.TWZ-TtoLNu-WtoLNu-Zto2L-DR1

The default MC is only WZto3LNu because the first central SF validation should
measure the prompt-lepton efficiency ratio Data/WZ prompt MC.  Optional background
MC can be produced later for subtraction or purity studies.
"""

BASE_DIR = os.environ.get("BASE_DIR", "NanoAOD/2024")
SUBMIT_DIR = os.environ.get("SUBMIT_DIR", "submit_2024_leptonSF")
PROD_SCRIPT = os.environ.get("PROD_SCRIPT", "01_prod_ntuple_leptonSF.sh")

NFILE_DATA = int(os.environ.get("NFILE_DATA", os.environ.get("NFILE", "10")))
NFILE_MC = int(os.environ.get("NFILE_MC", os.environ.get("NFILE", "10")))

MODES = [x.strip() for x in os.environ.get("LEPSF_MODES", "LepSFEl,LepSFMu").split(",") if x.strip()]
VALID_MODES = {"LepSFEl", "LepSFMu"}
for mode in MODES:
    if mode not in VALID_MODES:
        raise RuntimeError("Invalid LEPSF mode: %s" % mode)

DATA_GROUPS = [
    x.strip() for x in os.environ.get(
        "LEPSF_DATA_GROUPS",
        "Run2024.EGamma,Run2024.Muon,Run2024.MuonEG",
    ).split(",") if x.strip()
]
MC_GROUPS = [
    x.strip() for x in os.environ.get(
        "LEPSF_MC_GROUPS",
        "MC2024.WZto3LNu",
    ).split(",") if x.strip()
]
INCLUDE_BKG_MC = os.environ.get("INCLUDE_LEPSF_BKG_MC", "1") == "1"
BKG_GROUPS = [
    x.strip() for x in os.environ.get(
        "LEPSF_BKG_GROUPS",
        "MC2024.DYGto2LG-1Jets_Bin-MLL-4to50,MC2024.DYGto2LG-1Jets_Bin-MLL-50,"
        "MC2024.WGtoLNuG-1Jets,MC2024.ZZto4L,"
        "MC2024.TTLL_Bin-MLL-50,MC2024.TTLL_Bin-MLL-4to50,"
        "MC2024.TZQB-Zto2L-4FS_Bin-MLL-30,MC2024.TTto2L2Nu,"
        "MC2024.WWZ,MC2024.WZZ,MC2024.ZZZ,MC2024.TTH-HtoNon2B_Par-M-125,"
        "MC2024.TWZ-Tto2Q-WtoLNu-Zto2L-DR1,"
        "MC2024.TWZ-TtoLNu-Wto2Q-Zto2L-DR1,"
        "MC2024.TWZ-TtoLNu-WtoLNu-Zto2L-DR1",
    ).split(",") if x.strip()
]

os.makedirs(SUBMIT_DIR, exist_ok=True)

# Read the NanoAOD/2024/*.yaml dataset dictionaries using the same convention as
# the existing 01.1_submit.py.  The key is the local dataset group, for example
# Run2024.EGamma or MC2024.WZto3LNu.
dataset_info = {}
for yaml_file in glob(os.path.join(BASE_DIR, "*.yaml")):
    with open(yaml_file) as f:
        content = yaml.safe_load(f)
    for dataset_group, datasets in content.get("dataset", {}).items():
        dataset_info[dataset_group] = list(datasets.keys())

wanted_groups = set(DATA_GROUPS + MC_GROUPS + (BKG_GROUPS if INCLUDE_BKG_MC else []))
missing_groups = sorted([g for g in wanted_groups if g not in dataset_info])
if missing_groups:
    print("[submit_leptonSF] WARNING: these requested groups were not found in %s/*.yaml:" % BASE_DIR)
    for g in missing_groups:
        print("  -", g)

jobs = []
summary = []
for dataset_group in sorted(dataset_info.keys()):
    is_data = dataset_group in DATA_GROUPS
    is_signal_mc = dataset_group in MC_GROUPS
    is_bkg_mc = INCLUDE_BKG_MC and dataset_group in BKG_GROUPS
    if not (is_data or is_signal_mc or is_bkg_mc):
        continue

    nfile = NFILE_DATA if is_data else NFILE_MC
    for dataset in dataset_info[dataset_group]:
        dataset_txt = dataset.replace("/", ".")[1:] + ".txt"
        pattern = os.path.join(BASE_DIR, "*", dataset_group, dataset_txt)
        file_lists = glob(pattern)
        if not file_lists:
            print("[submit_leptonSF] WARNING: no file list matched", pattern)
            continue

        for file_list in file_lists:
            with open(file_list) as f:
                nfiles = sum(1 for line in f if line.strip() and not line.lstrip().startswith("#"))
            if nfiles == 0:
                print("[submit_leptonSF] WARNING: empty file list", file_list)
                continue
            njobs = int(ceil(1.0 * nfiles / nfile))

            for mode in MODES:
                job_name = "%s.%s" % (mode, os.path.basename(file_list)[:-4])
                cmd = (
                    "cd {submit_dir}; "
                    "create-batch3_mine bash ../{prod_script} {mode} ../{file_list} {nfile} "
                    "--jobName {job_name} -T --nJobs {njobs}"
                ).format(
                    submit_dir=SUBMIT_DIR,
                    prod_script=PROD_SCRIPT,
                    mode=mode,
                    file_list=file_list,
                    nfile=nfile,
                    job_name=job_name,
                    njobs=njobs,
                )
                jobs.append(cmd)
                summary.append((dataset_group, os.path.basename(file_list), mode, nfiles, nfile, njobs))

print("[submit_leptonSF] BASE_DIR       =", BASE_DIR)
print("[submit_leptonSF] SUBMIT_DIR     =", SUBMIT_DIR)
print("[submit_leptonSF] PROD_SCRIPT    =", PROD_SCRIPT)
print("[submit_leptonSF] modes          =", ",".join(MODES))
print("[submit_leptonSF] data groups    =", ",".join(DATA_GROUPS))
print("[submit_leptonSF] signal MC      =", ",".join(MC_GROUPS))
print("[submit_leptonSF] background MC  =", ",".join(BKG_GROUPS) if INCLUDE_BKG_MC else "disabled")
print("[submit_leptonSF] batch submissions:", len(jobs))

for row in summary:
    dataset_group, file_name, mode, nfiles, nfile, njobs = row
    print("  {mode:8s} {group:40s} {file:90s} nFiles={nfiles:5d} split={nfile:3d} nJobs={njobs:4d}".format(
        mode=mode,
        group=dataset_group,
        file=file_name,
        nfiles=nfiles,
        nfile=nfile,
        njobs=njobs,
    ))

for cmd in jobs:
    print(cmd)
    os.system(cmd)
