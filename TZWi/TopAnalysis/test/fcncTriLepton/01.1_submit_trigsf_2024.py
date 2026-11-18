#!/usr/bin/env python3
import os
from glob import glob
from math import ceil

MODES = ["ElElEl", "ElMuMu", "MuElEl", "MuMuMu"]
BASE_DIR = os.environ.get("TRIGSF_BASEDIR", "NanoAOD/2024")
SUBMIT_DIR = os.environ.get("TRIGSF_SUBMIT_DIR", "submit_2024_trigSF")
PROD_SCRIPT = os.environ.get("TRIGSF_PROD_SCRIPT", "01_prod_trigsf_ntuple_2024.sh")
CREATE_BATCH = os.environ.get("TRIGSF_CREATE_BATCH", "create-batch3_mine")

NFILE_DATA = int(os.environ.get("NFILE_DATA", os.environ.get("NFILE", "5")))
NFILE_MC = int(os.environ.get("NFILE_MC", os.environ.get("NFILE", "5")))

# For an unbiased data trigger-SF measurement, prefer an independent reference
# PD such as MET/JetMET.  EGamma/Muon/MuonEG are the analysis trigger PDs and
# should not be used as the denominator unless you are doing only a validation.
#DEFAULT_DATA_GROUPS = ["Run2024.MET", "Run2024.JetMET"]
DEFAULT_DATA_GROUPS = ["Run2024.JetMET"]
if os.environ.get("TRIGSF_ALLOW_LEPTON_DATA", "0") == "1":
    DEFAULT_DATA_GROUPS += ["Run2024.EGamma", "Run2024.Muon", "Run2024.MuonEG"]

DEFAULT_MC_GROUPS = ["MC2024.WZto3LNu"]

DATA_GROUPS = [x.strip() for x in os.environ.get("TRIGSF_DATASET_GROUPS", ",".join(DEFAULT_DATA_GROUPS)).split(",") if x.strip()]
MC_GROUPS = [x.strip() for x in os.environ.get("TRIGSF_MC_GROUPS", ",".join(DEFAULT_MC_GROUPS)).split(",") if x.strip()]
DATASET_GROUPS = DATA_GROUPS + MC_GROUPS

def list_filelists(base_dir, dataset_group):
    # Works for both:
    #   NanoAOD/2024/RD.Run2024/Run2024.MET/*.txt
    #   NanoAOD/2024/MC.RunIISummer24.central/MC2024.WZto3LNu/*.txt
    pattern = os.path.join(base_dir, "*", dataset_group, "*.txt")
    return sorted(glob(pattern))


def non_comment_line_count(path):
    n = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith("#"):
                n += 1
    return n


def files_per_job_for(flist):
    dataset_group = os.path.basename(os.path.dirname(flist))
    return NFILE_MC if dataset_group.startswith("MC") else NFILE_DATA


def main():
    os.makedirs(SUBMIT_DIR, exist_ok=True)

    all_filelists = []
    missing_groups = []
    for group in DATASET_GROUPS:
        found = list_filelists(BASE_DIR, group)
        if not found:
            missing_groups.append(group)
        all_filelists.extend(found)

    if missing_groups:
        print("[WARN] No file lists found for these dataset groups:")
        for group in missing_groups:
            print("       - %s" % group)
        if any(g in ("Run2024.MET", "Run2024.JetMET") for g in missing_groups):
            print("[WARN] Add an independent reference-trigger primary dataset, e.g. Run2024.MET or Run2024.JetMET, before measuring an unbiased data trigger SF.")

    if not all_filelists:
        raise RuntimeError("No input file lists found. Check TRIGSF_BASEDIR and TRIGSF_DATASET_GROUPS/TRIGSF_MC_GROUPS.")

    for flist in all_filelists:
        n_files = non_comment_line_count(flist)
        if n_files <= 0:
            print("[WARN] Empty file list: %s" % flist)
            continue

        n_file_per_job = files_per_job_for(flist)
        n_jobs = int(ceil(float(n_files) / float(n_file_per_job)))
        dataset0 = os.path.basename(flist)[:-4]

        for mode in MODES:
            job_name = "%s.%s" % (mode, dataset0)
            cmd = (
                "cd {submit_dir}; "
                "{create_batch} bash ../{prod_script} {mode} ../{flist} {nfile} "
                "--jobName {job_name} -T --nJobs {njobs}"
            ).format(
                submit_dir=SUBMIT_DIR,
                create_batch=CREATE_BATCH,
                prod_script=PROD_SCRIPT,
                mode=mode,
                flist=flist,
                nfile=n_file_per_job,
                job_name=job_name,
                njobs=n_jobs,
            )
            print(cmd)
            os.system(cmd)


if __name__ == "__main__":
    main()
