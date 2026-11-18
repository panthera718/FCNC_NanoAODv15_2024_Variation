#!/usr/bin/env python3
import yaml
import sys, os
import subprocess
import shlex
from glob import glob
from multiprocessing import Pool, cpu_count

def runCmd(cmd):
    os.system(cmd)

def shell_join(args):
    return " ".join(shlex.quote(str(x)) for x in args)

def dataset_name_from_txt(txt_path):
    """Convert a NanoAOD/2024/.../*.txt file name to a CMS dataset full name.

    Example:
      A.B.C.txt -> /A/B/C
    """
    base = os.path.basename(txt_path)
    if not base.endswith(".txt"):
        return None
    return "/" + base[:-4].replace(".", "/")

def build_unskimmed_filelist_map():
    """Map CMS dataset full names to their unskimmed NanoAOD txt file lists."""
    ret = {}
    for txt in glob("NanoAOD/2024/*/*/*.txt"):
        dset = dataset_name_from_txt(txt)
        if dset is None:
            continue
        if dset in ret:
            print("[WARNING] Duplicate unskimmed file list for dataset:", dset)
            print("          old:", ret[dset])
            print("          new:", txt)
        ret[dset] = txt
    return ret

if __name__ == '__main__':
    pool = Pool(processes=min(cpu_count(), 20))

    ## Load all information
    histSetFile = "config/histogramming.yaml"
    datasetConfig = {}
    for f in glob("NanoAOD/2024/*.yaml"): datasetConfig.update(yaml.safe_load(open(f))['dataset'])
    procInfo = yaml.safe_load(open("config/grouping.yaml"))["processes"]
    #systConfig = yaml.load(open("config/systematics.yaml"))["systematics"] ## FIXME: To be implemented later

    unskimmedFileLists = build_unskimmed_filelist_map()
    print("[INFO] Found %d unskimmed NanoAOD txt file lists." % len(unskimmedFileLists))

    ## Build CMS official dataset full name to process group mapping
    dsetfullname2procs = {}
    for proc in procInfo:
        for dataset in procInfo[proc]['datasets']:
            if dataset not in datasetConfig: continue
            for dsetfullname in datasetConfig[dataset]:
                if dsetfullname not in dsetfullname2procs: dsetfullname2procs[dsetfullname] = []
                dsetfullname2procs[dsetfullname].append(proc)

    if os.path.exists("raw_hist_2024"):
      path = "raw_hist_2024"
      subprocess.call(["rm", "-rf", path])
    ress = []
    for din in glob("ntuple_2024/*/*/*"):
        mode, dataset = din.split('/')[2:]
        dataset = '/'+dataset.replace('.', '/')
        if dataset not in dsetfullname2procs: continue

        # Only MC NanoAODSIM needs the unskimmed file list for sum sign(genWeight).
        # Data has no genWeight and should keep the old behavior.
        unskimmedFileList = None
        if dataset.endswith("/NANOAODSIM"):
            unskimmedFileList = unskimmedFileLists.get(dataset)
            if unskimmedFileList is None:
                print("[WARNING] No unskimmed NanoAOD txt file list found for MC dataset:", dataset)

        for proc in dsetfullname2procs[dataset]:
            config = procInfo[proc]
            if 'modes' in config and mode not in config['modes']: continue

            dout = "raw_hist_2024/%s/%s" % (mode, proc)

            #cut = config['cut'] if 'cut' in config else '1'
            cut = '1'
            if mode[:3] == "NPL":
                weight = config['FRweight']
            else:
                weight = config['weight'] if 'weight' in config else '1'

            cmd = ["tzwi-makehistograms", cut, weight, histSetFile, din, dout]
            if unskimmedFileList is not None:
                cmd.append(unskimmedFileList)

            res = pool.apply_async(os.system, (shell_join(cmd),))
            ress.append(res)

    for r in ress: r.get()
