import os
import glob
import ROOT

prod_dir = "MC_2022EE_Test07Feb"

local_path = "/afs/cern.ch/user/i/iehle/polZZTo4l/CMSSW_13_0_16/src/ZZAnalysis/NanoAnalysis/test/prod/{}".format(prod_dir)
eos_dir   = "/eos/user/i/iehle/Analysis/{}".format(prod_dir)

failures   = [d for d in os.listdir(local_path) if "Chunk" in d]

paths = {"Success" : [],
         "Failure" : []}
for fail in failures:
    log_dir = os.path.join(local_path, fail, "log")
    outfile_path = glob.glob("{}/*.out".format(log_dir))[0]
    with open(outfile_path, "r") as outfile:
        if "Writing fake empty file." in outfile.read():
            eos_path = os.path.join(eos_dir, fail)
            if "ZZ4lAnalysis.root" in os.listdir(eos_path):
                paths["Success"].append(fail)
            elif "ZZ4lAnalysis.root.empty" in os.listdir(eos_path):
                paths["Failure"].append(fail)

def get_root_files(fail):
    bad_chars = ["[", "]", " ", "'"]

    cfg_path = os.path.join(local_path, fail, "run_cfg.py")
    with open(cfg_path, "r") as cfg:
        cfg_lines = cfg.readlines()
        ffnames   = [l for l in cfg_lines if l.startswith('setConf("fileNames"')][0]
        ffname    = ffnames.replace('setConf("fileNames",', "").replace(")\n", "")
        f_paths   = "".join([char for char in list(ffname) if char not in bad_chars]).split(",")
    
    return f_paths

def run_test(rfile_path):
    pass

for fail in paths["Failure"]:
    f_paths = get_root_files(fail)