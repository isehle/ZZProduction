import os

def fill_samples(path, combined_samples):
    """Takes a dictionary of combined processes
    and replaces each process list with the path
    to its respective file. Example:
    'EW': ['WWZ','WZZ','ZZZ'] --> 'EW': {'WWZ': <wwz_path>, 'WZZ': <wzz_path>, ...}"""
    
    full_path = lambda base_dir, proc_dir: os.path.join(base_dir, proc_dir, "ZZ4lAnalysis.root")
    
    for proc_info in combined_samples.values():
        if proc_info["procs"] is not None:
            samples = {key: full_path(path, key) for key in proc_info["procs"]}
        else:
            samples = None
        
        proc_info.update({"samples": samples})
        del proc_info["procs"]

    return combined_samples

#path_MC = "/eos/user/i/iehle/Analysis/secondSamples/PROD_inclZZTo4LSamples_2022EE_MC_176cc0be"
path_MC = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022EE/'

combined_samples = {
    "EW": dict(
        procs      = ["WWZ", "WZZ", "ZZZ", "TTWW", "TTZZ"],
        fill_color = "#0331B9",
        line_color = "#000099",
    ),
    "q#bar{q}#rightarrow ZZ,Z#gamma*": dict(
        procs = ["ZZTo4l"],
        fill_color = "#99ccff",
        line_color = "#000099",
    ),
    "H(125)": dict(
        procs = ["VBFH125", "ggH125", "WplusH125", "WHminus125", "ZH125", "ttH125", "bbH125"],
        fill_color = "#ff9b9b",
        line_color = "#cc0000",
    ),
    "Z+X": dict(
        procs = None,
        fill_color = "#669966",
        line_color = "003300"
    )
}

final_samples = fill_samples(path_MC, combined_samples)

ang_vars = ["eta", "cos", "phi"]
hist_info = dict(
        prop  = "mass",
        which = "ZZ",
        reg   = "SR",
        xlow  = 70.,
        xhigh = 300.,
        step  = 4.,
        x_title = "m_{#it{4l}} (GeV)",
        logx    = True,
        logy    = False,
        xlabels = [80, 100, 200, 300, 400, 500],
        legend_loc = (0.72, 0.70, 0.94, 0.92)
    )
hist_info.update(
    dict(
        y_title = "Events" if hist_info["prop"] in ang_vars else "Events / {} GeV".format(int(hist_info["step"]))
    )
)

all_info = dict(
    proc_info = combined_samples,
    hist_info = hist_info
)
