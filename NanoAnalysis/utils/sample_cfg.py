import os

sample_paths = lambda base_dir, proc_list: {proc: os.path.join(base_dir, proc, "ZZ4lAnalysis.root") for proc in proc_list}

#path_MC = "/eos/user/i/iehle/Analysis/secondSamples/PROD_inclZZTo4LSamples_2022EE_MC_176cc0be"
path_MC = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022EE/'
ang_vars = ["eta", "cos", "phi"]

sample_info = {
    "EW": dict(
        samples    = sample_paths(path_MC, ["WWZ", "WZZ", "ZZZ", "TTWW", "TTZZ"]),
        fill_color = "#0331B9",
        line_color = "#000099",
    ),
    "q#bar{q}#rightarrow ZZ,Z#gamma*": dict(
        samples    = sample_paths(path_MC, ["ZZTo4l"]),
        fill_color = "#99ccff",
        line_color = "#000099",
    ),
    "H(125)": dict(
        samples    = sample_paths(path_MC, ["VBFH125", "ggH125", "WplusH125", "WHminus125", "ZH125", "ttH125", "bbH125"]),
        fill_color = "#ff9b9b",
        line_color = "#cc0000",
    ),
    "Z+X": dict(
        samples    = None,
        fill_color = "#669966",
        line_color = "#003300"
    )
}

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
    proc_info = sample_info,
    hist_info = hist_info
)
