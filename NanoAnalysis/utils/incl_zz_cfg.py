import os
import ROOT

sample_paths = lambda base_dir, proc_list: {proc: os.path.join(base_dir, proc, "ZZ4lAnalysis.root") for proc in proc_list}

lumi_2022 = 35.08424e3 # 1/pb 2022 C-G  (= 35.181930231/fb of full 355100_362760 Golden json - 0.097685694 of eraB that we don't use
lumi_CD   = 8.077e3 # 1/pb
lumi_EFG  = 27.007e3 # 1/pb

path_data_2022EFG = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_EFG/"

path_MC = "/eos/user/i/iehle/Analysis/inclZZ_MC_2022EE_cosThetaIncluded"

ang_vars = ["eta", "cos", "phi"]

# "To first order L should always include taus but you should check." -- Andrew
sample_info = {
    "q#bar{q}#rightarrow ZZ,Z#gamma*": dict(
        fill_color = "#99ccff",
        line_color = "#000099",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC, ["ZZTo4l"]),
                lum     = lumi_EFG
            )
        }
    ),
    "gg#rightarrow ZZ,Z#gamma*": dict(
        fill_color = "#4b78ff",
        line_color = "#000099",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC, ['ggTo2e2mu_Contin_MCFM701']), # Missing 4e(mu)
                lum     = lumi_EFG
            )
        }
    ),
    "WZ": dict(
        fill_color = "#ba03d4",
        line_color = "#8b1ba2",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC, ["WZto3LNu"]),
                lum     = lumi_EFG
            )
        }
    ),
    "DY": dict(
        fill_color = "#89a223",
        line_color = "#344d0e",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC, ["DYJetsToLL"]),
                lum     = lumi_EFG
            )
        }
    ),
    "t#bar{t}": dict(
        fill_color = "#cc0000",
        line_color = "#610018",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC, ['TTto2L2Nu_TuneCP5_13p6TeV_powheg']),
                lum     = lumi_EFG
            )
        }
    ),
    "VVV": dict(
        fill_color = "#eeeb53",
        line_color = "#ffbf00",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC, ["WWZ", "WZZ", "ZZZ"]),
                lum     = lumi_EFG
            )
        }
    ),
    "H(125)": dict(
        fill_color = "#ff9b9b",
        line_color = "#cc0000",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC, ["ggH125"]),
                lum     = lumi_EFG
            )
        }              
    )
}

# Not that ZpX is modelled as a landau function
# normalized to current lumi, so it can only be
# currently filled for mass histograms
Zpx = {
    "Z+X": dict(
        fill_color = "#669966",
        line_color = "#003300",
        eras       = {
            "2022_EE": dict(
                samples = None,
                lum     = lumi_EFG
            )
        }
    )
}

# ZZ mass: (70, 1002, 4)
# Z  mass: (40, 120, 2)

hist_info = dict(
        prop  = "pt",
        which = "Z2",
        reg   = "SR",
        xlow  = 0.,
        xhigh = 1500.,
        step  = 10,
        weight= False # Current weights seem to be false, need to check # only affects single lepton distributions
    )

plot_info = dict(
    ratio      = True,
    x_range    = (0, 1500),
    y_min      = 0,
    x_title    = r"${p_T}^{2l}$ $\text{GeV}$",
    logx       = False,
    logy       = False,
    xlabels    = None,
    legend_loc = (0.72, 0.70, 0.94, 0.92),
    out_file   = "{}_{}_step_{}_withRatio_sqrt".format(hist_info["which"], hist_info["prop"], hist_info["step"]),
    format     = "png"
)
plot_info.update(
    dict(
        y_title = "Events" if hist_info["prop"] in ang_vars else "Events / {} GeV".format(int(hist_info["step"]))
    )
)

all_info = dict(
    proc_info = sample_info,
    hist_info = hist_info,
    plot_info = plot_info,
    lumin     = lumi_EFG
)