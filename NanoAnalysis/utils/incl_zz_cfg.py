import os
import ROOT

sample_paths = lambda base_dir, proc_list: {proc: os.path.join(base_dir, proc, "ZZ4lAnalysis.root") for proc in proc_list}

lumi_2022 = 35.08424e3 # 1/pb 2022 C-G  (= 35.181930231/fb of full 355100_362760 Golden json - 0.097685694 of eraB that we don't use
lumi_CD   = 8.077e3 # 1/pb
lumi_EFG  = 27.007e3 # 1/pb

path_hZZ_2022EE    = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022EE/"
path_data_2022EFG = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_EFG/"

path_MC = "/eos/user/i/iehle/Analysis/secondSamples/PROD_inclZZTo4LSamples_2022EE_MC_176cc0be"
path_VVV = "/eos/user/i/iehle/Analysis/firstSamples/PROD_inclZZTo4LSamples_2022EE_MC_68597165d9e/"

# APPROXIMATELY 25% OF ZZTO4L JOBS FAILED IN THESE SAMPLES
path_MC_third = "/eos/user/i/iehle/Analysis/thirdSamples/PROD_inclZZTo4LSamples_2022EE_MC_289047f8"

ang_vars = ["eta", "cos", "phi"]

sample_info = {
    "q#bar{q}#rightarrow ZZ,Z#gamma*": dict(
        fill_color = "#99ccff",
        line_color = "#000099",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC_third, ["ZZTo4l"]),
                lum     = lumi_EFG
            )
        }
    ),
    "gg#rightarrow ZZ,Z#gamma*": dict(
        fill_color = "#4b78ff",
        line_color = "#000099",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC_third, ['ggTo2e2mu_Contin_MCFM701', 'ggTo2e2tau_Contin_MCFM701', 'ggTo2mu2tau_Contin_MCFM701']), # need 4l,
                lum     = lumi_EFG
            )
        }
    ),
    "WZ": dict(
        fill_color = "#ba03d4",
        line_color = "#8b1ba2",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC_third, ["WZto3LNu"]),
                lum     = lumi_EFG
            )
        }
    ),
    "DY": dict(
        fill_color = "#89a223",
        line_color = "#344d0e",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC_third, ["DYJetsToLL"]),
                lum     = lumi_EFG
            )
        }
    ),
    "t#bar{t}": dict(
        fill_color = "#cc0000",
        line_color = "#610018",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC_third, ['TTto2L2Nu_TuneCP5_13p6TeV_powheg']),
                lum     = lumi_EFG
            )
        }
    ),
    "VVV": dict(
        fill_color = "#eeeb53",
        line_color = "#ffbf00",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC_third, ["WWZ", "WZZ", "ZZZ"]),
                lum     = lumi_EFG
            )
        }
    ),
    "H(125)": dict(
        fill_color = "#ff9b9b",
        line_color = "#cc0000",
        eras       = {
            "2022_EE": dict(
                samples = sample_paths(path_MC_third, ["ggH125"]),
                lum     = lumi_EFG
            )
        }              
    ),
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

hist_info = dict(
        prop  = "pt",
        which = "Z1",
        reg   = "SR",
        xlow  = 0.,
        xhigh = 1500.,
        step  = 4.,
    )

plot_info = dict(
    x_range    = (0., 1500.),
    y_min      = 0,
    x_title    = "p_T^{#it{Z1}} (GeV)",
    logx       = True,
    logy       = True,
    xlabels    = [80, 100, 200, 300, 400, 500],
    legend_loc = (0.72, 0.70, 0.94, 0.92),
    out_file   = "{}_{}_step_{}_thirdSamples_TESTING".format(hist_info["which"], hist_info["prop"], hist_info["step"]),
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