import os
import ROOT

sample_paths = lambda base_dir, proc_list: {proc: os.path.join(base_dir, proc, "ZZ4lAnalysis.root") for proc in proc_list}

path_data_2022CD  = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_CD/"
path_data_2022EFG = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_EFG/"

path_MC_2018      = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/MC2018/" # For ggZZ
path_MC_2022      = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022/"
path_MC_2022EE    = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022EE/"

lumi_2022 = 35.08424e3 # 1/pb 2022 C-G  (= 35.181930231/fb of full 355100_362760 Golden json - 0.097685694 of eraB that we don't use
lumi_CD   = 8.077e3 # 1/pb
lumi_EFG  = 27.007e3 # 1/pb

#path_MC = "/eos/user/i/iehle/Analysis/secondSamples/PROD_inclZZTo4LSamples_2022EE_MC_176cc0be"
#path_VVV = "/eos/user/i/iehle/Analysis/firstSamples/PROD_inclZZTo4LSamples_2022EE_MC_68597165d9e/"
ang_vars = ["eta", "cos", "phi"]

sample_info = {
    "EW": dict(
        fill_color = "#0331B9",
        line_color = "#000099",
        eras       = {
            "2022_CD": dict(
                samples = sample_paths(path_MC_2022, ["WWZ", "WZZ", "ZZZ", "TTWW", "TTZZ"]),
                lum     = lumi_CD
            ),
            "2022_EFG": dict(
                samples = sample_paths(path_MC_2022EE, ["WWZ", "WZZ", "ZZZ", "TTWW", "TTZZ"]),
                lum     = lumi_EFG
            )
        }
    ),
    "q#bar{q}#rightarrow ZZ,Z#gamma*": dict(
        fill_color = "#99ccff",
        line_color = "#000099",
        eras       = {
            "2022_CD": dict(
                samples = sample_paths(path_MC_2022, ["ZZTo4l"]),
                lum     = lumi_CD
            ),
            "2022_EFG": dict(
                samples = sample_paths(path_MC_2022EE, ["ZZTo4l"]),
                lum     = lumi_EFG
            )
        }
    ),
    "gg#rightarrow ZZ,Z#gamma*": dict(
        fill_color = "#4b78ff",
        line_color = "#000099",
        eras       = {
            "2018": dict(
                samples = sample_paths(path_MC_2018, ["ggTo4mu_Contin_MCFM701", "ggTo4e_Contin_MCFM701", "ggTo4tau_Contin_MCFM701","ggTo2e2mu_Contin_MCFM701", "ggTo2e2tau_Contin_MCFM701", "ggTo2mu2tau_Contin_MCFM701"]),
                lum     = lumi_2022
            ),
        }                
    ),
    "H(125)": dict(
        fill_color = "#ff9b9b",
        line_color = "#cc0000",
        eras       = {
            "2022_CD": dict(
                samples = sample_paths(path_MC_2022, ["VBFH125", "ggH125", "WplusH125", "WHminus125", "ZH125", "ttH125", "bbH125"]),
                lum     = lumi_CD
            ),
            "2022_EFG": dict(
                samples = sample_paths(path_MC_2022EE, ["VBFH125", "ggH125", "WplusH125", "WHminus125", "ZH125", "ttH125", "bbH125"]),
                lum     = lumi_EFG
            )
        }              
    ),
    "Z+X": dict(
        fill_color = "#669966",
        line_color = "#003300",
        eras       = {
            "2022": dict(
                samples = None,
                lum     = lumi_2022
            )
        }
    ),
    "Data": dict(
        line_color   = ROOT.kBlack,
        marker_style = 20,
        marker_size  = 0.9,
        eras = {
            "2022_CD": dict(
                samples = {"2022_CD": os.path.join(path_data_2022CD, "ZZ4lAnalysis.root")}
            ),
            "2022_EFG": dict(
                samples = {"2022_EFG": os.path.join(path_data_2022EFG, "ZZ4lAnalysis.root")}
            )
        }
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
        legend_loc = (0.72, 0.70, 0.94, 0.92),
        format  = "png"
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
