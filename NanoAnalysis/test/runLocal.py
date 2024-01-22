#!/usr/bin/env python3
###
# Example for running the analysis locally, after customizing variables.
# Run with: 
# python runLocal.py
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter

#SampleToRun = "MCsync_Rereco"
#SampleToRun = "MCsync_UL"
#SampleToRun = "Data2022"
#SampleToRun = "MC2022"
SampleToRun = "MC_ZZto4L"

### Customize processing variables
#setConf("runMELA", False)
#setConf("bestCandByMELA", False)
setConf("APPLYMUCORR", False) #NOTE: mu corrections removed for comparision with mini since they are not deterministic on nanoAOD.

## Select specific events to debug
#setConf("preselection","run==316239  && luminosityBlock==226 && event==284613817")

## Force filling K factors and weights (default: all off)
#setConf("APPLY_K_NNLOQCD_ZZGG", 1) # 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
#setConf("APPLY_K_NNLOQCD_ZZQQB", True)
#setConf("APPLY_K_NNLOEW_ZZQQB", True)
#setConf("APPLY_QCD_GGF_UNCERT", True)

setConf("PROCESS_CR", True)
setConf("DEBUG", False)
setConf("SYNCMODE", True) # Force muon resolution correction with fixed +1 sigma smearing

json = None #replace this if needed

mc_samples = {
    "ZZ"             : "/store/mc/Run3Summer22EENanoAODv11/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/126X_mcRun3_2022_realistic_postEE_v1-v2/2810000/06da713c-4abc-4b80-a50b-d38bc3974588.root", # nevents = 1083240,
    "DYJets"         : "/store/mc/Run3Summer22EENanoAODv12/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/NANOAODSIM/forPOG_130X_mcRun3_2022_realistic_postEE_v6-v2/2530000/c8912838-7ac7-4eaa-97a9-808a1d70379f.root", # nevents = 1153196
    "ZpJets"         : None,
    #"WZto3LNu_1Jets" : "/store/mc/Run3Summer22EENanoAODv12/WZto3LNu-1Jets-4FS_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2520000/4b3776db-35c3-4e56-9082-82e89a10efa5.root", #nevents = 199326 # in principal incorrect
    "WZto3LNu"       : "/store/mc/Run3Summer22EENanoAODv11/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/126X_mcRun3_2022_realistic_postEE_v1-v2/2820000/12669125-e7d7-4d14-95c0-18617c9aeaa2.root", # nevents = 195160
    #"WZtoLNu2Q_1Jets": "/store/mc/Run3Summer22EENanoAODv12/WZtoLNu2Q-1Jets-4FS_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/c760afa6-aa32-49da-838f-faf02828209c.root", #nevents = 622612
    "WZtoLNu2Q_1Jets": None, # In principal not needed
    "TTZ"            : None,
    "Higgs"          : "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2520000/17b1cba9-8e8b-480f-88ed-059b3dac2884.root", # nevents = 385237
    "TTTo2L2Nu"      : "/store/mc/Run3Summer22EENanoAODv12/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/30000/66b84e6a-483b-4d00-ad08-d345cc957db1.root", # nevents = 701910
    "WWZ"            : "/store/mc/Run3Summer22EENanoAODv12/WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/30000/8a240e41-ee8e-42cb-9bda-6450fa763326.root", # nevents = 698660
    "WZZ"            : "/store/mc/Run3Summer22EENanoAODv12/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/30000/08125750-33be-4e9c-a02b-97e13fbd2336.root", # nevents = 407049
    "ZZZ"            : "/store/mc/Run3Summer22EENanoAODv12/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/50000/bdb14d5b-e93f-4cfc-9e37-7ea86034852e.root" # nevents = 624960
}


################################################################################
if SampleToRun == "Data2022" :
    # 2022 data sample from /MuonEG/Run2022D-PromptNanoAODv10_v1-v1/NANOAOD
    setConf("IsMC", False)
    setConf("LEPTON_SETUP", 2022)
    setConf("PD", "any")
    setConf("SAMPLENAME", "test")
    setConf("TRIGPASSTHROUGH", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/data/Run2022D/MuonEG/NANOAOD/PromptNanoAODv10_v2-v1/50000/68f42f42-3274-46ec-b23d-bfadc13012c2.root",
        ])


################################################################################
elif SampleToRun == "ggh125_UL" : ### 2018 UL test sample
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2018)
    setConf("DATA_TAG", "UL")
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/RunIISummer20UL18NanoAODv2/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/3B6A5CB5-2B7C-924D-85B4-FC3B0C1F4909.root",
        ])

################################################################################
elif SampleToRun == "MCsync_UL" :
    # Custom-reprocessed Rereco nanoAOD file with updated FSR and electron MVA,
    # no packing for genparticle p3; 26000 events
    # corresponding to:/store/mc/RunIISummer20UL17MiniAODv2/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/130000/3E4E8D55-3993-2B43-AF3B-7AB45BBE0BDA.root
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2017)
    setConf("NANOVERSION", 10) # variable defined as per nanoAOD v10 (notably electron_mvaHZZIso)
    setConf("DATA_TAG", "UL")
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR.root"])
#    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR_nopacking.root"]) # with no packing of muon eta, phi, mass


################################################################################
elif SampleToRun == "MCsync_Rereco" :
     # Custom-reprocessed Rereco nanoAOD file with updated FSR,
     # corresponding to:/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root
    setConf("APPLYMUCORR", True)
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("NANOVERSION", 9)
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_fixedFSR.root"])


################################################################################
elif SampleToRun == "MC2022" :
    # 2022 MC sample; for sync purposes, we use a JSON to select events matching the nanoAOD file:
    #/store/mc/Run3Summer22EEMiniAODv3/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/2540000/0b3804ca-2ae7-46df-bfab-49f92f76b047.root
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 52.234*0.0002745)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/Run3Summer22EENanoAODv11/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/126X_mcRun3_2022_realistic_postEE_v1-v1/2810000/e5e2fe04-7fb1-43ed-a81a-3acded81f0e7.root",
        ])
    json = {"1": [[1245, 1245],[1306, 1306],[1410, 1410],[1692, 1692],[1903, 1903],[1910, 1910],[1915, 1915],[1927, 1927],[1939, 1939],[1940, 1940],[1944, 1944],[1945, 1945],[1956, 1956],[1960, 1960],[1965, 1965],[1967, 1967],[1968, 1968],[1969, 1969],[2104, 2104]]}

################################################################################
    
elif SampleToRun == "MC_ZZto4L" :
    setConf("SAMPLENAME", "ZZto4L")
    setConf("XSEC", 0.0349) ## pp-->ZZ-->4l [pb] (https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV) # NOT THE VALUE REPORTED IN samplesNano_2022EE_MC.csv
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    ######## Test no MELA ###########
    setConf("runMELA", False)
    setConf("bestCandByMELA", False)
    #################################
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames", [f for f in mc_samples.values() if f is not None])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Local ZZAnalysis")
    parser.add_argument("--out_ext", help="string to add to end of output file", default="_Skim")
    parser.add_argument("--max_entries", help="maximum number of entries to process. Use -1 for all entries.", default=1000, type=int)
    parser.add_argument("--debug", help="run with DEBUG True (verbose output)", default=False)

    FLAGS = parser.parse_args()

    if FLAGS.debug: setConf("DEBUG", True)

    #####################################################################
    ### This import should be done AFTER all customizations (setConf calls)
    from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
    ######################################################################

    ### Tweak postprocessor parameters as necessary
    p.prefetch=True # Prefetch remote files
    p.longTermCache=True # keep prefetched files (useful for rerunning tests several times)
    if len(p.inputFiles) == 1 :
        p.haddFileName = None # Skip final hadd
    p.postfix = FLAGS.out_ext
    if FLAGS.max_entries != -1: p.maxEntries = FLAGS.max_entries

    if len(mc_samples.keys()) > 1:
        irr_bkg = "_".join(mc_samples.keys()).replace("ZZ_TTC_", "")
        p.outputDir  = "Processed/MC"
        p.haddFileName = "ZZ4lAnalysis_with{}_nev_{}_COSTHETA.root".format(irr_bkg, FLAGS.max_entries)

    ### Print out detailed candidate information for debug purposes
    #p.cut = None # Remove preselction
    #from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
    #insertAfter(p.modules,"lepFiller",dumpEvents(level=-1),getConf("NANOVERSION", 11)) 

    #p.branchsel=None #Read all branches
    #p.outputbranchsel=None #Output all branches

    #replace JSON
    #p.json = json

    ### Run the postprocessor
    p.run()
