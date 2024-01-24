#!/bin/env python3
### Example macro for filling standard histograms from H4l nanoAODs.
### Histograms are stored on a file and can then be plotted with

from __future__ import print_function
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons


pathMC2018 = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/MC2018/" # FIXME: Use 2018 MC for the time being
pathDATA_CD = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_CD/" # data for 2022 period
pathDATA_EFG = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231209_nano/Data2022_EFG/" # data for 2022EE period
pathMC2022 = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022/'
pathMC2022EE = '/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII/231214_nano/MC2022EE/'


ZmassValue = 91.1876

maxEntriesPerSample = 1e12 # Use only up to this number of events in each MC sample, for quick tests.



ROOT.TH1.SetDefaultSumw2()

####################################
def fillHistos(samplename, filename) :

    ### ---------------------
    ## ZZMass
    h_ZZMass2 = ROOT.TH1F("ZZMass_2GeV_"+samplename,
                          "ZZMass_2GeV_"+samplename,65,70.,200.)
    h_ZZMass2.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_ZZMass2_4mu = ROOT.TH1F("ZZMass_2GeV_4mu_"+samplename,
                              "ZZMass_2GeV_4mu_"+samplename,65,70.,200.)
    h_ZZMass2_4mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_ZZMass2_4e = ROOT.TH1F("ZZMass_2GeV_4e_"+samplename,
                             "ZZMass_2GeV_4e_"+samplename,65,70.,200.)
    h_ZZMass2_4e.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_ZZMass2_2e2mu = ROOT.TH1F("ZZMass_2GeV_2e2mu_"+samplename,
                                "ZZMass_2GeV_2e2mu_"+samplename,65,70.,200.)
    h_ZZMass2_2e2mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    ### ---------------------    
    h_ZZMass4 = ROOT.TH1F("ZZMass_4GeV_"+samplename,
                          "ZZMass_4GeV_"+samplename,233,70.,1002.)
    h_ZZMass4.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4.GetYaxis().SetTitle("Events / 4 GeV")
    #4mu
    h_ZZMass4_4mu = ROOT.TH1F("ZZMass_4GeV_4mu_"+samplename,
                              "ZZMass_4GeV_4mu_"+samplename,233,70.,1002.)
    h_ZZMass4_4mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4_4mu.GetYaxis().SetTitle("Events / 4 GeV")
    #4e
    h_ZZMass4_4e = ROOT.TH1F("ZZMass_4GeV_4e_"+samplename,
                             "ZZMass_4GeV_4e_"+samplename,233,70.,1002.)
    h_ZZMass4_4e.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4_4e.GetYaxis().SetTitle("Events / 4 GeV")
    #2e2mu
    h_ZZMass4_2e2mu = ROOT.TH1F("ZZMass_4GeV_2e2mu_"+samplename,
                                "ZZMass_4GeV_2e2mu_"+samplename,233,70.,1002.)
    h_ZZMass4_2e2mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4_2e2mu.GetYaxis().SetTitle("Events / 4 GeV")

    # h_ZZMass10 = ROOT.TH1F("ZZMass_10GeV_"+samplename,
    #                        "ZZMass_10GeV_"+samplename,93,70.,1000.)
    # h_ZZMass10.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    # h_ZZMass10.GetYaxis().SetTitle("Events / 10 GeV")

    ### ---------------------    
    ## Z1 and Z2 masses
    # Z1
    h_Z1Mass = ROOT.TH1F("Z1Mass_"+samplename,
                         "Z1Mass_"+samplename,40,40.,120.)
    h_Z1Mass.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_Z1Mass_4mu = ROOT.TH1F("Z1Mass_4mu_"+samplename,
                             "Z1Mass_4mu_"+samplename,40,40.,120.)
    h_Z1Mass_4mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_Z1Mass_4e = ROOT.TH1F("Z1Mass_4e_"+samplename,
                            "Z1Mass_4e_"+samplename,40,40.,120.)
    h_Z1Mass_4e.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_Z1Mass_2e2mu = ROOT.TH1F("Z1Mass_2e2mu_"+samplename,
                               "Z1Mass_2e2mu_"+samplename,40,40.,120.)
    h_Z1Mass_2e2mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    ### ---------------------    
    # Z2
    h_Z2Mass = ROOT.TH1F("Z2Mass_"+samplename,
                         "Z2Mass_"+samplename,54,12.,120.)
    h_Z2Mass.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_Z2Mass_4mu = ROOT.TH1F("Z2Mass_4mu_"+samplename,
                             "Z2Mass_4mu_"+samplename,54,12.,120.)
    h_Z2Mass_4mu.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_Z2Mass_4e = ROOT.TH1F("Z2Mass_4e_"+samplename,
                            "Z2Mass_4e_"+samplename,54,12.,120.)
    h_Z2Mass_4e.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_Z2Mass_2e2mu = ROOT.TH1F("Z2Mass_2e2mu_"+samplename,
                               "Z2Mass_2e2mu_"+samplename,54,12.,120.)
    h_Z2Mass_2e2mu.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    ### ---------------------    
    ## KD
    h_KD = ROOT.TH1F("KD_"+samplename,
                     "KD_"+samplename,10,0.,1.)
    h_KD.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD.GetYaxis().SetTitle("Events / 0.1")
    #4mu
    h_KD_4mu = ROOT.TH1F("KD_4mu_"+samplename,
                         "KD_4mu_"+samplename,10,0.,1.)
    h_KD_4mu.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD_4mu.GetYaxis().SetTitle("Events / 0.1")
    #4e
    h_KD_4e = ROOT.TH1F("KD_4e_"+samplename,
                        "KD_4e_"+samplename,10,0.,1.)
    h_KD_4e.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD_4e.GetYaxis().SetTitle("Events / 0.1")
    #2e2mu
    h_KD_2e2mu = ROOT.TH1F("KD_2e2mu_"+samplename,
                           "KD_2e2mu_"+samplename,10,0.,1.)
    h_KD_2e2mu.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD_2e2mu.GetYaxis().SetTitle("Events / 0.1")

    ### ---------------------    
    # 2D plots: Z1mass vs Z2mass
    h2_Z1Mass_Z2Mass = ROOT.TH2F("Z1MassVsZ2Mass_"+samplename,
                                 "Z1MassVsZ2Mass_"+samplename,
                                 40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")
    #4mu
    h2_Z1Mass_Z2Mass_4mu = ROOT.TH2F("Z1MassVsZ2Mass_4mu_"+samplename,
                                     "Z1MassVsZ2Mass_4mu_"+samplename,
                                     40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass_4mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass_4mu.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")
    #4e
    h2_Z1Mass_Z2Mass_4e = ROOT.TH2F("Z1MassVsZ2Mass_4e_"+samplename,
                                    "Z1MassVsZ2Mass_4e_"+samplename,
                                    40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass_4e.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass_4e.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")
    #2e2mu
    h2_Z1Mass_Z2Mass_2e2mu = ROOT.TH2F("Z1MassVsZ2Mass_2e2mu_"+samplename,
                                       "Z1MassVsZ2Mass_2e2mu_"+samplename,
                                       40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass_2e2mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass_2e2mu.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")

    ### ---------------------    
    # 2D plots: ZZMass vs KD
    h2_ZZMass_KD = ROOT.TH2F("ZZMassVsKD_"+samplename,
                             "ZZMassVsKD_"+samplename,
                             65,70.,200.,10,0.,1.)
    h2_ZZMass_KD.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")
    #4mu
    h2_ZZMass_KD_4mu = ROOT.TH2F("ZZMassVsKD_4mu_"+samplename,
                                 "ZZMassVsKD_4mu_"+samplename,
                                 65,70.,200.,10,0.,1.)
    h2_ZZMass_KD_4mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD_4mu.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")
    #4e
    h2_ZZMass_KD_4e = ROOT.TH2F("ZZMassVsKD_4e_"+samplename,
                                "ZZMassVsKD_4e_"+samplename,
                                65,70.,200.,10,0.,1.)
    h2_ZZMass_KD_4e.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD_4e.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")
    #2e2mu
    h2_ZZMass_KD_2e2mu = ROOT.TH2F("ZZMassVsKD_2e2mu_"+samplename,
                                   "ZZMassVsKD_2e2mu_"+samplename,
                                   65,70.,200.,10,0.,1.)
    h2_ZZMass_KD_2e2mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD_2e2mu.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")
    
    ### ---------------------    
    ### BLIND plots
    ## ZZMass
    h_ZZMass2_blind = ROOT.TH1F("ZZMass_2GeV_blind_"+samplename,
                                "ZZMass_2GeV_blind_"+samplename,65,70.,200.)
    h_ZZMass2_blind.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2_blind.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_ZZMass2_blind_4mu = ROOT.TH1F("ZZMass_2GeV_blind_4mu_"+samplename,
                                    "ZZMass_2GeV_blind_4mu_"+samplename,65,70.,200.)
    h_ZZMass2_blind_4mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2_blind_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_ZZMass2_blind_4e = ROOT.TH1F("ZZMass_2GeV_blind_4e_"+samplename,
                                   "ZZMass_2GeV_blind_4e_"+samplename,65,70.,200.)
    h_ZZMass2_blind_4e.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2_blind_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_ZZMass2_blind_2e2mu = ROOT.TH1F("ZZMass_2GeV_blind_2e2mu_"+samplename,
                                      "ZZMass_2GeV_blind_2e2mu_"+samplename,65,70.,200.)
    h_ZZMass2_blind_2e2mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2_blind_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    ### ---------------------   
    h_ZZMass4_blind = ROOT.TH1F("ZZMass_4GeV_blind_"+samplename,
                                "ZZMass_4GeV_blind_"+samplename,233,70.,1002.)
    h_ZZMass4_blind.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4_blind.GetYaxis().SetTitle("Events / 4 GeV")
    #4mu
    h_ZZMass4_blind_4mu = ROOT.TH1F("ZZMass_4GeV_blind_4mu_"+samplename,
                                    "ZZMass_4GeV_blind_4mu_"+samplename,233,70.,1002.)
    h_ZZMass4_blind_4mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4_blind_4mu.GetYaxis().SetTitle("Events / 4 GeV")
    #4e
    h_ZZMass4_blind_4e = ROOT.TH1F("ZZMass_4GeV_blind_4e_"+samplename,
                                   "ZZMass_4GeV_blind_4e_"+samplename,233,70.,1002.)
    h_ZZMass4_blind_4e.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4_blind_4e.GetYaxis().SetTitle("Events / 4 GeV")
    #2e2mu
    h_ZZMass4_blind_2e2mu = ROOT.TH1F("ZZMass_4GeV_blind_2e2mu_"+samplename,
                                      "ZZMass_4GeV_blind_2e2mu_"+samplename,233,70.,1002.)
    h_ZZMass4_blind_2e2mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4_blind_2e2mu.GetYaxis().SetTitle("Events / 4 GeV")

    ### ---------------------   
    ## Z1 and Z2 masses
    # Z1
    h_Z1Mass_blind = ROOT.TH1F("Z1Mass_blind_"+samplename,
                               "Z1Mass_blind_"+samplename,40,40.,120.)
    h_Z1Mass_blind.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_blind.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_Z1Mass_blind_4mu = ROOT.TH1F("Z1Mass_blind_4mu_"+samplename,
                                   "Z1Mass_blind_4mu_"+samplename,40,40.,120.)
    h_Z1Mass_blind_4mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_blind_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_Z1Mass_blind_4e = ROOT.TH1F("Z1Mass_blind_4e_"+samplename,
                                  "Z1Mass_blind_4e_"+samplename,40,40.,120.)
    h_Z1Mass_blind_4e.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_blind_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_Z1Mass_blind_2e2mu = ROOT.TH1F("Z1Mass_blind_2e2mu_"+samplename,
                                     "Z1Mass_blind_2e2mu_"+samplename,40,40.,120.)
    h_Z1Mass_blind_2e2mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h_Z1Mass_blind_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    ### ---------------------   
    # Z2
    h_Z2Mass_blind = ROOT.TH1F("Z2Mass_blind_"+samplename,
                               "Z2Mass_blind_"+samplename,54,12.,120.)
    h_Z2Mass_blind.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_blind.GetYaxis().SetTitle("Events / 2 GeV")
    #4mu
    h_Z2Mass_blind_4mu = ROOT.TH1F("Z2Mass_blind_4mu_"+samplename,
                                   "Z2Mass_blind_4mu_"+samplename,54,12.,120.)
    h_Z2Mass_blind_4mu.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_blind_4mu.GetYaxis().SetTitle("Events / 2 GeV")
    #4e
    h_Z2Mass_blind_4e = ROOT.TH1F("Z2Mass_blind_4e_"+samplename,
                                  "Z2Mass_blind_4e_"+samplename,54,12.,120.)
    h_Z2Mass_blind_4e.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_blind_4e.GetYaxis().SetTitle("Events / 2 GeV")
    #2e2mu
    h_Z2Mass_blind_2e2mu = ROOT.TH1F("Z2Mass_blind_2e2mu_"+samplename,
                                     "Z2Mass_blind_2e2mu_"+samplename,54,12.,120.)
    h_Z2Mass_blind_2e2mu.GetXaxis().SetTitle("m_{#it{Z2}} (GeV)")
    h_Z2Mass_blind_2e2mu.GetYaxis().SetTitle("Events / 2 GeV")

    ### ---------------------   
    ## KD
    h_KD_blind = ROOT.TH1F("KD_blind_"+samplename,
                           "KD_blind_"+samplename,10,0.,1.)
    h_KD_blind.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD_blind.GetYaxis().SetTitle("Events / 0.1")
    #4mu
    h_KD_blind_4mu = ROOT.TH1F("KD_blind_4mu_"+samplename,
                               "KD_blind_4mu_"+samplename,10,0.,1.)
    h_KD_blind_4mu.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD_blind_4mu.GetYaxis().SetTitle("Events / 0.1")
    #4e
    h_KD_blind_4e = ROOT.TH1F("KD_blind_4e_"+samplename,
                              "KD_blind_4e_"+samplename,10,0.,1.)
    h_KD_blind_4e.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD_blind_4e.GetYaxis().SetTitle("Events / 0.1")
    #2e2mu
    h_KD_blind_2e2mu = ROOT.TH1F("KD_blind_2e2mu_"+samplename,
                                 "KD_blind_2e2mu_"+samplename,10,0.,1.)
    h_KD_blind_2e2mu.GetXaxis().SetTitle("#it{D}_{bkg}^{kin}")
    h_KD_blind_2e2mu.GetYaxis().SetTitle("Events / 0.1")

    ### ---------------------   
    # 2D plots: Z1mass vs Z2mass
    h2_Z1Mass_Z2Mass_blind = ROOT.TH2F("Z1MassVsZ2Mass_blind_"+samplename,
                                       "Z1MassVsZ2Mass_blind_"+samplename,
                                       40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass_blind.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass_blind.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")
    #4mu
    h2_Z1Mass_Z2Mass_blind_4mu = ROOT.TH2F("Z1MassVsZ2Mass_blind_4mu_"+samplename,
                                           "Z1MassVsZ2Mass_blind_4mu_"+samplename,
                                           40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass_blind_4mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass_blind_4mu.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")
    #4e
    h2_Z1Mass_Z2Mass_blind_4e = ROOT.TH2F("Z1MassVsZ2Mass_blind_4e_"+samplename,
                                          "Z1MassVsZ2Mass_blind_4e_"+samplename,
                                          40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass_blind_4e.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass_blind_4e.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")
    #2e2mu
    h2_Z1Mass_Z2Mass_blind_2e2mu = ROOT.TH2F("Z1MassVsZ2Mass_blind_2e2mu_"+samplename,
                                             "Z1MassVsZ2Mass_blind_2e2mu_"+samplename,
                                             40,40.,120.,54,12.,120.)
    h2_Z1Mass_Z2Mass_blind_2e2mu.GetXaxis().SetTitle("m_{#it{Z1}} (GeV)")
    h2_Z1Mass_Z2Mass_blind_2e2mu.GetYaxis().SetTitle("m_{#it{Z2}} (GeV)")

    ### ---------------------   
    # 2D plots: ZZMass vs KD
    h2_ZZMass_KD_blind = ROOT.TH2F("ZZMassVsKD_blind_"+samplename,
                                   "ZZMassVsKD_blind_"+samplename,
                                   65,70.,200.,10,0.,1.)
    h2_ZZMass_KD_blind.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD_blind.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")
    #4mu
    h2_ZZMass_KD_blind_4mu = ROOT.TH2F("ZZMassVsKD_blind_4mu_"+samplename,
                                       "ZZMassVsKD_blind_4mu_"+samplename,
                                       65,70.,200.,10,0.,1.)
    h2_ZZMass_KD_blind_4mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD_blind_4mu.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")
    #4e
    h2_ZZMass_KD_blind_4e = ROOT.TH2F("ZZMassVsKD_blind_4e_"+samplename,
                                      "ZZMassVsKD_blind_4e_"+samplename,
                                      65,70.,200.,10,0.,1.)
    h2_ZZMass_KD_blind_4e.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD_blind_4e.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")
    #2e2mu
    h2_ZZMass_KD_blind_2e2mu = ROOT.TH2F("ZZMassVsKD_blind_2e2mu_"+samplename,
                                         "ZZMassVsKD_blind_2e2mu_"+samplename,
                                         65,70.,200.,10,0.,1.)
    h2_ZZMass_KD_blind_2e2mu.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h2_ZZMass_KD_blind_2e2mu.GetYaxis().SetTitle("#it{D}_{bkg}^{kin}")



    f = ROOT.TFile.Open(filename)

    event = f.Events
    event.SetBranchStatus("*", 0)
    event.SetBranchStatus("run", 1)
    event.SetBranchStatus("luminosityBlock", 1)
    event.SetBranchStatus("*Muon*", 1)
    event.SetBranchStatus("*Electron*", 1)
    event.SetBranchStatus("*ZZCand*", 1)
    event.SetBranchStatus("bestCandIdx", 1)
    event.SetBranchStatus("HLT_passZZ4l", 1)
    nEntries = event.GetEntries() 

    isMC = False
    if(samplename == "Data"):
        print("Data: sel=", nEntries)
    else:
        isMC = True
        event.SetBranchStatus("overallEventWeight",1)

        # Get sum of weights
        runs = f.Runs
        nRuns = runs.GetEntries()
        iRun = 0
        genEventCount = 0
        genEventSumw = 0.
        while iRun < nRuns and runs.GetEntry(iRun) :
            genEventCount += runs.genEventCount
            genEventSumw += runs.genEventSumw
            iRun +=1
        print (samplename, ": gen=", genEventCount, "sel=",nEntries, "sumw=", genEventSumw)
        # run only up to a certain number of evt
        # if nEntries>maxEntriesPerSample :
        #    genEventSumw = genEventSumw*maxEntriesPerSample/nEntries
        #    nEntries=maxEntriesPerSample
        #    print("   scaling to:", nEntries, "sumw=", genEventSumw )


    iEntry=0
    printEntries=max(5000,nEntries/10)
    while iEntry<nEntries and event.GetEntry(iEntry):
        iEntry+=1
        if iEntry%printEntries == 0 : print("Processing", iEntry)

        bestCandIdx = event.bestCandIdx

        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)
        if(bestCandIdx != -1 and event.HLT_passZZ4l): 
            weight = 1.
            ZZs = Collection(event, 'ZZCand')
            theZZ = ZZs[bestCandIdx]        
            if isMC : 
                weight = (event.overallEventWeight*theZZ.dataMCWeight/genEventSumw)
            ## ZZmass
            m4l=theZZ.mass
            h_ZZMass2.Fill(m4l,weight)
            h_ZZMass4.Fill(m4l,weight)
            # h_ZZMass10.Fill(m4l,weight)
            # ## Z1Mass
            mZ1=theZZ.Z1mass
            h_Z1Mass.Fill(mZ1,weight)
            # ## Z2Mass
            mZ2=theZZ.Z2mass
            h_Z2Mass.Fill(mZ2,weight)
            # ## KD
            KD=theZZ.KD
            h_KD.Fill(KD,weight)

            # # 2D histo
            h2_Z1Mass_Z2Mass.Fill(mZ1,mZ2,weight)
            h2_ZZMass_KD.Fill(m4l,KD,weight)

            # plots per final state
            Z1flav = theZZ.Z1flav
            Z2flav = theZZ.Z2flav
            if(Z1flav==-169 and Z2flav==-169):
                h_ZZMass2_4mu.Fill(m4l,weight)
                h_ZZMass4_4mu.Fill(m4l,weight)
                h_Z1Mass_4mu.Fill(mZ1,weight)
                h_Z2Mass_4mu.Fill(mZ2,weight)
                h_KD_4mu.Fill(KD,weight)
                h2_Z1Mass_Z2Mass_4mu.Fill(mZ1,mZ2,weight)
                h2_ZZMass_KD_4mu.Fill(m4l,KD,weight)
            elif(Z1flav==-121 and Z2flav==-121):
                h_ZZMass2_4e.Fill(m4l,weight)
                h_ZZMass4_4e.Fill(m4l,weight)
                h_Z1Mass_4e.Fill(mZ1,weight)
                h_Z2Mass_4e.Fill(mZ2,weight)
                h_KD_4e.Fill(KD,weight)
                h2_Z1Mass_Z2Mass_4e.Fill(mZ1,mZ2,weight)
                h2_ZZMass_KD_4e.Fill(m4l,KD,weight) 
            elif((Z1flav==-169 and Z2flav==-121) or 
                 (Z1flav==-121 and Z2flav==-169)):
                h_ZZMass2_2e2mu.Fill(m4l,weight)
                h_ZZMass4_2e2mu.Fill(m4l,weight)
                h_Z1Mass_2e2mu.Fill(mZ1,weight)
                h_Z2Mass_2e2mu.Fill(mZ2,weight)
                h_KD_2e2mu.Fill(KD,weight)
                h2_Z1Mass_Z2Mass_2e2mu.Fill(mZ1,mZ2,weight)
                h2_ZZMass_KD_2e2mu.Fill(m4l,KD,weight) 
            else:
                print('error in Zflav ',Z1flav,Z2flav)


            # ### BLIND plots
            # if m4l < 105. or m4l > 140.:
            #     h_Z1Mass_blind.Fill(mZ1,weight)
            #     h_Z2Mass_blind.Fill(mZ2,weight)
            #     h_KD_blind.Fill(KD,weight)
            #     h2_Z1Mass_Z2Mass_blind.Fill(mZ1,mZ2,weight)
            #     h2_ZZMass_KD_blind.Fill(m4l,KD,weight)

            #     # plots per final state
            #     if(Z1flav==-169 and Z2flav==-169):
            #         h_Z1Mass_blind_4mu.Fill(mZ1,weight)
            #         h_Z2Mass_blind_4mu.Fill(mZ2,weight)
            #         h_KD_blind_4mu.Fill(KD,weight)
            #         h2_Z1Mass_Z2Mass_blind_4mu.Fill(mZ1,mZ2,weight)
            #         h2_ZZMass_KD_blind_4mu.Fill(m4l,KD,weight)
            #     elif(Z1flav==-121 and Z2flav==-121):
            #         h_Z1Mass_blind_4e.Fill(mZ1,weight)
            #         h_Z2Mass_blind_4e.Fill(mZ2,weight)
            #         h_KD_blind_4e.Fill(KD,weight)
            #         h2_Z1Mass_Z2Mass_blind_4e.Fill(mZ1,mZ2,weight)
            #         h2_ZZMass_KD_blind_4e.Fill(m4l,KD,weight) 
            #     elif((Z1flav==-169 and Z2flav==-121) or 
            #          (Z1flav==-121 and Z2flav==-169)):
            #         h_Z1Mass_blind_2e2mu.Fill(mZ1,weight)
            #         h_Z2Mass_blind_2e2mu.Fill(mZ2,weight)
            #         h_KD_blind_2e2mu.Fill(KD,weight)
            #         h2_Z1Mass_Z2Mass_blind_2e2mu.Fill(mZ1,mZ2,weight)
            #         h2_ZZMass_KD_blind_2e2mu.Fill(m4l,KD,weight) 
            #     else:
            #         print('error in Zflav ',Z1flav,Z2flav)

            # Example on how to get the four leptons of the candidates, ordered as
            # [Z1l1, Z2l2, Z2l1, Z2l2]
            #leps = getLeptons(theZZ, event)
            #print(leps[3].pt)
        
    f.Close()
    
    histos = [h_ZZMass2, h_ZZMass2_4mu, h_ZZMass2_4e, h_ZZMass2_2e2mu,
              h_ZZMass4, h_ZZMass4_4mu, h_ZZMass4_4e, h_ZZMass4_2e2mu,
              h_Z1Mass, h_Z1Mass_4mu, h_Z1Mass_4e, h_Z1Mass_2e2mu,
              h_Z2Mass, h_Z2Mass_4mu, h_Z2Mass_4e, h_Z2Mass_2e2mu,  
              h_KD, h_KD_4mu, h_KD_4e, h_KD_2e2mu, 
              h2_Z1Mass_Z2Mass, h2_Z1Mass_Z2Mass_4mu, h2_Z1Mass_Z2Mass_4e, h2_Z1Mass_Z2Mass_2e2mu,
              h2_ZZMass_KD, h2_ZZMass_KD_4mu, h2_ZZMass_KD_4e, h2_ZZMass_KD_2e2mu,
              # h_Z1Mass_blind, h_Z1Mass_blind_4mu, h_Z1Mass_blind_4e, h_Z1Mass_blind_2e2mu,
              # h_Z2Mass_blind, h_Z2Mass_blind_4mu, h_Z2Mass_blind_4e, h_Z2Mass_blind_2e2mu,
              # h_KD_blind, h_KD_blind_4mu, h_KD_blind_4e, h_KD_blind_2e2mu,
              # h2_Z1Mass_Z2Mass_blind, h2_Z1Mass_Z2Mass_blind_4mu, h2_Z1Mass_Z2Mass_blind_4e, h2_Z1Mass_Z2Mass_blind_2e2mu,
              # h2_ZZMass_KD_blind, h2_ZZMass_KD_blind_4mu, h2_ZZMass_KD_blind_4e, h2_ZZMass_KD_blind_2e2mu
              ]

    return histos


def runMC(outFile): 

    if '2018' in outFile:
        samples = [
            # ggZZ from 2018
            dict(name = "ggTo4mu",filename = pathMC2018+
                        "ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo4e",filename = pathMC2018+
                        "ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo4tau",filename = pathMC2018+
                        "ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo2e2mu",filename = pathMC2018+
                        "ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root"),       
            dict(name = "ggTo2e2tau",filename = pathMC2018+
                        "ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root"),
            dict(name = "ggTo2mu2tau",filename = pathMC2018+
                        "ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root"),
        ]
    elif '2022EE' in outFile:
        samples = [
            dict(name = "ggH125",filename = pathMC2022EE+
                        "ggH125/ZZ4lAnalysis.root"),
            dict(name = "VBF125",filename = pathMC2022EE+
                        "VBFH125/ZZ4lAnalysis.root"),
            dict(name = "WplusH125",filename = pathMC2022EE+
                        "WplusH125/ZZ4lAnalysis.root"),
            dict(name = "WHminus125",filename = pathMC2022EE+
                        "WHminus125/ZZ4lAnalysis.root"),
            dict(name = "ZH125",filename = pathMC2022EE+
                        "ZH125/ZZ4lAnalysis.root"),
            dict(name = "ttH125",filename = pathMC2022EE+
                        "ttH125/ZZ4lAnalysis.root"),
            dict(name = "bbH125",filename = pathMC2022EE+
                        "bbH125/ZZ4lAnalysis.root"),
            
            dict(name = "ZZTo4l",filename = pathMC2022EE+
                        "ZZTo4l/ZZ4lAnalysis.root"),

            dict(name = "WWZ",filename = pathMC2022EE+
                        "WWZ/ZZ4lAnalysis.root"),
            dict(name = "WZZ",filename = pathMC2022EE+
                        "WZZ/ZZ4lAnalysis.root"),
            dict(name = "ZZZ",filename = pathMC2022EE+
                        "ZZZ/ZZ4lAnalysis.root"),
            dict(name = "TTWW",filename = pathMC2022EE+
                        "TTWW/ZZ4lAnalysis.root"),
            dict(name = "TTZZ",filename = pathMC2022EE+
                        "TTZZ/ZZ4lAnalysis.root"),
        ]
    elif '2022' in outFile:
        samples = [
            dict(name = "ggH125",filename = pathMC2022+
                        "ggH125/ZZ4lAnalysis.root"),
            dict(name = "VBF125",filename = pathMC2022+
                        "VBFH125/ZZ4lAnalysis.root"),
            dict(name = "WplusH125",filename = pathMC2022+
                        "WplusH125/ZZ4lAnalysis.root"),
            dict(name = "WHminus125",filename = pathMC2022+
                        "WHminus125/ZZ4lAnalysis.root"),
            dict(name = "ZH125",filename = pathMC2022+
                        "ZH125/ZZ4lAnalysis.root"),
            dict(name = "ttH125",filename = pathMC2022+
                        "ttH125/ZZ4lAnalysis.root"),
            dict(name = "bbH125",filename = pathMC2022+
                        "bbH125/ZZ4lAnalysis.root"),

            dict(name = "ZZTo4l",filename = pathMC2022+
                        "ZZTo4l/ZZ4lAnalysis.root"),

            dict(name = "WWZ",filename = pathMC2022+
                        "WWZ/ZZ4lAnalysis.root"),
            dict(name = "WZZ",filename = pathMC2022+
                        "WZZ/ZZ4lAnalysis.root"),
            dict(name = "ZZZ",filename = pathMC2022+
                        "ZZZ/ZZ4lAnalysis.root"),
            dict(name = "TTWW",filename = pathMC2022+
                        "TTWW/ZZ4lAnalysis.root"),
            dict(name = "TTZZ",filename = pathMC2022+
                        "TTZZ/ZZ4lAnalysis.root"),
        ]


    of = ROOT.TFile.Open(outFile,"recreate") 
    
    for s in samples:
         histos = fillHistos(s["name"], s["filename"])
         for h in histos:
             of.WriteObject(h,h.GetName())

    of.Close()

def runData(outFile):

    of = ROOT.TFile.Open(outFile,"recreate") 

    if 'CD' in outFile:
        hs_data = fillHistos("Data", pathDATA_CD + "ZZ4lAnalysis.root")
    elif 'EFG' in outFile:
        hs_data = fillHistos("Data", pathDATA_EFG + "ZZ4lAnalysis.root")

    for h in hs_data:
        h.SetBinErrorOption(ROOT.TH1.kPoisson)
        of.WriteObject(h,h.GetName())

    of.Close()

if __name__ == "__main__" :

    print('Running 2018')
    runMC('origCJLST_histos/H4l_MC2018.root')
    print('Running 2022')
    runMC('origCJLST_histos/H4l_MC2022.root')
    print('Running 2022EE')
    runMC('origCJLST_histos/H4l_MC2022EE.root')
    print('Running C-D data')
    runData('origCJLST_histos/H4l_Data_CD.root')
    print('Running E-F-G data')
    runData('origCJLST_histos/H4l_Data_EFG.root')