### Draw and decorate plots produced with H4l_fill.py.
# This is just a quick example, for illustration purposes only!!!
# It lacks all the plots and features of the full miniAOD plotter.
#
# run 
# python3 H4l_draw_mZZ_full2022.py

from __future__ import print_function
import glob
import optparse
import os
import os.path as osp
import sys
from datetime import date
import math
import ctypes
import ROOT
import CMSGraphics, CMS_lumi
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True

inFilenameMC2018     = "origCJLST_histos/H4l_MC2018.root"
inFilenameMC2022     = 'origCJLST_histos/H4l_MC2022.root'
inFilenameMC2022EE   = 'origCJLST_histos/H4l_MC2022EE.root'
inFilenameData2022   = "origCJLST_histos/H4l_Data_CD.root"
inFilenameData2022EE = "origCJLST_histos/H4l_Data_EFG.root"
outFilename = "Plots_inclusive_TESTING.root"

## output directory
today = date.today()
print('Creating output dir...')
out_dir = str(today)+'_plots_mZZ_inclusive'
os.makedirs(out_dir, exist_ok=True) #check if output dir exist

### 2018 plots
#Lum = 59.74 # 1/fb
#pathMC = "/eos/user/n/namapane/H4lnano/220420/"
#pathDATA = "/eos/user/n/namapane/H4lnano/220420/Data2018/"

# lumi
lumi_2022 = 35.08424 # 1/fb 2022 C-G  (= 35.181930231/fb of full 355100_362760 Golden json - 0.097685694 of eraB that we don't use
lumi_CD   = 8.077 # 1/fb
lumi_EFG  = 27.007 # 1/fb


# plots options
blindPlots = True
blindHLow = 105.
blindHHi  = 140.
blindHM   = 500.
epsilon=0.1
addEmptyBins = True


# Set style matching the one used for HZZ plots
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetPadTopMargin(0.05)  
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.16) 
ROOT.gStyle.SetPadRightMargin(0.03)
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetStripDecimals(True)
ROOT.gStyle.SetTickLength(0.03, "XYZ")
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetTitleOffset(1.00, "X")
ROOT.gStyle.SetTitleOffset(1.25, "Y")
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")

canvasSizeX=910
canvasSizeY=700



#ZX estaimation parameters - taken from 2018 data - approx. normalization, just for visualization purposes
# CURRENTLY JUST SCALING 2018 DATA TO CURRENT LUMI
def getZX(h_model) :
    n_entries = 10000
    bin_down  = 70.
    bin_up    = 3000.
    lumi2018  = 59.7*1000. # to normalize
   
    f_4e_comb    = ROOT.TF1("f_4e_comb", "TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_4mu_comb   = ROOT.TF1("f_4mu_comb","TMath::Landau(x, [0], [1])", bin_down, bin_up)
    f_2e2mu_comb = ROOT.TF1("f_2e2mu_comb","[0]*TMath::Landau(x, [1], [2]) + [3]*TMath::Landau(x, [4], [5])", bin_down, bin_up)

    f_4e_comb.SetParameters(141.9, 21.3)
    f_4mu_comb.SetParameters(130.4, 15.6)
    f_2e2mu_comb.SetParameters(0.45,131.1,18.1, 0.55,133.8,18.9)

    yield_Comb_4e_2018    = 19.42/lumi2018
    yield_Comb_4mu_2018   = 50.72/lumi2018
    yield_Comb_2e2mu_2018 = 63.87/lumi2018

    h_4e=h_model.Clone("ZX_4e")
    h_4e.Reset()
#    h_4e.SetFillColor(ROOT.TColor.GetColor("#0331B9"))
    h_4mu=h_4e.Clone("ZX_4mu")
    h_2e2mu=h_4e.Clone("ZX_2e2mu")
    
    h_4e.FillRandom("f_4e_comb"   , n_entries)
    h_4mu.FillRandom("f_4mu_comb"  , n_entries)
    h_2e2mu.FillRandom("f_2e2mu_comb", n_entries)

    h_4e.Scale(yield_Comb_4e_2018/h_4e.Integral())
    h_4mu.Scale(yield_Comb_4mu_2018/h_4mu.Integral())
    h_2e2mu.Scale(yield_Comb_2e2mu_2018/h_2e2mu.Integral())
    

    h_total=h_4e.Clone("ZX_tot")
    h_total.Add(h_4mu)
    h_total.Add(h_2e2mu)
    print("Z+X integral", h_total.Integral())
    return h_total


#####################
def printCanvases(type="png", path=".") :
    canvases = ROOT.gROOT.GetListOfCanvases()
    for c in canvases :
        c.Print(path+"/"+c.GetTitle()+"."+type)

def printCanvas(c, type="png", name=None, path="." ) :
    if name == None : name = c.GetTitle()
    name=name.replace(">","")
    name=name.replace("<","")
    name=name.replace(" ","_")
    c.Print(path+"/"+name+"."+type)


######################
def Stack_full2022(f2018, f2022, f2022EE, version = "_4GeV_"):
    name = "ZZMass" + version

    #------------EW------------------#
    # 2022 (C-D)
    WWZ  = f2022.Get(name+"WWZ")
    WZZ  = f2022.Get(name+"WZZ")
    ZZZ  = f2022.Get(name+"ZZZ")
    TTWW = f2022.Get(name+"TTWW")
    TTZZ = f2022.Get(name+"TTZZ")
    EWSamples = [WZZ, ZZZ, TTWW, TTZZ]
    EW = WWZ.Clone("h_EW")
    for i in EWSamples:
        EW.Add(i,1.)
    EW.Scale(lumi_CD*1000.)

    # 2022EE (E-G)
    WWZee  = f2022EE.Get(name+"WWZ")
    WZZee  = f2022EE.Get(name+"WZZ")
    ZZZee  = f2022EE.Get(name+"ZZZ")
    TTWWee = f2022EE.Get(name+"TTWW")
    TTZZee = f2022EE.Get(name+"TTZZ")
    EWSamplesee = [WZZee, ZZZee, TTWWee, TTZZee]
    EWee = WWZee.Clone("h_EWee")
    for i in EWSamplesee:
        EWee.Add(i,1.)
    EWee.Scale(lumi_EFG*1000.)

    # full 2022 histo
    EW.Add(EWee,1.) #full 2022
            
    EW.SetLineColor(ROOT.TColor.GetColor("#000099"))
    EW.SetFillColor(ROOT.TColor.GetColor("#0331B9"))

    
    #-----------qqZZ---------------#
    # 2022 (C-D)
    ZZTo4laa = f2022.Get(name+"ZZTo4l")
    ZZTo4laa.Scale(lumi_CD*1000.) 
    ZZTo4l = ZZTo4laa.Clone("h_ZZTo4l")
    # 2022EE (E-G)
    ZZTo4lee = f2022EE.Get(name+"ZZTo4l")
    ZZTo4lee.Scale(lumi_EFG*1000.) 
    # full 2022 histo
    ZZTo4l.Add(ZZTo4lee,1.) #full 2022
       
    ZZTo4l.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ZZTo4l.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    
    #-----------signal------------#
    # 2022 (C-D)
    VBF125     = f2022.Get(name+"VBF125")
    ggH125     = f2022.Get(name+"ggH125")
    WplusH125  = f2022.Get(name+"WplusH125")
    WminusH125 = f2022.Get(name+"WHminus125")
    ZH125      = f2022.Get(name+"ZH125")
    ttH125     = f2022.Get(name+"ttH125")
    bbH125     = f2022.Get(name+"bbH125")

    signalSamples = [ggH125, WplusH125, WminusH125, ZH125, ttH125, bbH125]
    signal = VBF125.Clone("h_signal")
    for i in signalSamples:
        signal.Add(i, 1.)
    signal.Scale(lumi_CD*1000.) 

    # 2022EE (E-G)
    VBF125ee     = f2022EE.Get(name+"VBF125")
    ggH125ee     = f2022EE.Get(name+"ggH125")
    WplusH125ee  = f2022EE.Get(name+"WplusH125")
    WminusH125ee = f2022EE.Get(name+"WHminus125")
    ZH125ee      = f2022EE.Get(name+"ZH125")
    ttH125ee     = f2022EE.Get(name+"ttH125")
    bbH125ee     = f2022EE.Get(name+"ttH125")
    
    signalSamplesee = [ggH125ee, WplusH125ee, WminusH125ee, ZH125ee, ttH125ee, bbH125ee]
    signalee = VBF125ee.Clone("h_signalee")
    for i in signalSamplesee:
        signalee.Add(i, 1.)
    signalee.Scale(lumi_EFG*1000.)

    # full 2022 histo
    signal.Add(signalee,1.) # full 2022
      
    signal.SetLineColor(ROOT.TColor.GetColor("#cc0000"))
    signal.SetFillColor(ROOT.TColor.GetColor("#ff9b9b"))

    
    #------------ggTo-----------------#
    # from 2018 for now
    ggTo4mu     = f2018.Get(name+"ggTo4mu") 
    ggTo4e      = f2018.Get(name+"ggTo4e")
    ggTo4tau    = f2018.Get(name+"ggTo4tau")
    ggTo2e2mu   = f2018.Get(name+"ggTo2e2mu")
    ggTo2e2tau  = f2018.Get(name+"ggTo2e2tau")
    ggTo2mu2tau = f2018.Get(name+"ggTo2mu2tau")

    ggZZSamples = [ ggTo4e, ggTo4tau, ggTo2e2mu, ggTo2e2tau, ggTo2mu2tau]
    ggToZZ = ggTo4mu.Clone("h_ggTo")
    for i in ggZZSamples:
        ggToZZ.Add(i,1.)
    ggToZZ.Scale(lumi_2022*1000.) #full 2022   
    ggToZZ.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ggToZZ.SetFillColor(ROOT.TColor.GetColor("#4b78ff"))  
    
    ### ZX
    # from 2018 for now
    hzx=getZX(signal)
    hzx.Scale(lumi_2022*1000.) #full 2022
    hzx.SetLineColor(ROOT.TColor.GetColor("#003300"))
    hzx.SetFillColor(ROOT.TColor.GetColor("#669966"))
    
    
    #------------------Stack----------#
    if version=="_4GeV_" :
        hs = ROOT.THStack("Stack_4GeV", "; m_{#it{4l}} (GeV) ; Events / 4 GeV" )
    elif version=="_10GeV_" :
        hs = ROOT.THStack("Stack_10GeV", "; m_{#it{4l}} (GeV) ; Events / 10 GeV" )
    else:
        hs = ROOT.THStack("Stack_2GeV", "; m_{#it{4l}} (GeV) ; Events / 2 GeV" )

    hs.Add(hzx,"HISTO")
    hs.Add(EW,"HISTO")
    hs.Add(ggToZZ,"HISTO")
    hs.Add(ZZTo4l,"HISTO")
    hs.Add(signal,"HISTO")
    
    return hs, [hzx, EW, ggToZZ, ZZTo4l, signal]


### Get a TGraph for data, blinded if required
def dataGraph (f1, f2, version = "_4GeV_", blind = True):
    name = "ZZMass"+ version
    hd1 = f1.Get(name+"Data")
    hd2 = f2.Get(name+"Data")
    hd = hd1.Clone('h_data') # full 2022
    hd.Add(hd2,1.)
    
    nbinsIn = hd.GetNbinsX()
    nbins = 0
    
    x = np.array([0.]*nbinsIn, dtype='double')
    y = np.array([0.]*nbinsIn, dtype='double')
    errX = np.array([0.]*nbinsIn, dtype='double')  
    UpErr = np.array([0.]*nbinsIn, dtype='double')
    LowErr = np.array([0.]*nbinsIn, dtype='double')


    for i in range (1, nbinsIn):
        if blind and ((hd.GetBinCenter(i)>=blindHLow and hd.GetBinCenter(i)<=blindHHi) or hd.GetBinCenter(i)>=blindHM) : continue
        if not addEmptyBins and hd.GetBinContent(i) == 0 : continue 
        x[nbins]      = hd.GetBinCenter(i)
        y[nbins]      = hd.GetBinContent(i)
        UpErr[nbins]  = hd.GetBinErrorUp(i)
        LowErr[nbins] = hd.GetBinErrorLow(i)
        nbins += 1
    
    breakpoint()

    Data = ROOT.TGraphAsymmErrors(nbins,x,y,errX,errX,LowErr,UpErr)
    Data.SetMarkerStyle(20)
    Data.SetLineColor(ROOT.kBlack)
    Data.SetMarkerSize(0.9)
    return Data



## --------------------------------
if __name__ == "__main__" :
    

    fMC2018     = ROOT.TFile.Open(inFilenameMC2018,"READ")
    fMC2022     = ROOT.TFile.Open(inFilenameMC2022,"READ")
    fMC2022EE   = ROOT.TFile.Open(inFilenameMC2022EE,"READ")
    fData2022   = ROOT.TFile.Open(inFilenameData2022,"READ")
    fData2022EE = ROOT.TFile.Open(inFilenameData2022EE,"READ")
    of = ROOT.TFile.Open(outFilename,"recreate")


    # Labels for log plots
    xlabelsv = [80, 100, 200, 300, 400, 500]
    label_margin = -0.1
    xlabels=[None]*len(xlabelsv)
    for i, label in enumerate(xlabelsv): 
        xlabels[i] = ROOT.TLatex(label, label_margin , str(label));
        xlabels[i].SetTextAlign(23)
        xlabels[i].SetTextFont(42)
        xlabels[i].SetTextSize(0.04)


    ## --- full 2022 m4l plot
    HStack, h_list = Stack_full2022(fMC2018, fMC2022, fMC2022EE)
    HData = dataGraph(fData2022, fData2022EE, blind=blindPlots)
    HStack_hm = HStack.Clone()
    HData_hm = HData.Clone()

    Canvas = ROOT.TCanvas("M4l_full2022","M4l_full2022",canvasSizeX,canvasSizeY)
    Canvas.SetTicks()
    Canvas.SetLogx()
    #ymaxd=HData.GetMaximum()
    xmin=ctypes.c_double(0.)
    ymin=ctypes.c_double(0.)
    xmax=ctypes.c_double(0.)
    ymax=ctypes.c_double(0.)
    HData.ComputeRange(xmin,ymin,xmax,ymax)
    yhmax=math.ceil(max(HStack.GetMaximum(), ymax.value))
    HStack.SetMaximum(yhmax)
    HStack.Draw("histo")
    HStack.GetXaxis().SetRangeUser(70., 300.)
    if blindPlots:
         ROOT.gPad.GetRangeAxis(xmin,ymin,xmax,ymax)
         bblind = ROOT.TBox(blindHLow, 0, blindHHi, ymax.value-epsilon)
         bblind.SetFillColor(ROOT.kGray)
         bblind.SetFillStyle(3002)
         bblind.Draw()
    HData.Draw("samePE1")
    # Hide labels and rewrite them
    HStack.GetXaxis().SetLabelSize(0)
    for label in xlabels :
        label.Draw()
    ROOT.gPad.RedrawAxis()

    legend = ROOT.TLegend(0.72,0.70,0.94,0.92)
    legend.AddEntry(h_list[4],"H(125)","f")
    legend.AddEntry(h_list[3],"q#bar{q}#rightarrow ZZ,Z#gamma*","f")
    legend.AddEntry(h_list[2],"gg#rightarrow ZZ,Z#gamma*","f")
    legend.AddEntry(h_list[1],"EW","f")
    legend.AddEntry(h_list[0],"Z+X","f")
    legend.AddEntry(HData,"Data", "p")
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetTextFont(43)
    legend.SetTextSize(20)
    legend.Draw()

    #draw CMS and lumi text
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText      = "Preliminary"
    CMS_lumi.lumi_sqrtS     = "35.1 fb-1 (13.6 TeV)"
    CMS_lumi.cmsTextSize    = 1 #0.6
    CMS_lumi.lumiTextSize   = 0.7 #0.46
    CMS_lumi.extraOverCmsTextSize = 0.75
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(Canvas, 0, 0)
    
    Canvas.Update() #very important!!!
    #Canvas.Write()

    
    ### Zoomed m4l
    HStack_z, h_list = Stack_full2022(fMC2018, fMC2022, fMC2022EE, "_2GeV_")
    HData_z = dataGraph(fData2022, fData2022EE, "_2GeV_", blind=blindPlots)
    Canvas_z = ROOT.TCanvas("M4l_full2022_z","M4l_full2022_z",canvasSizeX,canvasSizeY)
    Canvas_z.SetTicks()
    HData_z.ComputeRange(xmin,ymin,xmax,ymax)
    yhmax=math.ceil(max(HStack_z.GetMaximum(), ymax.value))
    HStack_z.SetMaximum(yhmax)
    HStack_z.Draw("histo")
    HStack_z.GetXaxis().SetRangeUser(70., 170.)
    if blindPlots:
         ROOT.gPad.GetRangeAxis(xmin,ymin,xmax,ymax)
         bblind_z = ROOT.TBox(blindHLow, 0, blindHHi, ymax.value-epsilon)
         bblind_z.SetFillColor(ROOT.kGray)
         bblind_z.SetFillStyle(3002)
         bblind_z.Draw()
    HData_z.Draw("samePE1")
    ROOT.gPad.RedrawAxis()
    
    legend_z = ROOT.TLegend(0.72,0.70,0.94,0.92)
    legend_z.AddEntry(h_list[4],"H(125)","f")
    legend_z.AddEntry(h_list[3],"q#bar{q}#rightarrow ZZ,Z#gamma*","f")
    legend_z.AddEntry(h_list[2],"gg#rightarrow ZZ,Z#gamma*","f")
    legend_z.AddEntry(h_list[1],"EW","f")
    legend_z.AddEntry(h_list[0],"Z+X","f")
    legend_z.AddEntry(HData,"Data", "p")
    legend_z.SetFillColor(ROOT.kWhite)
    legend_z.SetLineColor(ROOT.kWhite)
    legend_z.SetTextFont(43)
    legend_z.SetTextSize(20)
    legend_z.Draw()
    
    #draw CMS and lumi text
    CMS_lumi.writeExtraText = True
    CMS_lumi.extraText      = "Preliminary"
    CMS_lumi.lumi_sqrtS     = "35.1 fb-1 (13.6 TeV)"
    CMS_lumi.cmsTextSize    = 1 #0.6
    CMS_lumi.lumiTextSize   = 0.7 #0.46
    CMS_lumi.extraOverCmsTextSize = 0.75
    CMS_lumi.relPosX = 0.12
    CMS_lumi.CMS_lumi(Canvas_z, 0, 0)
    
    Canvas_z.Update()
    
    # ### Zoomed high mass
    # Canvas_hm = ROOT.TCanvas("M4l_hm","M4l_hm",canvasSizeX,canvasSizeY)
    # Canvas_hm.SetTicks()
    # Canvas_hm.SetLogy()
    # HStack_hm.Draw("histo")
    # HStack_hm.GetXaxis().SetRangeUser(170.,1000.)
    # HData_hm.Draw("samePE1")
    # Canvas_hm.Update()
    
    # ### High mass, 10 GeV
    # HStack10, h_list = Stack(fMC2018, fMC2022, "_10GeV_")
    # HData10 = dataGraph(fData, "_10GeV_", blind=blindPlots)
    # Canvas10 = ROOT.TCanvas("M4l_hm10","M4l_hm10",canvasSizeX,canvasSizeY)
    # Canvas10.SetTicks()
    # Canvas10.SetLogy()
    # HStack10.Draw("histo")
    # HStack10.GetXaxis().SetRangeUser(170.,1000.)
    # HData10.Draw("samePE0E1")
    # Canvas10.Update()

    
    printCanvases(path=out_dir)