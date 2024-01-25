### Example macro for filling standard histograms from H4l nanoAODs.
### Histograms are stored on a file and can then be plotted with

from __future__ import print_function

import os
import datetime

import math
import ctypes
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons

import CMSGraphics, CMS_lumi
from H4l_draw_mZZ_periods2022 import getZX

# 1/fb--> 1/pb 2022EFG MC - Run3Summer22EE from https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
Lum = 27.007e3  # w/ normtag, golden JSON

ang_vars = ["eta", "cos", "phi"]

class ZZHists:
    def __init__(self):
        self._samples = None
        self.hists = None
        self.props = None
        self.ZmassValue = 91.1876
        self._cand_maps = dict(
                            ZZ = dict(
                                SR = dict(
                                    idx = "bestCandIdx",
                                    col = lambda prop: "ZZCand_{}".format(prop)
                                ),
                                CR = dict(
                                    idx = lambda reg: "ZLLbest{}Idx".format(reg),
                                    col = lambda prop: "ZLLCand_{}".format(prop)
                                )
                            ),
                            Z1 = dict(
                                SR = dict(
                                    idx = "bestZIdx",
                                    col = lambda prop: "ZCand_{}".format(prop)
                                )
                            )
                        )
        
    def get_props(self, samplename, info):
        """Return dictionary of histogram properties based off
        of the histogram info and the samplename."""
        props = dict(
            nbins = int((info["xhigh"] - info["xlow"])/info["step"]),
            xlow  = info["xlow"],
            xhigh = info["xhigh"],
        )
        if samplename is not None:
            unit = "GeV" if info["prop"] not in ang_vars else ""
            if info["prop"] == "cos":
                name = "cosTheta_{}".format(1*(info["which"]=="Z1") + 3*(info["which"]=="Z2"))
                name += "_per{}_{}".format(info["step"], samplename)
            else:
                name = "{}{}_per{}{}_".format(info["which"], info["prop"], info["step"], unit)+samplename
            props.update(dict(
                            name  = name,
                            title = name,
            ))
        return props
    
    def select_cand(self, info):
        """Uses the info dictionary to return the correct column, index,
        and which cand (i.e. ZZ, Z1, Z2...) from the property (i.e. "mass")
        and region (i.e. "SR"), returning which directly."""
        prop, which, reg = info["prop"], info["which"], info["reg"]
        if reg == "SR":
            idx  = self._cand_maps[which][reg]["idx"]
            col  = self._cand_maps[which][reg]["col"](prop)
        else:
            idx = self._cand_maps[which]["CR"]["idx"](reg)
            col = self._cand_maps[which]["CR"]["col"](prop)
        return col, idx, which
    
    def _fill_hist(self, fname, proc, info):
        """Fills and weighs 1D histograms of the
        intereseted property (info["prop"]) for the
        process whose path is given by fname."""
        col, idx, which = self.select_cand(info)
        props = self.get_props(proc, info)

        df = ROOT.RDataFrame("Events", fname)
        if proc != "2022_CD" and proc != "2022_EFG": #data
            Runs = ROOT.RDataFrame("Runs", fname)
            df = df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue()))

        filt = df.Filter("({}!=-1) && HLT_passZZ4l".format(idx)).Define(which, "{}[{}]".format(col, idx))


        if proc != "2022_CD" and proc != "2022_EFG": #data
            wgt = filt.Define("weight", "overallEventWeight*ZZCand_dataMCWeight/genEventSumw")

            #Need to broadcast "which" (a scalar float) to an RVec in order to multiply it by the RVec weight
            vec = wgt.Define(which+"_vec", "return ROOT::VecOps::RVec<Float_t>(weight.size(), {});".format(which))
            
            #hist = vec.Histo1D((props["name"], props["title"], props["nbins"], props["xlow"], props["xhigh"]), which+"_vec", "weight").GetValue()
            hist = vec.Histo1D((props["name"], props["title"], 233, 70., 1002.), which+"_vec", "weight").GetValue()
        else:
            #hist = filt.Histo1D((props["name"], props["title"], props["nbins"], props["xlow"], props["xhigh"]), which).GetValue()
            hist = filt.Histo1D((props["name"], props["title"], 233, 70., 1002.), which).GetValue()

        return hist

    def fillHistos(self, **kwargs):
        """Combines individual process histograms into
        each category (ex. EW: WWZ+WZZ+ZZZ+...), and
        scales them by the luminosity."""
        hists_weighted = {}
        for cat, cat_dict in kwargs["proc_info"].items():
            era_hists = []
            
            for proc_dict in cat_dict["eras"].values():
                if proc_dict["samples"] is None:
                    continue
                
                samples = list(proc_dict["samples"].items())

                proc_0, fname_0 = samples[0][0], samples[0][1]
                hist = self._fill_hist(fname_0, proc_0, kwargs["hist_info"])
                Hist = hist.Clone()
                
                for sample in samples[1:]:
                    proc, fname = sample[0], sample[1]
                    Hist.Add(self._fill_hist(fname, proc, kwargs["hist_info"]))
                
                if cat != "Data":
                    Hist.Scale(proc_dict["lum"])
                else:
                    Hist.SetBinErrorOption(ROOT.TH1.kPoisson)
                era_hists.append(Hist)
            
            if len(era_hists) > 0:
                final_hist = era_hists[0]
                for hi in era_hists[1:]:
                    final_hist.Add(hi, 1.)
                hists_weighted[cat] = final_hist

        return hists_weighted
    
    def runZZ(self, **kwargs):
        """Calling function for the class,
        acting as a wrapper for the primary class functions
        and setting the hists property."""
        hists = self.fillHistos(**kwargs)
        if "Z+X" in kwargs["proc_info"]:
            hists.update({
                "Z+X": getZX(hists["H(125)"], "fs_4l")
            })
            hists["Z+X"].Scale(kwargs["proc_info"]["Z+X"]["eras"]["2022_EE"]["lum"])

        self.hists = dict(sorted(hists.items(), key = lambda item: item[1].Integral()))

class ZZPlotter:
    """Plotting class for histograms created with the
    ZZHists class."""
    def __init__(self, zz):
        self.zz = zz

    def _set_gStyle(self):
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

    def _get_xlabels(self, xlabelsv):
        # Labels for log plots
        label_margin = -0.1
        xlabels=[None]*len(xlabelsv)
        for i, label in enumerate(xlabelsv): 
            xlabels[i] = ROOT.TLatex(label, label_margin , str(label));
            xlabels[i].SetTextAlign(23)
            xlabels[i].SetTextFont(42)
            xlabels[i].SetTextSize(0.04)

        return xlabels
    
    def _get_hStacks(self, **kwargs):
        """Create HStack from individual histograms."""
        hist_info = kwargs["hist_info"]

        HStack = ROOT.THStack("Stack_{}{}".format(hist_info["which"], hist_info["prop"]),
                          "; {} ; {}".format(hist_info["x_title"], hist_info["y_title"]))

        for proc, h in self.zz.hists.items():
            h.GetXaxis().SetTitle(hist_info["x_title"])
            h.GetYaxis().SetTitle(hist_info["y_title"])

            if proc != "Data":
                line_color = kwargs["proc_info"][proc]["line_color"]
                fill_color = kwargs["proc_info"][proc]["fill_color"]
                
                h.SetLineColor(ROOT.TColor.GetColor(line_color))
                h.SetFillColor(ROOT.TColor.GetColor(fill_color))
                HStack.Add(h, "HISTO")

        return HStack
    
    def _get_data(self, **kwargs):
        data_info = kwargs["proc_info"]["Data"]

        hist_data = self.zz.hists["Data"]
        
        nbinsIn = hist_data.GetNbinsX()
        nbins = 0
        
        x = np.array([0.]*nbinsIn, dtype='double')
        y = np.array([0.]*nbinsIn, dtype='double')
        errX = np.array([0.]*nbinsIn, dtype='double')  
        UpErr = np.array([0.]*nbinsIn, dtype='double')
        LowErr = np.array([0.]*nbinsIn, dtype='double')

        for i in range (1, nbinsIn):
            x[nbins]      = hist_data.GetBinCenter(i)
            y[nbins]      = hist_data.GetBinContent(i)
            UpErr[nbins]  = hist_data.GetBinErrorUp(i)
            LowErr[nbins] = hist_data.GetBinErrorLow(i)
            nbins += 1

        Data = ROOT.TGraphAsymmErrors(nbins, x, y, errX, errX, LowErr, UpErr)
        Data.SetMarkerStyle(data_info["marker_style"])
        Data.SetLineColor(data_info["line_color"])
        Data.SetMarkerSize(data_info["marker_size"])
        
        return Data
    
    def _get_legend(self, loc=(0.72, 0.70, 0.94, 0.92)):
        Legend = ROOT.TLegend(*loc)
        for proc, h in self.zz.hists.items():
            style = "f" if proc != "Data" else "p"
            Legend.AddEntry(h, proc, style)

        Legend.SetFillColor(ROOT.kWhite)
        Legend.SetLineColor(ROOT.kWhite)
        Legend.SetTextFont(43)
        Legend.SetTextSize(20)
        Legend.Draw()

        return Legend
    
    def _set_lumi(self, Canvas, **kwargs):
        """Set up CMS lumi info--this info should be set
        in the configuration file and saved in kwargs."""
        CMS_lumi.writeExtraText       = True
        CMS_lumi.extraText            = "Preliminary"
        #CMS_lumi.lumi_sqrtS           = "{} fb-1 (13.6 TeV)".format(round(Lum*1e-3, 1))
        CMS_lumi.lumi_sqrtS           = "27.01 fb-1 (13.6 TeV)"
        CMS_lumi.cmsTextSize          = 1
        CMS_lumi.lumiTextSize         = 0.7
        CMS_lumi.extraOverCmsTextSize = 0.75
        CMS_lumi.relPosX              = 0.12
        CMS_lumi.CMS_lumi(Canvas, 0, 0)
        Canvas.Update()


    def _get_axes(self, Canvas, **kwargs):
        x_ax, y_ax = "linX", "linY"
        if kwargs["logx"]:
            Canvas.SetLogx()
            x_ax = "logX"
        if kwargs["logy"]:
            Canvas.SetLogy()
            y_ax = "logY"
        axes = "{}{}".format(x_ax, y_ax)

        return axes
    
    def _save_plot(self, Canvas, **kwargs):
        """Set up plot directory (based on the date),
        the plot title (based on the plot info) and saves
        the plot."""
        unit = "GeV" if kwargs["prop"] not in ang_vars else ""
        
        axes = self._get_axes(Canvas, **kwargs)
        axes_labs = "per{}{}_{}to{}{}_{}".format(kwargs["step"], unit, kwargs["xlow"], kwargs["xhigh"], unit, axes)

        date_direc = "plots/{}".format(datetime.date.today())
        os.makedirs(date_direc, exist_ok=True)

        reg_direc = os.path.join(date_direc, kwargs["reg"])
        os.makedirs(reg_direc, exist_ok=True)

        plot_path = "{}_{}_my_2022EE_samples_sorted.{}".format(kwargs["which"]+"_"+kwargs["prop"],
                                      axes_labs,
                                      kwargs["format"]
                                    )

        Canvas.SaveAs(os.path.join(reg_direc, plot_path))

    def plot(self, **kwargs):
        """Main plotting function, primarily acting
        as a wrapper for the class's primary functions."""

        self._set_gStyle()

        Canvas = ROOT.TCanvas(kwargs["hist_info"]["prop"], kwargs["hist_info"]["prop"],900, 700)
        Canvas.SetTicks()

        HStack = self._get_hStacks(**kwargs)
        yhmax  = math.ceil(max(HStack.GetMaximum(), 0.)) + 15.

        HStack.Draw("histo")
        HStack.GetXaxis().SetRangeUser(kwargs["hist_info"]["xlow"], kwargs["hist_info"]["xhigh"])
        HStack.GetYaxis().SetRangeUser(0.0, yhmax)
        if "Data" in kwargs["proc_info"]:
            Data = self._get_data(**kwargs)
            Data.Draw("samePE1")

        Legend = self._get_legend(kwargs["hist_info"]["legend_loc"])

        self._set_lumi(Canvas, **kwargs)

        if kwargs["hist_info"]["xlabels"] is not None:
            xlabels = self._get_xlabels(kwargs["hist_info"]["xlabels"])
            for label in xlabels:
                label.Draw()
            ROOT.gPad.RedrawAxis()

        self._save_plot(Canvas, **kwargs["hist_info"])

def get_rdfs(fname):
    Runs = ROOT.RDataFrame("Runs", fname)
    df = ROOT.RDataFrame("Events", fname)
    df = df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue()))

    filt = df.Filter("(bestCandIdx!=-1)&&(HLT_passZZ4l)")
    wgt = filt.Define("weight", "overallEventWeight*ZZCand_dataMCWeight/genEventSumw")

    rdf = wgt.Define("ZZ", "ZZCand_mass[bestCandIdx]")
    rdf = rdf.Define("ZZ_vec", "return ROOT::VecOps::RVec<Float_t>(weight.size(), ZZ);") 

    return rdf

def get_zx_rdf(info):
    sample_dict = info["proc_info"]["H(125)"]["samples"]
    samples = list(sample_dict.items())
    proc_0, fname_0 = samples[0][0], samples[0][1]

    hist = get_rdfs(fname_0).Histo1D(("mass_"+proc_0, "mass_"+proc_0, 233, 70., 1002.), "ZZ_vec", "weight").GetValue()
    Hist = hist.Clone()
    for sample in samples[1:]:
        proc, fname = sample
        hist = get_rdfs(fname).Histo1D(("mass_"+proc, "mass_"+proc, 233, 70., 1002.), "ZZ_vec", "weight").GetValue()
        Hist.Add(hist)
    Hist.Scale(Lum)

    return getZX(Hist, "fs_4l")

def validation(info):
    f2022 = ROOT.TFile.Open('H4l_MC2022_Test23Jan_SigOnly.root', "READ")
    name = "ZZMass_4GeV_"

    #-----------signal------------#
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
    signal.Scale(Lum)

    zpx     = getZX(signal, "fs_4l")
    zpx_rdf = get_zx_rdf(info)
    zpx.Add(zpx_rdf, -1) # Subtract

    zpx.SetTitle("Z+X: Original - RDataFrame")

    Canvas = ROOT.TCanvas("Z+X: Original - RDataFrame", "Z+X: Original - RDataFrame", 900, 700)
    zpx.DrawCopy()
    Canvas.SaveAs("plots/{}/zpx_validation.png".format(datetime.date.today()))


if __name__ == "__main__":
    import importlib.util
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Configuration")
    parser.add_argument("cfg_path")
    args = vars(parser.parse_args())

    spec = importlib.util.spec_from_file_location("cfg", args["cfg_path"])
    cfg_script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg_script)

    info = cfg_script.all_info

    #validation(info)

    ZZ = ZZHists()
    ZZ.runZZ(**info)

    ZZPlot = ZZPlotter(ZZ)
    ZZPlot.plot(**info)