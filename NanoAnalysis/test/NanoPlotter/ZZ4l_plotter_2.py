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
        self._filt_sr   = lambda df: df.Filter("bestCandIdx!=-1 && HLT_passZZ4l")
        self._filt_cr   = lambda df, reg: df.Filter(f"ZLLbest{reg}Idx!=-1 && HLT_passZZ4l")
        self._sig_col   = lambda which, prop: "ZZCand_{}{}".format("" if which=="ZZ" else which, prop)
        self._cr_col    = lambda which, prop: "ZLLCand_{}{}".format("" if which=="ZZ" else which, prop)

    def _define_sr_cand(self, df, prop, which):
        col = "ZZCand_{}{}".format("" if which=="ZZ" else which, prop)
        return df.Define(which, "{}[bestCandIdx]".format(col))
            
    def _fill_hist(self, proc_info, hist_info):
        """Fills and weighs 1D histograms of the
        intereseted property (info["prop"]) for the
        process whose path is given by fname."""

        proc = list(proc_info["samples"].keys())[0]
        samples = proc_info["samples"].values()
        
        name = hist_info["which"] + "_" + proc
        title = hist_info["which"] + "_" + proc
        xhigh, xlow = hist_info["xhigh"], hist_info["xlow"]
        nbins = int((xhigh - xlow)/hist_info["step"])

        df = ROOT.RDataFrame("Events", samples)

        breakpoint()

        if proc != "Data" and proc != "Data":
            Runs = ROOT.RDataFrame("Runs", samples)
            df = df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue()))

        if hist_info["reg"]=="SR":
            col = self._sig_col(hist_info["which"], hist_info["prop"])
            filt = self._filt_sr(df).Define(hist_info["which"], "{}[bestCandIdx]".format(col))
        else:
            col = self._cr_col(hist_info["which"], hist_info["prop"])
            filt = self._filt_cr(df, hist_info["reg"]).Define(hist_info["which"], "{}[ZLLbest{}Idx]".format(col, hist_info["reg"]))

        if proc != "Data" and proc != "Data":
            # In CRs ZZCand_dataMCWeight is 0.0 ...
            if hist_info["reg"]=="SR":
                wgt = filt.Define("weight", "overallEventWeight*ZZCand_dataMCWeight/genEventSumw")
                #Need to broadcast "which" (a scalar float) to an RVec in order to multiply it by the RVec weight
                vec = wgt.Define(hist_info["which"]+"_vec", "return ROOT::VecOps::RVec<Float_t>(weight.size(), {});".format(hist_info["which"]))
                
                hist = vec.Histo1D((name, title, nbins, xlow, xhigh), hist_info["which"]+"_vec", "weight").GetValue()
            else:
                wgt = filt.Define("weight", "overallEventWeight/genEventSumw")
                hist = wgt.Histo1D((name, title, nbins, xlow, xhigh), hist_info["which"], "weight").GetValue()

            hist.Scale(proc_info["lum"])
        else:
            hist = filt.Histo1D((name, title, nbins, xlow, xhigh), hist_info["which"]).GetValue()
            hist.SetBinErrorOption(ROOT.TH1.kPoisson)

        return hist

    def fillHistos(self, **kwargs):
        """Combines individual process histograms into
        each category (ex. EW: WWZ+WZZ+ZZZ+...), and
        scales them by the luminosity."""
        hists_weighted = {}
        for cat, cat_dict in kwargs["proc_info"].items():
            for proc_dict in cat_dict["eras"].values():
                if proc_dict["samples"] is None:
                    continue

                hists_weighted[cat] = self._fill_hist(proc_dict, kwargs["hist_info"])

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
        plot_info = kwargs["plot_info"]

        HStack = ROOT.THStack("Stack_{}{}".format(hist_info["which"], hist_info["prop"]),
                          "; {} ; {}".format(plot_info["x_title"], plot_info["y_title"]))

        for proc, h in self.zz.hists.items():
            h.GetXaxis().SetTitle(plot_info["x_title"])
            h.GetYaxis().SetTitle(plot_info["y_title"])

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
    
    def _save_plot(self, Canvas, reg, **kwargs):
        """Set up plot directory (based on the date),
        the plot title (based on the plot info) and saves
        the plot."""

        date_direc = "plots/{}".format(datetime.date.today())
        os.makedirs(date_direc, exist_ok=True)

        reg_direc = os.path.join(date_direc, reg)
        os.makedirs(reg_direc, exist_ok=True)

        plot_path = "{}_{}_{}.{}".format(kwargs["out_file"],
                                      "logX" if kwargs["logx"] else "linX",
                                      "logY" if kwargs["logy"] else "linY",
                                      kwargs["format"])

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
        HStack.GetXaxis().SetRangeUser(*kwargs["plot_info"]["x_range"])

        if kwargs["plot_info"]["logy"] and kwargs["plot_info"]["y_min"]==0:
            yhmin = 1e-3
        else:
            yhmin = kwargs["plot_info"]["y_min"]

        HStack.SetMinimum(yhmin)
        HStack.SetMaximum(yhmax)

        if "Data" in kwargs["proc_info"]:
            Data = self._get_data(**kwargs)
            Data.Draw("samePE1")

        Legend = self._get_legend(kwargs["plot_info"]["legend_loc"])

        self._set_lumi(Canvas, **kwargs)

        if kwargs["plot_info"]["xlabels"] is not None:
            xlabels = self._get_xlabels(kwargs["plot_info"]["xlabels"])
            for label in xlabels:
                label.Draw()
            ROOT.gPad.RedrawAxis()

        if kwargs["plot_info"]["logx"]:
            Canvas.SetLogx()
        if kwargs["plot_info"]["logy"]:
            Canvas.SetLogy()

        self._save_plot(Canvas, kwargs["hist_info"]["reg"], **kwargs["plot_info"])

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

    ZZ = ZZHists()
    ZZ.runZZ(**info)

    ZZPlot = ZZPlotter(ZZ)
    ZZPlot.plot(**info)