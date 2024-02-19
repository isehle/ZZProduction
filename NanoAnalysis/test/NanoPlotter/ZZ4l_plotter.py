### Example macro for filling standard histograms from H4l nanoAODs.
### Histograms are stored on a file and can then be plotted with

from __future__ import print_function

import os
import datetime

import time

import math
import ctypes
import numpy as np
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons

import CMSGraphics, CMS_lumi
from H4l_draw_mZZ_periods2022 import getZX
from RDF_Helpers import *

# 1/fb--> 1/pb 2022EFG MC - Run3Summer22EE from https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis
Lum = 27.007e3  # w/ normtag, golden JSON

ang_vars = ["eta", "cos", "phi"]

class ZZHists:
    def __init__(self):
        self._samples = None
        self.hists = None
        self.props = None
        self.ZmassValue = 91.1876

        self.new_col    = None

        self._idx       = lambda reg: "bestCandIdx" if reg == "SR" else "ZLLbest{}Idx".format(reg)
        self._filt_df   = lambda df, reg: df.Filter("{}!=-1 && HLT_passZZ4l".format(self._idx(reg)))

        self._sig_col   = lambda which, prop: "ZZCand_{}{}".format("" if which=="ZZ" else which, prop)
        self._cr_col    = lambda which, prop: "ZLLCand_{}{}".format("" if which=="ZZ" else which, prop)
        self._sel_col   = lambda which, prop, reg: self._sig_col(which, prop) if reg=="SR" else self._cr_col(which, prop)

        self._lep_col    = lambda lep, prop: f"{lep}_{prop}"
        self._lep_arrs   = lambda df, Z, prop: df.AsNumpy(["Electron_charge",
                                                           "Muon_charge",
                                                           "ZZCand_{}flav".format(Z),
                                                           self._lep_col("Electron", prop),
                                                           self._lep_col("Muon", prop)])

        self._weight_arrs = lambda df: df.AsNumpy(["overallEventWeight",
                                                   "ZZCand_dataMCWeight",
                                                   "genEventSumw"])
        
        self._get_lep_arr = lambda dict, key: np.array(
            list(
                map(
                    lambda x: x[0], dict[key]
                )
            )
        )

        self._charge_mask = lambda charges, lep: charges==-1 if lep=="lm" else charges==1

    def _fill_lep_var(self, filt, hist_info, prop):
        theZed, theLep = hist_info["which"].split("_")

        charge = -1 if theLep=="lm" else 1

        z_flav = "ZZCand_{}flav".format(theZed)
        z_l1idx = "ZZCand_{}l1Idx".format(theZed)
        z_l2idx = "ZZCand_{}l2Idx".format(theZed)

        el_prop = "Electron_{}".format(prop)
        mu_prop = "Muon_{}".format(prop)

        self.new_col = "Lepton_{}".format(prop)

        filt = filt.Define(self.new_col, "LepFromZ({}, Electron_charge, {}, Muon_charge, {}, {}, {}, {})".format(\
            el_prop, mu_prop, z_flav, z_l1idx, z_l2idx, charge
            ))

        return filt

    def _cand_from_df(self, df, hist_info, prop=None):
        reg_filt = self._filt_df(df, hist_info["reg"]) if self.new_col is None else df

        true_prop = hist_info["prop"] if prop is None else prop

        if "_" in hist_info["which"]:
            return self._fill_lep_var(reg_filt, hist_info, true_prop)
        else:
            self.new_col = "{}_{}".format(hist_info["which"], true_prop)
            sel_col = self._sel_col(hist_info["which"], true_prop, hist_info["reg"])
            return reg_filt.Define(self.new_col, "{}[{}]".format(sel_col, self._idx(hist_info["reg"])))
    
    def _defCosTheta(self, df, hist_info):
        eta_cand = self.new_col
        self.new_col = "{}_cos".format(hist_info["which"])
        return df.Define(self.new_col, "cos(2*atan(exp(-{})))".format(eta_cand))        

    def _getHist(self, filt, hist_info, proc):
        name = hist_info["which"] + "_" + proc
        title = hist_info["which"] + "_" + proc
        xhigh, xlow = hist_info["xhigh"], hist_info["xlow"]
        nbins = int((xhigh - xlow)/hist_info["step"])

        if proc != "Data":
            if (hist_info["reg"] == "SR") and ("_" not in hist_info["which"]):
                wgt = filt.Define("weight", "overallEventWeight*ZZCand_dataMCWeight/genEventSumw")
                vec = wgt.Define(self.new_col+"_vec", "return ROOT::VecOps::RVec<Float_t>(weight.size(), {});".format(self.new_col))
                hist = vec.Histo1D((name, title, nbins, xlow, xhigh), self.new_col+"_vec", "weight").GetValue()
            else:
                wgt = filt.Define("weight", "overallEventWeight/genEventSumw")
                hist = wgt.Histo1D((name, title, nbins, xlow, xhigh), self.new_col, "weight").GetValue()
        else:
            hist = filt.Histo1D((name, title, nbins, xlow, xhigh), self.new_col).GetValue()
            hist.SetBinErrorOption(ROOT.TH1.kPoisson)

        return hist
    
    def _fill_hist(self, proc_info, hist_info):
        """Fills and weighs 1D histograms of the
        intereseted property (info["prop"]) for the
        process whose path is given by fname."""

        proc = list(proc_info["samples"].keys())[0]
        samples = proc_info["samples"].values()

        df = ROOT.RDataFrame("Events", samples)
        if proc != "Data":
            Runs = ROOT.RDataFrame("Runs", samples)
            df = df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue()))

        if hist_info["prop"] == "cos":
            eta_df = self._cand_from_df(df, hist_info, prop = "eta") # theta is defined by the eta of selected cands
            filt   = self._defCosTheta(eta_df, hist_info)
        else:
            filt = self._cand_from_df(df, hist_info)

        hist = self._getHist(filt, hist_info, proc)

        self.new_col = None

        if proc != "Data": hist.Scale(proc_info["lum"])

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
    
    def _get_sqrt(self, hist):
        hist_sqrt = hist.Clone("hist_sqrt")
        for i in range(1, hist.GetNbinsX()+1):
            bin_content = hist.GetBinContent(i)
            bin_error   = hist.GetBinError(i)

            # Temporary fix, why do we have some bins with negative values??
            if bin_content < 0: bin_content = 0
            
            new_content = ROOT.TMath.Sqrt(bin_content)
            new_error   = bin_error / (2*new_content) if new_content != 0 else 0

            hist_sqrt.SetBinContent(i, new_content)
            hist_sqrt.SetBinError(i, new_error)

        return hist_sqrt

    def _get_ratio(self, plot_info):
        Sig = self.zz.hists['q#bar{q}#rightarrow ZZ,Z#gamma*'].Clone("Sig")
        Sig.Add(self.zz.hists["gg#rightarrow ZZ,Z#gamma*"])

        Bkg = self.zz.hists["DY"].Clone("Bkg")
        for proc in ['WZ', 't#bar{t}', 'VVV', 'H(125)']:
            Bkg.Add(self.zz.hists[proc])
        
        Sqrt_Bkg = self._get_sqrt(Bkg)

        Ratio = Sig.Clone("Ratio")
        Ratio.SetLineColor(1)
        Ratio.SetMarkerStyle(21)
        Ratio.SetTitle("")
        Ratio.Sumw2()
        Ratio.SetStats(0)
        Ratio.Divide(Sqrt_Bkg)

        y = Ratio.GetYaxis()
        y.SetTitle(r"$Sig/\sqrt{Bkg} $")
        y.SetNdivisions(505)
        y.SetTitleSize(20)
        y.SetTitleFont(43)
        y.SetTitleOffset(1.55)
        y.SetLabelFont(43)
        y.SetLabelSize(15)

        # Adjust x-axis settings
        x = Ratio.GetXaxis()
        x.SetTitleSize(30)
        x.SetTitleFont(43)
        x.SetTitleOffset(1.0)
        x.SetLabelFont(43)
        x.SetLabelSize(15)
        x.SetLabelSize(30)
        x.SetTitle(plot_info["x_title"])

        return Ratio
    
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

        yhmax = math.ceil(max(HStack.GetMaximum(), 0.)) + 15.
        if plot_info["logy"] and plot_info["y_min"]==0:
            yhmin = 1e-3
        else:
            yhmin = kwargs["plot_info"]["y_min"]

        HStack.SetMinimum(yhmin)
        HStack.SetMaximum(yhmax)

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

        return Legend
    
    def _set_lumi(self, pad, **kwargs):
        """Set up CMS lumi info--this info should be set
        in the configuration file and saved in kwargs."""
        CMS_lumi.writeExtraText       = True
        CMS_lumi.extraText            = "Preliminary"
        CMS_lumi.lumi_sqrtS           = "{} fb-1 (13.6 TeV)".format(round(kwargs["lumin"], 3)*1e-3)
        CMS_lumi.cmsTextSize          = 1
        CMS_lumi.lumiTextSize         = 0.7
        CMS_lumi.extraOverCmsTextSize = 0.75
        CMS_lumi.relPosX              = 0.12
        CMS_lumi.CMS_lumi(pad, 0, 0)
    
    def _out_path(self, reg, **kwargs):
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
        
        return os.path.join(reg_direc, plot_path)
    
    def _write_pads(self, Canvas, **kwargs):
        hist_pad = ROOT.TPad("hist_pad", "hist_pad", 0, 0.3, 1, 0.95)

        hist_pad.SetBottomMargin(0)
        hist_pad.SetGridx()

        if kwargs["logx"]: hist_pad.SetLogx()
        if kwargs["logy"]: hist_pad.SetLogy()

        Canvas.cd()

        ratio_pad = ROOT.TPad("ratio_pad", "ratio_pad", 0, 0.1, 1, 0.3)

        ratio_pad.SetTopMargin(0)
        ratio_pad.SetBottomMargin(0.35)
        ratio_pad.SetGridx()

        if kwargs["logx"]: ratio_pad.SetLogx()
        if kwargs["ratio_diff"] > 1e2: ratio_pad.SetLogy()

        return hist_pad, ratio_pad
    
    def canvas_setup(self, **kwargs):

        self._set_gStyle()

        plot_info = kwargs["plot_info"]

        c_size = (1200, 1200) if "c_size" not in plot_info else plot_info["c_size"]
        Canvas = ROOT.TCanvas(kwargs["hist_info"]["prop"], kwargs["hist_info"]["prop"], c_size[0], c_size[1])
        
        return Canvas

    def plot(self, **kwargs):
        """Main plotting function, primarily acting
        as a wrapper for the class's primary functions."""

        HStack = self._get_hStacks(**kwargs)

        data = False
        if "Data" in kwargs["proc_info"]:
            Data = self._get_data(**kwargs)
            data = True

        plot_path = self._out_path(kwargs["hist_info"]["reg"], **kwargs["plot_info"])

        Canvas = self.canvas_setup(**kwargs)
        Legend = self._get_legend(kwargs["plot_info"]["legend_loc"])

        if kwargs["plot_info"]["ratio"]:
            Ratio = self._get_ratio(kwargs["plot_info"])

            kwargs["plot_info"]["ratio_diff"] = Ratio.GetMaximum() - Ratio.GetMinimum()

            hist_pad, ratio_pad = self._write_pads(Canvas, **kwargs["plot_info"])

            hist_pad.Draw()
            hist_pad.cd()
            HStack.Draw("histo")
            HStack.GetXaxis().SetRangeUser(*kwargs["plot_info"]["x_range"])

            if data:
                Data.Draw("samePE1")
            
            Legend.Draw()

            self._set_lumi(hist_pad, **kwargs)
            hist_pad.Update()

            Canvas.cd()

            ratio_pad.Draw()

            ratio_pad.cd()

            Ratio.Draw("ep")
        else:
            HStack.Draw("histo")
            if data:
                Data.Draw("samePE1")
            Legend.Draw()

            self._set_lumi(Canvas, **kwargs)
            Canvas.Update()

        Canvas.SaveAs(plot_path)


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