### Example macro for filling standard histograms from H4l nanoAODs.
### Histograms are stored on a file and can then be plotted with

from __future__ import print_function

import os
import datetime

import math
import ctypes
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
        prop, which, reg = info["prop"], info["which"], info["reg"]
        if reg == "SR":
            idx  = self._cand_maps[which][reg]["idx"]
            col  = self._cand_maps[which][reg]["col"](prop)
        else:
            idx = self._cand_maps[which]["CR"]["idx"](reg)
            col = self._cand_maps[which]["CR"]["col"](prop)
        return col, idx, which
    
    def _fill_hist(self, fname, proc, info):
            col, idx, which = self.select_cand(info)
            props = self.get_props(proc, info)

            Runs = ROOT.RDataFrame("Runs", fname)
            df = ROOT.RDataFrame("Events", fname)
            df = df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue()))

            filt = df.Filter("({}!=-1) && HLT_passZZ4l".format(idx))
            wgt = filt.Define("weight", "overallEventWeight*ZZCand_dataMCWeight/genEventSumw")

            hist = wgt.Define(which, "{}[{}]".format(col, idx)) # --> Need to broadcast which (a scalar float) to an RVec in order to multiply it by the RVec weight
            hist = hist.Define(which+"_vec", "return ROOT::VecOps::RVec<Float_t>(weight.size(), {});".format(which))
            
            hist_weighted = hist.Histo1D((props["name"], props["title"], props["nbins"], props["xlow"], props["xhigh"]), which+"_vec", "weight").GetValue()

            return hist_weighted

    def fillHistos(self, **kwargs):
        hists_weighted = {}
        for key, proc_dict in kwargs["proc_info"].items():
            if proc_dict["samples"] is None:
                continue
            samples = list(proc_dict["samples"].items())
            proc_0, fname_0 = samples[0][0], samples[0][1]
            hist = self._fill_hist(fname_0, proc_0, kwargs["hist_info"])
            Hist = hist.Clone()
            for sample in samples[1:]:
                proc, fname = sample[0], sample[1]
                Hist.Add(self._fill_hist(fname, proc, kwargs["hist_info"]))
            Hist.Scale(Lum)
            hists_weighted[key] = Hist

        return hists_weighted
    
    def runZZ(self, **kwargs):
        hists_weighted = self.fillHistos(**kwargs)
        hists_weighted.update({
            "Z+X": getZX(hists_weighted["H(125)"], "fs_4l")
        })
        hists_weighted["Z+X"].Scale(Lum)

        # Explicit ordering to best compare to Alessandra's plot
        self.hists = {
            "Z+X":  hists_weighted["Z+X"],
            "EW" :  hists_weighted["EW"],
            "q#bar{q}#rightarrow ZZ,Z#gamma*": hists_weighted["q#bar{q}#rightarrow ZZ,Z#gamma*"],
            "H(125)": hists_weighted["H(125)"]
        }

class ZZPlotter:
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
        hist_info = kwargs["hist_info"]


        HStack = ROOT.THStack("Stack_{}{}".format(hist_info["which"], hist_info["prop"]),
                          "; {} ; {}".format(hist_info["x_title"], hist_info["y_title"]))

        for proc, h in self.zz.hists.items():
            h.GetXaxis().SetTitle(hist_info["x_title"])
            h.GetYaxis().SetTitle(hist_info["y_title"])

            line_color = kwargs["proc_info"][proc]["line_color"]
            fill_color = kwargs["proc_info"][proc]["fill_color"]

            h.SetLineColor(ROOT.TColor.GetColor(line_color))
            h.SetFillColor(ROOT.TColor.GetColor(fill_color))

            HStack.Add(h, "HISTO")

        yhmax  = math.ceil(max(HStack.GetMaximum(), 0.))
        HStack.SetMinimum(1.0)
        HStack.SetMaximum(yhmax)
        HStack.Draw("histo")
        HStack.GetXaxis().SetRangeUser(hist_info["xlow"], hist_info["xhigh"])

        return HStack
    
    def _get_legend(self, loc=(0.72, 0.70, 0.94, 0.92)):
        Legend = ROOT.TLegend(*loc)
        for proc, h in self.zz.hists.items():
            Legend.AddEntry(h, proc, "f")

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
        CMS_lumi.lumi_sqrtS           = "{} fb-1 (13.6 TeV)".format(round(Lum*1e-3, 1))
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
        unit = "GeV" if kwargs["prop"] not in ang_vars else ""
        
        axes = self._get_axes(Canvas, **kwargs)
        axes_labs = "per{}{}_{}to{}{}_{}".format(kwargs["step"], unit, kwargs["xlow"], kwargs["xhigh"], unit, axes)

        direc = "plots/{}".format(datetime.date.today())
        os.makedirs(direc, exist_ok=True)

        plot_path = "{}_{}_CONFIG_TEST.png".format(kwargs["prop"],
                                                  axes_labs,
                                                )

        Canvas.SaveAs(os.path.join(direc, plot_path))

    def plot(self, **kwargs):

        self._set_gStyle()

        Canvas = ROOT.TCanvas(kwargs["hist_info"]["prop"], kwargs["hist_info"]["prop"],900, 700)
        Canvas.SetTicks()

        HStack = self._get_hStacks(**kwargs)
        Legend = self._get_legend(kwargs["hist_info"]["legend_loc"])

        self._set_lumi(Canvas, **kwargs)

        if kwargs["hist_info"]["xlabels"] is not None:
            xlabels = self._get_xlabels(kwargs["hist_info"]["xlabels"])
            for label in xlabels:
                label.Draw()
            ROOT.gPad.RedrawAxis()

        self._save_plot(Canvas, **kwargs["hist_info"])

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