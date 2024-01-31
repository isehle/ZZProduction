import ROOT
from pprint import pprint
import json

class fakeRates:
    def __init__(self, reg, cat, v, cfg):
        self.reg      = reg
        self.cat      = cat
        self.verbose  = v
        self.cfg      = cfg
        self.filt_cr  = lambda df: df.Filter(f"ZLLbest{reg}Idx!=-1")
        self.filt_sr  = lambda df: df.Filter("bestCandIdx!=-1")
        
        self.tot_el   = lambda df: df.Sum("nElectron").GetValue()
        self.tot_mu   = lambda df: df.Sum("nMuon").GetValue()

        self.pass_el  = lambda df: df.Sum("Electron_ZZFullId").GetValue()
        self.pass_mu  = lambda df: df.Sum("Muon_ZZFullId").GetValue()
        
        self.fail_el = lambda df, pass_el: df.Sum("Electron_ZZRelaxedId").GetValue() - pass_el
        self.fail_mu = lambda df, pass_mu: df.Sum("Muon_ZZRelaxedId").GetValue() - pass_mu

    def _get_rdf(self, samples):
        return ROOT.RDataFrame("Events", samples)
    
    def _rates_by_reg(self, rdf):
        pass_el = self.pass_el(rdf)
        pass_mu = self.pass_mu(rdf)

        fail_el = self.fail_el(rdf, pass_el)
        fail_mu = self.fail_mu(rdf, pass_mu)

        return dict(
            Electron = dict(
                pass_perc = pass_el/self.n_el,
                fail_perc = fail_el/self.n_el
            ),
            Muon = dict(
                pass_perc = pass_mu/self.n_mu,
                fail_perc = fail_mu/self.n_mu
            )
        )  
    
    def run_cat(self, category):
        proc_dict = self.cfg["proc_info"][category]
        fake_rates = {}
        for era in proc_dict["eras"]:
            samples = proc_dict["eras"][era]["samples"].values()

            tot_rdf = self._get_rdf(samples)
            self.n_el = self.tot_el(tot_rdf)
            self.n_mu = self.tot_mu(tot_rdf)

            sr_rdf  = self.filt_sr(tot_rdf)

            fake_rates[era] = {}
            if self.reg == None:
                cols     = map(str, list(tot_rdf.GetColumnNames()))
                reg_cols = [col for col in cols if "ZLLbest" in col]
                for col in reg_cols:
                    cr_rdf               = tot_rdf.Filter(f"{col}!=-1")
                    reg                  = col.replace("ZLLbest", "").replace("Idx","")
                    fake_rates[era][reg] = self._rates_by_reg(cr_rdf)
            else:
                cr_rdf                    = self.filt_cr(tot_rdf)
                fake_rates[era][self.reg] = self._rates_by_reg(cr_rdf)

            fake_rates[era]["SR"]    = self._rates_by_reg(sr_rdf)
            fake_rates[era]["Total"] = self._rates_by_reg(tot_rdf)

        if self.verbose:
            print(f"\n{category}")
            print("=="*30)
            pprint(fake_rates, width=50, compact=True)

        return fake_rates
    
    def run(self, outfile=None):
        fake_rates = {}
        if self.cat is None:
            for cat in self.cfg["proc_info"].keys():
                if cat == "H(125)":
                    break
                else:
                    fake_rates[cat] = self.run_cat(cat)
        else:
            fake_rates[self.cat] = self.run_cat(self.cat)

        if outfile is not None:
            with open(outfile, "w") as outfile:
                json.dump(fake_rates, outfile, indent=4)

if __name__=="__main__":
    import importlib.util
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Configuration")
    parser.add_argument("cfg_path")
    parser.add_argument("--reg", help="Control region", choices=("2P2F", "3P1F", "SS", "SIP"), required=False)
    parser.add_argument("--cat", help="Category", type=str, required=False)
    parser.add_argument("--outfile", help="Name for output json file (if desired)", type=str)
    parser.add_argument("--v", help="verbosity: print output if true", action="store_true")
    args = vars(parser.parse_args())

    spec = importlib.util.spec_from_file_location("cfg", args["cfg_path"])
    cfg_script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg_script)

    info = cfg_script.all_info

    fr = fakeRates(args["reg"], args["cat"], args["v"], info)
    fr.run(outfile=args["outfile"])