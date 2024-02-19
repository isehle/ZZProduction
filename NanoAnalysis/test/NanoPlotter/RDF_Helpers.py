import ROOT
import numpy as np
                          
ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
using ROOT::RVecI;
RVecF getLepProp(RVecF &Lep_prop, RVecI &Lep_charge, short l1Idx, short l2Idx, int charge){
    RVecF goodLepProp;

    if ( (Lep_charge.size() > l1Idx) && (Lep_charge.size() > l2Idx) ){
        int the_l1_charge = Lep_charge.at(l1Idx);
        int the_l2_charge = Lep_charge.at(l2Idx);
        
        float the_prop;

        if ( the_l1_charge == charge ){
            the_prop = Lep_prop.at(l1Idx);
        }
        else if ( the_l2_charge == charge ){
            the_prop = Lep_prop.at(l2Idx);
        }

        goodLepProp.push_back(the_prop);
    }

    return goodLepProp;
}
""")
                          
ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
using ROOT::RVecI;
RVecF LepFromZ(RVecF &El_prop, RVecI &El_charge, RVecF &Mu_prop, RVecI &Mu_charge, RVecI &Z_flav, ROOT::RVec<Short_t> &Z_l1Idx, ROOT::RVec<Short_t> &Z_l2Idx, int charge){
    RVecF goodLepProp;
    
    short the_l1_idx = Z_l1Idx.at(0);
    short the_l2_idx = Z_l2Idx.at(0);

    int the_flav = Z_flav.at(0);

    if ( the_flav == -121 ){
        goodLepProp = getLepProp(El_prop, El_charge, the_l1_idx, the_l2_idx, charge);
    }
    else if ( the_flav == -169 ){
        goodLepProp = getLepProp(Mu_prop, Mu_charge, the_l1_idx, the_l2_idx, charge);
    }

    return goodLepProp;
}
""")


