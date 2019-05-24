#include "H2MuAnalyzer/MakeHistos/interface/SampleID.h"


// assign an ID to the sample based on its parent particle PDG ID   -- XWZ 20.11.2018
// generic notations: 00 for quark, 10 for lepton, 20 for boson
int getSampleID(TString name, const bool verbose) {
    if (verbose) std::cout << "\nInside getSampleID for sample " << name << std::endl;
    // 0 for data
    if ( name.Contains("SingleMu") ) 	return 0;
    // positive value for signal
    else if ( name.Contains("H2Mu_ggH_120") )       return 120000025;
    else if ( name.Contains("H2Mu_VBF_120") )       return 120010225;  // represented by one up and one down quark
    else if ( name.Contains("H2Mu_WH_pos_120") )    return 1202425;
    else if ( name.Contains("H2Mu_WH_neg_120") )    return 1202425;
    else if ( name.Contains("H2Mu_ZH_120") )        return 1202325;
    else if ( name.Contains("H2Mu_ttH_120") )       return 120060625;
  
    else if ( name.Contains("H2Mu_ggH_130") )       return 130000025;
    else if ( name.Contains("H2Mu_VBF_130") )       return 130010225;  // represented by one up and one down quark
    else if ( name.Contains("H2Mu_WH_pos_130") )    return 1302425;
    else if ( name.Contains("H2Mu_WH_neg_130") )    return 1302425;
    else if ( name.Contains("H2Mu_ZH_130") )        return 1302325;
    else if ( name.Contains("H2Mu_ttH_130") )       return 130060625;

    else if ( name.Contains("H2Mu_ggH") ) 	return 25;
    else if ( name.Contains("H2Mu_VBF") ) 	return 10225;  // represented by one up and one down quark
    else if ( name.Contains("H2Mu_WH") ) 	return 2425;
    else if ( name.Contains("H2Mu_ZH") ) 	return 2325;
    else if ( name.Contains("H2Mu_ttH") ) 	return 60625;
    // negative value for bkg
    else if ( name.Contains("ZJets") ) 		return -23; // can add more entries if need to distinguish between AMC/MG, 0j/1j/2j samples
    else if ( name.Contains("tt_ll") or name == "tt" ) 	return -606;
    else if ( name.Contains("tW_pos") or name.Contains("tW_neg") ) return -624;
    else if ( name == "tHq" )			return -62500; // 00 for quark in general
    else if ( name == "tHW" )                   return -62524; 
    else if ( name == "tZq" ) 			return -62300; // 00 for quark in general
    else if ( name == "tZW" ) 			return -62324;
    else if ( name == "ttZ" ) 			return -60623;
    else if ( name == "ttH" ) 			return -60625;
    else if ( name == "ttWW")                   return -6062424;
    else if ( name.Contains("ttW") )            return -60624;

    else if ( name.Contains("WWW") )            return -242424;
    else if ( name.Contains("WWZ") )            return -242423;
    else if ( name.Contains("WZZ") )            return -242323;
    else if ( name.Contains("ZZZ") )            return -232323;
    else if ( name.Contains("WW") ) 		return -2424;
    else if ( name.Contains("WZ") ) 		return -2423;
    else if ( name.Contains("ZZ_4l_gg"))	return -23230000;
    else if ( name.Contains("ZZ") ) 		return -2323;
    return -999;  //just so that this function can be compiled. do not expect -999 to show in any case
}




