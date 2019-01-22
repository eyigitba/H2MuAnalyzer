#include "H2MuAnalyzer/MakeHistos/interface/SampleID.h"


// assign an ID to the sample based on its parent particle PDG ID   -- XWZ 20.11.2018
// generic notations: 00 for quark, 10 for lepton, 20 for boson
int getSampleID(TString name, const bool verbose) {
    if (verbose) std::cout << "\nInside getSampleID for sample " << name << std::endl;
    // 0 for data
    if ( name.Contains("SingleMu") ) 	return 0;
    // positive value for signal
    else if ( name.Contains("H2Mu_ggH") ) 	return 25;
    else if ( name.Contains("H2Mu_VBF") ) 	return 010225;  // represented by one up and one down quark
    else if ( name.Contains("H2Mu_WH") ) 	return 2425;
    else if ( name.Contains("H2Mu_ZH") ) 	return 2325;
    else if ( name.Contains("H2Mu_ttH") ) 	return 060625;
    // negative value for bkg
    else if ( name.Contains("ZJets") ) 		return -23; // can add more entries if need to distinguish between AMC/MG, 0j/1j/2j samples
    else if ( name.Contains("tt_ll") or name == "tt" ) 	return -0606;
    else if ( name.Contains("tW_pos") or name.Contains("tW_neg") ) return -0624;
    else if ( name == "tZq" ) 			return -062300; // 00 for quark in general
    else if ( name == "tZW" ) 			return -062324;
    else if ( name.Contains("ttW") ) 		return -060624;
    else if ( name == "ttZ" ) 			return -060623;
    else if ( name == "ttH" ) 			return -060625;
    else if ( name.Contains("WW") ) 		return -2424;
    else if ( name.Contains("WZ") ) 		return -2423;
    else if ( name.Contains("ZZ") ) 		return -2323;
    else if ( name.Contains("WWW") ) 		return -242424;
    else if ( name.Contains("WWZ") ) 		return -242423;
    else if ( name.Contains("WZZ") ) 		return -242323;
    else if ( name.Contains("ZZZ") ) 		return -232323;

    return -999;  //just so that this function can be compiled. do not expect -999 to show in any case
}




