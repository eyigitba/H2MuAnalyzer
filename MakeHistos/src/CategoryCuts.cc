
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"


bool InCategory( const ObjectSelectionConfig & cfg, const NTupleBranches & br,
		 const std::string sel, const bool verbose) {
  
  if (verbose) std::cout << "\nInside InCategory, checking if event passes " << sel << std::endl;
  
  bool PASS = false;

  if (sel == "NONE") PASS = true;
  
  else if (sel == "Mu1Pos") {
    if (verbose) std::cout << "  * Applying Mu1Pos cuts" << std::endl;

    for (int i = 0; i < br.nMuPairs; i++) {
      int iMu1 = br.muPairs->at(i).iMu1;

      if (br.muons->at(iMu1).charge == 1) PASS = true;
      if (verbose) std::cout << "    - Mu1 charge = " << br.muons->at(iMu1).charge << ", PASS = " << PASS << std::endl;
      break; // Only look at the first muon pair
    } // End loop for (int i = 0; i < nMuPairs; i++)
  } // End if (sel == "Mu1Pos")
  
  else if (sel == "Mu1Neg") {
    if (verbose) std::cout << "  * Applying Mu1Neg cuts" << std::endl;
    
    for (int i = 0; i < br.nMuPairs; i++) {
      int iMu1 = br.muPairs->at(i).iMu1;
      
      if (br.muons->at(iMu1).charge == -1) PASS = true;
      if (verbose) std::cout << "    - Mu1 charge = " << br.muons->at(iMu1).charge << ", PASS = " << PASS << std::endl;
      break; // Only look at the first muon pair
    } // End loop for (int i = 0; i < nMuPairs; i++)
  } // End if (sel == "Mu1Neg")

  // else if (sel == "WHlep") {  // real toy category, needs more study on cuts
  //   if (verbose) std::cout << "  * Applying WHlep cuts" << std::endl;
    
  //   if (br.nMuons == 3) {
  // 	PASS = true;
  // 	for (int i = 0; i < br.nMuons; i++) {
  // 	  if ( !MuonPass(br.muons->at(i),20) ) PASS = false;
  // 	}
  //   }
  //   if (br.nMuons == 2 && br.nEles == 1) {
  // 	if (br.eles->at(0).pt > 20) PASS = true;
  //   }
  // } // End if (sel == "WHlep")

  else if (sel == "ZHmu") {  // real toy category, needs more study on cuts
    if (verbose) std::cout << "  * Applying ZHmu cuts" << std::endl;

    if (br.nMuons == 4) {
	for (int i = 0; i < br.nMuPairs; i ++) {
	  if (br.muPairs->at(i).iMu1 != 0 && br.muPairs->at(i).iMu1 != 1 && br.muPairs->at(i).iMu2 != 0 && br.muPairs->at(i).iMu2 != 1) {
	    if (br.muPairs->at(i).mass > 80 && br.muPairs->at(i).mass < 100) PASS = true;
	  }
	} 
    }
  } // End if (sel == "ZHmu")

  else if (sel == "ZHele") {   // real toy category, needs more study on cuts
    if (verbose) std::cout << "  * Applying ZHele cuts" << std::endl;

    if (br.nMuons == 2 && br.nEles == 2) {
	if (br.eles->at(0).pt > 20 && br.eles->at(1).pt > 20) PASS = true;
    }
  } // End if (sel == "ZHele")


  else if (sel == "WH_3l_mu") {  // WH --> 3 muon category
    if (verbose) std::cout << "  * Applying WH_3l_mu cuts" << std::endl;

    // Require MET > 30, no medium b-tags, exactly 2 selected opposite-charge muon pairs, and 1 pair close to the Higgs mass
    if (br.met->pt < 30) return false;
    if (SelectedJets(cfg, br, "BTagMedium").size() > 0) return false;

    // EleInfos eles = SelectedEles(cfg, br);
    MuPairInfos muPairs = SelectedMuPairs(cfg, br);

    // if (eles.size()    != 0) return false; // Need to add in requirement of 0 selected electrons - AWB 03.10.2018
    if (muPairs.size() != 2) return false; // Exactly 2 opposite-charge muon pairs in event

    bool H_mass = false;
    for (const auto & pair : muPairs) {
      if ( MuPairMass(pair, cfg.mu_pt_corr) > 105 &&
	   MuPairMass(pair, cfg.mu_pt_corr) < 160 ) H_mass = true;
    }
    if (not H_mass) return false;

    PASS = true;
  } // End if (sel == "WH_3l_mu")


  else if (sel == "WZ_3l_val_mu") {  // WZ --> 3 muon validation category
    if (verbose) std::cout << "  * Applying WZ_3l_val_mu cuts" << std::endl;

    // Require MET > 40, no medium b-tags, exactly 2 selected opposite-charge muon pairs, and 1 pair close to the Z mass
    if (br.met->pt < 40) return false;
    if (SelectedJets(cfg, br, "BTagMedium").size() > 0) return false;

    // EleInfos eles = SelectedEles(cfg, br);
    MuPairInfos muPairs = SelectedMuPairs(cfg, br);

    // if (eles.size()    != 0) return false; // Need to add in requirement of 0 selected electrons - AWB 03.10.2018
    if (muPairs.size() != 2) return false; // Exactly 2 opposite-charge muon pairs in event

    bool Z_mass = false;
    for (const auto & pair : muPairs) {
      if ( abs(MuPairMass(pair, cfg.mu_pt_corr) - 91) < 10 ) Z_mass = true;
    }
    if (not Z_mass) return false;

    PASS = true;
  } // End if (sel == "WZ_3l_val_mu")


  else if (sel == "ttZ_3l_val_mu") {  // ttZ --> 3 muon validation category
    if (verbose) std::cout << "  * Applying ttZ_3l_val_mu cuts" << std::endl;

    // Require MET > 40, >= 3 jets, >= 1 medium or >= 1 tight b-tags, exactly 2 selected opposite-charge muon pairs, and 1 pair close to the Z mass
    if (br.met->pt < 40) return false;
    if (SelectedJets(cfg, br).size() < 3) return false;
    if (SelectedJets(cfg, br, "BTagMedium").size() < 2 &&
	SelectedJets(cfg, br, "BTagTight").size() < 1) return false;

    // EleInfos eles = SelectedEles(cfg, br);
    MuPairInfos muPairs = SelectedMuPairs(cfg, br);

    // if (eles.size()    != 0) return false; // Need to add in requirement of 0 selected electrons - AWB 03.10.2018
    if (muPairs.size() != 2) return false; // Exactly 2 opposite-charge muon pairs in event

    bool Z_mass = false;
    for (const auto & pair : muPairs) {
      if ( abs(MuPairMass(pair, cfg.mu_pt_corr) - 91) < 10 ) Z_mass = true;
    }
    if (not Z_mass) return false;

    PASS = true;
  } // End if (sel == "ttZ_3l_val_mu")


  else if (sel == "ttW_3l_val_mu") {  // ttW --> 3 muon validation category
    if (verbose) std::cout << "  * Applying ttW_3l_val_mu cuts" << std::endl;

    // Require MET > 50, >= 3 jets, >= 1 medium or >= 1 tight b-tags, exactly 2 selected opposite-charge muon pairs, and 0 pairs close to the Z mass
    if (br.met->pt < 50) return false;
    if (SelectedJets(cfg, br).size() < 3) return false;
    if (SelectedJets(cfg, br, "BTagMedium").size() < 2 &&
	SelectedJets(cfg, br, "BTagTight").size() < 1) return false;

    // EleInfos eles = SelectedEles(cfg, br);
    MuPairInfos muPairs = SelectedMuPairs(cfg, br);

    // if (eles.size()    != 0) return false; // Need to add in requirement of 0 selected electrons - AWB 03.10.2018
    if (muPairs.size() != 2) return false; // Exactly 2 opposite-charge muon pairs in event

    bool Z_mass = false;
    for (const auto & pair : muPairs) {
      if ( abs(MuPairMass(pair, cfg.mu_pt_corr) - 91) < 10 ) Z_mass = true;
    }
    if (Z_mass) return false;

    PASS = true;
  } // End if (sel == "ttW_3l_val_mu")


  else if (sel == "Z_3l_val_mu") {  // Z --> 3 muon validation category (non-prompt)
    if (verbose) std::cout << "  * Applying Z_3l_val_mu cuts" << std::endl;

    // Require MET < 30, no medium b-tags, exactly 2 selected opposite-charge muon pairs, and 1 pair close to the Z mass
    if (br.met->pt > 30) return false;
    if (SelectedJets(cfg, br, "BTagMedium").size() > 0) return false;

    // EleInfos eles = SelectedEles(cfg, br);
    MuPairInfos muPairs = SelectedMuPairs(cfg, br);

    // if (eles.size()    != 0) return false; // Need to add in requirement of 0 selected electrons - AWB 03.10.2018
    if (muPairs.size() != 2) return false; // Exactly 2 opposite-charge muon pairs in event

    bool Z_mass = false;
    for (const auto & pair : muPairs) {
      if ( abs(MuPairMass(pair, cfg.mu_pt_corr) - 91) < 10 ) Z_mass = true;
    }
    if (not Z_mass) return false;

    PASS = true;
  } // End if (sel == "Z_3l_val_mu")


  else if (sel == "ttbar_3l_val_mu") {  // ttbar --> 3 muon validation category (non-prompt)
    if (verbose) std::cout << "  * Applying ttbar_3l_val_mu cuts" << std::endl;

    // Require MET > 40, exactly 1 medium b-tag, exactly 2 selected opposite-charge muon pairs, and 0 pairs close to the Z mass
    if (br.met->pt < 40) return false;
    if (SelectedJets(cfg, br, "BTagMedium").size() != 1) return false;

    // EleInfos eles = SelectedEles(cfg, br);
    MuPairInfos muPairs = SelectedMuPairs(cfg, br);

    // if (eles.size()    != 0) return false; // Need to add in requirement of 0 selected electrons - AWB 03.10.2018
    if (muPairs.size() != 2) return false; // Exactly 2 opposite-charge muon pairs in event

    bool Z_mass = false;
    for (const auto & pair : muPairs) {
      if ( abs(MuPairMass(pair, cfg.mu_pt_corr) - 91) < 10 ) Z_mass = true;
    }
    if (Z_mass) return false;

    PASS = true;
  } // End if (sel == "ttbar_3l_val_mu")


  else {
    std::cout << "\nInside CategoryCuts.cc, don't recognize category " << sel << std::endl;
    assert(false);
  }

  return PASS;

}
