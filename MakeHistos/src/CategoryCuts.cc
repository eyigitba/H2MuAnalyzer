
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"
#include "TLorentzVector.h"


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

//  else if (sel.compare("WHlepM") == 0) {  // real toy category, needs more study on cuts
//    if (verbose) std::cout << "  * Applying WHlepM cuts" << std::endl;
//    
//    if (br.nBTight == 0 and (br.met)->pt > 20) {
//      if (br.nMuons == 3) {
//	pass = true;
//	for (int i = 0; i < br.nMuons; i++) {
//	  if ( !MuonPass(br.muons->at(i),20) or abs(br.muons->at(i).d0_PV)>0.02 ) pass = false;
// 	}
//      }//end if (br.nMuons == 3)
//    }//end if (br.nBTight == 0 and (br.met)->pt > 30)
//  }//end if (sel.compare("WHlep") == 0)
//
//  else if (sel.compare("WHlepE") == 0) {
//    if (verbose) std::cout << "  * Applying WHlepE cuts" << std::endl;
//    if (br.nBTight == 0 and (br.met)->pt > 20) {
//      if (br.nMuons == 2 and br.nEles == 1) {
//        if (br.eles->at(0).pt > 20 and abs(br.eles->at(0).d0_PV)<0.02 ) pass = true;
//      }
//    }//end if (br.nBTight == 0 and (br.met)->pt > 20)
//  }//end if (sel.compare("WHlepE") == 0)
//
//
//  else if (sel.compare("WHhad2J") == 0) {
//    if (verbose) std::cout << "  * Applying WHhad2J cuts" << std::endl;
//    
//    if (not InCategory(br, "WHlep", verbose) and br.nMuons == 2 and br.nBTight == 0 ) {
//	for (int i = 0; i < br.nJetPairs; i++) {
//	  JetPairInfo & dijet = br.jetPairs->at(i);
//	  if ( DijetPass(br, dijet) and abs(dijet.mass-80)<20 and abs(dijet.eta)<4 and dijet.dR<3 ) pass = true;
//	}//end for (int i = 0; i < br.nJetPairs; i++)
//    }
//  }//end if (sel.compare("WHhad2J") == 0) 
//
//  else if (sel.compare("WHhad1J") == 0) {
//    if (verbose) std::cout << "  * Applying WHhad1J cuts" << std::endl;
//
//    if (not InCategory(br, "WHlep", verbose) and br.nMuons == 2 and br.nBTight == 0 ) {
//      if (not InCategory(br, "WHhad2J", verbose) ) {
//	for (int i = 0; i < br.nJets; i++) {
//	   JetInfo & jet = br.jets->at(i);
//	   if ( JetPass(jet) and abs(jet.eta)<2 and not BJetPass(jet) ) pass = true;
//	}// end for (int i = 0; i < br.nJets; i++) {
//      }
//    }
//  }//end if (sel.compare("WHhad1J") == 0)

//  else if (sel == "ZHmu") {  // real toy category, needs more study on cuts
//    if (verbose) std::cout << "  * Applying ZHmu cuts" << std::endl;
//
//    if (br.nMuons == 4) {
//	for (int i = 0; i < br.nMuPairs; i ++) {
//	  if (br.muPairs->at(i).iMu1 != 0 && br.muPairs->at(i).iMu1 != 1 && br.muPairs->at(i).iMu2 != 0 && br.muPairs->at(i).iMu2 != 1) {
//	    if (br.muPairs->at(i).mass > 80 && br.muPairs->at(i).mass < 100) PASS = true;
//	  }
//	} 
//    }
//  } // End if (sel == "ZHmu")
//
//  else if (sel == "ZHele") {   // real toy category, needs more study on cuts
//    if (verbose) std::cout << "  * Applying ZHele cuts" << std::endl;
//
//    if (br.nMuons == 2 and br.nEles ==2) {
//      EleInfo & ele1 = br.eles->at(0);
//      EleInfo & ele2 = br.eles->at(1);
//      if(ele1.pt > 20 and ele2.pt > 20 and ele1.charge + ele2.charge == 0) {
//    	TLorentzVector v1, v2, vz;
//    	v1.SetPtEtaPhiM( ele1.pt, ele1.eta, ele1.phi, 0.0005);
//    	v2.SetPtEtaPhiM( ele2.pt, ele2.eta, ele2.phi, 0.0005);
//    	vz = v1 + v2;
//        std::cout << vz.M() << std::endl;
//    	if( vz.M()>80 and vz.M()<100 ) pass = true;
//      }
//    }
//  } // end if (sel.compare("ZHele") == 0)

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
