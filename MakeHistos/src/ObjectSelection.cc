#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"

// Configure constants related to object selection
void ConfigureObjectSelection( ObjectSelectionConfig & cfg, const std::string _year, const std::string _opt ) {

  if (_year == "Legacy2016") {
    cfg.year = _year;

    // Muon selection
    cfg.mu_pt_corr = "KaMu";   // Muon pT correction: "PF", "Roch", or "KaMu"
    cfg.mu_pt_min  = 20.0;     // Minimum muon pT
    cfg.mu_eta_max =  2.4;     // Maximum muon |eta|
    cfg.mu_ID_cut  = "medium"; // Muon ID: "loose", "medium", or "tight"
    cfg.mu_iso_max = 0.25;     // Maximum muon relative isolation

    // Electron selection
    cfg.ele_pt_min  = 10.0;     // Minimum electron pT
    cfg.ele_eta_max =  2.5;     // Maximum electron |eta|
    cfg.ele_ID_cut  = "medium"; // Electron ID: "loose", "medium", or "tight"
    cfg.ele_iso_max = 0.25;     // Maximum electron relative isolation

    // Jet selection
    cfg.jet_pt_min     = 30.0;   // Minimum jet pT
    cfg.jet_eta_max    =  4.7;   // Maximum jet |eta|
    cfg.jet_mu_dR_min  =  0.4;   // Minimum dR(jet, muon)
    cfg.jet_ele_dR_min =  0.4;   // Minimum dR(jet, electron)
    cfg.jet_btag_cuts  = {0.5426, 0.8484, 0.9535}; // CSVv2 recommendation from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    cfg.jet_PU_ID_cut  = "NONE"; // Jet passes PU ID cut

    // Higgs candidate selection
    cfg.muPair_Higgs = "sort_OS_sum_muon_pt";
  } // End if (_year == "Legacy2016")

  else if (_year == "2016" || _year == "2017" || _year == "2018") {
    cfg.year = _year;

    // Muon selection
    if (_year == "2016") cfg.mu_pt_corr  = "KaMu";   // Muon pT correction: "PF", "Roch", or "KaMu"
    else 		 cfg.mu_pt_corr  = "Roch";   // Muon pT correction: "PF", "Roch", or "KaMu"
    cfg.mu_pt_min   = 20.0;     // Minimum muon pT
    cfg.mu_eta_max  =  2.4;     // Maximum muon |eta|
    cfg.mu_ID_cut   = "medium"; // Muon ID: "loose", "medium", or "tight"
    cfg.mu_mIso_max =  0.4;     // Maximum muon relative miniIsolation
    cfg.mu_d0_max   =  0.05;    // Maximum muon |dXY| from vertex
    cfg.mu_dZ_max   =  0.1;     // Maximum muon |dZ| from vertex
    cfg.mu_SIP_max  =  8.0;     // Maximum muon impact parameter significance
    cfg.mu_seg_min  =  -99;     // Minimum muon segment compatibility
    cfg.mu_MVA_min  =  -99;     // Minimum prompt muon lepton MVA (lepMVA) BDT score
    cfg.mu_CSV_max  = "NONE";   // Veto muons with pT < 20 GeV with matching jet passing b-tag threshold

    // Electron selection
    cfg.ele_pt_min  = 10.0;     // Minimum electron pT
    cfg.ele_eta_max =  2.5;     // Maximum electron |eta|
    cfg.ele_ID_cut  = "loose";  // Electron ID: "loose", "medium", or "tight"
    cfg.ele_mIso_max =  0.4;    // Maximum electron relative miniIsolation
    cfg.ele_d0_max  =   0.05;   // Maximum electron |dXY| from vertex
    cfg.ele_dZ_max  =   0.1;    // Maximum electron |dZ| from vertex
    cfg.ele_SIP_max =   8.0;    // Maximum electron impact parameter significance
    cfg.ele_MVA_min =   -99;    // Minimum prompt electron lepton MVA (lepMVA) BDT score
    cfg.ele_CSV_max = "NONE";   // Veto electrons with pT < 20 GeV with matching jet passing b-tag threshold

    // Jet selection
    cfg.jet_pt_min     = 30.0;    // Minimum jet pT
    cfg.jet_eta_max    =  4.7;    // Maximum jet |eta|
    cfg.jet_mu_dR_min  =  0.4;    // Minimum dR(jet, muon)
    cfg.jet_ele_dR_min =  0.4;    // Minimum dR(jet, electron)
    cfg.jet_PU_ID_cut  = "loose"; // Jet passes PU ID cut
    if (_year == "2016") cfg.jet_btag_cuts  = {0.2217, 0.6321, 0.8953}; // DeepCSV recommendation from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
    if (_year == "2017") cfg.jet_btag_cuts  = {0.1522, 0.4941, 0.8001}; // DeepCSV recommendation from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    if (_year == "2018") cfg.jet_btag_cuts  = {0.1241, 0.4184, 0.7527}; // DeepCSV recommendation from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X

    // Higgs candidate selection
    cfg.muPair_Higgs = "sort_OS_sum_muon_pt";

    if (_opt == "lepMVA") {       // LepMVA pre-selection from TOP-18-008, for 3-lepton and 4-lepton channels
      cfg.mu_pt_min   = 10.0;     // Lower minimum muon pT for higher acceptance
      cfg.mu_seg_min  = 0.30;     // Minimum muon segment compatibility
      cfg.mu_MVA_min  = -0.4;     // Minimum prompt muon lepton MVA (lepMVA) BDT score
      cfg.mu_CSV_max  = "loose";  // Veto muons with pT < 20 GeV with matching jet passing b-tag threshold

      cfg.ele_pt_min  = 10.0;     // Minimum electron pT
      cfg.ele_MVA_min = -0.4;     // Minimum prompt electron lepton MVA (lepMVA) BDT score
      cfg.ele_CSV_max = "loose";  // Veto electrons with pT < 20 GeV with matching jet passing b-tag threshold

      cfg.jet_pt_min  = 20.0;     // Lower minimum jet pT for higher acceptance
    }

  } // End if (_year == "2016" || _year == "2017" || _year == "2018")

  else {
    std::cout << "Inside ConfigureEventWeight, invalid year = " << _year << std::endl;
    assert(false);
  }

} // End function: ConfigureObjectSelection()


// Select muons passing ID and kinematic cuts
bool MuonPass ( const ObjectSelectionConfig & cfg, const MuonInfo & muon, const NTupleBranches & br, const bool verbose ) {

  if ( cfg.mu_pt_min   != -99 )
    if ( MuonPt(muon, cfg.mu_pt_corr) < cfg.mu_pt_min   ) return false;
  if ( cfg.mu_eta_max  != -99 )
    if ( fabs(muon.eta)               > cfg.mu_eta_max  ) return false;
  if ( cfg.mu_ID_cut   != "NONE" )
    if ( MuonID(muon, cfg.mu_ID_cut) != true            ) return false;
  if ( cfg.mu_iso_max  != -99 )
    if ( muon.relIso                  > cfg.mu_iso_max  ) return false;
  if ( cfg.mu_mIso_max != -99 )
    if ( muon.miniIso                 > cfg.mu_mIso_max ) return false;
  if ( cfg.mu_d0_max   != -99 )
    if ( fabs(muon.d0_PV)             > cfg.mu_d0_max   ) return false;
  if ( cfg.mu_dZ_max   != -99 )
    if ( fabs(muon.dz_PV)             > cfg.mu_dZ_max   ) return false;
  if ( cfg.mu_SIP_max  != -99 )
    if ( muon.SIP_3D                  > cfg.mu_SIP_max  ) return false;
  if ( cfg.mu_seg_min  != -99 )
    if ( muon.segCompat               < cfg.mu_seg_min  ) return false;
  if ( cfg.mu_MVA_min  != -99 )
    if ( muon.lepMVA                  < cfg.mu_MVA_min  ) return false;
  if ( cfg.mu_CSV_max  != "NONE" ) {
    // Drop low-pT muons which are too close to b-tagged jets
    for ( const auto & jet : (*br.jets) ) {
      if ( JetMuonClean(cfg, jet, muon, br) < 2 ) return false;
    }
  }

  return true;
} // End function: bool MuonPass()


// Select electrons passing ID and kinematic cuts
bool ElePass ( const ObjectSelectionConfig & cfg, const EleInfo & ele, const NTupleBranches & br, const bool verbose ) {

  if ( cfg.ele_pt_min  != -99 )
    if ( ele.pt                      < cfg.ele_pt_min   ) return false;
  if ( cfg.ele_eta_max != -99 )
    if ( fabs(ele.eta)               > cfg.ele_eta_max  ) return false;
  if ( cfg.ele_ID_cut  != "NONE" )
    if ( EleID(ele, cfg.ele_ID_cut) != true             ) return false;
  if ( cfg.ele_iso_max != -99 )
    if ( ele.relIso                  > cfg.ele_iso_max  ) return false;
  if ( cfg.ele_mIso_max != -99 )
    if ( ele.miniIso                 > cfg.ele_mIso_max ) return false;
  if ( cfg.ele_d0_max   != -99 )
    if ( fabs(ele.d0_PV)             > cfg.ele_d0_max   ) return false;
  if ( cfg.ele_dZ_max   != -99 )
    if ( fabs(ele.dz_PV)             > cfg.ele_dZ_max   ) return false;
  if ( cfg.ele_SIP_max != -99 )
    if ( ele.SIP_3D                  > cfg.ele_SIP_max  ) return false;
  if ( cfg.ele_MVA_min  != -99 )
    if ( ele.lepMVA                  < cfg.ele_MVA_min  ) return false;
  if ( cfg.ele_CSV_max  != "NONE" ) {
    // Drop low-pT electrons which are too close to b-tagged jets
    for ( const auto & jet : (*br.jets) ) {
      if ( JetEleClean(cfg, jet, ele, br) < 2 ) return false;
    }
  }

  return true;
} // End function: bool ElectronPass()


// Select jets passing ID and kinematic cuts
bool JetPass( const ObjectSelectionConfig & cfg, const JetInfo & jet, const NTupleBranches & br, const std::string sel, const bool verbose ) {

  if (verbose) std::cout << "  * Inside JetPass, selection = " << sel << std::endl;
  if (verbose) std::cout << "  * Jet pT = " << jet.pt << ", eta = " << jet.eta << ", phi = " << jet.phi
                         << ", puID = " << jet.puID << ", CSV = " << JetCSV(jet, "CSVv2") << ", deepCSV = " << JetCSV(jet) << std::endl;

  float jet_CSV = JetCSV(jet, (cfg.year == "Legacy2016" ? "CSVv2" : "deepCSV"));

  if ( jet.pt         < cfg.jet_pt_min     )         return false;
  if ( fabs(jet.eta)  > cfg.jet_eta_max    )         return false;
  if ( !JetPUID(jet, cfg.jet_PU_ID_cut, cfg.year ) ) return false;

  if (sel == "Central"    &&  abs(jet.eta) >  2.4) return false;
  if (sel == "Forward"    &&  abs(jet.eta) <= 2.4) return false;
  if (sel == "BTagLoose"  && (abs(jet.eta) > 2.4 || jet_CSV < cfg.jet_btag_cuts.at(0))) return false;
  if (sel == "BTagMedium" && (abs(jet.eta) > 2.4 || jet_CSV < cfg.jet_btag_cuts.at(1))) return false;
  if (sel == "BTagTight"  && (abs(jet.eta) > 2.4 || jet_CSV < cfg.jet_btag_cuts.at(2))) return false;

  // Clean jet against selected muons and electrons
  if ( cfg.jet_mu_dR_min > 0 )
    for ( const auto & mu : (*br.muons) )
      if ( (JetMuonClean(cfg, jet, mu, br) % 2) != 1 ) return false;
  if ( cfg.jet_ele_dR_min > 0 )
    for ( const auto & ele : (*br.eles) )
      if ( (JetEleClean(cfg, jet, ele, br) % 2) != 1 ) return false;

  if (verbose) std::cout << "    - PASSED selection!!!" << std::endl;

  return true;
} // End function: bool JetPass()


// Decide whether jet or lepton passes jet-lepton cleaning algorithm
// Returns a bit mask: 0 for both fail, 1 for jet passes, 2 for lepton passes, 3 for both pass
int JetMuonClean( const ObjectSelectionConfig & cfg, const JetInfo & jet, const MuonInfo & mu, const NTupleBranches & br, const bool verbose ) {

  TLorentzVector vLep = FourVec(mu, cfg.mu_pt_corr);
  TLorentzVector vJet = FourVec(jet);

  // Don't clean if separation is large enough
  if (vLep.DeltaR(vJet) > cfg.jet_mu_dR_min) return 3;

  // Temporarily modify config to ignore jet-lepton cleaning
  ObjectSelectionConfig cfg_mod = cfg;
  cfg_mod.mu_CSV_max = "NONE";  // Don't drop low-pT leptons too close to b-tagged jets
  cfg_mod.jet_mu_dR_min = -99;  // Don't drop jets too close to leptons

  // No cleaning to do if either the lepton or the jet fail the basic cuts
  bool lepPass = MuonPass(cfg_mod,  mu, br);
  bool jetPass = JetPass (cfg_mod, jet, br);
  if (!jetPass || !lepPass) return (1*jetPass + 2*lepPass);

  if (cfg.mu_CSV_max == "NONE") {  // Always remove jet if it is too close to the lepton
    jetPass = false;
  } else {                         // Remove leptons with pT < 20 matched to b-tagged jets
    assert(cfg.mu_CSV_max == "loose");
    float jet_CSV = JetCSV(jet, (cfg.year == "Legacy2016" ? "CSVv2" : "deepCSV"));  // Will switch to deepCSV in all years eventually - AWB 2019.01.19
    if (vLep.Pt() > 20 || jet_CSV < cfg.jet_btag_cuts.at(0)) {
      jetPass = false;
    } else {
      lepPass = false;
      std::cout << "Jet with pT " << jet.pt << ", eta " << jet.eta << ", phi " << jet.phi << ", CSV " << jet_CSV
		<< " has dR(mu) " << vLep.DeltaR(vJet)  << " where lepMVA(mu) = " << mu.lepMVA
		<< " and muon pT = " << vLep.Pt() << ", eta = " << mu.eta << ", phi = " << mu.phi
		<< ", matched jet-with-muon pT = " << mu.pt/mu.jet_ptRatio << ", CSV = " << mu.jet_deepCSV << std::endl;
    }
  }

  return (1*jetPass + 2*lepPass);
} // End function: int JetMuonClean()

// Decide whether jet or lepton passes jet-lepton cleaning algorithm
// Returns a bit mask: 0 for both fail, 1 for jet passes, 2 for lepton passes, 3 for both pass
int JetEleClean( const ObjectSelectionConfig & cfg, const JetInfo & jet, const EleInfo & ele, const NTupleBranches & br, const bool verbose ) {

  TLorentzVector vLep = FourVec(ele);
  TLorentzVector vJet = FourVec(jet);

  // Don't clean if separation is large enough
  if (vLep.DeltaR(vJet) > cfg.jet_ele_dR_min) return 3;

  // Temporarily modify config to ignore jet-lepton cleaning
  ObjectSelectionConfig cfg_mod = cfg;
  cfg_mod.ele_CSV_max = "NONE";  // Don't drop low-pT leptons too close to b-tagged jets
  cfg_mod.jet_ele_dR_min = -99;  // Don't drop jets too close to leptons

  // No cleaning to do if either the lepton or the jet fail the basic cuts
  bool lepPass = ElePass(cfg_mod, ele, br);
  bool jetPass = JetPass(cfg_mod, jet, br);
  if (!jetPass || !lepPass) return (1*jetPass + 2*lepPass);

  if (cfg.ele_CSV_max == "NONE") {  // Always remove jet if it is too close to the lepton
    jetPass = false;
  } else {                         // Remove leptons with pT < 20 matched to b-tagged jets
    assert(cfg.ele_CSV_max == "loose");
    float jet_CSV = JetCSV(jet, (cfg.year == "Legacy2016" ? "CSVv2" : "deepCSV"));  // Will switch to deepCSV in all years eventually - AWB 2019.01.19
    if (vLep.Pt() > 20 || jet_CSV < cfg.jet_btag_cuts.at(0)) {
      jetPass = false;
    } else {
      lepPass = false;
      if (verbose) std::cout << "Jet with pT " << jet.pt << ", eta " << jet.eta << ", phi " << jet.phi << ", CSV " << jet_CSV
			     << " has dR(ele) " << vLep.DeltaR(vJet)  << " where lepMVA(ele) = " << ele.lepMVA
			     << " and electron pT = " << vLep.Pt() << ", eta = " << ele.eta << ", phi = " << ele.phi
			     << ", matched jet-with-ele pT = " << ele.pt/ele.jet_ptRatio << ", CSV = " << ele.jet_deepCSV << std::endl;
    }
  }

  return (1*jetPass + 2*lepPass);
} // End function: int JetEleClean()


// Return selected muons in event
MuonInfos SelectedMuons ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const bool verbose ) {

  MuonInfos selMuons;
  for (const auto & muon : (*br.muons)) {
    if ( MuonPass(cfg, muon, br) ) {
	selMuons.push_back(muon);
    }
  }
  return selMuons;
}


// Return selected muons in event
MuPairInfos SelectedMuPairs ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel, const bool verbose ) {

  assert(sel == "OS" || sel == "SS" || sel.length() == 0);

  MuPairInfos selMuPairs;
  for (const auto & muPair : (*br.muPairs)) {
    if ( MuonPass(cfg, br.muons->at(muPair.iMu1), br) &&
	 MuonPass(cfg, br.muons->at(muPair.iMu2), br) ) {
      if (sel == "OS" && muPair.charge != 0) continue;
      if (sel == "SS" && muPair.charge == 0) continue;
      selMuPairs.push_back(muPair);
    }
  }
  return selMuPairs;
}


// Return selected dimuon Higgs candidate pair
MuPairInfo SelectedCandPair ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const bool verbose ) {

  MuPairInfo candPair;

  if (cfg.muPair_Higgs == "sort_OS_sum_muon_pt") {
    float sum_mu_pt = -99;
    for (const auto & muPair : SelectedMuPairs(cfg, br)) {
      if ( MuonPt(br.muons->at(muPair.iMu1), cfg.mu_pt_corr) +
	   MuonPt(br.muons->at(muPair.iMu2), cfg.mu_pt_corr) <= sum_mu_pt) continue;
      sum_mu_pt = MuonPt(br.muons->at(muPair.iMu1), cfg.mu_pt_corr) + MuonPt(br.muons->at(muPair.iMu2), cfg.mu_pt_corr);
      candPair  = muPair;
    }
  } // End if (cfg.muPair_Higgs == "sort_OS_sum_muon_pt")
  
  else if (cfg.muPair_Higgs == "sort_OS_dimuon_pt") {
    float dimu_pt = -99;
    for (const auto & muPair : SelectedMuPairs(cfg, br)) {
      if ( MuPairPt(muPair, cfg.mu_pt_corr) <= dimu_pt ) continue;
      dimu_pt  = MuPairPt(muPair, cfg.mu_pt_corr);
      candPair = muPair;
    }
  } // End if (cfg.muPair_Higgs == "sort_OS_dimuon_pt")

  else if (cfg.muPair_Higgs == "sort_OS_dimuon_mass") {
    float dimu_mass = -99;
    for (const auto & muPair : SelectedMuPairs(cfg, br)) {
      if ( MuPairMass(muPair, cfg.mu_pt_corr) <= dimu_mass ) continue;
      dimu_mass = MuPairMass(muPair, cfg.mu_pt_corr);
      candPair = muPair;
    }
  } //End if (cfg.muPair_Higgs == "sort_OS_dimuon_mass")

  else if (cfg.muPair_Higgs == "sort_WH_3_mu_v1") {
    assert( SelectedMuPairs(cfg, br).size() == 2 ); // Only viable for events with 3 selected muons
    int iMuW = -99; // Index of muon from W
    int jMuW = -99; // Index of alternate muon from W

    // Loop over selected di-muon pairs
    for (const auto & muPair : SelectedMuPairs(cfg, br)) {
      // Require mass to fall inside signal window
      if ( MuPairMass(muPair, cfg.mu_pt_corr) < 110 ) continue;
      if ( MuPairMass(muPair, cfg.mu_pt_corr) > 150 ) continue;
      // If we have no other candidate yet, use this pair
      if ( iMuW < 0 ) {
	candPair = muPair;
	// Muon not in Higgs candidate pair presumably from W
	for (int i = 0; i < int((br.muons)->size()); i++) {
	  if ( i != candPair.iMu1 && i != candPair.iMu2 &&
	       MuonPass(cfg, br.muons->at(i), br) ) {
	    assert(iMuW < 0); // There should be only one option
	    iMuW = i;
	  }
	}
	continue;
      }
      // If we already have another candidate pair, find the W muon for this pair
      for (int j = 0; j < int((br.muons)->size()); j++) {
	if ( j != muPair.iMu1 && j != muPair.iMu2 &&
	     MuonPass(cfg, br.muons->at(j), br) ) jMuW = j;
      }
      // Expect MT(W muon, MET) < 150 GeV for higher efficiency, large smeared tail
      if ( ( FourVec(br.muons->at(iMuW), cfg.mu_pt_corr, "T") + FourVec(*br.met) ).M() > 150 &&
	   ( FourVec(br.muons->at(jMuW), cfg.mu_pt_corr, "T") + FourVec(*br.met) ).M() < 150 ) {
	candPair = muPair;
	iMuW     = jMuW;
	continue;
      }
      // If both are > 150 GeV, or both are < 150 GeV, pick the pair with the higher di-muon pT
      if ( (( FourVec(br.muons->at(iMuW), cfg.mu_pt_corr, "T") + FourVec(*br.met) ).M() > 150) ==
	   (( FourVec(br.muons->at(jMuW), cfg.mu_pt_corr, "T") + FourVec(*br.met) ).M() > 150) &&
	   MuPairPt(muPair, cfg.mu_pt_corr) > MuPairPt(candPair, cfg.mu_pt_corr) ) {
	candPair = muPair;
	iMuW     = jMuW;
	continue;
      }
    } // End loop: for (const auto & muPair : SelectedMuPairs(cfg, br))

    if (iMuW < 0) {
      candPair = SelectedMuPairs(cfg, br).at(0);
    }

  } // End if (cfg.muPair_Higgs == "sort_WH_3_mu_v1")

  else {
    std::cout << "Invalid configuration for muPair_Higgs = " << cfg.muPair_Higgs << std::endl;
    assert(false);
  }
  
  return candPair;
} // End function: MuPairInfo CandPair()


// Return selected electrons in event
EleInfos SelectedEles ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const bool verbose ) {

  EleInfos selEles;
  for (const auto & ele : (*br.eles)) {
    if ( ElePass(cfg, ele, br) ) {
	selEles.push_back(ele);
    }
  }
  return selEles;
}


// Return selected jets in event
JetInfos SelectedJets ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel, const bool verbose ) {

  JetInfos selJets;
  for (const auto & jet : (*br.jets)) {
    if ( JetPass(cfg, jet, br, sel) ) {
      selJets.push_back(jet);
    }
  }
  return selJets;
}


// Return selected dijet pairs in event
JetPairInfos SelectedJetPairs ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel, const bool verbose ) {

  JetPairInfos selJetPairs;
  for (const auto & jetPair : (*br.jetPairs)) {
    JetInfo jet1 = br.jets->at(jetPair.iJet1);
    JetInfo jet2 = br.jets->at(jetPair.iJet2);
    if ( JetPass(cfg, jet1, br, sel) &&
	 JetPass(cfg, jet2, br, sel) ) {
      selJetPairs.push_back(jetPair);
    }
  }
  return selJetPairs;
}


