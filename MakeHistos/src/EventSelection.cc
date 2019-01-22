#include "H2MuAnalyzer/MakeHistos/interface/EventSelection.h"

// Configure constants related to event selection
void ConfigureEventSelection( EventSelectionConfig & cfg, std::string _year ) {

  if (_year == "2016") {
    cfg.year = _year;

    // Muon selection requirements
    cfg.mu_trig_pt_min    = 26.0;  // Minmum triggering muon pT (or lead muon if no trigger required)
    cfg.mu_trig_HLT_match = true;  // Require muon to match HLT trigger
    cfg.muPair_mass_min   = 60.0;  // Minimum invariant mass
    cfg.muPair_OS         = true;  // Require opposite-sign pair
  } // End if (_year == "2016")

  else if (_year == "2017") {
    cfg.year = _year;

    // Muon selection requirements
    cfg.mu_trig_pt_min    = 30.0;  // Minmum triggering muon pT (or lead muon if no trigger required)
    cfg.mu_trig_HLT_match = true;  // Require muon to match HLT trigger
    cfg.muPair_mass_min   = 60.0;  // Minimum invariant mass
    cfg.muPair_OS         = true;  // Require opposite-sign pair
  } // End if (_year == "2017")

  else {
    std::cout << "Inside ConfigureEventSelection, invalid year = "<< _year << std::endl;
    assert(false);
  }

} // End function: ConfigureEventSelection()


bool PassSelection( const NTupleBranches & br, const EventSelectionConfig & evt,
		    const ObjectSelectionConfig & obj, const std::string selection,
		    const bool verbose ) {

  if (verbose) std::cout << "\nInside PassSelection, checking if event passes " << selection << std::endl;

  bool pass = false;

  //////////////////////////
  ///  No selection cut  ///
  //////////////////////////  
  if (selection == "NONE") pass = true;


  ///////////////////////////////
  ///  Pre-selection in 2016  ///
  ///////////////////////////////
  else if (selection == "Presel2016") {
    if (verbose) std::cout << "\n  * Applying Presel2016 cuts" << std::endl;
    if (verbose) std::cout << "    - Using " << obj.mu_pt_corr << " muon pT corrections"         << std::endl;
    if (verbose) std::cout << "    - Lead (triggering) muon pT > "      << evt.mu_trig_pt_min    << std::endl;
    if (verbose) std::cout << "    - Triggering muon must match HLT = " << evt.mu_trig_HLT_match << std::endl;
    if (verbose) std::cout << "    - Minimum dimuon invariant mass > "  << evt.muPair_mass_min   << std::endl;
    if (verbose) std::cout << "    - Dimuon require opposite sign = "   << evt.muPair_OS         << std::endl;

    assert(evt.year == "2016" && obj.year == "2016");

    // Loop over all muon pairs
    for (int i = 0; i < br.nMuPairs; i++) {
      MuPairInfo & Pair = br.muPairs->at(i);
      MuonInfo   & Mu1  = br.muons->at(Pair.iMu1);
      MuonInfo   & Mu2  = br.muons->at(Pair.iMu2);

      // Require both muons to pass ID cuts
      if ( !MuonPass(obj, Mu1) || !MuonPass(obj, Mu2) ) continue;
      // Require one muon to pass trigger cut
      if ( evt.mu_trig_HLT_match ) {
	if ( !(MuonTrig(Mu1, obj.year) && MuonPt(Mu1, obj.mu_pt_corr) > evt.mu_trig_pt_min) ||
	     !(MuonTrig(Mu2, obj.year) && MuonPt(Mu2, obj.mu_pt_corr) > evt.mu_trig_pt_min) ) continue;
      }
      // Require a minimum dimuon invariant mass
      if (MuPairMass(Pair, obj.mu_pt_corr) < evt.muPair_mass_min) continue;
      // Apply dimuon opposite charge requirement
      if (evt.muPair_OS && Pair.charge != 0) continue;

      // If all conditions pass, break out of the loop
      pass = true;
      break;
    } // End loop: for (int i = 0; i < br.nMuPairs; i++)

  } // End if (selection == "Presel2016")


  ///////////////////////////////
  ///  Pre-selection in 2017  ///
  ///////////////////////////////
  else if (selection == "Presel2017") {
    if (verbose) std::cout << "\n  * Applying Presel2017 cuts" << std::endl;
    if (verbose) std::cout << "    - Using " << obj.mu_pt_corr << " muon pT corrections"         << std::endl;
    if (verbose) std::cout << "    - Lead (triggering) muon pT > "      << evt.mu_trig_pt_min    << std::endl;
    if (verbose) std::cout << "    - Triggering muon must match HLT = " << evt.mu_trig_HLT_match << std::endl;
    if (verbose) std::cout << "    - Minimum dimuon invariant mass > "  << evt.muPair_mass_min   << std::endl;
    if (verbose) std::cout << "    - Dimuon require opposite sign = "   << evt.muPair_OS         << std::endl;
    
    assert(evt.year == "2017" && obj.year == "2017");

    // Loop over all muon pairs
    for (int i = 0; i < br.nMuPairs; i++) {
      MuPairInfo & Pair = br.muPairs->at(i);
      MuonInfo   & Mu1  = br.muons->at(Pair.iMu1);
      MuonInfo   & Mu2  = br.muons->at(Pair.iMu2);

      // Require both muons to pass ID cuts
      if ( !MuonPass(obj, Mu1) || !MuonPass(obj, Mu2) ) continue;
      // Require one muon to pass trigger cut
      if ( evt.mu_trig_HLT_match ) {
	if ( !(MuonTrig(Mu1, obj.year) && MuonPt(Mu1, obj.mu_pt_corr) > evt.mu_trig_pt_min) ||
	     !(MuonTrig(Mu2, obj.year) && MuonPt(Mu2, obj.mu_pt_corr) > evt.mu_trig_pt_min) ) continue;
      }
      // Require a minimum dimuon invariant mass
      if (MuPairMass(Pair, obj.mu_pt_corr) < evt.muPair_mass_min) continue;
      // Apply dimuon opposite charge requirement
      if (evt.muPair_OS && Pair.charge != 0) continue;

      // If all conditions pass, break out of the loop
      pass = true;
      break;
    } // End loop: for (int i = 0; i < br.nMuPairs; i++)

  } // End else if (selection == "Presel2017")

 
  else {
    std::cout << "\nInside EventSelection.cc, don't recognize selection " << selection << std::endl;
    assert(false);
  }

  return pass;
  
}
