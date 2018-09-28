#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"

// Configure constants related to object selection
void ConfigureObjectSelection( ObjectSelectionConfig & cfg, const std::string _year ) {

  if (_year == "2016") {
    cfg.year = _year;

    // Muon selection
    cfg.mu_pt_corr = "KaMu";   // Muon pT correction: "PF", "Roch", or "KaMu"
    cfg.mu_pt_min  = 20.0;     // Minimum muon pT
    cfg.mu_eta_max =  2.4;     // Maximum muon |eta|
    cfg.mu_ID_cut  = "medium"; // Muon ID: "loose", "medium", or "tight"
    cfg.mu_iso_max = 0.25;     // Maximum muon relative isolation

    // Jet selection
    cfg.jet_pt_min    = 30.0;   // Minimum jet pT
    cfg.jet_eta_max   =  4.7;   // Maximum jet |eta|
    cfg.jet_PU_ID_cut = "NONE"; // Jet passes PU ID cut
  } // End if (_year == "2016")

  else if (_year == "2017") {
    cfg.year = _year;

    // Muon selection
    cfg.mu_pt_corr = "KaMu";   // Muon pT correction: "PF", "Roch", or "KaMu"
    cfg.mu_pt_min  = 20.0;     // Minimum muon pT
    cfg.mu_eta_max =  2.4;     // Maximum muon |eta|
    cfg.mu_ID_cut  = "medium"; // Muon ID: "loose", "medium", or "tight"
    cfg.mu_iso_max = 0.25;     // Maximum muon relative isolation

    // Jet selection
    cfg.jet_pt_min    = 30.0;    // Minimum jet pT
    cfg.jet_eta_max   =  4.7;    // Maximum jet |eta|
    cfg.jet_PU_ID_cut = "loose"; // Jet passes PU ID cut
  } // End if (_year == "2017")

  else {
    std::cout << "Inside ConfigureEventWeight, invalid year = " << _year << std::endl;
    assert(false);
  }

} // End function: ConfigureObjectSelection()


// Select muons passing ID and kinematic cuts
bool MuonPass ( const ObjectSelectionConfig & cfg, const MuonInfo & muon, const bool verbose ) {

  if ( MuonPt(muon, cfg.mu_pt_corr) < cfg.mu_pt_min  ) return false;
  if ( fabs(muon.eta)               > cfg.mu_eta_max ) return false;
  if ( MuonID(muon, cfg.mu_ID_cut) != true           ) return false;
  if ( muon.relIso                  > cfg.mu_iso_max ) return false;

  return true;
} // End function: bool MuonPass()


// Select jets passing ID and kinematic cuts
bool JetPass( const ObjectSelectionConfig & cfg, const JetInfo & jet, const bool verbose ) {

  if ( jet.pt        < cfg.jet_pt_min  )                    return false;
  if ( fabs(jet.eta) > cfg.jet_eta_max )                    return false;
  if ( JetPUID(jet, cfg.jet_PU_ID_cut, cfg.year ) != true ) return false;

  return true;
} // End function: bool JetPass()

