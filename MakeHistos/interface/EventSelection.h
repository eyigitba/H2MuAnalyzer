
#ifndef EVENT_SELECTION
#define EVENT_SELECTION

#include <iostream>

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"


struct EventSelectionConfig {

  // General settings
  std::string year = "NONE";

  // Muon selection requirements
  float mu_trig_pt_min    = -99.0;  // Minmum triggering muon pT (or lead muon if no trigger required)
  bool  mu_trig_HLT_match = false;  // Require muon to match HLT trigger
  float muPair_mass_min   = -99.0;  // Minimum invariant mass
  bool  muPair_OS         = false;  // Require opposite-sign pair


  void Print() {
    std::cout << "\n*** EventSelectionConfig for year = " << year << " ***" << std::endl;
    std::cout << "Muons: trig_pt_min = " << mu_trig_pt_min << ", trig_HLT_match = " << mu_trig_HLT_match << std::endl;
    std::cout << "MuPairs: mass_min = " << muPair_mass_min << ", OS = " << muPair_OS << "\n" << std::endl;
  } // End function: void Print()

}; // End struct EventSelectionConfig

void ConfigureEventSelection( EventSelectionConfig & cfg, const std::string _year );

bool PassSelection( const NTupleBranches & br, const EventSelectionConfig & evt,
		    const ObjectSelectionConfig & obj, const std::string selection,
		    const bool verbose = false);


#endif  // #ifndef EVENT_SELECTION
