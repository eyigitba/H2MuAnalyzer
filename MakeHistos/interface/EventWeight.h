
#ifndef EVENT_WEIGHT
#define EVENT_WEIGHT

#include <iostream>
#include <assert.h>

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"


struct EventWeightConfig {  // Default values taken from 2016

  // General settings
  std::string year = "NONE";

  // Weights to apply
  bool PU         = false; // PU_wgt   (Pileup reweighting SF)
  bool muon_ID    = false; // MuID_SF  (RECO muon ID SF)
  bool muon_Iso   = false; // MuIso_SF (RECO muon isolation SF)
  bool trig_IsoMu = false; // IsoMu_SF (HLT IsoMu trigger SF)
  bool GEN        = false; // GEN_wgt  (NLO MC GEN weight, ±1)

  // Systematics to use
  std::string SYS = "noSys"; 

  void Print() {
    std::cout << "\n*** EventWeightConfig for year = " << year << " ***" << std::endl;
    std::cout << "Weights to use: PU = " << PU << ", muon_ID = " << muon_ID << ", muon_Iso = " << muon_Iso
	      << ", trig_IsoMu = " << trig_IsoMu << ", GEN = " << GEN << std::endl;
    std::cout << "Systematics to use: " << SYS << std::endl;
  } // End function: void Print()

}; // End struct EventWeightConfig


void ConfigureEventWeight( EventWeightConfig & cfg, const std::string _year, const std::string _SYS = "noSys");

// TODO: implement systematic up/down weights for PU, muon_ID, muon_Iso, and trig_IsoMu - AWB 27.09.2018
float MuonWeight (const NTupleBranches & br, const EventWeightConfig & cfg, const bool verbose = false);
float EventWeight(const NTupleBranches & br, const EventWeightConfig & cfg, const bool verbose = false);


#endif  // #ifndef EVENT_WEIGHT
