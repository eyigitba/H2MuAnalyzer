#ifndef OBJECT_SELECTION
#define OBJECT_SELECTION

#include <iostream>
#include <cmath>

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"
#include "H2MuAnalyzer/MakeHistos/interface/ObjectHelper.h"

struct ObjectSelectionConfig {  // Default values taken from 2016

  // General settings
  std::string year = "NONE";
  
  // Muon selection
  std::string mu_pt_corr = "NONE";  // Muon pT correction: "PF", "Roch", or "KaMu"
  float       mu_pt_min  =  -99.0;  // Minimum muon pT
  float       mu_eta_max =  -99.0;  // Maximum muon |eta|
  std::string mu_ID_cut  = "NONE";  // Muon ID: "loose", "medium", or "tight"
  float       mu_iso_max =  -99.0;  // Maximum muon relative isolation

  // Jet selection
  float jet_pt_min          =  -99.0;  // Minimum jet pT
  float jet_eta_max         =  -99.0;  // Maximum jet |eta|
  std::string jet_PU_ID_cut = "NONE";  // Jet passes PU ID cut


  void Print() {
    std::cout << "\n*** ObjectSelectionConfig for year = " << year << " ***" << std::endl;
    std::cout << "Muons: pt_corr = " << mu_pt_corr << ", pt_min = " << mu_pt_min << ", eta_max = " << mu_eta_max << ", ID_cut = " << mu_ID_cut << ", iso_max = " << mu_iso_max << std::endl;
    std::cout << "Jets: pt_min = " << jet_pt_min << ", eta_max = " << jet_eta_max << ", PU_ID_cut = " << jet_PU_ID_cut << "\n" << std::endl;
  } // End function: void Print()

}; // End struct ObjectSelectionConfig

void ConfigureObjectSelection( ObjectSelectionConfig & cfg, const std::string _year );

bool MuonPass ( const ObjectSelectionConfig & cfg, const MuonInfo & muon, const bool verbose = false );
bool JetPass  ( const ObjectSelectionConfig & cfg, const JetInfo & jet,   const bool verbose = false );


#endif  // #ifndef OBJECT_SELECTION
