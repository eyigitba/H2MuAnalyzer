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

  // Electron selection
  float       ele_pt_min  =  -99.0;  // Minimum electron pT
  float       ele_eta_max =  -99.0;  // Maximum electron |eta|
  std::string ele_ID_cut  = "NONE";  // Electron ID: "loose", "medium", or "tight"
  float       ele_iso_max =  -99.0;  // Maximum electron relative isolation

  // Jet selection
  float jet_pt_min       =      -99.0;  // Minimum jet pT
  float jet_eta_max      =      -99.0;  // Maximum jet |eta|
  float jet_mu_dR_min    =      -99.0;  // Minimum dR(muon, jet)
  std::string jet_PU_ID_cut = "NONE";   // Jet passes PU ID cut
  std::vector<float> jet_btag_cuts = {99,99,99}; // Loose, medium, and tight b-tag thresholds 

  // Higgs candidate pair selection
  std::string muPair_Higgs = "NONE"; // How to choose dimuon candidate pair
                                     // e.g. sort_OS_sum_muon_pt, sort_OS_dimuon_pt, etc.

  void Print() {
    std::cout << "\n*** ObjectSelectionConfig for year = " << year << " ***" << std::endl;
    std::cout << "Muons: pt_corr = " << mu_pt_corr << ", pt_min = " << mu_pt_min << ", eta_max = " << mu_eta_max << ", ID_cut = " << mu_ID_cut << ", iso_max = " << mu_iso_max << std::endl;
    std::cout << "Electrons: pt_min = " << ele_pt_min << ", eta_max = " << ele_eta_max << ", ID_cut = " << ele_ID_cut << ", iso_max = " << ele_iso_max << std::endl;
    std::cout << "Jets: pt_min = " << jet_pt_min << ", eta_max = " << jet_eta_max << ", PU_ID_cut = " << jet_PU_ID_cut << "\n" << std::endl;
  } // End function: void Print()

}; // End struct ObjectSelectionConfig

void ConfigureObjectSelection( ObjectSelectionConfig & cfg, const std::string _year );

bool MuonPass ( const ObjectSelectionConfig & cfg, const MuonInfo & muon, const bool verbose = false );
bool ElePass  ( const ObjectSelectionConfig & cfg, const EleInfo & ele, const bool verbose = false );
bool JetPass  ( const ObjectSelectionConfig & cfg, const JetInfo & jet, const MuonInfos & muons, const std::string sel = "", const bool verbose = false );

MuonInfos    SelectedMuons    ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const bool verbose = false );  // Return selected muons in event
MuPairInfos  SelectedMuPairs  ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel = "OS", const bool verbose = false );  // Return selected dimuon pairs in event
MuPairInfo   SelectedCandPair ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const bool verbose = false );  // Return selected dimuon Higgs candidate pair
EleInfos     SelectedEles     ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const bool verbose = false );  // Return selected electrons in event
JetInfos     SelectedJets     ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel = "", const bool verbose = false );  // Return selected jets in event
JetPairInfos SelectedJetPairs ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel = "", const bool verbose = false );  // Return selected dijet pairs in event

#endif  // #ifndef OBJECT_SELECTION
