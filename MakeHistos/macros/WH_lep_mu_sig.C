
//////////////////////////////////////////////
///   Macro to study WH, W-->lv category   ///
///                                        ///
/// Specifically, how to correctly choose  ///
///    2 muons coming from Higgs decay     ///
///                                        ///
///     Andrew Brinkerhoff  25.09.2018     ///
//////////////////////////////////////////////

// Basic ROOT includes to read and write files
#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h" // List of branches in the NTuple tree
#include "H2MuAnalyzer/MakeHistos/interface/HistoHelper.h"        // Use to book and fill output histograms
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"    // Common object selections
#include "H2MuAnalyzer/MakeHistos/interface/EventSelection.h"     // Common event selections
#include "H2MuAnalyzer/MakeHistos/interface/EventWeight.h"        // Common event weights
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"       // Common category definitions

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 1;     // Maximum index of input files to process
const int MAX_EVT  = 10000; // Maximum number of events to process
const int PRT_EVT  = 1000;  // Print every N events
const float SAMP_WGT = 1.0;
const bool verbose = false; // Print extra information

const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos/170315_105045/0000";
const TString SAMPLE   = "H2Mu_WH_pos";

const std::string YEAR = "2016";
const std::string SLIM = "Slim";
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {"Presel2016"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE", "3GenMu"}; // Event selection categories, also applied in parallel

// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void WH_lep_mu_sig( TString sample = "", TString in_dir = "", TString out_dir = "",
		    std::vector<TString> in_files = {}, TString out_file_str = "",
		    int max_evt = 0, int prt_evt = 0, float samp_weight = 1.0) {

  // Set variables to hard-coded values if they are not initialized
  if (sample.Length()  == 0) sample  	 = SAMPLE;
  if (in_dir.Length()  == 0) in_dir  	 = IN_DIR;
  if (out_dir.Length() == 0) out_dir 	 = OUT_DIR;
  if (max_evt          == 0) max_evt 	 = MAX_EVT;
  if (prt_evt          == 0) prt_evt 	 = PRT_EVT;
  if (samp_weight      == 0) samp_weight = SAMP_WGT;


  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  TString in_file_name;
  for (int i = 0; i < in_files.size(); i++) {
    in_file_name.Form("%s/%s", in_dir.Data(), in_files.at(i).Data());
    std::cout << "Adding file " << in_file_name.Data() << std::endl;
    in_file_names.push_back(in_file_name);
  }
  if (in_files.size() == 0) {
    for (int i = MIN_FILE; i <= MAX_FILE; i++) {
      in_file_name.Form("%s/tuple_%d.root", in_dir.Data(), i);
      std::cout << "Adding file " << in_file_name.Data() << std::endl;
      in_file_names.push_back(in_file_name.Data());
    }
  }

  // Open all input files
  for (int i = 0; i < in_file_names.size(); i++) {
    if ( !gSystem->AccessPathName(in_file_names.at(i)) )
      file_tmp = TFile::Open( in_file_names.at(i) ); // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
      return;
    }
  }

  // Initialize set of pointers to all branches in tree
  NTupleBranches br;

  // Initialize empty map of histogram names to histograms
  std::map<TString, TH1*> h_map_1D;
  std::map<TString, TH2*> h_map_2D;

  // Add trees from the input files to the TChain
  TChain * in_chain = new TChain("dimuons/tree");
  for (int i = 0; i < in_file_names.size(); i++) {
    in_chain->Add( in_file_names.at(i) );

    // Set branch addresses, from interface/LoadNTupleBranches.h
    if (sample.Contains("SingleMu"))
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM}, false); // Options in {} include "JES", "Flags", and "SFs"
    else
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM, "GEN", "Wgts"}, false); // Options in {} include "JES", "Flags", and "SFs"
  }

  // Configuration for object selection, event selection, and object weighting
  ObjectSelectionConfig obj_sel;
  EventSelectionConfig  evt_sel;
  EventWeightConfig     evt_wgt;
  ConfigureObjectSelection(obj_sel, YEAR);
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);
  if (verbose) obj_sel.Print();
  if (verbose) evt_sel.Print();
  if (verbose) evt_wgt.Print();

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later


  std::cout << "\n******* About to enter the event loop *******" << std::endl;
  for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++) {
    
    if (iEvt > max_evt && max_evt > 0) break;
    if ( (iEvt % prt_evt) == 0 ) {
      std::cout << "\n*********************" << std::endl;
      std::cout << "Looking at event " << iEvt <<  std::endl;
      std::cout << "*********************" << std::endl;
    }
    
    if (verbose) std::cout << "Before running GetEntry, event = " << br.event;
    
    in_chain->GetEntry(iEvt);
    
    if (verbose) std::cout << "... after, event = " << br.event << std::endl;

    // For original 2016 and some 2017 NTuples, convert "SlimJets" collection into regular jets
    JetInfos jets_tmp;
    if (YEAR == "2016") {
      jets_tmp = ConvertSlimJets(*(br.slimJets));
      br.jets  = &jets_tmp;
    }


    ///////////////////////////////////////////
    ///  Apply selection and category cuts  ///
    ///////////////////////////////////////////

    // Check if event passes basic selection cuts defined in src/SelectionCuts.cc
    bool pass_sel_cuts = true;
    for (int i = 0; i < SEL_CUTS.size(); i++) {
      if (not PassSelection(br, evt_sel, obj_sel, SEL_CUTS.at(i), verbose)) pass_sel_cuts = false;
    }

    if (pass_sel_cuts) {
      
      // Get event weight for MC, defined in src/EventWeight.cc
      float event_wgt = ( sample.Contains("SingleMu") ? 1.0 : EventWeight(br, evt_wgt, verbose) );

      // Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc
      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	if (not PassSelection(br, evt_sel, obj_sel, OPT_CUTS.at(iOpt), verbose)) continue;
	
	// Loop through category cuts defined in src/CategoryCuts.cc
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	  
	  if (CAT_CUTS.at(iCat) == "3GenMu") {
	    bool hasWToMu = false;
	    for (const auto & genMu : (*br.genMuons)) {
	      if (abs(genMu.mother_ID) == 24) hasWToMu = true;
	    }
	    if (not hasWToMu) continue;
	  }
	  else if (not InCategory(obj_sel, br, CAT_CUTS.at(iCat), verbose)) continue;
	  
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUTS.at(iOpt) << ", in category " << CAT_CUTS.at(iCat) << std::endl;
	  std::string h_pre = (std::string)sample + "_"+OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)+"_";
	  
	  /////////////////////////////////
	  ///  Generate and fill plots  ///
	  /////////////////////////////////

	  // Store event weights
	  BookAndFill(h_map_1D, h_pre+"PU_wgt",    400, -2, 2, br.PU_wgt  );
	  BookAndFill(h_map_1D, h_pre+"muon_wgt",  400, -2, 2, event_wgt / (br.PU_wgt * br.GEN_wgt) );
	  BookAndFill(h_map_1D, h_pre+"GEN_wgt",   400, -2, 2, br.GEN_wgt );
	  BookAndFill(h_map_1D, h_pre+"event_wgt", 400, -2, 2, event_wgt  );

	  TLorentzVector MET_vec = FourVec(*br.met);
	  TLorentzVector mu_H_pos_vec;
	  TLorentzVector mu_H_neg_vec;
	  TLorentzVector mu_W_vec;
	  TLorentzVector GEN_W_vec;
	  int W_charge = 0;

	  // Plot GEN muon kinematics from Higgs
	  for (const auto & genMu : (*br.genMuons)) {
	    
	    if (     genMu.mother_ID  == 25) BookAndFill(h_map_1D, h_pre+"GEN_mu_H_pt", 100, 0, 200, genMu.pt, event_wgt );
	    if ( abs(genMu.mother_ID) == 24) BookAndFill(h_map_1D, h_pre+"GEN_mu_W_pt", 100, 0, 200, genMu.pt, event_wgt );
	  
	    if (     genMu.mother_ID  == 25) BookAndFill(h_map_1D, h_pre+"GEN_mu_H_eta", 50, -2.5, 2.5, genMu.eta, event_wgt );
	    if ( abs(genMu.mother_ID) == 24) BookAndFill(h_map_1D, h_pre+"GEN_mu_W_eta", 50, -2.5, 2.5, genMu.eta, event_wgt );

	    if (genMu.mother_ID == 25 && genMu.charge ==  1) mu_H_pos_vec = FourVec(genMu);
	    if (genMu.mother_ID == 25 && genMu.charge == -1) mu_H_neg_vec = FourVec(genMu);
	    if ( abs(genMu.mother_ID) == 24) {
	      mu_W_vec  = FourVec(genMu);
	      GEN_W_vec = FourVec(br.genParents->at(genMu.mother_idx));
	      W_charge = genMu.charge;
	    }

	  } // End loop: for (const auto & genMu : (*br.genMuons))

	  TLorentzVector mu_H_OS_vec = (W_charge == 1 ? mu_H_neg_vec : mu_H_pos_vec);
	  TLorentzVector mu_H_SS_vec = (W_charge == 1 ? mu_H_pos_vec : mu_H_neg_vec);
	  TLorentzVector nu_W_vec    = GEN_W_vec - mu_W_vec;
	  
	  TLorentzVector dimu_H_vec = mu_H_pos_vec + mu_H_neg_vec;
	  TLorentzVector dimu_W_vec = mu_H_OS_vec + mu_W_vec;
	  TLorentzVector munu_W_vec = mu_W_vec + nu_W_vec;

	  // Get transverse component of muon 4-vectors
	  TLorentzVector mu_H_SS_vecT;
	  TLorentzVector mu_W_vecT;
	  TLorentzVector nu_W_vecT;
	  mu_H_SS_vecT.SetPtEtaPhiM(mu_H_SS_vec.Pt(), 0, mu_H_SS_vec.Phi(), 0);
	  mu_W_vecT   .SetPtEtaPhiM(mu_W_vec.Pt(),    0, mu_W_vec.Phi(),    0);
	  nu_W_vecT   .SetPtEtaPhiM(nu_W_vec.Pt(),    0, nu_W_vec.Phi(),    0);

	  // GEN-level dimuon mass from true Higgs and OS muon pair from Higgs and W
	  BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_mass", 100, 0, 500, dimu_H_vec.M(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_mass", 100, 0, 500, dimu_W_vec.M(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_munu_mass", 100, 0, 500, munu_W_vec.M(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_munu_MT",   100, 0, 500, (mu_W_vecT + nu_W_vecT).M(), event_wgt );

	  // GEN-level muon and dimuon pT from true Higgs and OS muon pair from Higgs and W
	  BookAndFill(h_map_1D, h_pre+"GEN_H_SS_mu_pt", 100, 0, 500, mu_H_SS_vec.Pt(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_H_OS_mu_pt", 100, 0, 500, mu_H_OS_vec.Pt(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_mu_pt",    100, 0, 500, mu_W_vec.Pt(),    event_wgt );

	  BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_pt",  100, 0, 500, dimu_H_vec.Pt(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_pt",  100, 0, 500, dimu_W_vec.Pt(),  event_wgt );

	  BookAndFill(h_map_1D, h_pre+"GEN_H_SS_minus_W_mu_pt", 100, -250, 250, mu_H_SS_vec.Pt() - mu_W_vec.Pt(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_H_SS_over_W_mu_pt",  100,  0.0, 4.0, mu_H_SS_vec.Pt() / mu_W_vec.Pt(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_dimu_pt",  100, -250, 250, dimu_H_vec.Pt() - dimu_W_vec.Pt(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_dimu_pt",   100,  0.0, 4.0, dimu_H_vec.Pt() / dimu_W_vec.Pt(), event_wgt );

	  // GEN-level dimuon dR from true Higgs and OS muon pair from Higgs and W
	  BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_dR",  32, 0, 6.4, mu_H_pos_vec.DeltaR(mu_H_neg_vec), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_dR", 32, 0, 6.4, mu_H_OS_vec.DeltaR(mu_W_vec), event_wgt );

	  // GEN-level dimuon dEta from true Higgs and OS muon pair from Higgs and W
	  BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_dEta",  50, -5, 5, mu_H_pos_vec.Eta() - mu_H_neg_vec.Eta(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_dEta", 50, -5, 5, mu_H_OS_vec.Eta() - mu_W_vec.Eta(), event_wgt );

	  // GEN-level dimuon dPhi from true Higgs and OS muon pair from Higgs and W
	  BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_dPhi",  32, -3.2, 3.2, mu_H_pos_vec.DeltaPhi(mu_H_neg_vec), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_dPhi", 32, -3.2, 3.2, mu_H_OS_vec.DeltaPhi(mu_W_vec), event_wgt );

	  // Transverse mass ("MT") between muon from Higgs or W and MET vector
	  BookAndFill(h_map_1D, h_pre+"GEN_H_mu_MET_MT", 100, 0, 500, (MET_vec + mu_H_SS_vecT).M() );
	  BookAndFill(h_map_1D, h_pre+"GEN_W_mu_MET_MT", 100, 0, 500, (MET_vec + mu_W_vecT).M() );
	  BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_mu_MET_MT", 100, -250, 250, (MET_vec + mu_H_SS_vecT).M() - (MET_vec + mu_W_vecT).M() );
	  BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_mu_MET_MT",  100,  0.0, 4.0, (MET_vec + mu_H_SS_vecT).M() / (MET_vec + mu_W_vecT).M() );

	  // Plot specific cases where OS dimuon pair from W falls inside Higgs candidate mass window
	  if (dimu_W_vec.M() > 110 && dimu_W_vec.M() < 150) {
	    BookAndFill(h_map_1D, h_pre+"GEN_H_SS_mu_pt_mass_110_150",  100, 0, 500, mu_H_SS_vec.Pt(), event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_W_mu_pt_mass_110_150",     100, 0, 500, mu_W_vec.Pt(),    event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_pt_mass_110_150",   100, 0, 500, dimu_H_vec.Pt(),  event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_pt_mass_110_150",   100, 0, 500, dimu_W_vec.Pt(),  event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_H_mu_MET_MT_mass_110_150", 100, 0, 500, (MET_vec + mu_H_SS_vecT).M(), event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_W_mu_MET_MT_mass_110_150", 100, 0, 500, (MET_vec + mu_W_vecT).M(), event_wgt );
	    
	    BookAndFill(h_map_1D, h_pre+"GEN_H_SS_minus_W_mu_pt_mass_110_150_MT_150",  100, -250, 250, mu_H_SS_vec.Pt() - mu_W_vec.Pt(), event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_H_SS_over_W_mu_pt_mass_110_150_MT_150",   100,  0.0, 4.0, mu_H_SS_vec.Pt() / mu_W_vec.Pt(), event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_dimu_pt_mass_110_150_MT_150",   100, -250, 250, dimu_H_vec.Pt() - dimu_W_vec.Pt(), event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_dimu_pt_mass_110_150_MT_150",    100,  0.0, 4.0, dimu_H_vec.Pt() / dimu_W_vec.Pt(), event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_mu_MET_MT_mass_110_150_MT_150", 100, -250, 250, (MET_vec + mu_H_SS_vecT).M() - (MET_vec + mu_W_vecT).M(), event_wgt );
	    BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_mu_MET_MT_mass_110_150_MT_150",  100,  0.0, 4.0, (MET_vec + mu_H_SS_vecT).M() / (MET_vec + mu_W_vecT).M(), event_wgt );

	    // Plot specific cases where MT(muon from Higgs, MET) < 150
	    if ((MET_vec + mu_H_SS_vecT).M() < 150) {
	      BookAndFill(h_map_1D, h_pre+"GEN_H_SS_mu_pt_mass_110_150_MT_150",  100, 0, 500, mu_H_SS_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_W_mu_pt_mass_110_150_MT_150",     100, 0, 500, mu_W_vec.Pt(),    event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_pt_mass_110_150_MT_150",   100, 0, 500, dimu_H_vec.Pt(),  event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_pt_mass_110_150_MT_150",   100, 0, 500, dimu_W_vec.Pt(),  event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_mu_MET_MT_mass_110_150_MT_150", 100, 0, 500, (MET_vec + mu_H_SS_vecT).M(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_W_mu_MET_MT_mass_110_150_MT_150", 100, 0, 500, (MET_vec + mu_W_vecT).M(), event_wgt );

	      BookAndFill(h_map_1D, h_pre+"GEN_H_SS_minus_W_mu_pt_mass_110_150_MT_150",  100, -250, 250, mu_H_SS_vec.Pt() - mu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_SS_over_W_mu_pt_mass_110_150_MT_150",   100,  0.0, 4.0, mu_H_SS_vec.Pt() / mu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_dimu_pt_mass_110_150_MT_150",   100, -250, 250, dimu_H_vec.Pt() - dimu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_dimu_pt_mass_110_150_MT_150",    100,  0.0, 4.0, dimu_H_vec.Pt() / dimu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_mu_MET_MT_mass_110_150_MT_150", 100, -250, 250, (MET_vec + mu_H_SS_vecT).M() - (MET_vec + mu_W_vecT).M(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_mu_MET_MT_mass_110_150_MT_150",  100,  0.0, 4.0, (MET_vec + mu_H_SS_vecT).M() / (MET_vec + mu_W_vecT).M(), event_wgt );
	    } // End if ((MET_vec + mu_H_SS_vecT).M() < 150)

	    // Plot specific cases where MT(muon from Higgs, MET) < 110
	    if ((MET_vec + mu_H_SS_vecT).M() < 110) {
	      BookAndFill(h_map_1D, h_pre+"GEN_H_SS_mu_pt_mass_110_150_MT_110",  100, 0, 500, mu_H_SS_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_W_mu_pt_mass_110_150_MT_110",     100, 0, 500, mu_W_vec.Pt(),    event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_dimu_pt_mass_110_150_MT_110",   100, 0, 500, dimu_H_vec.Pt(),  event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_W_dimu_pt_mass_110_150_MT_110",   100, 0, 500, dimu_W_vec.Pt(),  event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_mu_MET_MT_mass_110_150_MT_110", 100, 0, 500, (MET_vec + mu_H_SS_vecT).M(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_W_mu_MET_MT_mass_110_150_MT_110", 100, 0, 500, (MET_vec + mu_W_vecT).M(), event_wgt );

	      BookAndFill(h_map_1D, h_pre+"GEN_H_SS_minus_W_mu_pt_mass_110_150_MT_110",  100, -250, 250, mu_H_SS_vec.Pt() - mu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_SS_over_W_mu_pt_mass_110_150_MT_110",   100,  0.0, 4.0, mu_H_SS_vec.Pt() / mu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_dimu_pt_mass_110_150_MT_110",   100, -250, 250, dimu_H_vec.Pt() - dimu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_dimu_pt_mass_110_150_MT_110",    100,  0.0, 4.0, dimu_H_vec.Pt() / dimu_W_vec.Pt(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_minus_W_mu_MET_MT_mass_110_150_MT_110", 100, -250, 250, (MET_vec + mu_H_SS_vecT).M() - (MET_vec + mu_W_vecT).M(), event_wgt );
	      BookAndFill(h_map_1D, h_pre+"GEN_H_over_W_mu_MET_MT_mass_110_150_MT_110",  100,  0.0, 4.0, (MET_vec + mu_H_SS_vecT).M() / (MET_vec + mu_W_vecT).M(), event_wgt );
	    } // End if ((MET_vec + mu_H_SS_vecT).M() < 150)

	  } // End if (dimu_W_vec.M() > 110 && dimu_W_vec.M() < 150)

	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
    } // End conditional: if (pass_sel_cuts)


  } // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;

  std::cout << "\n******* normalizing histos, weight =  " << samp_weight << " *******" << std::endl;
  if (h_map_1D.empty()) std::cout << "h_map_1D is empty" << std::endl;
  else {
    for (std::map<TString, TH1*>::iterator it_term = h_map_1D.begin() ; it_term != h_map_1D.end() ; it_term++) {
        it_term->second ->Scale(samp_weight);
    }
  }
  if (h_map_2D.empty()) std::cout << "h_map_2D is empty" << std::endl;
  else {
    for (std::map<TString, TH2*>::iterator it_term = h_map_2D.begin() ; it_term != h_map_2D.end() ; it_term++) {
        it_term->second ->Scale(samp_weight);
    }
  }

  // Concatenate basic selection cuts into string for output file name
  std::string selStr = "";
  for (int i = 0; i < SEL_CUTS.size(); i++) selStr = selStr+SEL_CUTS.at(i)+"_";
  selStr.pop_back();

  // Place histograms made with separate selection and cateogry cuts into separate files
  for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
    for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
      std::string optCatStr = OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat);
      
      // Create output file
      TString out_file_name;
      if (out_file_str.Length() > 0) out_file_name.Form( "%s/histos_%s_%s_%s_%s.root",    out_dir.Data(), sample.Data(), selStr.c_str(), optCatStr.c_str(), out_file_str.Data() );
      else                           out_file_name.Form( "%s/histos_%s_%s_%s_%d_%d.root", out_dir.Data(), sample.Data(), selStr.c_str(), optCatStr.c_str(), MIN_FILE, MAX_FILE );
      std::cout << "\nCreating output file " << out_file_name.Data() << std::endl;
      TFile* out_file = TFile::Open( out_file_name, "RECREATE" );
      
      // Write output file
      if (verbose) std::cout << "\nWriting output file " << out_file_name.Data() << std::endl;
      out_file->cd();

      // Write out 1D histograms
      for (std::map<TString, TH1*>::iterator it = h_map_1D.begin(); it != h_map_1D.end(); ++it) {
	std::string h_name = it->second->GetName();
	if (h_name.find(optCatStr+"_") != std::string::npos) {
	  // Remove optional selection and category cuts from histogram names
// cat name may be useful
	  h_name.erase( h_name.find(optCatStr+"_"), optCatStr.length() + 1 );
	  it->second->SetName(h_name.c_str());
	  it->second->SetTitle(h_name.c_str());
	  std::cout << "  * Writing 1D histogram " << it->second->GetName() << std::endl;
	  it->second->Write();
	}
      }
      // Write out 2D histograms
      for (std::map<TString, TH2*>::iterator it = h_map_2D.begin(); it != h_map_2D.end(); ++it) {
	std::string h_name = it->second->GetName();
	if (h_name.find(optCatStr+"_") != std::string::npos) {
	  // Remove optional selection and category cuts from histogram names
	  h_name.erase( h_name.find(optCatStr+"_"), optCatStr.length() + 1 );
	  it->second->SetName(h_name.c_str());
	  it->second->SetTitle(h_name.c_str());
	  std::cout << "  * Writing 2D histogram " << it->second->GetName() << std::endl;
	  it->second->Write();
	}
      }
  
      out_file->Write();
      std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
      
    } // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
  } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
      
  std::cout << "\nExiting WH_lep_mu_sig()\n";
  
} // End void WH_lep_mu_sig()
