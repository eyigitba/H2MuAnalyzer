
/////////////////////////////////////////////////////////////
///   Macro to study backgrounds in WH, W-->lv category   ///
///          (WZ, ttW, ttZ, Z+fake, ttbar+fake)           ///
///            Andrew Brinkerhoff  03.09.2018             ///
/////////////////////////////////////////////////////////////

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
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
// const TString SAMPLE   = "SingleMu";

const std::string YEAR = "2016";
const std::string SLIM = "Slim";  // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {"Presel2016"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"3LooseMu", "3TightMu"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE", "WZ_3l_val_mu", "ttZ_3l_val_mu", "ttW_3l_val_mu",
                                           "Z_3l_val_mu", "ttbar_3l_val_mu", "WH_3l_mu"}; // Event selection categories, also applied in parallel

// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void WH_lep_bkg_val( TString sample = "", TString in_dir = "", TString out_dir = "",
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
      // in_file_name.Form("%s/NTuple_0.root", in_dir.Data());
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

  evt_sel.muPair_mass_min = 12; // Not looking exclusively at Z --> mu-mu categories
  obj_sel.mu_pt_min       = 10; // Lower muon pT threshold for muons not from Higgs
  obj_sel.muPair_Higgs = "sort_WH_3_mu_v1"; // Choose Higgs candidate based on MT(W muon, MET) 

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
    if (SLIM == "Slim") {
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
	std::string OPT_CUT = OPT_CUTS.at(iOpt);
	if (OPT_CUT == "3LooseMu") {
	  // Require exactly 3 selected loose muons whose charge sums to +/-1
	  obj_sel.mu_pt_min  = 10.0;
	  obj_sel.mu_ID_cut  = "medium";
	  obj_sel.mu_iso_max = 0.25;
	  if ( SelectedMuPairs(obj_sel, br).size() != 2 ) continue;
	  if ( MuPairMass( SelectedMuPairs(obj_sel, br).at(0), PTC ) < 12 ) continue;
	  if ( MuPairMass( SelectedMuPairs(obj_sel, br).at(1), PTC ) < 12 ) continue;
	}
	else if (OPT_CUT == "3TightMu") {
	  // Require exactly 3 selected loose muons whose charge sums to +/-1
	  obj_sel.mu_pt_min  = 10.0;
	  obj_sel.mu_ID_cut  = "medium";
	  obj_sel.mu_iso_max = 0.25;
	  if ( SelectedMuPairs(obj_sel, br).size() != 2 ) continue;
	  // Require that all 3 of those muons pass the tight selection as well
	  obj_sel.mu_pt_min  = 20.0;
	  obj_sel.mu_ID_cut  = "tight";
	  obj_sel.mu_iso_max = 0.12;
	  if ( SelectedMuPairs(obj_sel, br).size() != 2 ) continue;
	  if ( MuPairMass( SelectedMuPairs(obj_sel, br).at(0), PTC ) < 12 ) continue;
	  if ( MuPairMass( SelectedMuPairs(obj_sel, br).at(1), PTC ) < 12 ) continue;
	}
	else if (not PassSelection(br, evt_sel, obj_sel, OPT_CUT, verbose)) continue;
	
	// Loop through category cuts defined in src/CategoryCuts.cc
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {

	  //////////////////////////////////////////////////////
	  ///  Compute variables relevent for category cuts  ///
	  //////////////////////////////////////////////////////

	  // Get selected muons in event
	  MuonInfos muons = SelectedMuons(obj_sel, br);
	  assert(muons.size() == 3);
	  int sum_mu_charge = muons.at(0).charge + muons.at(1).charge + muons.at(2).charge;
	  assert(abs(sum_mu_charge) == 1);
	  
	  MuPairInfo Z_true;  // The true Z boson dimuon pair, if it exists
	  MuPairInfo Z_pair;  // The selected OS dimuon pair closest to the Z mass
	  MuonInfo   nonZ_mu;
	  MuPairInfo nonZ_pair;
	  MuonInfo   SS_mu1;
	  MuonInfo   SS_mu2;
	  MuonInfo   OS_mu;
	  
	  // Loop over selected dimuon pairs
	  float Z_mass = 9999;
	  for (const auto & muPair : SelectedMuPairs(obj_sel, br)) {
	    // Check if the pair is closest to the Z mass
	    if ( abs(MuPairMass(muPair, PTC) - 91) < abs(Z_mass - 91) ) {
	      Z_mass    = MuPairMass(muPair, PTC);
	      nonZ_pair = Z_pair;
	      Z_pair    = muPair;
	    }
	    // Check if the pair matches a GEN muon pair from Z
	    if (not sample.Contains("SingleMu"))
	      if ( IsGenMatched( muPair, *br.muons, *br.genMuons, "Z") )
		Z_true = muPair;
	  }
	  assert(Z_pair.mass > 0 && nonZ_pair.mass > 0 && Z_mass < 9998); // We should have found a Z and non-Z candidate

	  // Find the muon not in the Z pair, and the same-sign and opposite-sign muons
	  for (const auto & mu : muons) {
	    if ( mu.pt != br.muons->at(Z_pair.iMu1).pt && mu.eta != br.muons->at(Z_pair.iMu1).eta &&
		 mu.pt != br.muons->at(Z_pair.iMu2).pt && mu.eta != br.muons->at(Z_pair.iMu2).eta ) {
	      assert(nonZ_mu.pt <= 0); // We should not have found a W candidate before
	      nonZ_mu = mu;
	    }
	    if (mu.charge != sum_mu_charge) OS_mu = mu;
	    for (const auto & mu2 : muons) {
	      if (mu.charge == mu2.charge && (mu.pt != mu2.pt || mu.eta != mu2.eta)) {
		SS_mu1 = ( MuonPt(mu, PTC) > MuonPt(mu2, PTC) ? mu : mu2 );
		SS_mu2 = ( MuonPt(mu, PTC) > MuonPt(mu2, PTC) ? mu2 : mu );
	      }
	    }
	  }
	  assert(nonZ_mu.pt >= obj_sel.mu_pt_min); // We should always find a W candidate
	  assert(SS_mu1.charge == SS_mu2.charge && OS_mu.charge == -1*sum_mu_charge); // We should always find two SS and one OS muon
	  float nonZ_mu_MT = ( FourVec(nonZ_mu, PTC, "T", *br.muons) + FourVec(*br.met) ).M();
	  
	  ///////////////////////////////////////////
	  ///  Apply the category selection cuts  ///
	  ///////////////////////////////////////////

	  std::string CAT_CUT = CAT_CUTS.at(iCat);
	  if ( not InCategory(obj_sel, br, CAT_CUT, verbose) ) continue;
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUT << ", in category " << CAT_CUT << std::endl;
	  std::string h_pre = (std::string)sample + "_"+OPT_CUT+"_"+CAT_CUT+"_";

	  
	  /////////////////////////////////
	  ///  Generate and fill plots  ///
	  /////////////////////////////////

	  // Store event weights
	  if (not sample.Contains("SingleMu")) {
	    BookAndFill(h_map_1D, h_pre+"PU_wgt",    40, -2, 2, br.PU_wgt  );
	    BookAndFill(h_map_1D, h_pre+"muon_wgt",  40, -2, 2, event_wgt / (br.PU_wgt * br.GEN_wgt) );
	    BookAndFill(h_map_1D, h_pre+"GEN_wgt",   40, -2, 2, br.GEN_wgt );
	    BookAndFill(h_map_1D, h_pre+"event_wgt", 40, -2, 2, event_wgt  );
	  }

	  // Plot kinematic histograms
	  BookAndFill(h_map_1D, h_pre+"nJets",     8, -0.5, 7.5, SelectedJets(obj_sel, br).size(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nBJets",    4, -0.5, 3.5, SelectedJets(obj_sel, br, "BTagMedium").size(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nJetsCent", 5, -0.5, 4.5, SelectedJets(obj_sel, br, "Central").size(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nJetsFwd",  5, -0.5, 4.5, SelectedJets(obj_sel, br, "Forward").size(), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"MET", 20, 0, 200, br.met->pt, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"MHT", 20, 0, 200, br.mht->pt, event_wgt );

	  BookAndFill(h_map_1D, h_pre+"mu1_pt", 30, 0, 300, MuonPt(muons.at(0), PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"mu2_pt", 20, 0, 200, MuonPt(muons.at(1), PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"mu3_pt", 20, 0, 100, MuonPt(muons.at(2), PTC), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"Z_mu1_pt",   30, 0, 300, MuonPt(br.muons->at(Z_pair.iMu1), PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"Z_mu2_pt",   20, 0, 200, MuonPt(br.muons->at(Z_pair.iMu2), PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonZ_mu_pt", 20, 0, 200, MuonPt(                  nonZ_mu, PTC), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_pt",  20, 0, 200, MuonPt(OS_mu,  PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_mu1_pt", 20, 0, 200, MuonPt(SS_mu1, PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_mu2_pt", 20, 0, 100, MuonPt(SS_mu2, PTC), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"mu1_eta", 24, -2.4, 2.4, muons.at(0).eta, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"mu2_eta", 24, -2.4, 2.4, muons.at(1).eta, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"mu3_eta", 24, -2.4, 2.4, muons.at(2).eta, event_wgt );

	  BookAndFill(h_map_1D, h_pre+"Z_mu1_eta",   24, -2.4, 2.4, br.muons->at(Z_pair.iMu1).eta, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"Z_mu2_eta",   24, -2.4, 2.4, br.muons->at(Z_pair.iMu2).eta, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonZ_mu_eta", 24, -2.4, 2.4,                   nonZ_mu.eta, event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_iso",  15, 0, 0.3, OS_mu.relIso,  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_mu1_iso", 15, 0, 0.3, SS_mu1.relIso, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_mu2_iso", 15, 0, 0.3, SS_mu2.relIso, event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_tightID",  2, -0.5, 1.5, OS_mu.isTightID,  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_mu1_tightID", 2, -0.5, 1.5, SS_mu1.isTightID, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_mu2_tightID", 2, -0.5, 1.5, SS_mu2.isTightID, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"sum_mu_charge",  2, -0.5, 1.5, sum_mu_charge,    event_wgt );

	  if (not sample.Contains("SingleMu"))
	    BookAndFill(h_map_1D, h_pre+"Z_mass_true", 40, 71, 111, MuPairMass(Z_true, PTC), event_wgt, false ); // Don't include overflow
	  BookAndFill(h_map_1D, h_pre+"Z_mass_zoom",   40, 71, 111, MuPairMass(Z_pair, PTC), event_wgt, false ); // Don't include overflow
	  BookAndFill(h_map_1D, h_pre+"Z_mass_on",      8, 71, 111, MuPairMass(Z_pair, PTC), event_wgt, false ); // Don't include overflow
	  BookAndFill(h_map_1D, h_pre+"Z_mass_off",    40,  0, 400, MuPairMass(Z_pair, PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"Z_pt",          30,  0, 300, MuPairPt  (Z_pair, PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonZ_mass",     40,  0, 400, MuPairMass(nonZ_pair, PTC), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonZ_mu_MT",    40,  0, 200, nonZ_mu_MT, event_wgt );

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
      
  std::cout << "\nExiting WH_lep_bkg_val()\n";
  
} // End void WH_lep_bkg_val()
