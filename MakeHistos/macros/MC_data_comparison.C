
/////////////////////////////////////////////////////////
///   Macro to make data/MC comparison for variables  ///
///                                                   ///
///              Xunwu Zuo  12.09.2018                ///
/////////////////////////////////////////////////////////

// Basic ROOT includes to read and write files
#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h" // List of branches in the NTuple tree
#include "H2MuAnalyzer/MakeHistos/interface/HistoHelper.h"        // Use to book and fill output histograms
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"    // Common object selection 
#include "H2MuAnalyzer/MakeHistos/interface/EventSelection.h"     // Common event selection
#include "H2MuAnalyzer/MakeHistos/interface/EventWeight.h"       // Common event weight
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"       // Common category definitions

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 3;     // Maximum index of input files to process
const int MAX_EVT  = 10000; // Maximum number of events to process
const int PRT_EVT  = 1000;   // Print every N events
const float SAMP_WGT = 1.0;
const bool verbose = false; // Print extra information

// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/SingleMuon/SingleMu_2017F/180802_164117/0000";
// const TString SAMPLE   = "SingleMu_2017F";
const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/SingleMuon/SingleMu_2018A/190426_124916/0000";
const TString SAMPLE   = "SingleMu_2018A";
const std::string YEAR = "2017";
const std::string SLIM = "notSlim"; // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR  = "plots";
const TString HIST_TREE = "Hist"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both. Not in use in this macro

const std::vector<std::string> SEL_CUTS = {"Presel2017"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel

// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void MC_data_comparison( TString sample = "", TString in_dir = "", TString out_dir = "",
		         std::vector<TString> in_files = {}, TString out_file_str = "",
		         int max_evt = 0, int prt_evt = 0, float samp_weight = 1.0,
			 TString hist_tree = "" ) {

  // Set variables to hard-coded values if they are not initialized
  if (sample.Length()    == 0) sample  	   = SAMPLE;
  if (in_dir.Length()    == 0) in_dir  	   = IN_DIR;
  if (out_dir.Length()   == 0) out_dir 	   = OUT_DIR;
  if (max_evt            == 0) max_evt 	   = MAX_EVT;
  if (prt_evt            == 0) prt_evt 	   = PRT_EVT;
  if (samp_weight        == 0) samp_weight = SAMP_WGT;
  if (hist_tree.Length() == 0) hist_tree   = HIST_TREE;

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
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
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
//  BookForMCvsData(h_map_1D, (std::string) sample, OPT_CUTS, CAT_CUTS); // initialize histomap 1D, all histos should be claimed here - XWZ 20.09.2018

  // Configuration for object selection, event selection, and object weighting
  ObjectSelectionConfig obj_sel;
  EventSelectionConfig  evt_sel;
  EventWeightConfig     evt_wgt;
  ConfigureObjectSelection(obj_sel, YEAR);
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

//  evt_sel.muPair_mass_min = 105; // Require at least one Higgs candidate pair, default is 60
  // obj_sel.muPair_Higgs = "sort_WH_3_mu_v1"; // Choose Higgs candidate based on MT(W muon, MET)

  if (verbose) obj_sel.Print();
  if (verbose) evt_sel.Print();
  if (verbose) evt_wgt.Print();

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later


  // Add trees from the input files to the TChain
  TChain * in_chain = new TChain("dimuons/tree");
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    in_chain->Add( in_file_names.at(i) );

    // Set branch addresses, from interface/LoadNTupleBranches.h
    if (sample.Contains("SingleMu"))
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM}, false); // Options in {} include "JES", "Flags", and "SFs"
    else
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM, "GEN", "Wgts"}, false); // Options in {} include "JES", "Flags", and "SFs"
  }

  std::cout << "\n******* About to enter the event loop *******" << std::endl;
  for (UInt_t iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++) {
    
    if (iEvt > max_evt && max_evt > 0) break;
    if ( (iEvt % prt_evt) == 0 ) {
      std::cout << "\n*********************" << std::endl;
      std::cout << "Looking at event " << iEvt <<  std::endl;
      std::cout << "*********************" << std::endl;
    }
    
    if (verbose) std::cout << "Before running GetEntry, event = " << br.event;
    
    in_chain->GetEntry(iEvt);
    
    if (verbose) std::cout << "... after, event = " << br.event << std::endl;

    // For 2016 NTuples, convert "SlimJets" collection into regular jets
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

    //////////////////////////
    ///  Get event weight  ///
    //////////////////////////
    

    if (pass_sel_cuts) {
      // Get event weight for MC, defined in src/EventWeight.cc  (samp_weight applied as scale in the end)
      float event_wgt = ( sample.Contains("SingleMu") ? 1.0 : EventWeight(br, evt_wgt, verbose) );    
  
      // Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc
      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	// if (not PassSelection(br, OPT_CUTS.at(iOpt), verbose)) continue;
	
	// Loop through category cuts defined in src/CategoryCuts.cc
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	  //if (not InCategory(br, CAT_CUTS.at(iCat), verbose)) continue;
	  
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUTS.at(iOpt) << ", in category " << CAT_CUTS.at(iCat) << std::endl;
	  std::string h_pre = (std::string)sample + "_"+OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)+"_";
	  
	  /////////////////////////////////
	  ///  Generate and fill plots  ///
	  /////////////////////////////////
	  

	  MuPairInfo dimu  = SelectedCandPair(obj_sel, br);
	  if ( dimu.mass < 70 ||
               dimu.mass > 170 ) continue;  // 70-110 for Z mass validation, 105-160 for analysis window


	  MuonInfos  muons = SelectedMuons(obj_sel, br);
          EleInfos   eles  = SelectedEles(obj_sel, br);
          JetInfos   jets  = SelectedJets(obj_sel, br);

	  MuonInfo   mu_1 = br.muons->at(dimu.iMu1);
          MuonInfo   mu_2 = br.muons->at(dimu.iMu2);

	  TLorentzVector dimu_vec = FourVec( dimu, PTC);
  	  TLorentzVector mu1_vec  = FourVec(br.muons->at(dimu.iMu1), PTC);
  	  TLorentzVector mu2_vec  = FourVec(br.muons->at(dimu.iMu2), PTC);
	  TLorentzVector met_vec;

	  // Function from interface/HistoHelper.h that books a histogram (if it has not already been booked), then fills it
	  if (not sample.Contains("SingleMu") or dimu_vec.M() < 120 or dimu_vec.M() > 130) {
	    BookAndFill(h_map_1D, h_pre+"dimuon_mass", 50, 110, 160, dimu_vec.M(), event_wgt);
	  }
	  BookAndFill(h_map_1D, h_pre+"dimuon_pt", 50, 20, 520, dimu_vec.Pt(), event_wgt);
	  BookAndFill(h_map_1D, h_pre+"dimuon_eta", 50, -5, 5, dimu.eta, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"dimuon_delta_eta", 50, -5, 5, dimu.dEta, event_wgt);
          BookAndFill(h_map_1D, h_pre+"dimuon_delta_phi", 50, -4, 4, dimu.dPhi, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"dimuon_dR", 50, 0, 5, dimu.dR, event_wgt);
	  float d0_diff = 2 * (mu_1.d0_PV - mu_2.d0_PV) / (mu_1.charge - mu_2.charge);
	  BookAndFill(h_map_1D, h_pre+"dimuon_d0_diff", 50, -0.1, 0.1, d0_diff, event_wgt);

          BookAndFill(h_map_1D, h_pre+"leading_muon_pt", 50, 20, 270, mu1_vec.Pt(), event_wgt);
          BookAndFill(h_map_1D, h_pre+"leading_muon_eta", 50, -2.5, 2.5, mu_1.eta, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"leading_muon_d0", 50, -0.1, 0.1, mu_1.d0_PV, event_wgt);
          BookAndFill(h_map_1D, h_pre+"subleading_muon_pt", 50, 20, 270, mu2_vec.Pt(), event_wgt);
          BookAndFill(h_map_1D, h_pre+"subleading_muon_eta", 50, -2.5, 2.5, mu_2.eta, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"subleading_muon_d0", 50, -0.1, 0.1, mu_2.d0_PV, event_wgt);
//          BookAndFill(h_map_1D, h_pre+, event_wgt);
//          BookAndFill(h_map_1D, h_pre+, event_wgt);

	  if (SelectedJetPairs(obj_sel, br).size() > 0) {
	    JetPairInfo dijet = SelectedJetPairs(obj_sel, br).at(0);
	    JetInfo      jet1 = br.jets->at(dijet.iJet1);
  	    JetInfo      jet2 = br.jets->at(dijet.iJet2);
	    BookAndFill(h_map_1D, h_pre+"dijet_mass_1000", 50, 0, 1000, dijet.mass, event_wgt);
	    BookAndFill(h_map_1D, h_pre+"dijet_mass_200", 50, 0, 200, dijet.mass, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_pt_800", 50, 0, 800, dijet.pt, event_wgt);
	    BookAndFill(h_map_1D, h_pre+"dijet_pt_200", 50, 0, 200, dijet.pt, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_eta", 50, -10, 10, dijet.eta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_delta_eta", 50, -10, 10, dijet.dEta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_delta_phi", 50, -4, 4, dijet.dPhi, event_wgt);
	    BookAndFill(h_map_1D, h_pre+"dijet_dR", 50, 0, 10, dijet.dR, event_wgt);
// add jet selection, leading in the event or in the pair
            BookAndFill(h_map_1D, h_pre+"jet1_pt", 50, 30, 530, jet1.pt, event_wgt);
            BookAndFill(h_map_1D, h_pre+"jet1_eta", 50, -5, 5, jet1.eta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"jet2_pt", 50, 30, 530, jet2.pt, event_wgt);
            BookAndFill(h_map_1D, h_pre+"jet2_eta", 50, -5, 5, jet2.eta, event_wgt);
	  } // end of if (SelectedJetPairs(obj_sel, br).size() > 0) 
	  if (SelectedJets(obj_sel, br).size() > 0) {
	    JetInfo      jet0 = br.jets->at(0);
	    BookAndFill(h_map_1D, h_pre+"jet0_pt", 50, 30, 530, jet0.pt, event_wgt);
            BookAndFill(h_map_1D, h_pre+"jet0_eta", 50, -5, 5, jet0.eta, event_wgt);
	  } // end of if (SelectedJets(obj_sel, br).size() > 0)

	  BookAndFill(h_map_1D, h_pre+"MET", 50, 0, 300, (br.met)->pt, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"mht_pt", 50, 0, 300, (br.mht)->pt, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"nJets", 6, -0.5, 5.5, br.nJets, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"nBJets", 5, 0, 5, br.nBMed, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"nMuons", 5, 0, 5, br.nMuons, event_wgt);
          BookAndFill(h_map_1D, h_pre+"nElectrons", 4, 0, 4, br.nEles, event_wgt); 
          BookAndFill(h_map_1D, h_pre+"nVertices", 50, 0, 100, br.nVertices, event_wgt); 
          BookAndFill(h_map_1D, h_pre+"nPU", 50, 0, 100, br.nPU, event_wgt); 	  
	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
    } // End conditional: if (pass_sel_cuts)

    // Before we exit the loop over all events, copy maps of histograms to local memory
    // *** Super-hacky!!! Should instead pass h_map_1D to BookAndFill function! *** - AWB 20.08.2018
//    if (iEvt+1 == in_chain->GetEntries() || (iEvt == max_evt && max_evt > 0)) {
//      std::cout << "\nAt event " << iEvt << ", copying maps to local area" << std::endl;
//      RetrieveMap(h_map_1D);
//      RetrieveMap(h_map_2D);
//    }
    
  } // End loop: for (UInt_t iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
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
      
  std::cout << "\nExiting MC_data_comparison()\n";
  
} // End void MC_data_comparison()
