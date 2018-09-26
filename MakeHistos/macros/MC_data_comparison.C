
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
#include "H2MuAnalyzer/MakeHistos/interface/SelectionCuts.h"      // Common event selections
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"       // Common category definitions
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelections.h"   // Common object selections 

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

const TString IN_DIR   = "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/SingleMuon/SingleMu_2017F/180802_164117/0000";
const TString SAMPLE   = "SingleMu_2017F";
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {"Presel2017"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE", "WHlep", "ZHmu", "ZHele"}; // Event selection categories, also applied in parallel

// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void MC_data_comparison( TString sample = "", TString in_dir = "", TString out_dir = "",
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
    in_file_name.Form("root://eoscms.cern.ch/%s/%s", in_dir.Data(), in_files.at(i).Data());
    std::cout << "Adding file " << in_file_name.Data() << std::endl;
    in_file_names.push_back(in_file_name);
  }
  if (in_files.size() == 0) {
    for (int i = MIN_FILE; i <= MAX_FILE; i++) {
      in_file_name.Form("root://eoscms.cern.ch/%s/tuple_%d.root", in_dir.Data(), i);
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
  BookForMCvsData(h_map_1D, (std::string) sample, OPT_CUTS, CAT_CUTS); // initialize histomap 1D, all histos should be claimed here - XWZ 20.09.2018

  // Add trees from the input files to the TChain
  TChain * in_chain = new TChain("dimuons/tree");
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    in_chain->Add( in_file_names.at(i) );

    // Set branch addresses, from interface/LoadNTupleBranches.h
    SetBranchAddresses(*in_chain, br, {"Wgts","Flags"}, false); // Options in {} include "JES", "Flags", and "SFs"
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

    ///////////////////////////////////////////
    ///  Apply selection and category cuts  ///
    ///////////////////////////////////////////

    // Check if event passes basic selection cuts defined in src/SelectionCuts.cc
    bool pass_sel_cuts = true;
    for (int i = 0; i < SEL_CUTS.size(); i++) {
      if (not PassSelection(br, SEL_CUTS.at(i), verbose)) pass_sel_cuts = false;
    }

    //////////////////////////
    ///  Get event weight  ///
    //////////////////////////
    float PU_wgt = 1;
    float muon_wgt = 1;
    float GEN_wgt = 1;
    if (not sample.Contains("SingleMu")) {
      PU_wgt 	= br.PU_wgt;
      muon_wgt 	= br.IsoMu_SF_3 * br.MuIso_SF_3 * br.MuID_SF_3;
      GEN_wgt 	= br.GEN_wgt;
    }
    float event_wgt = PU_wgt * muon_wgt * GEN_wgt;


    if (pass_sel_cuts) {
      
      // Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc
      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	if (not PassSelection(br, OPT_CUTS.at(iOpt), verbose)) continue;
	
	// Loop through category cuts defined in src/CategoryCuts.cc
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	  if (not InCategory(br, CAT_CUTS.at(iCat), verbose)) continue;
	  
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUTS.at(iOpt) << ", in category " << CAT_CUTS.at(iCat) << std::endl;
	  std::string h_pre = (std::string)sample + "_"+OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)+"_";
	  
	  /////////////////////////////////
	  ///  Generate and fill plots  ///
	  /////////////////////////////////
	  
	  // Loop over opposite-charge dimuon pairs
	  bool found_good_dimu = false;
	  for (UInt_t iPair = 0; iPair < br.nMuPairs; iPair++) { // added dimuon selection -- XWZ 25.09.2018
	    
	    if (verbose) std::cout << "  * Pair " << iPair+1 << " has mass = " << br.muPairs->at(iPair).mass << std::endl;
	    if (found_good_dimu) continue;
	    found_good_dimu = DimuPass(br, br.muPairs->at(iPair));
	    if (!found_good_dimu) continue;
	    // Function from interface/HistoHelper.h that books a histogram (if it has not already been booked), then fills it
//	    if (br.muPairs->at(iPair).mass_Roch < 120 or br.muPairs->at(iPair).mass_Roch > 130) {
	    BookAndFill(h_map_1D, h_pre+"dimuon_mass", 50, 110, 160, br.muPairs->at(iPair).mass_Roch, event_wgt);
//	    }
	    BookAndFill(h_map_1D, h_pre+"dimuon_pt", 50, 20, 520, br.muPairs->at(iPair).pt_Roch, event_wgt);
	    BookAndFill(h_map_1D, h_pre+"dimuon_eta", 50, -5, 5, br.muPairs->at(iPair).eta, event_wgt);
	    BookAndFill(h_map_1D, h_pre+"dimuon_delta_eta", 50, -5, 5, br.muPairs->at(iPair).dEta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dimuon_delta_phi", 50, -4, 4, br.muPairs->at(iPair).dPhi, event_wgt);

            BookAndFill(h_map_1D, h_pre+"leading_muon_pt", 50, 20, 270, br.muons->at( br.muPairs->at(iPair).iMu1 ).pt_Roch, event_wgt);
            BookAndFill(h_map_1D, h_pre+"leading_muon_eta", 50, -2.5, 2.5, br.muons->at( br.muPairs->at(iPair).iMu1 ).eta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"subleading_muon_pt", 50, 20, 270, br.muons->at( br.muPairs->at(iPair).iMu2 ).pt_Roch, event_wgt);
            BookAndFill(h_map_1D, h_pre+"subleading_muon_eta", 50, -2.5, 2.5, br.muons->at( br.muPairs->at(iPair).iMu2 ).eta, event_wgt);
//            BookAndFill(h_map_1D, h_pre+, event_wgt);
//            BookAndFill(h_map_1D, h_pre+, event_wgt);
	  } // End loop: for (UInt_t iPair = 0; iPair < nPairs; iPair++)

	  if (!found_good_dimu) continue;
	  bool found_good_dijet = false;
	  for (UInt_t iDijet = 0; iDijet < br.nJetPairs; iDijet++) { // added dijet selection -- XWZ 25.09.2018
	    if (found_good_dijet) continue;
	    found_good_dijet = DijetPass(br, br.jetPairs->at(iDijet));
	    if (!found_good_dijet) continue;

	    BookAndFill(h_map_1D, h_pre+"dijet_mass", 100, 0, 1000, br.jetPairs->at(iDijet).mass, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_pt", 60, 0, 600, br.jetPairs->at(iDijet).pt, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_eta", 50, -10, 10, br.jetPairs->at(iDijet).eta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_delta_eta", 50, -10, 10, br.jetPairs->at(iDijet).dEta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"dijet_delta_phi", 50, -4, 4, br.jetPairs->at(iDijet).dPhi, event_wgt);
// add jet selection, leading in the event or in the pair
            BookAndFill(h_map_1D, h_pre+"leading_jet_pt", 50, 30, 530, br.jets->at( br.jetPairs->at(iDijet).iJet1 ).pt, event_wgt);
            BookAndFill(h_map_1D, h_pre+"leading_jet_eta", 50, -5, 5, br.jets->at( br.jetPairs->at(iDijet).iJet1 ).eta, event_wgt);
            BookAndFill(h_map_1D, h_pre+"subleading_jet_pt", 50, 30, 530, br.jets->at( br.jetPairs->at(iDijet).iJet2 ).pt, event_wgt);
            BookAndFill(h_map_1D, h_pre+"subleading_jet_eta", 50, -5, 5, br.jets->at( br.jetPairs->at(iDijet).iJet2 ).eta, event_wgt);
	  } // End loop: for (UInt_t iDijet = 0; iDijet < br.nJetPairs; iDijet++)

	  BookAndFill(h_map_1D, h_pre+"MET", 40, 0, 200, (br.met)->pt, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"nJets", 6, -0.5, 5.5, br.nJets, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"nBJets", 5, 0, 5, br.nBMed, event_wgt);
	  BookAndFill(h_map_1D, h_pre+"nMuons", 5, 0, 5, br.nMuons, event_wgt);
          BookAndFill(h_map_1D, h_pre+"nElectrons", 4, 0, 4, br.nEles, event_wgt); 
          BookAndFill(h_map_1D, h_pre+"nVertices", 100, 0, 100, br.nVertices, event_wgt); 
          BookAndFill(h_map_1D, h_pre+"nPU", 100, 0, 100, br.nPU, event_wgt); 	  
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
