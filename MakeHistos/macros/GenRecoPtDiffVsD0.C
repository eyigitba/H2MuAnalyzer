
///////////////////////////////////////////////////////////
///  Macro to study dependence of RECO - GEN pT vs. d0  ///
///                                                     ///
///            Andrew Brinkerhoff 24.08.2018            ///
///////////////////////////////////////////////////////////

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

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to GenRecoPtDiffVsD0 when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 1;     // Maximum index of input files to process
const int MAX_EVT  = 5000;  // Maximum number of events to process
const int PRT_EVT  = 100;   // Print every N events
const bool verbose = false; // Print extra information

const TString IN_DIR   = "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/180802_165055/0000";
const TString SAMPLE   = "ZJets_AMC";
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"d0_m9","d0_m8","d0_m7","d0_m6","d0_m5","d0_m4","d0_m3","d0_m2","d0_m1","d0_0", // Event selection categories
					   "d0_p1","d0_p2","d0_p3","d0_p4","d0_p5","d0_p6","d0_p7","d0_p8","d0_p9"};       // Also applied in parallel

// Command-line options for running in batch.  Running "root -b -l -q macros/GenRecoPtDiffVsD0.C" will use hard-coded options above.
void GenRecoPtDiffVsD0( TString sample = "", TString in_dir = "", TString out_dir = "",
			std::vector<TString> in_files = {}, TString out_file_str = "",
			int max_evt = 0, int prt_evt = 0) {

  // Set variables to hard-coded values if they are not initialized
  if (sample.Length()  == 0) sample  = SAMPLE;
  if (in_dir.Length()  == 0) in_dir  = IN_DIR;
  if (out_dir.Length() == 0) out_dir = OUT_DIR;
  if (max_evt          == 0) max_evt = MAX_EVT;
  if (prt_evt          == 0) prt_evt = PRT_EVT;

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
  std::map<TString, TH1*> h_map_1D_loc;
  std::map<TString, TH2*> h_map_2D_loc;

  // Add trees from the input files to the TChain
  TChain * in_chain = new TChain("dimuons/tree");
  if (verbose) std::cout << "Adding branches to tree" << std::endl;
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    in_chain->Add( in_file_names.at(i) );

    // Set branch addresses, from interface/LoadNTupleBranches.h
    SetBranchAddresses(*in_chain, br, {"GEN"}, false); // Options in {} include "GEN", "JES", "Flags", and "SFs"
  }
  if (verbose) std::cout << "Finished with SetBranchAddress calls" << std::endl;

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


    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Custom selection: require exactly one GEN Z --> mu-mu pair, matched to a pair of RECO muons //
    /////////////////////////////////////////////////////////////////////////////////////////////////

    pass_sel_cuts = false;

    int iGen1 = -99; // Positive GEN muon
    int iGen2 = -99; // Negative GEN muon
    int iMu1  = -99; // Positive RECO muon
    int iMu2  = -99; // Negative RECO muon

    for (UInt_t iGen = 0; iGen < br.nGenMuPairs; iGen++) {
      if (br.genMuPairs->at(iGen).mother_ID != 23) continue;  // Require a Z boson
      if (br.genMuPairs->at(iGen).postFSR   !=  0) continue;  // Require pre-FSR pair
      int iGen1_tmp = br.genMuPairs->at(iGen).iMu1;
      int iGen2_tmp = br.genMuPairs->at(iGen).iMu2;

      iGen1 = ( br.genMuons->at(iGen1_tmp).pt > br.genMuons->at(iGen2_tmp).pt ? iGen1_tmp : iGen2_tmp );
      iGen2 = ( br.genMuons->at(iGen1_tmp).pt < br.genMuons->at(iGen2_tmp).pt ? iGen1_tmp : iGen2_tmp );

      for (UInt_t iMu = 0; iMu < br.nMuPairs; iMu++) {
	int iMu1_tmp = br.muPairs->at(iMu).iMu1;
	int iMu2_tmp = br.muPairs->at(iMu).iMu2;

	iMu1 = ( br.muons->at(iMu1_tmp).pt > br.muons->at(iMu2_tmp).pt ? iMu1_tmp : iMu2_tmp );
	iMu2 = ( br.muons->at(iMu1_tmp).pt < br.muons->at(iMu2_tmp).pt ? iMu1_tmp : iMu2_tmp );

	// Require that both RECO muons have a GEN match
	if ( abs(br.muons->at(iMu1).eta - br.genMuons->at(iGen1).eta) > 0.02 ||
	     abs(br.muons->at(iMu1).phi - br.genMuons->at(iGen1).phi) > 0.02 ||
	     br.muons->at(iMu1).charge != br.genMuons->at(iGen1).charge ) continue;
	if ( abs(br.muons->at(iMu2).eta - br.genMuons->at(iGen2).eta) > 0.02 ||
	     abs(br.muons->at(iMu2).phi - br.genMuons->at(iGen2).phi) > 0.02 ||
	     br.muons->at(iMu2).charge != br.genMuons->at(iGen2).charge ) continue;

	pass_sel_cuts = true;
	break;
      } // End loop: for (UInt_t iMu = 0; iMu < br.nMuPairs; iMu++) 
      if (pass_sel_cuts) break;
    } // End loop: for (UInt_t iGen = 0; iGen < br.nGenMuPairs; iGen++)

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // End custom selection: require exactly one GEN Z --> mu-mu pair, matched to a pair of RECO muons //
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    if (pass_sel_cuts) {

      GenMuonInfo gen1 = br.genMuons->at(iGen1);
      GenMuonInfo gen2 = br.genMuons->at(iGen2);
      MuonInfo mu1 = br.muons->at(iMu1);
      MuonInfo mu2 = br.muons->at(iMu2);

      // Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc
      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	// if (not PassSelection(br, OPT_CUTS.at(iOpt), verbose)) continue;  // Not using centrally-defined selection cuts

	BookAndFill( "h_"+OPT_CUTS.at(iOpt)+"_d0_mu1", 1000, -0.05, 0.05, mu1.d0_PV*mu1.charge );
	BookAndFill( "h_"+OPT_CUTS.at(iOpt)+"_d0_mu2", 1000, -0.05, 0.05, mu2.d0_PV*mu2.charge );
	
	// Loop through category cuts defined in src/CategoryCuts.cc
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	  // if (not InCategory(br, CAT_CUTS.at(iCat), verbose)) continue;  // Not using centrally-defined category cuts

	  // See if mu1 or mu2 falls inside the d0 bin
	  bool mu1_inCat = false;
	  bool mu2_inCat = false;
	  // Loop over possible central d0 values
	  for (int iD0 = 0; iD0 < ((CAT_CUTS.size() + 1) / 2); iD0++) {
	    double d0 = 0.001 * iD0; // Central d0 value in the bin
	    if ( (iD0 == 0 && CAT_CUTS.at(iCat).find("d0_0") != string::npos) || 
		 (iD0 != 0 && CAT_CUTS.at(iCat).find( std::to_string(iD0) ) != string::npos) ) {
	      // Switch central d0 value to negative if needed
	      if (CAT_CUTS.at(iCat).find("d0_m") != string::npos) d0 *= -1;
	      // See if mu1 or mu2 passes the d0 cuts
	      if ( (mu1.d0_PV * mu1.charge > d0 - 0.0005) && (mu1.d0_PV * mu1.charge < d0 + 0.0005) )
		mu1_inCat = true;
	      if ( (mu2.d0_PV * mu2.charge > d0 - 0.0005) && (mu2.d0_PV * mu2.charge < d0 + 0.0005) )
		mu2_inCat = true;
	    }
	  } // End loop: for (int iD0 = 0; iD0 < 7; iD0++)
	  
	  if (not (mu1_inCat || mu2_inCat) ) continue;
	  
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUTS.at(iOpt) << ", in category " << CAT_CUTS.at(iCat) << std::endl;
	  std::string h_pre = "h_"+OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)+"_";
	  
	  /////////////////////////////////
	  ///  Generate and fill plots  ///
	  /////////////////////////////////

	  if (mu1_inCat) {
	    if (verbose) std::cout << "Gen1 pT = " << gen1.pt << ", eta = " << gen1.eta << ", phi = " << gen1.phi << std::endl;
	    if (verbose) std::cout << "Mu1 pT = " << mu1.pt << ", eta = " << mu1.eta << ", phi = " << mu1.phi << ", d0 = " << mu1.d0_PV << ", charge = " << mu1.charge << std::endl;
	  
	    BookAndFill( h_pre+"dPt_Roch_mu1",     200,  -5,  5,        mu1.pt_Roch - gen1.pt );
	    BookAndFill( h_pre+"dPt_PF_mu1",       200,  -5,  5,        mu1.pt      - gen1.pt );
	    BookAndFill( h_pre+"dRelPt_Roch_mu1",  200, -10, 10,   100*(mu1.pt_Roch - gen1.pt) / gen1.pt );
	    BookAndFill( h_pre+"dRelPt_PF_mu1",    200, -10, 10,   100*(mu1.pt      - gen1.pt) / gen1.pt );
	    BookAndFill( h_pre+"dRelPt2_Roch_mu1", 200, -20, 20, 10000*(mu1.pt_Roch - gen1.pt) / pow(gen1.pt, 2) );
	    BookAndFill( h_pre+"dRelPt2_PF_mu1",   200, -20, 20, 10000*(mu1.pt      - gen1.pt) / pow(gen1.pt, 2) );
	    BookAndFill( h_pre+"dPhi_mu1", 200, -2, 2, 1000*(mu1.phi - gen1.phi)*mu1.charge );
	    BookAndFill( h_pre+"dEta_mu1", 200, -2, 2, 1000*(mu1.eta - gen1.eta) );
	    BookAndFill( h_pre+"d0_mu1", 200, -0.01, 0.01, mu1.d0_PV*mu1.charge );
	  }

	  if (mu2_inCat) {
	    if (verbose) std::cout << "Gen2 pT = " << gen2.pt << ", eta = " << gen2.eta << ", phi = " << gen2.phi << std::endl;
	    if (verbose) std::cout << "Mu2 pT = " << mu2.pt << ", eta = " << mu2.eta << ", phi = " << mu2.phi << ", d0 = " << mu2.d0_PV << ", charge = " << mu2.charge << std::endl;
	    
	    BookAndFill( h_pre+"dPt_Roch_mu2",     200,  -5,  5,        mu2.pt_Roch - gen2.pt );
	    BookAndFill( h_pre+"dPt_PF_mu2",       200,  -5,  5,        mu2.pt      - gen2.pt );
	    BookAndFill( h_pre+"dRelPt_Roch_mu2",  200, -10, 10,   100*(mu2.pt_Roch - gen2.pt) / gen2.pt );
	    BookAndFill( h_pre+"dRelPt_PF_mu2",    200, -10, 10,   100*(mu2.pt      - gen2.pt) / gen2.pt );
	    BookAndFill( h_pre+"dRelPt2_Roch_mu2", 200, -20, 20, 10000*(mu2.pt_Roch - gen2.pt) / pow(gen2.pt, 2) );
	    BookAndFill( h_pre+"dRelPt2_PF_mu2",   200, -20, 20, 10000*(mu2.pt      - gen2.pt) / pow(gen2.pt, 2) );
	    BookAndFill( h_pre+"dPhi_mu2", 200, -2, 2, 1000*(mu2.phi - gen2.phi)*mu2.charge );
	    BookAndFill( h_pre+"dEta_mu2", 200, -2, 2, 1000*(mu2.eta - gen2.eta) );
	    BookAndFill( h_pre+"d0_mu2", 200, -0.01, 0.01, mu2.d0_PV*mu2.charge );
	  }
	  
	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
    } // End conditional: if (pass_sel_cuts)

    // Before we exit the loop over all events, copy maps of histograms to local memory
    // *** Super-hacky!!! Should instead pass h_map_1D_loc to BookAndFill function! *** - AWB 20.08.2018
    if (iEvt+1 == in_chain->GetEntries() || (iEvt == max_evt && max_evt > 0)) {
      std::cout << "\nAt event " << iEvt << ", copying maps to local area" << std::endl;
      RetreiveMap(h_map_1D_loc);
      RetreiveMap(h_map_2D_loc);
    }
    
  } // End loop: for (UInt_t iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;

  // Concatenate basic selection cuts into string for output file name
  std::string selStr = "";
  for (int i = 0; i < SEL_CUTS.size(); i++) selStr = selStr+SEL_CUTS.at(i)+"_";
  selStr.pop_back();

  // Place histograms made with separate selection cuts into separate files
  for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
    std::string optStr = OPT_CUTS.at(iOpt);

    // Create output file
    TString out_file_name;
    if (out_file_str.Length() > 0) out_file_name.Form( "%s/histos_%s_%s_%s_%s.root",    out_dir.Data(), sample.Data(), selStr.c_str(), optStr.c_str(), out_file_str.Data() );
    else                           out_file_name.Form( "%s/histos_%s_%s_%s_%d_%d.root", out_dir.Data(), sample.Data(), selStr.c_str(), optStr.c_str(), MIN_FILE, MAX_FILE );
    std::cout << "\nCreating output file " << out_file_name.Data() << std::endl;
    TFile* out_file = TFile::Open( out_file_name, "RECREATE" );
    
    // Write output file
    if (verbose) std::cout << "\nWriting output file " << out_file_name.Data() << std::endl;
    out_file->cd();
    
    // Write out 1D histograms
    for (std::map<TString, TH1*>::iterator it = h_map_1D_loc.begin(); it != h_map_1D_loc.end(); ++it) {
      std::string h_name = it->second->GetName();
      if (h_name.find(optStr+"_") != std::string::npos) {
	// Remove optional selection and category cuts from histogram names
	h_name.erase( h_name.find(optStr+"_"), optStr.length() + 1 );
	it->second->SetName(h_name.c_str());
	it->second->SetTitle(h_name.c_str());
	std::cout << "  * Writing 1D histogram " << it->second->GetName() << std::endl;
	it->second->Write();
      }
    }
    // Write out 2D histograms
    for (std::map<TString, TH2*>::iterator it = h_map_2D_loc.begin(); it != h_map_2D_loc.end(); ++it) {
      std::string h_name = it->second->GetName();
      if (h_name.find(optStr+"_") != std::string::npos) {
	// Remove optional selection and category cuts from histogram names
	h_name.erase( h_name.find(optStr+"_"), optStr.length() + 1 );
	it->second->SetName(h_name.c_str());
	it->second->SetTitle(h_name.c_str());
	std::cout << "  * Writing 2D histogram " << it->second->GetName() << std::endl;
	it->second->Write();
      }
    }
    
    out_file->Write();
    std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
    
  } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
      
  std::cout << "\nExiting GenRecoPtDiffVsD0()\n";
  
} // End void GenRecoPtDiffVsD0()
