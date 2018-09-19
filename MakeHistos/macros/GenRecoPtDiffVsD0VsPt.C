
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
// Options passed in as arguments to GenRecoPtDiffVsD0VsPt when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 1;     // Maximum index of input files to process
const int MAX_EVT  = 5000;  // Maximum number of events to process
const int PRT_EVT  = 100;   // Print every N events
const float SAMP_WGT = 1.0;
const bool verbose = false; // Print extra information

const TString IN_DIR   = "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ZJets_MG/170313_225726/0000";
const TString SAMPLE   = "ZJets_MG";
// const TString IN_DIR   = "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/180802_165055/0000";
// const TString SAMPLE   = "ZJets_AMC";
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"pt_20_35", "pt_35_42", "pt_42_50", "pt_50_inf"};  // 4 muon pT bins: 20 - 35, 35 - 42.5, 42.5 - 50, > 50 GeV

// Command-line options for running in batch.  Running "root -b -l -q macros/GenRecoPtDiffVsD0VsPt.C" will use hard-coded options above.
void GenRecoPtDiffVsD0VsPt( TString sample = "", TString in_dir = "", TString out_dir = "",
			    std::vector<TString> in_files = {}, TString out_file_str = "",
			    int max_evt = 0, int prt_evt = 0, float samp_weight = 1.0) {
  
  // Set variables to hard-coded values if they are not initialized
  if (sample.Length()  == 0) sample  = SAMPLE;
  if (in_dir.Length()  == 0) in_dir  = IN_DIR;
  if (out_dir.Length() == 0) out_dir = OUT_DIR;
  if (max_evt          == 0) max_evt = MAX_EVT;
  if (prt_evt          == 0) prt_evt = PRT_EVT;
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
    
    int iGen1 = -99; // Higher-pT GEN muon
    int iGen2 = -99; // Lower-pT GEN muon
    int iMu1  = -99; // Higher-pT RECO muon
    int iMu2  = -99; // Lower-pT RECO muon
    
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
      
      for (int iMuon = 0; iMuon < 2; iMuon++) {
	
	GenMuonInfo gen = (iMuon == 0 ? br.genMuons->at(iGen1) : br.genMuons->at(iGen2) );
	MuonInfo mu     = (iMuon == 0 ? br.muons->at(iMu1)     : br.muons->at(iMu2)     );
	
	// Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc
	for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	  // if (not PassSelection(br, OPT_CUTS.at(iOpt), verbose)) continue;  // Not using centrally-defined selection cuts
	  
	  BookAndFill(h_map_1D, "h_"+OPT_CUTS.at(iOpt)+"_d0", 10000, -500,  500, 10000*mu.d0_PV*mu.charge );
	  BookAndFill(h_map_1D, "h_"+OPT_CUTS.at(iOpt)+"_pt", 10000,    0, 1000, mu.pt );
	
	  // Loop through category cuts defined in src/CategoryCuts.cc
	  for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	    // if (not InCategory(br, CAT_CUTS.at(iCat), verbose)) continue;  // Not using centrally-defined category cuts

	    // See if mu1 or mu2 falls inside the d0 bin
	    bool mu_inCat = false;
	    // Loop over pT ranges (each covers ~1/4 of the Z --> mu-mu phase space)
	    for (int iPt = 0; iPt < CAT_CUTS.size(); iPt++) {
	      double pt[2] = {20, 999};
	      if        (CAT_CUTS.at(iCat).find("pt_20_35") != string::npos) {
		pt[0] = 20;
		pt[1] = 35;
	      } else if (CAT_CUTS.at(iCat).find("pt_35_42") != string::npos) {
		pt[0] = 35;
		pt[1] = 42.5;
	      } else if (CAT_CUTS.at(iCat).find("pt_42_50") != string::npos) {
		pt[0] = 42.5;
		pt[1] = 50;
	      } else if (CAT_CUTS.at(iCat).find("pt_50_inf") != string::npos) {
		pt[0] = 50;
		pt[1] = 999;
	      }

	      // See if muon passes the pT cuts
	      if (pt[0] < mu.pt && mu.pt < pt[1])
		mu_inCat = true;
	    } // End loop: for (int iPt = 0; iPt < CAT_CUTS.size(); iPt++)
	  
	    if (not mu_inCat) continue;
	  
	    if (verbose) std::cout << "\nPassed cut " << OPT_CUTS.at(iOpt) << ", in category " << CAT_CUTS.at(iCat) << std::endl;
	    std::string h_pre = "h_"+OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)+"_";
	  
	    /////////////////////////////////
	    ///  Generate and fill plots  ///
	    /////////////////////////////////

	    // d0 range of +/-50 microns covers ~99% of Z --> mu-mu events

	    if (verbose) std::cout << "Gen pT = " << gen.pt << ", eta = " << gen.eta << ", phi = " << gen.phi << std::endl;
	    if (verbose) std::cout << "Mu pT = " << mu.pt << ", eta = " << mu.eta << ", phi = " << mu.phi << ", d0 = " << mu.d0_PV << ", charge = " << mu.charge << std::endl;

	    BookAndFill( h_map_1D, h_pre+"d0", 10000, -500,  500, 10000*mu.d0_PV*mu.charge );
	    BookAndFill( h_map_1D, h_pre+"pt", 10000,    0, 1000, mu.pt );
	    
	    BookAndFill( h_map_2D, h_pre+"dPt_Roch_vs_d0",       1000, -50, 50, 400,  -4,  4, 10000*mu.d0_PV*mu.charge,        mu.pt_Roch - gen.pt );
	    BookAndFill( h_map_2D, h_pre+"dPt_PF_vs_d0",         1000, -50, 50, 400,  -4,  4, 10000*mu.d0_PV*mu.charge,        mu.pt      - gen.pt );
	    BookAndFill( h_map_2D, h_pre+"dRelPt_Roch_vs_d0",    1000, -50, 50, 400,  -4,  4, 10000*mu.d0_PV*mu.charge,   100*(mu.pt_Roch - gen.pt) / gen.pt );
	    BookAndFill( h_map_2D, h_pre+"dRelPt_PF_vs_d0",      1000, -50, 50, 400,  -4,  4, 10000*mu.d0_PV*mu.charge,   100*(mu.pt      - gen.pt) / gen.pt );

	    BookAndFill( h_map_2D, h_pre+"dRelPt1p4_Roch_vs_d0", 1000, -50, 50, 400, -64, 64, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt_Roch - gen.pt) / pow(gen.pt, 1.4) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt1p4_PF_vs_d0",   1000, -50, 50, 400, -64, 64, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt      - gen.pt) / pow(gen.pt, 1.4) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt1p6_Roch_vs_d0", 1000, -50, 50, 400, -32, 32, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt_Roch - gen.pt) / pow(gen.pt, 1.6) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt1p6_PF_vs_d0",   1000, -50, 50, 400, -32, 32, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt      - gen.pt) / pow(gen.pt, 1.6) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt1p8_Roch_vs_d0", 1000, -50, 50, 400, -16, 16, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt_Roch - gen.pt) / pow(gen.pt, 1.8) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt1p8_PF_vs_d0",   1000, -50, 50, 400, -16, 16, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt      - gen.pt) / pow(gen.pt, 1.8) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p0_Roch_vs_d0", 1000, -50, 50, 400,  -8,  8, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt_Roch - gen.pt) / pow(gen.pt, 2.0) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p0_PF_vs_d0",   1000, -50, 50, 400,  -8,  8, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt      - gen.pt) / pow(gen.pt, 2.0) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p2_Roch_vs_d0", 1000, -50, 50, 400,  -4,  4, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt_Roch - gen.pt) / pow(gen.pt, 2.2) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p2_PF_vs_d0",   1000, -50, 50, 400,  -4,  4, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt      - gen.pt) / pow(gen.pt, 2.2) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p4_Roch_vs_d0", 1000, -50, 50, 400,  -2,  2, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt_Roch - gen.pt) / pow(gen.pt, 2.4) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p4_PF_vs_d0",   1000, -50, 50, 400,  -2,  2, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt      - gen.pt) / pow(gen.pt, 2.4) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p6_Roch_vs_d0", 1000, -50, 50, 400,  -1,  1, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt_Roch - gen.pt) / pow(gen.pt, 2.6) );
	    BookAndFill( h_map_2D, h_pre+"dRelPt2p6_PF_vs_d0",   1000, -50, 50, 400,  -1,  1, 10000*mu.d0_PV*mu.charge, 10000*(mu.pt      - gen.pt) / pow(gen.pt, 2.6) );

	    int eta_sign = gen.eta > 0 ? 1 : -1;
	    BookAndFill( h_map_2D, h_pre+"dPhi_vs_d0", 1000, -50, 50, 400, -1, 1, 10000*mu.d0_PV*mu.charge, 1000*(mu.phi - gen.phi)*mu.charge );
	    BookAndFill( h_map_2D, h_pre+"dEta_vs_d0", 1000, -50, 50, 400, -1, 1, 10000*mu.d0_PV*mu.charge, 1000*(mu.eta - gen.eta)*mu.charge*eta_sign );

	  } // End loop: for (int iMuon = 0; iMuon < 2; iMuon++)
	    
	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
    } // End conditional: if (pass_sel_cuts)

    // Before we exit the loop over all events, copy maps of histograms to local memory
    // *** Super-hacky!!! Should instead pass h_map_1D_loc to BookAndFill function! *** - AWB 20.08.2018
//    if (iEvt+1 == in_chain->GetEntries() || (iEvt == max_evt && max_evt > 0)) {
//      std::cout << "\nAt event " << iEvt << ", copying maps to local area" << std::endl;
//      RetreiveMap(h_map_1D_loc);
//      RetreiveMap(h_map_2D_loc);
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
    for (std::map<TString, TH1*>::iterator it = h_map_1D.begin(); it != h_map_1D.end(); ++it) {
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
    for (std::map<TString, TH2*>::iterator it = h_map_2D.begin(); it != h_map_2D.end(); ++it) {
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
      
  std::cout << "\nExiting GenRecoPtDiffVsD0VsPt()\n";
  
} // End void GenRecoPtDiffVsD0VsPt()
