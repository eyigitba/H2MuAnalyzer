
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
#include "H2MuAnalyzer/MakeHistos/interface/MassCalPlots.h"  	  // class and config for mass calibration study
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"    // Common object selections
#include "H2MuAnalyzer/MakeHistos/interface/EventSelection.h"     // Common event selections
#include "H2MuAnalyzer/MakeHistos/interface/EventWeight.h"        // Common event weights
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"       // Common category definitions
#include "H2MuAnalyzer/MakeHistos/interface/MiniNTupleHelper.h"   // "PlantTree" and "BookBranch" functions
#include "H2MuAnalyzer/MakeHistos/interface/ReadMVA.h"            // Read and evaluate XMLs for MVA

// #include "H2MuAnalyzer/MakeHistos/interface/SampleDatabase2016.h" // Input data and MC samples

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 10;     // Maximum index of input files to process
const int MAX_EVT  = 10000; // Maximum number of events to process
const int PRT_EVT  = 1000;  // Print every N events
const float SAMP_WGT = 1.0;
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information


const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/prod-v18-pre-tag/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZJets_MG_1/190521_174140/0000";
const TString SAMPLE   = "ZJets_MG_1";
//const TString IN_DIR   = "/eos/cms/store/user/bortigno/h2mm/ntuples/2016/94X_v3/STR/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/190625_204843/0000";
//const TString SAMPLE   = "ZJets_AMC";
//const TString IN_DIR   = "/eos/cms/store/user/bortigno/h2mm/ntuples/2016/94X_v3/STR/SingleMuon/SingleMu_2016B/190625_203851/0000";
//const TString SAMPLE   = "SingleMu_2016B";

const std::string YEAR  = "2017";
const std::string SLIM  = "notSlim";  // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR   = "plots";
const TString HIST_TREE = "HistTree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both

const std::vector<std::string> SEL_CUTS = {"Presel2017"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE", "inclusive_01jet", "inclusive_2jets"}; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void MassCalibration( TString sample = "", TString in_dir = "", TString out_dir = "",
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

  // Add trees from the input files to the TChain
  TChain * in_chain = new TChain("dimuons/tree");
  for (int i = 0; i < in_file_names.size(); i++) {
    in_chain->Add( in_file_names.at(i) );
  // for (int i = 0; i < samp->filenames.size(); i++) {
  //   in_chain->Add( samp->filenames.at(i) );
    // Set branch addresses, from interface/LoadNTupleBranches.h
    if (sample.Contains("SingleMu"))
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM}, false); // Options in {} include "JES", "Flags", and "SFs"
    else
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM, "GEN", "Wgts"}, false); // Options in {} include "JES", "Flags", and "SFs"
  }

  // Initialize empty map of histogram names to output histograms
  std::map<TString, TH1*> h_map_1D;
  std::map<TString, TH2*> h_map_2D;

  gROOT->cd(); // Navigate to "local" memory, so all histograms are not saved in out_tuple

  // Configuration for object selection, event selection, and object weighting
  ObjectSelectionConfig obj_sel;
  EventSelectionConfig  evt_sel;
  EventWeightConfig     evt_wgt;
  ConfigureObjectSelection(obj_sel, YEAR);
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

  evt_sel.muPair_mass_min = 12; // Allow masses down to 12 GeV (instead of 60 GeV) for background studies

  if (verbose) obj_sel.Print();
  if (verbose) evt_sel.Print();
  if (verbose) evt_wgt.Print();

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later


  /////////////////////////////////////////
  ////   Config for mass calibration   ////
  /////////////////////////////////////////
  MassCalConfig 	mc_cfg;
  ConfigureMassCal( mc_cfg, "Z", "d0_diff");

  std::map<TString, MassCalPlots*> mc_map_PF;
  std::map<TString, MassCalPlots*> mc_map_Roch;
  std::map<TString, MassCalPlots*> mc_map_Kinfit;
  std::map<TString, MassCalPlots*> mc_map_KinRoch;
  std::map<TString, MassCalPlots*> mc_map_Kingood;
  std::map<TString, MassCalPlots*> mc_map_Kind0;

  std::map<TString, MassCalPlots*> mc_map_Kind0_BB;
  std::map<TString, MassCalPlots*> mc_map_Kind0_BE;
  std::map<TString, MassCalPlots*> mc_map_Kind0_EE;

//  std::map<TString, MassCalPlots*> mc_map_Kind0_N50_N15;
//  std::map<TString, MassCalPlots*> mc_map_Kind0_N15_N05;
//  std::map<TString, MassCalPlots*> mc_map_Kind0_P05_P15;
//  std::map<TString, MassCalPlots*> mc_map_Kind0_P15_P50;

  for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
    for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	TString optCatStr = OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat);
	mc_map_PF     [ optCatStr ]   = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_PF", mc_cfg );
	mc_map_Roch   [ optCatStr ]   = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Roch", mc_cfg );
	mc_map_Kinfit [ optCatStr ]   = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kinfit", mc_cfg );
        mc_map_KinRoch[ optCatStr ]   = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_KinRoch", mc_cfg );
	mc_map_Kingood[ optCatStr ]   = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_good_Kinfit", mc_cfg );
	mc_map_Kind0  [ optCatStr ]   = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin", mc_cfg );

	mc_map_Kind0_BB[ optCatStr ]  = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin_BB", mc_cfg );
	mc_map_Kind0_BE[ optCatStr ]  = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin_BE", mc_cfg );
	mc_map_Kind0_EE[ optCatStr ]  = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin_EE", mc_cfg );

//	mc_map_Kind0_N50_N15[ optCatStr ]  = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_N50_N15", mc_cfg );
//        mc_map_Kind0_N15_N05[ optCatStr ]  = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_N15_N05", mc_cfg );
//        mc_map_Kind0_P05_P15[ optCatStr ]  = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_P05_P15", mc_cfg );
//	mc_map_Kind0_P15_P50[ optCatStr ]  = new MassCalPlots( (sample.Contains("SingleMu") ? "data" :  sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_P15_P50", mc_cfg );

    }
  }

  std::cout << "\n******* About to enter the loop over " << in_chain->GetEntries() << " events *******" << std::endl;
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
    if (not pass_sel_cuts) continue;

    // Get event weight for MC, defined in src/EventWeight.cc
    bool isData = sample.Contains("SingleMu");
    float event_wgt = ( isData ? 1.0 : EventWeight(br, evt_wgt, verbose) );

      /////////////////////////////////////////////////////////////////////////////////////////
      ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
      /////////////////////////////////////////////////////////////////////////////////////////

    for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
      std::string OPT_CUT = OPT_CUTS.at(iOpt);
      
      for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	TString optCatStr = OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat);
	if ( not InCategory(obj_sel, br, CAT_CUTS.at(iCat), verbose) ) continue;

	if ( SelectedMuPairs(obj_sel, br).size()!=1 or SelectedMuons(obj_sel, br).size()!=2 ) continue;
	MuPairInfo    dimu;
	dimu = SelectedCandPair(obj_sel, br);

	MuonInfo mu1, mu2, muP, muN;
	mu1 = br.muons->at(dimu.iMu1);
	mu2 = br.muons->at(dimu.iMu2);
	float dimu_mass_KinRoch = dimu.mass_kinfit * dimu.mass_Roch / dimu.mass;
	float dimu_pt_KinRoch	= dimu.pt_kinfit * dimu.pt_Roch / dimu.pt;
 	float mu1_pt_KinRoch	= mu1.pt_kinfit * mu1.pt_Roch / mu1.pt;
	float mu2_pt_KinRoch    = mu2.pt_kinfit * mu2.pt_Roch / mu2.pt;

	if (mu1.charge == 1)  muP = mu1, muN = mu2;
 	else muP = mu2, muN = mu1;	
 	float d0_diff = (muP.d0_PV - muN.d0_PV) / 2.0;
	float d0_mean = (mu1.d0_PV + mu2.d0_PV) / 2.0;

	mc_map_PF     [optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_PF",      dimu.mass, 	 d0_diff, event_wgt, false);
	mc_map_Roch   [optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Roch",    dimu.mass_Roch,    d0_diff, event_wgt, false);
  	mc_map_Kinfit [optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kinfit",  dimu.mass_kinfit,  d0_diff, event_wgt, false);
	mc_map_KinRoch[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_KinRoch", dimu_mass_KinRoch, d0_diff, event_wgt, false);
	if (mu1.d0_PV_kinfit > -1 and mu2.d0_PV_kinfit > -1) {
	  mc_map_Kingood[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_good_Kinfit", dimu.mass_kinfit, d0_diff, event_wgt, false);
	  float d0kin_diff = (muP.d0_PV_kinfit - muN.d0_PV_kinfit) / 2.0;
	  float d0kin_mean = (muP.d0_PV_kinfit + muN.d0_PV_kinfit) / 2.0;
	  mc_map_Kind0[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin", dimu.mass_kinfit, d0kin_diff, event_wgt, false);

	  if (abs(mu1.eta) < 0.9 and abs(mu2.eta) < 0.9) mc_map_Kind0_BB[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_BB", dimu_mass_KinRoch, d0kin_diff, event_wgt, false);
	  if (abs(mu1.eta) > 1.2 and abs(mu2.eta) > 1.2) mc_map_Kind0_EE[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_EE", dimu_mass_KinRoch, d0kin_diff, event_wgt, false);
	  if (abs(mu1.eta) < 0.9 and abs(mu2.eta) > 1.2) mc_map_Kind0_BE[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_BE", dimu_mass_KinRoch, d0kin_diff, event_wgt, false);
	  if (abs(mu1.eta) > 1.2 and abs(mu2.eta) < 0.9) mc_map_Kind0_BE[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_BE", dimu_mass_KinRoch, d0kin_diff, event_wgt, false);

//	  if ( d0_diff > -0.005 and d0_diff < -0.0015 )       mc_map_Kind0_N50_N15[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_N50_N15", dimu.mass_kinfit, d0kin_diff, event_wgt, false);
//	  else if ( d0_diff > -0.0015 and d0_diff < -0.0005 ) mc_map_Kind0_N15_N05[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_N15_N05", dimu.mass_kinfit, d0kin_diff, event_wgt, false);
//	  else if ( d0_diff > 0.0005 and d0_diff < 0.0015 )   mc_map_Kind0_P05_P15[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_P05_P15", dimu.mass_kinfit, d0kin_diff, event_wgt, false);
//	  else if ( d0_diff > 0.0015 and d0_diff < 0.005 )    mc_map_Kind0_P15_P50[optCatStr]->FillEvent( (sample.Contains("SingleMu") ? "data" : sample) + "_" + optCatStr + "_Kin_vs_d0kin_d0PV_P15_P50", dimu.mass_kinfit, d0kin_diff, event_wgt, false);
	}

      } // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
    } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)

  } // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;


  // Place histograms made with separate selection and cateogry cuts into separate files
  for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
    for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
      TString optCatStr = OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat);
      
      // Create output file
      TString out_file_name;
      if (out_file_str.Length() > 0) out_file_name.Form( "%s/MassCal_%s_%s_%s.root",    out_dir.Data(), sample.Data(), optCatStr.Data(), out_file_str.Data() );
      else                           out_file_name.Form( "%s/MassCal_%s_%s_%d_%d.root", out_dir.Data(), sample.Data(), optCatStr.Data(), MIN_FILE, MAX_FILE );
      std::cout << "\nCreating output file " << out_file_name.Data() << std::endl;
      TFile* out_file = TFile::Open( out_file_name, "RECREATE" );
      
      // Write output file
      if (verbose) std::cout << "\nWriting output file " << out_file_name.Data() << std::endl;
      out_file->cd();

      TDirectory* PF_dir = out_file->mkdir("PF");
      TDirectory* Roch_dir = out_file->mkdir("Roch");
      TDirectory* Kinfit_dir = out_file->mkdir("Kinfit");
      TDirectory* KinRoch_dir = out_file->mkdir("KinRoch");
      TDirectory* Kingood_dir = out_file->mkdir("good_Kinfit");    
      TDirectory* Kind0_dir = out_file->mkdir("Kin_vs_d0kin");

      TDirectory* Kind0_BB_dir = out_file->mkdir("Kin_vs_d0kin_BB");
      TDirectory* Kind0_BE_dir = out_file->mkdir("Kin_vs_d0kin_BE");
      TDirectory* Kind0_EE_dir = out_file->mkdir("Kin_vs_d0kin_EE");

//      TDirectory* Kind0_N50_N15_dir  = out_file->mkdir("Kin_vs_d0kin_d0PV_N50_N15");
//      TDirectory* Kind0_N15_N05_dir  = out_file->mkdir("Kin_vs_d0kin_d0PV_N15_N05");
//      TDirectory* Kind0_P05_P15_dir  = out_file->mkdir("Kin_vs_d0kin_d0PV_P05_P15");
//      TDirectory* Kind0_P15_P50_dir  = out_file->mkdir("Kin_vs_d0kin_d0PV_P15_P50");
 
      MassCalPlots* mc_PF	 = mc_map_PF[optCatStr];
      MassCalPlots* mc_Roch	 = mc_map_Roch[optCatStr];
      MassCalPlots* mc_Kinfit	 = mc_map_Kinfit[optCatStr];
      MassCalPlots* mc_KinRoch   = mc_map_KinRoch[optCatStr]; 
      MassCalPlots* mc_Kingood   = mc_map_Kingood[optCatStr]; 
      MassCalPlots* mc_Kind0     = mc_map_Kind0[optCatStr];

      MassCalPlots* mc_Kind0_BB     = mc_map_Kind0_BB[optCatStr];
      MassCalPlots* mc_Kind0_BE     = mc_map_Kind0_BE[optCatStr];
      MassCalPlots* mc_Kind0_EE     = mc_map_Kind0_EE[optCatStr];

//      MassCalPlots* mc_Kind0_N50_N15    = mc_map_Kind0_N50_N15[optCatStr];
//      MassCalPlots* mc_Kind0_N15_N05    = mc_map_Kind0_N15_N05[optCatStr];
//      MassCalPlots* mc_Kind0_P05_P15    = mc_map_Kind0_P05_P15[optCatStr];
//      MassCalPlots* mc_Kind0_P15_P50    = mc_map_Kind0_P15_P50[optCatStr];

      PF_dir->cd();
      mc_PF->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_PF->mass_plots.begin(); it != mc_PF->mass_plots.end(); ++it) {
	std::string h_name = it->second->GetName();
	h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
	it->second->SetName(h_name.c_str());
	it->second->SetTitle(h_name.c_str());
  	it->second->Write();
      }

      Roch_dir->cd();
      mc_Roch->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_Roch->mass_plots.begin(); it != mc_Roch->mass_plots.end(); ++it) {
	std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
  	it->second->Write(); 
      }

      Kinfit_dir->cd();
      mc_Kinfit->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_Kinfit->mass_plots.begin(); it != mc_Kinfit->mass_plots.end(); ++it) {
	std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
  	it->second->Write();
      }

      KinRoch_dir->cd();
      mc_KinRoch->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_KinRoch->mass_plots.begin(); it != mc_KinRoch->mass_plots.end(); ++it) {
	std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
   	it->second->Write();
      }

      Kingood_dir->cd();
      mc_Kingood->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_Kingood->mass_plots.begin(); it != mc_Kingood->mass_plots.end(); ++it) {
        std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
        it->second->Write();
      }

      Kind0_dir->cd();
      mc_Kind0->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_Kind0->mass_plots.begin(); it != mc_Kind0->mass_plots.end(); ++it) {
        std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
        it->second->Write();
      }



      Kind0_BB_dir->cd();
      mc_Kind0_BB->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_Kind0_BB->mass_plots.begin(); it != mc_Kind0_BB->mass_plots.end(); ++it) {
        std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
        it->second->Write();
      }

      Kind0_BE_dir->cd();
      mc_Kind0_BE->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_Kind0_BE->mass_plots.begin(); it != mc_Kind0_BE->mass_plots.end(); ++it) {
        std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
        it->second->Write();
      }

      Kind0_EE_dir->cd();
      mc_Kind0_EE->summary_plot_1D->Write();
      for (std::map<TString, TH1*>::iterator it = mc_Kind0_EE->mass_plots.begin(); it != mc_Kind0_EE->mass_plots.end(); ++it) {
        std::string h_name = it->second->GetName();
        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
        it->second->SetName(h_name.c_str());
        it->second->SetTitle(h_name.c_str());
        it->second->Write();
      }

      
//      Kind0_N50_N15_dir->cd();
//      mc_Kind0_N50_N15->summary_plot_1D->Write();
//      for (std::map<TString, TH1*>::iterator it = mc_Kind0_N50_N15->mass_plots.begin(); it != mc_Kind0_N50_N15->mass_plots.end(); ++it) {
//        std::string h_name = it->second->GetName();
//        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
//        it->second->SetName(h_name.c_str());
//        it->second->SetTitle(h_name.c_str());
//        it->second->Write();
//      }
//
//      Kind0_N15_N05_dir->cd();
//      mc_Kind0_N15_N05->summary_plot_1D->Write();
//      for (std::map<TString, TH1*>::iterator it = mc_Kind0_N15_N05->mass_plots.begin(); it != mc_Kind0_N15_N05->mass_plots.end(); ++it) {
//        std::string h_name = it->second->GetName();
//        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
//        it->second->SetName(h_name.c_str());
//        it->second->SetTitle(h_name.c_str());
//        it->second->Write();
//      }
//
//      Kind0_P05_P15_dir->cd();
//      mc_Kind0_P05_P15->summary_plot_1D->Write();
//      for (std::map<TString, TH1*>::iterator it = mc_Kind0_P05_P15->mass_plots.begin(); it != mc_Kind0_P05_P15->mass_plots.end(); ++it) {
//        std::string h_name = it->second->GetName();
//        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
//        it->second->SetName(h_name.c_str());
//        it->second->SetTitle(h_name.c_str());
//        it->second->Write();
//      }
//
//      Kind0_P15_P50_dir->cd();
//      mc_Kind0_P15_P50->summary_plot_1D->Write();
//      for (std::map<TString, TH1*>::iterator it = mc_Kind0_P15_P50->mass_plots.begin(); it != mc_Kind0_P15_P50->mass_plots.end(); ++it) {
//        std::string h_name = it->second->GetName();
//        h_name.erase( h_name.find(optCatStr+"_"), optCatStr.Length() + 1 );
//        it->second->SetName(h_name.c_str());
//        it->second->SetTitle(h_name.c_str());
//        it->second->Write();
//      }


      out_file->Write();
      out_file->Close();
      std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
      
    } // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
  } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)


  std::cout << "\nExiting MassCalibration()\n";
  
} // End void MassCalibration()
