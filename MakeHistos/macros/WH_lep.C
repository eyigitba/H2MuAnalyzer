
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
#include "H2MuAnalyzer/MakeHistos/interface/MiniNTupleHelper.h"   // "PlantTree" and "BookBranch" functions
#include "H2MuAnalyzer/MakeHistos/interface/ReadMVA.h"            // Read and evaluate XMLs for MVA

// #include "H2MuAnalyzer/MakeHistos/interface/SampleDatabase2016.h" // Input data and MC samples

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 1;     // Maximum index of input files to process
const int MAX_EVT  = 1000; // Maximum number of events to process
const int PRT_EVT  = 100;  // Print every N events
const float SAMP_WGT = 1.0;
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information

const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/WplusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8/H2Mu_WH_pos_125";
const TString SAMPLE   = "H2Mu_WH_pos_125";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZJets_MG_1";
// const TString SAMPLE   = "ZJets_MG_1";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/WZ_3l";
// const TString SAMPLE   = "WZ_3l";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/SingleMuon/SingleMu_2017D";
// const TString SAMPLE   = "SingleMu";

const std::string YEAR  = "2017";
const std::string SLIM  = "notSlim";  // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR   = "plots";
const TString HIST_TREE = "HistTree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both

// Cuts which every event must pass, applied in sequence
const std::vector<std::string> SEL_CUTS = {"Presel2017"};
// Multiple selection cuts, applied independently in parallel
// const std::vector<std::string> OPT_CUTS = {"3mu", "3mu_allMass", "3mu_hiPt", "e2mu", "e2mu_allMass"};
const std::vector<std::string> OPT_CUTS = {"3lep", "3mu", "e2mu"};
// Category selection cuts, also applied in parallel
// *** IMPORTANT!!! No category name may start with a sub-string which is identical to another entire category name! ***
const std::vector<std::string> CAT_CUTS = { "looseLepMVA_noZ_noBtag_mass12",
					    "medLepMVA_noZ_noBtag_mass12" };


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void WH_lep( TString sample = "", TString in_dir = "", TString out_dir = "",
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
    // for (int i = MIN_FILE; i <= MAX_FILE; i++) {
    //   in_file_name.Form("%s/tuple_%d.root", in_dir.Data(), i);
    //   std::cout << "Adding file " << in_file_name.Data() << std::endl;
    //   in_file_names.push_back(in_file_name.Data());
    // }
    in_file_name.Form("%s/NTuple_0.root", in_dir.Data());
    std::cout << "Adding file " << in_file_name.Data() << std::endl;
    in_file_names.push_back(in_file_name.Data());
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

  // ///////////////////////////////////////////////////////
  // ///  Input samples: MC signal, MC background, data  ///
  // ///////////////////////////////////////////////////////
  
  // std::map<TString, Sample*> samples;
  // std::cout << "\nAbout to get the samples" << std::endl;

  // // Load samples into our map. "UF" if you are on the UF servers
  // // or "CERN" if you at CERN. "ALL" specifies that we want to load the Data
  // // and all of the MC samples. Can loop through and remove the ones you don't want
  // // to use if you desire or just grab the ones you care about from the map.
  // GetSamples2016(samples, "CERN_hiM", sample );
  // Sample * samp = samples[sample];

  // // Open all input files
  // for (int i = 0; i < samp->filenames.size(); i++) {
  //   if ( !gSystem->AccessPathName(samp->filenames.at(i)) )
  //     file_tmp = TFile::Open( samp->filenames.at(i) ); // Check if file exists
  //   if (!file_tmp) {
  //     std::cout << "ERROR: could not open data file " << samp->filenames.at(i) << std::endl;
  //     return;
  //   }
  // }


  // Initialize set of pointers to all branches in input tree
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
  // float lumi_SF = samp->getLumiScaleFactor(LUMI);
  // std::cout << "For LUMI = " << LUMI << ", lumi_SF = " << lumi_SF << std::endl;


  // Concatenate basic selection cuts into string for output file name
  std::string selStr = "";
  for (int i = 0; i < SEL_CUTS.size(); i++) selStr = selStr+SEL_CUTS.at(i)+"_";
  selStr.pop_back();

  // Create output file for ntuple
  TString out_tuple_name;
  if (out_file_str.Length() > 0)   out_tuple_name.Form( "%s/tuple_%s_%s_%s.root",    out_dir.Data(), sample.Data(), selStr.c_str(), out_file_str.Data() );
  else                             out_tuple_name.Form( "%s/tuple_%s_%s_%d_%d.root", out_dir.Data(), sample.Data(), selStr.c_str(), MIN_FILE, MAX_FILE );
  TFile * out_tuple = TFile::Open( out_tuple_name, "RECREATE" );
  TTree * out_tree  = PlantTree("tree", "tree");

  // Initialize empty map of histogram names to output histograms
  std::map<TString, TH1*> h_map_1D;
  std::map<TString, TH2*> h_map_2D;

  // Initialize empty map of branch names to ouput branches
  std::map<TString, float>       b_map_flt;
  std::map<TString, int>         b_map_int;
  std::map<TString, std::string> b_map_str;
  // Track the number of elements in each map (should be the same in each event)
  int n_map_flt = -99;
  int n_map_int = -99;
  int n_map_str = -99;

  gROOT->cd(); // Navigate to "local" memory, so all histograms are not saved in out_tuple


  // Configuration for object selection, event selection, and object weighting
  ObjectSelectionConfig obj_sel;
  EventSelectionConfig  evt_sel;
  EventWeightConfig     evt_wgt;
  ConfigureObjectSelection(obj_sel, YEAR, "lepMVA");
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

  evt_sel.muPair_mass_min = 12; // Allow masses down to 12 GeV (instead of 60 GeV) for background studies

  if (verbose) obj_sel.Print();
  if (verbose) evt_sel.Print();
  if (verbose) evt_wgt.Print();

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later

  std::cout << "\n******* About to load 2D LepMVA efficiency scale factor histograms *******" << std::endl;
  TH2F * lepSF_mu_T  = LoadSFsLepMVA(YEAR,  "mu", "T");
  TH2F * lepSF_mu_M  = LoadSFsLepMVA(YEAR,  "mu", "M");
  TH2F * lepSF_mu_L  = LoadSFsLepMVA(YEAR,  "mu", "L");
  TH2F * lepSF_ele_T = LoadSFsLepMVA(YEAR, "ele", "T");
  TH2F * lepSF_ele_M = LoadSFsLepMVA(YEAR, "ele", "M");
  TH2F * lepSF_ele_L = LoadSFsLepMVA(YEAR, "ele", "L");

  std::cout << "\n******* About to load XML files for signal-background BDTs *******" << std::endl;
  MVA::MVA BDT_XWZ_noMass( "data/XMLs/WH_3l/Xunwu/2019_05_12/2017_WH_ele_against_inclu_trimvar_all_sig_all_bkg_ge0j/",
		       "weights/2017_WH_ele_against_inclu_trimvar_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml",
		       "BDTG_UF_v2" );
  MVA::MVA BDT_XWZ_withMass  ( "data/XMLs/WH_3l/Xunwu/2019_05_12/2017_WH_ele_against_inclu_trimvar_with_mass_all_sig_all_bkg_ge0j/",
			       "weights/2017_WH_ele_against_inclu_trimvar_with_mass_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml",
			       "BDTG_UF_v2" );

  MVA::MVA BDT_AWB_2_noMass( "data/XMLs/WH_3l/Andrew/2019_05_15/f_Opt_AWB_noMass_v2_all_sig_all_bkg_ge0j/",
			     "weights/f_Opt_AWB_noMass_v2_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml",
			     "BDTG_UF_v1" );
  MVA::MVA BDT_AWB_2_withMass( "data/XMLs/WH_3l/Andrew/2019_05_15/f_Opt_AWB_withMass_v2_all_sig_all_bkg_ge0j/",
			       "weights/f_Opt_AWB_withMass_v2_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml",
			       "BDTG_UF_v1" );
  MVA::MVA BDT_AWB_2_retrain( "data/XMLs/WH_3l/Andrew/2019_05_15/f_Opt_noMassBDT_mass_v2_all_sig_all_bkg_ge0j/",
			      "weights/f_Opt_noMassBDT_mass_v2_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml",
			      "BDTG_UF_v1" );

  MVA::MVA BDT_AWB_3_noMass( "data/XMLs/WH_3l/Andrew/2019_05_20/f_Opt_AWB_noMass_v3_resWgt_all_sig_all_bkg_resWgt/",
			     "weights/f_Opt_AWB_noMass_v3_resWgt_all_sig_all_bkg_resWgt_BDTG_UF_v1.weights.xml",
			     "BDTG_UF_v1" );
  MVA::MVA BDT_AWB_3_withMass( "data/XMLs/WH_3l/Andrew/2019_05_20/f_Opt_AWB_withMass_v3_all_sig_all_bkg_/",
			       "weights/f_Opt_AWB_withMass_v3_all_sig_all_bkg__BDTG_UF_v1.weights.xml",
			       "BDTG_UF_v1" );
  MVA::MVA BDT_AWB_3_retrain( "data/XMLs/WH_3l/Andrew/2019_05_20/f_Opt_noMassBDT_mass_v3_all_sig_all_bkg_ge0j/",
			      "weights/f_Opt_noMassBDT_mass_v3_all_sig_all_bkg_ge0j_BDTG_UF_v1.weights.xml",
			      "BDTG_UF_v1" );

  std::cout << "\n******* About to enter the loop over " << in_chain->GetEntries() << " events *******" << std::endl;
  for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++) {
    
    if (iEvt > max_evt && max_evt > 0) break;
    if ( (iEvt % prt_evt) == 0 ) {
      std::cout << "\n*********************" << std::endl;
      std::cout << "Looking at event " << iEvt <<  std::endl;
      std::cout << "*********************" << std::endl;
    }
    
    if (verbose) std::cout << "Before running GetEntry, event = " << br.event->event;
    
    in_chain->GetEntry(iEvt);
    
    if (verbose) std::cout << "... after, event = " << br.event->event << std::endl;

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

    // Fill a map of LepMVA weights corresponding to each OPT+CAT category
    std::map<TString, float> lepMVA_wgts;
    for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
      for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	lepMVA_wgts[OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)] = 1.0;
      }
    }

    // Initialize the selected Higgs candidate dimuon pair
    MuPairInfo        H_pair;     // When filling histograms, have only one candidate pair at a time
    TLorentzVector    H_pair_vec;
      
    bool MU = false;  // Says whether event has 3 muons or 1 electron + 2 muons

    ///////////////////////////////////////////////////////////////
    ///  Loop through possible Higgs candidate pairs (up to 2)  ///
    ///////////////////////////////////////////////////////////////

    for (int iPair = 0; iPair < 2; iPair++) {

      // Each candidate pair defines a unique "event" - i.e. one collision "event" can sometimes enter the tree / histograms twice
      // For each "event", we set all the branch values to their default value (-99)
      for (std::map<TString, float>::iterator       it = b_map_flt.begin(); it != b_map_flt.end(); ++it) b_map_flt[it->first] = -99.0;
      for (std::map<TString, int>::iterator         it = b_map_int.begin(); it != b_map_int.end(); ++it) b_map_int[it->first] = -99;
      for (std::map<TString, std::string>::iterator it = b_map_str.begin(); it != b_map_str.end(); ++it) b_map_str[it->first] = "-99";


      /////////////////////////////////////////////////////////////////////////////////////////
      ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
      /////////////////////////////////////////////////////////////////////////////////////////

      bool pass_opt_cuts = false; // Check if event passes at least one optional cut
      bool pass_cat_cuts = false; // Check if event passes at least one category cut

      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	std::string OPT_CUT = OPT_CUTS.at(iOpt);
	ASSERT( OPT_CUT == "3lep" || OPT_CUT == "3mu" || OPT_CUT == "e2mu" || OPT_CUT == "3lep_allMass" || OPT_CUT == "3mu_allMass" || OPT_CUT == "e2mu_allMass" || OPT_CUT == "3lep_hiPt" || OPT_CUT == "3mu_hiPt",
	       "OPT_CUT == '3lep' || OPT_CUT == '3mu' || OPT_CUT == 'e2mu' || OPT_CUT == '3lep_allMass' || OPT_CUT == '3mu_allMass' || OPT_CUT == 'e2mu_allMass' || OPT_CUT == '3lep_hiPt' || OPT_CUT == '3mu_hiPt'" );

	// Choose Higgs mass (110 - 160 GeV) or all mass (> 12 GeV) requirements
	bool allMass = ( OPT_CUT.find("allMass") != std::string::npos );
	// Require Higgs candidate pair to have the vector pT value
	bool hiPt    = ( OPT_CUT.find("hiPt")    != std::string::npos );

	// Assume event fails the optional cut
	BookAndFill(b_map_int, out_tree, sample, "OPT_"+OPT_CUT, 0 );
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
	  // Assume event fails each optional x category cut as well
	  BookAndFill(b_map_int, out_tree, sample+"_"+OPT_CUT, "OPT_"+OPT_CUT+"_CAT_"+CAT_CUTS.at(iCat), 0 );
	  // Fill default weights for LepMVA scale factors
	  BookAndFill(b_map_flt, out_tree, sample+"_"+OPT_CUT, "lepMVA_wgt_OPT_"+OPT_CUT+"_CAT_"+CAT_CUTS.at(iCat), -99.0 );
	}

	// All categories require >=2 opposite-charge muons, and == 3 leptons with charge summing to +/-1
	if (true) {
	  if ( SelectedEles(obj_sel, br).size() == 0 ) {
	    if ( SelectedMuPairs(obj_sel, br).size() != 2 ) continue;
	    MU = true;
	  } else if ( SelectedEles(obj_sel, br).size() == 1 ) {
	    if ( SelectedMuPairs(obj_sel, br).size() != 1 ) continue;
	    if ( iPair != 0 ) continue;  // Only one Higgs pair possible, hence iPair must be 0
	    MU = false;
	  } else continue;

	  // Choose a unique Higgs candidate pair
	  H_pair     = SelectedMuPairs(obj_sel, br).at(iPair);
	  H_pair_vec = FourVec(H_pair, PTC);
	  if ( !allMass ) {
	    // Require selected Higgs candidate dimuon pair to fall inside mass window
	    if ( H_pair_vec.M() < 110 ||
		 H_pair_vec.M() > 160 ) continue;
	  } else if (MU) {
	    // Pick the candidate closer to the Z mass as the "H_pair"
	    int iZ = ( abs(FourVec( SelectedMuPairs(obj_sel, br).at(0), PTC ).M() - 91) <
		       abs(FourVec( SelectedMuPairs(obj_sel, br).at(1), PTC ).M() - 91) ? 0 : 1 );
	    if (iZ != iPair) continue;
	  }
	  if ( MU && hiPt ) {
	    // Require Higgs candidate pair to have the higher vector pT value
	    if ( H_pair_vec.Pt() < FourVec(SelectedMuPairs(obj_sel, br).at((iPair+1) % 2), PTC).Pt() ) continue;
	  }
	}

	// Special cuts for 3 muon category
	if (OPT_CUT == "3mu" || OPT_CUT == "3mu_allMass" || OPT_CUT == "3mu_hiPt") {
	  if (!MU) continue;
	}
	else if (OPT_CUT == "e2mu" || OPT_CUT == "e2mu_allMass") {
	  if (MU) continue;
	}
	pass_opt_cuts = true;

	// Now we know the event passes the optional cut
	b_map_int["OPT_"+OPT_CUT] = 1;

	
	//////////////////////////////////////////////////////////////////
	/// Loop through category cuts defined in src/CategoryCuts.cc  ///
	//////////////////////////////////////////////////////////////////
	
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {

	  //////////////////////////////////////////////////////
	  ///  Compute variables relevent for category cuts  ///
	  //////////////////////////////////////////////////////

	  // Get selected muons and electrons in event
	  MuonInfos muons = SelectedMuons(obj_sel, br);
	  EleInfos  eles  = SelectedEles(obj_sel, br);
	  JetInfos  jets  = SelectedJets(obj_sel, br);
	  int sum_lep_charge = 0;
	  if (MU) {
	    ASSERT(muons.size() == 3 && eles.size() == 0, "muons.size() == 3 && eles.size() == 0");
	    sum_lep_charge = muons.at(0).charge + muons.at(1).charge + muons.at(2).charge;
	  } else {
	    ASSERT(muons.size() == 2 && eles.size() == 1, "muons.size() == 2 && eles.size() == 1");
	    sum_lep_charge = muons.at(0).charge + muons.at(1).charge + eles.at(0).charge;
	  }
	  ASSERT(abs(sum_lep_charge) == 1, "abs(sum_lep_charge) == 1");
	  
	  MuPairInfo H_true;  // Dimuon pair matched to GEN-level Higgs decay
	  MuPairInfo Z_true;  // Dimuon pair matched to GEN-level Z boson decay
	  MuonInfo   muH1;    // Higher pT muon from chosen Higgs candidate pair
	  MuonInfo   muH2;    // Lower pT muon from chosen Higgs candidate pair
	  MuonInfo   muW;     // Muon not from chosen Higgs candidate pair
	  EleInfo    ele;     // Electron
	  MuonInfo   muSS1;   // Same-sign muon with higher pT
	  MuonInfo   muSS2;   // Same-sign muon with lower pT
	  MuonInfo   muOS;    // Opposite-sign muon

	  TLorentzVector H_GEN_vec;
	  TLorentzVector Z_GEN_vec;
	  TLorentzVector H_true_vec;
	  TLorentzVector Z_true_vec;
	  TLorentzVector muH1_vec;
	  TLorentzVector muH2_vec;
	  TLorentzVector lep_vec;      // Lepton not chosen to be in Higgs candidate pair
	  TLorentzVector lep_vecT;     // Transverse component only (eta = 0, pZ = 0, mass = 0)
	  TLorentzVector OS_pair_vec;  // OS lepton pair not chosen as Higgs candidate (mu-mu or ele-mu)
	  TLorentzVector lepSS1_vec;   // SS lepton with higher pT
	  TLorentzVector lepSS2_vec;   // SS lepton with lower pT
	  TLorentzVector muSS_vec;     // SS muon from the chosen Higgs candidate pair
	  TLorentzVector muSS_vecT;    // Transverse component only (eta = 0, pZ = 0, mass = 0)
	  TLorentzVector muOS_vec;     // OS muon
	  TLorentzVector muOS_vecT;    // Transverse component only (eta = 0, pZ = 0, mass = 0)

	  muH1     = br.muons->at(H_pair.iMu1);
	  muH2     = br.muons->at(H_pair.iMu2);
	  muH1_vec = FourVec(muH1, PTC);
	  muH2_vec = FourVec(muH2, PTC);

	  // Look for a GEN Higgs or Z boson
	  if (not isData) {
	    for (int i = 0; i < br.nGenParents; i++) {
	      if (abs((br.genParents->at(i)).ID) == 25) {
		if (FourVec(br.genParents->at(i)).M() > H_GEN_vec.M()) {
		  H_GEN_vec = FourVec(br.genParents->at(i));
		}
	      }
	      if (abs((br.genParents->at(i)).ID) == 23) {
		if (FourVec(br.genParents->at(i)).M() > Z_GEN_vec.M() &&
		    FourVec(br.genParents->at(i)).M() > 0) {
		  Z_GEN_vec = FourVec(br.genParents->at(i));
		}
	      }
	    }
	  }

	  // Loop over selected dimuon pairs
	  for (const auto & muPair : SelectedMuPairs(obj_sel, br)) {
	    // Check if the pair matches a GEN muon pair from H
	    if (not isData) {
	      if ( IsGenMatched( muPair, *br.muons, *br.genMuons, *br.genParents, "H") ) {
		H_true     = muPair;
		H_true_vec = FourVec(muPair, PTC);
	      }
	      if ( IsGenMatched( muPair, *br.muons, *br.genMuons, *br.genParents, "Z") ) {
		Z_true     = muPair;
		Z_true_vec = FourVec(muPair, PTC);
	      }
	    }
	    // Check if the pair does not match the selected Higgs candidate pair
	    if ( muPair.iMu1 != H_pair.iMu1 || muPair.iMu2 != H_pair.iMu2 ) {
	      OS_pair_vec = FourVec(muPair, PTC);
	    }
	  }

	  // For 3 muon events
	  if (MU) {
	    // Find the muon not in the H pair, and the same-sign and opposite-sign muons
	    for (const auto & mu : muons) {
	      // Find the W muon candidate and the same-sign muon from the Higgs
	      if ( mu.pt != muH1.pt && mu.eta != muH1.eta &&
		   mu.pt != muH2.pt && mu.eta != muH2.eta ) {
		ASSERT(muW.pt <= 0, "muW.pt <= 0"); // We should not have found a W candidate before
		muW = mu;
		lep_vec  = FourVec(mu, PTC);
		lep_vecT = FourVec(mu, PTC, "T");
	      } else if ( mu.charge == sum_lep_charge ) {
		muSS_vec  = FourVec(mu, PTC);
		muSS_vecT = FourVec(mu, PTC, "T");
	      }
	      // Find the opposte-charge muon
	      if (mu.charge != sum_lep_charge) {
		muOS      = mu;
		muOS_vec  = FourVec(mu, PTC);
		muOS_vecT = FourVec(mu, PTC, "T");
	      }

	      // Loop through all other muons
	      for (const auto & mu2 : muons) {
		// Skip if mu2 is the same as mu
		if (mu.charge == mu2.charge && (mu.pt != mu2.pt || mu.eta != mu2.eta)) {
		  muSS1 = ( MuonPt(mu, PTC) > MuonPt(mu2, PTC) ? mu : mu2 );
		  muSS2 = ( MuonPt(mu, PTC) > MuonPt(mu2, PTC) ? mu2 : mu );
		  lepSS1_vec  = FourVec(muSS1, PTC);
		  lepSS2_vec  = FourVec(muSS2, PTC);
		}
	      }
	    }
	    ASSERT( muW.pt >= obj_sel.mu_pt_min, "muW.pt >= obj_sel.mu_pt_min" ); // We should always find a W candidate
	    ASSERT( muSS1.charge == muSS2.charge && muOS.charge == -1*sum_lep_charge,
		   "muSS1.charge == muSS2.charge && muOS.charge == -1*sum_lep_charge" ); // We should always find two SS and one OS muon
	  }
	  // For ele + 2 muon events
	  else {
	    ele      = eles.at(0);
	    lep_vec  = FourVec(ele);
	    lep_vecT = FourVec(ele, "T");
	    for (const auto & mu : muons) {
	      if (mu.charge == ele.charge) { muSS1 = mu;  muSS_vec = FourVec(mu, PTC); muSS_vecT = FourVec(mu, PTC, "T"); }
	      else                         { muOS  = mu;  muOS_vec = FourVec(mu, PTC); muOS_vecT = FourVec(mu, PTC, "T"); }
	    }
	    lepSS1_vec  = ( muSS_vec.Pt() >= lep_vec.Pt() ? muSS_vec : lep_vec );
	    lepSS2_vec  = ( muSS_vec.Pt() <  lep_vec.Pt() ? muSS_vec : lep_vec );
	    OS_pair_vec = lep_vec + muOS_vec;

	    ASSERT( muSS1.charge == ele.charge && muOS.charge == -1*sum_lep_charge,
		   "muSS1.charge == ele.charge && muOS.charge == -1*sum_lep_charge"); // We should always find one SS and one OS muon
	  }
	  
	  TLorentzVector MET_vec     = FourVec(*br.met);
	  TLorentzVector lep_MET_vec = lep_vecT + MET_vec;
	  TLorentzVector evt_vec     = FourVec(muons, PTC, eles, jets);
	  TLorentzVector MHT_vec     = FourVec(muons, PTC, eles, jets, "T");
	  MHT_vec.RotateZ(TMath::Pi());
	  TLorentzVector lep_MHT_vec = lep_vecT + MHT_vec;
	  
	  ///////////////////////////////////////////
	  ///  Apply the category selection cuts  ///
	  ///////////////////////////////////////////

	  std::string CAT_CUT   = CAT_CUTS.at(iCat);
	  std::string CAT_UNCUT = CAT_CUT; // Track what sub-strings don't match any known cuts
	  bool pass_cat_cut = true;
	  if ( CAT_CUT.find("mass12") != std::string::npos ) {
	    if ( OS_pair_vec.M() < 12 )                               { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("mass12"), std::string("mass12").length() );
	  }
	  if ( CAT_CUT.find("noZ") != std::string::npos ) {
	    if ( (abs(OS_pair_vec.M() - 91) < 5 && MU) ||
		  abs(H_pair_vec.M()  - 91) < 5        )             { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("noZ"), std::string("noZ").length() );
	  }
	  if ( CAT_CUT.find("onZ") != std::string::npos ) {
	    if ( (abs(OS_pair_vec.M() - 91) > 10 || !MU) &&
		  abs(H_pair_vec.M()  - 91) > 10          )           { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("onZ"), std::string("onZ").length() );
	  }
	  if ( CAT_CUT.find("mt150") != std::string::npos ) {
	    if ( lep_MET_vec.M() > 150 )                              { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("mt150"), std::string("mt150").length() );
	  }
	  if ( CAT_CUT.find("noBtag") != std::string::npos ) {
	    if ( SelectedJets(obj_sel, br, "BTagMedium").size() > 0 ) { pass_cat_cut = false; continue; }
	    if ( SelectedJets(obj_sel, br, "BTagLoose").size()  > 1 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("noBtag"), std::string("noBtag").length() );
	  }
	  if ( CAT_CUT.find("looseLepMVA") != std::string::npos ) {
	    if (!LepMVA(muH1, YEAR, "L") || !LepMVA(muH2, YEAR, "L") ) { pass_cat_cut = false; continue; }
	    if ( MU && !LepMVA(muW, YEAR, "L") )                       { pass_cat_cut = false; continue; }
	    if (!MU && !LepMVA(ele, YEAR, "L") )                       { pass_cat_cut = false; continue; }
	    if (MU) lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_L, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_L, muH2.pt, muH2.eta) * LepMVASF(lepSF_mu_L,  muW.pt, muW.eta);
	    else    lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_L, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_L, muH2.pt, muH2.eta) * LepMVASF(lepSF_ele_L, ele.pt, ele.eta);
	    CAT_UNCUT.erase( CAT_UNCUT.find("looseLepMVA"), std::string("looseLepMVA").length() );
	  }
	  if ( CAT_CUT.find("medLepMVA") != std::string::npos ) {
	    if ( !LepMVA(muH1, YEAR, "M") || !LepMVA(muH2, YEAR, "M") ) { pass_cat_cut = false; continue; }
	    if ( MU && !LepMVA(muW, YEAR, "M") )                        { pass_cat_cut = false; continue; }
	    if (!MU && !LepMVA(ele, YEAR, "M") )                        { pass_cat_cut = false; continue; }
	    if (MU) lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_M, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_M, muH2.pt, muH2.eta) * LepMVASF(lepSF_mu_M,  muW.pt, muW.eta);
	    else    lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_M, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_M, muH2.pt, muH2.eta) * LepMVASF(lepSF_ele_M, ele.pt, ele.eta);
	    CAT_UNCUT.erase( CAT_UNCUT.find("medLepMVA"), std::string("medLepMVA").length() );
	  }
	  if ( CAT_CUT.find("muMeleT_lepMVA") != std::string::npos ) {
	    if ( !LepMVA(muH1, YEAR, "M") || !LepMVA(muH2, YEAR, "M") ) { pass_cat_cut = false; continue; }
	    if ( MU && !LepMVA(muW, YEAR, "M") )                        { pass_cat_cut = false; continue; }
	    if (!MU && !LepMVA(ele, YEAR, "T") )                        { pass_cat_cut = false; continue; }
	    if (MU) lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_M, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_M, muH2.pt, muH2.eta) * LepMVASF(lepSF_mu_M,  muW.pt, muW.eta);
	    else    lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_M, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_M, muH2.pt, muH2.eta) * LepMVASF(lepSF_ele_T, ele.pt, ele.eta);
	    CAT_UNCUT.erase( CAT_UNCUT.find("muMeleT_lepMVA"), std::string("muMeleT_lepMVA").length() );
	  }
	  if ( CAT_CUT.find("tightLepMVA") != std::string::npos ) {
	    if ( !LepMVA(muH1, YEAR, "T") || !LepMVA(muH2, YEAR, "T") ) { pass_cat_cut = false; continue; }
	    if ( MU && !LepMVA(muW, YEAR, "T") )                        { pass_cat_cut = false; continue; }
	    if (!MU && !LepMVA(ele, YEAR, "Y") )                        { pass_cat_cut = false; continue; }
	    if (MU) lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_T, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_T, muH2.pt, muH2.eta) * LepMVASF(lepSF_mu_T,  muW.pt, muW.eta);
	    else    lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] = LepMVASF(lepSF_mu_T, muH1.pt, muH1.eta) * LepMVASF(lepSF_mu_T, muH2.pt, muH2.eta) * LepMVASF(lepSF_ele_T, ele.pt, ele.eta);
	    CAT_UNCUT.erase( CAT_UNCUT.find("tightLepMVA"), std::string("tightLepMVA").length() );
	  }
	  if ( CAT_CUT.find("tightLepCut") != std::string::npos ) {
	    if (        muH1_vec.Pt() < 20 || !muH1.isTightID || muH1.relIso > 0.12 ) { pass_cat_cut = false; continue; };
	    if (        muH2_vec.Pt() < 20 || !muH2.isTightID || muH2.relIso > 0.12 ) { pass_cat_cut = false; continue; };
	    if (  MU && (lep_vec.Pt() < 20 || !muW.isTightID  || muW.relIso > 0.12) ) { pass_cat_cut = false; continue; };
	    if ( !MU && (lep_vec.Pt() < 20 || !ele.isTightID  || ele.relIso > 0.12) ) { pass_cat_cut = false; continue; };
	    CAT_UNCUT.erase( CAT_UNCUT.find("tightLepCut"), std::string("tightLepCut").length() );
	  }
	  // Remove "_" characters left over after all known category sub-strings have been removed
	  while ( CAT_UNCUT.find("_") != std::string::npos ) {
	    CAT_UNCUT.erase( CAT_UNCUT.find("_"), std::string("_").length() );
	  }
	  // If CAT_CUT contains unknown elements, look for it in the central repository
	  if ( CAT_UNCUT.length() > 0 ) pass_cat_cut = InCategory(obj_sel, br, CAT_CUT, verbose);

	  ///  *** SPECIAL CUTS ***  ///
	  // // Only keep signal events if the chosen Higgs candidate pair really comes from the Higgs
	  // if (  sample.Contains("H2Mu") && H_true_vec.M() != H_pair_vec.M() ) pass_cat_cut = false;
	  // Throw away background events if there is a real di-muon pair from a Higgs
	  if ( !sample.Contains("H2Mu") && H_true_vec.M() > 0 ) pass_cat_cut = false;


	  if (not pass_cat_cut) continue;
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUT << ", in category " << CAT_CUT << std::endl;
	  std::string h_pre = (std::string) sample + "_"+OPT_CUT+"_"+CAT_CUT+"_";

	  // Now we know the event passes the optional + category cut
	  b_map_int["OPT_"+OPT_CUT+"_CAT_"+CAT_CUT] = 1;
	  pass_cat_cuts = true;

	  // Store the lepMVA weight separately (varies depending on category selection)
	  float lepMVA_wgt = (isData ? 1.0 : lepMVA_wgts[OPT_CUT+"_"+CAT_CUT]);
	  b_map_flt["lepMVA_wgt_OPT_"+OPT_CUT+"_CAT_"+CAT_CUT] = lepMVA_wgt;


	  //////////////////////////////////////
	  ///  Compute a few more variables  ///
	  //////////////////////////////////////

	  // Leptons ordered by pT
	  TLorentzVector lep1_vec;
	  TLorentzVector lep2_vec;
	  TLorentzVector lep3_vec;

	  if ( lepSS1_vec.Pt() > muOS_vec.Pt() ) {
	    lep1_vec = lepSS1_vec;
	    lep2_vec = (muOS_vec.Pt() > lepSS2_vec.Pt() ? muOS_vec : lepSS2_vec);
	    lep3_vec = (muOS_vec.Pt() > lepSS2_vec.Pt() ? lepSS2_vec : muOS_vec);
	  } else {
	    lep1_vec = muOS_vec;
	    lep2_vec = lepSS1_vec;
	    lep3_vec = lepSS2_vec;
	  }

	  TLorentzVector SS_pair_vec = lepSS1_vec + lepSS2_vec;
	  TLorentzVector trilep_vec  = lep_vec + H_pair_vec;

	  float muSS_MVA = -9;
	  float ele_MVA  = -9;
	  float lepSS1_iso;
	  float lepSS2_iso;
	  float lepSS1_MVA;
	  float lepSS2_MVA;
	  float lepSS1_SIP;
	  float lepSS2_SIP;
	  float lepSS1_seg;
	  float lepSS2_seg;
	  bool  lepSS1_tID;
	  bool  lepSS2_tID;
	  if (MU) {
	    lepSS1_iso = muSS1.relIso;
	    lepSS2_iso = muSS2.relIso;
	    lepSS1_MVA = muSS1.lepMVA;
	    lepSS2_MVA = muSS2.lepMVA;
	    lepSS1_SIP = muSS1.SIP_3D;
	    lepSS2_SIP = muSS2.SIP_3D;
	    lepSS1_seg = muSS1.segCompat;
	    lepSS2_seg = muSS2.segCompat;
	    lepSS1_tID = muSS1.isTightID;
	    lepSS2_tID = muSS2.isTightID;
	  } else {
	    muSS_MVA   = muSS1.lepMVA;
	    ele_MVA    = ele.lepMVA;
	    lepSS1_iso = (muSS_vec.Pt() > ele.pt ? muSS1.relIso : ele.relIso);
	    lepSS2_iso = (muSS_vec.Pt() > ele.pt ? ele.relIso : muSS1.relIso);
	    lepSS1_MVA = (muSS_vec.Pt() > ele.pt ? muSS1.lepMVA : ele.lepMVA);
	    lepSS2_MVA = (muSS_vec.Pt() > ele.pt ? ele.lepMVA : muSS1.lepMVA);
	    lepSS1_SIP = (muSS_vec.Pt() > ele.pt ? muSS1.SIP_3D : ele.SIP_3D);
	    lepSS2_SIP = (muSS_vec.Pt() > ele.pt ? ele.SIP_3D : muSS1.SIP_3D);
	    lepSS1_seg = (muSS_vec.Pt() > ele.pt ? muSS1.segCompat : -0.01);
	    lepSS2_seg = (muSS_vec.Pt() > ele.pt ? -0.01 : muSS1.segCompat);
	    lepSS1_tID = (muSS_vec.Pt() > ele.pt ? muSS1.isTightID : ele.isTightID);
	    lepSS2_tID = (muSS_vec.Pt() > ele.pt ? ele.isTightID : muSS1.isTightID);
	  }


	  ///////////////////////////////////////////////////
	  ///  Generate and fill branches and histograms  ///
	  ///////////////////////////////////////////////////

	  // Tuple containing maps and common options for "BookAndFill"
	  std::tuple< const TString, std::map<TString, float> &, TTree *, std::map<TString, TH1*> &, const TString >
	    tupF{ hist_tree, b_map_flt, out_tree, h_map_1D, h_pre };
	  std::tuple< const TString, std::map<TString, int> &,   TTree *, std::map<TString, TH1*> &, const TString >
	    tupI{ hist_tree, b_map_int, out_tree, h_map_1D, h_pre };

	  // Store sample and event information
	  BookAndFill(b_map_str, out_tree, h_pre, "sample",   sample );
	  BookAndFill(b_map_int, out_tree, h_pre, "event",    br.event->event );
	  BookAndFill(b_map_flt, out_tree, h_pre, "samp_wgt", samp_weight );  // Sample cross section x luminosity weight

	  // Store event weights
	  BookAndFill(tupF, "PU_wgt",    40, -2, 2, isData ? 1.0 : br.PU_wgt  );
	  BookAndFill(tupF, "muon_wgt",  40, -2, 2, isData ? 1.0 : MuonWeight(br, evt_wgt, verbose) );
	  BookAndFill(tupF, "GEN_wgt",   40, -2, 2, isData ? 1.0 : br.GEN_wgt );
	  BookAndFill(tupF, "event_wgt", 40, -2, 2, isData ? 1.0 : event_wgt  );  // Event-by-event PU, NLO, eff SF, etc. weights

	  // Scale event weight by lepMVA weight (varies by category)
	  float cat_evt_wgt = event_wgt * lepMVA_wgt;
	  ASSERT(lepMVA_wgt > 0, "lepMVA_wgt > 0");

	  // Plot kinematic histograms
	  BookAndFill(tupI, "nPV",   12, -0.5, 59.5, br.nVertices, cat_evt_wgt );
	  BookAndFill(tupI, "nEles",  2,  0.5,  1.5, eles.size(),  cat_evt_wgt );

	  // BookAndFill(tupI, "nJets",       8, -0.5, 7.5, jets.size(),                                    cat_evt_wgt );
	  // BookAndFill(tupI, "nBJetsLoose", 6, -0.5, 5.5, SelectedJets(obj_sel, br, "BTagLoose").size(),  cat_evt_wgt );
	  // BookAndFill(tupI, "nBJetsMed",   4, -0.5, 3.5, SelectedJets(obj_sel, br, "BTagMedium").size(), cat_evt_wgt );
	  // BookAndFill(tupI, "nBJetsTight", 4, -0.5, 3.5, SelectedJets(obj_sel, br, "BTagTight").size(),  cat_evt_wgt );
	  BookAndFill(tupI, "nJetsCent",   5, -0.5, 4.5, SelectedJets(obj_sel, br, "Central").size(),    cat_evt_wgt );
	  BookAndFill(tupI, "nJetsFwd",    5, -0.5, 4.5, SelectedJets(obj_sel, br, "Forward").size(),    cat_evt_wgt );

	  BookAndFill(tupF, "jet1_pt",  30,    0, 150, (jets.size() > 0 ? jets.at(0).pt      : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet2_pt",  20,    0, 100, (jets.size() > 1 ? jets.at(1).pt      : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet1_eta", 20, -5.0, 5.0, (jets.size() > 0 ? jets.at(0).eta     : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet2_eta", 20, -5.0, 5.0, (jets.size() > 1 ? jets.at(1).eta     : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet1_CSV", 50, -1.0, 1.0, (jets.size() > 0 ? jets.at(0).deepCSV : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet2_CSV", 50, -1.0, 1.0, (jets.size() > 1 ? jets.at(1).deepCSV : -9), cat_evt_wgt, false );

	  BookAndFill(tupF, "MET",  20, 0,  200, MET_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "MHT",  20, 0,  200, MHT_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "MASS", 30, 0, 1500, evt_vec.M(),  cat_evt_wgt );

	  BookAndFill(tupF, "lep1_pt", 30, 0, 300, lep1_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "lep2_pt", 20, 0, 200, lep2_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "lep3_pt", 20, 0, 100, lep3_vec.Pt(), cat_evt_wgt );

	  BookAndFill(tupF, "muH1_pt", 30, 0, 300, muH1_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "muH2_pt", 20, 0, 200, muH2_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "lep_pt",  20, 0, 200, lep_vec.Pt(),  cat_evt_wgt );

	  BookAndFill(tupF, "muOS_pt",   20, 0, 200, muOS_vec.Pt(),   cat_evt_wgt );
	  BookAndFill(tupF, "muSS_pt",   20, 0, 200, muSS_vec.Pt(),   cat_evt_wgt );  // "Wrong combo" of lep_pt
	  BookAndFill(tupF, "lepSS1_pt", 20, 0, 200, lepSS1_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "lepSS2_pt", 20, 0, 100, lepSS2_vec.Pt(), cat_evt_wgt );

	  // BookAndFill(tupF, "lep1_eta", 24, -2.4, 2.4, lep1_vec.Eta(), cat_evt_wgt );
	  // BookAndFill(tupF, "lep2_eta", 24, -2.4, 2.4, lep2_vec.Eta(), cat_evt_wgt );
	  // BookAndFill(tupF, "lep3_eta", 24, -2.4, 2.4, lep3_vec.Eta(), cat_evt_wgt );

	  BookAndFill(tupF, "muH1_eta", 12, -2.4, 2.4, muH1_vec.Eta(), cat_evt_wgt );
	  BookAndFill(tupF, "muH2_eta", 12, -2.4, 2.4, muH2_vec.Eta(), cat_evt_wgt );
	  BookAndFill(tupF, "lep_eta",  12, -2.4, 2.4, lep_vec.Eta(),  cat_evt_wgt );

	  // BookAndFill(tupF, "muOS_iso",   15, 0, 0.3, muOS.relIso, cat_evt_wgt );
	  // BookAndFill(tupF, "lepSS1_iso", 15, 0, 0.3, lepSS1_iso,  cat_evt_wgt );
	  // BookAndFill(tupF, "lepSS2_iso", 15, 0, 0.3, lepSS2_iso,  cat_evt_wgt );

	  // BookAndFill(tupI, "muOS_tightID",   2, -0.5, 1.5, muOS.isTightID, cat_evt_wgt );
	  // BookAndFill(tupI, "lepSS1_tightID", 2, -0.5, 1.5, lepSS1_tID,     cat_evt_wgt );
	  // BookAndFill(tupI, "lepSS2_tightID", 2, -0.5, 1.5, lepSS2_tID,     cat_evt_wgt );

	  BookAndFill(tupF, "muOS_lepMVA",   40, -1, 1, muOS.lepMVA, cat_evt_wgt );
	  BookAndFill(tupF, "lepSS1_lepMVA", 40, -1, 1, lepSS1_MVA,  cat_evt_wgt );
	  BookAndFill(tupF, "lepSS2_lepMVA", 40, -1, 1, lepSS2_MVA,  cat_evt_wgt );
	  // BookAndFill(h_map_2D, h_pre+"lepSS2_vs_lepSS1_lepMVA", 40, -1, 1, 40, -1, 1, lepSS1_MVA, lepSS2_MVA,  cat_evt_wgt );

	  BookAndFill(tupF, "muSS_lepMVA",  40, -1, 1, muSS_MVA,  cat_evt_wgt, false );
	  BookAndFill(tupF, "ele_lepMVA",   40, -1, 1, ele_MVA,   cat_evt_wgt, false );
	  // BookAndFill(h_map_2D, h_pre+"ele_vs_muSS_lepMVA", 40, -1, 1, 40, -1, 1, muSS_MVA, ele_MVA, cat_evt_wgt, false );

	  // BookAndFill(tupF, "muOS_SIP",   50, 0, 10, muOS.SIP_3D, cat_evt_wgt );
	  // BookAndFill(tupF, "lepSS1_SIP", 50, 0, 10, lepSS1_SIP,  cat_evt_wgt );
	  // BookAndFill(tupF, "lepSS2_SIP", 50, 0, 10, lepSS2_SIP,  cat_evt_wgt );

	  // BookAndFill(tupF, "muOS_segCompat",   50, 0, 1, muOS.segCompat, cat_evt_wgt );
	  // BookAndFill(tupF, "lepSS1_segCompat", 50, 0, 1, lepSS1_seg,     cat_evt_wgt, false );
	  // BookAndFill(tupF, "lepSS2_segCompat", 50, 0, 1, lepSS2_seg,     cat_evt_wgt, false );

	  BookAndFill(tupF, "H_mass_GEN",          50, 110, 160, isData ? -99 : H_GEN_vec.M(),           cat_evt_wgt, false );  // Don't include overflow
	  BookAndFill(tupF, "Z_mass_GEN",          40,   1, 201, isData ? -99 : Z_GEN_vec.M(),           cat_evt_wgt, false );  // Don't include overflow
	  BookAndFill(tupF, "H_mass_true",         50, 110, 160, isData ? -99 : MuPairMass(H_true, PTC), cat_evt_wgt, false );  // Don't include overflow
	  BookAndFill(tupF, "Z_mass_true",         40,   1, 201, isData ? -99 : MuPairMass(Z_true, PTC), cat_evt_wgt, false );  // Don't include overflow
	  BookAndFill(tupF, "H_mass_GEN_or_true",  50, 110, 160, ((H_GEN_vec.M() > 0) || (MuPairMass(H_true, PTC) > 0)) ? H_pair_vec.M() : -99, cat_evt_wgt, false );

	  BookAndFill(b_map_flt, out_tree, h_pre, "H_pair_mass",                      H_pair_vec.M() );
	  BookAndFill(h_map_1D,            h_pre+ "H_pair_mass_zoomH", 100, 110, 160, H_pair_vec.M(),  cat_evt_wgt, false );  // Don't include overflow
	  BookAndFill(h_map_1D,            h_pre+ "H_pair_mass_zoomZ",  40,  81, 101, H_pair_vec.M(),  cat_evt_wgt, false );  // Don't include overflow
	  BookAndFill(h_map_1D,            h_pre+ "H_pair_mass_window", 10, 110, 160, H_pair_vec.M(),  cat_evt_wgt, false );  // Don't include overflow
	  BookAndFill(h_map_1D,            h_pre+ "H_pair_mass_wide",   40,   0, 400, H_pair_vec.M(),  cat_evt_wgt );
	  BookAndFill(h_map_1D,            h_pre+ "OS_pair_mass_zoomZ", 40,  81, 101, OS_pair_vec.M(), cat_evt_wgt, false );  // Don't include overflow

	  BookAndFill(tupF, "H_pair_mass_err",   40,   0,   4, H_pair.massErr,   cat_evt_wgt );
	  BookAndFill(tupF, "H_pair_pt",         30,   0, 300, H_pair_vec.Pt(),  cat_evt_wgt );
	  BookAndFill(tupF, "OS_pair_mass",      40,   0, 400, OS_pair_vec.M(),  cat_evt_wgt );  // "Wrong combo" of H_pair_mass
	  BookAndFill(tupF, "OS_pair_pt",        30,   0, 300, OS_pair_vec.Pt(), cat_evt_wgt );  // "Wrong combo" of H_pair_pt
	  BookAndFill(tupF, "SS_pair_mass",      40,   0, 400, SS_pair_vec.M(),  cat_evt_wgt );
	  BookAndFill(tupF, "SS_pair_pt",        30,   0, 300, SS_pair_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "trilep_mass",       30,   0, 600, trilep_vec.M(),   cat_evt_wgt );


	  // More useful variables for WH vs. Z+jets and WH vs. WZ
	  BookAndFill(tupI, "lep_charge",           3, -1.5, 1.5, sum_lep_charge,                    cat_evt_wgt );
	  BookAndFill(tupF, "lep_H_pair_dEta",     20, -5.0, 5.0, SignedDEta( lep_vec, H_pair_vec ), cat_evt_wgt );
	  BookAndFill(tupF, "muSS_OS_pair_dEta",   20, -5.0, 5.0, SignedDEta(muSS_vec, OS_pair_vec), cat_evt_wgt );  // "Wrong combo" of lep_H_pair_dEta
	  BookAndFill(tupF, "lep_H_pair_dR",       20,    0, 6.0, lep_vec .DeltaR(H_pair_vec),       cat_evt_wgt );
	  BookAndFill(tupF, "muSS_OS_pair_dR",     20,    0, 6.0, muSS_vec.DeltaR(OS_pair_vec),      cat_evt_wgt );  // "Wrong combo" of lep_H_pair_dR
	  BookAndFill(tupF, "lep_muSS_cosThStar",  10, -1.0, 1.0, CosThetaStar(lep_vec, muSS_vec),   cat_evt_wgt );
	  BookAndFill(tupF, "muSS_lep_cosThStar",  10, -1.0, 1.0, CosThetaStar(muSS_vec, lep_vec),   cat_evt_wgt );  // "Wrong combo" of lep_muSS_cosThStar
	  BookAndFill(tupF, "lep_muOS_cosThStar",  10, -1.0, 1.0, CosThetaStar( lep_vec, muOS_vec),  cat_evt_wgt );
	  BookAndFill(tupF, "muSS_muOS_cosThStar", 10, -1.0, 1.0, CosThetaStar(muSS_vec, muOS_vec),  cat_evt_wgt );  // "Wrong combo" of lep_muOS_cosThStar
	  BookAndFill(tupF, "lep_muSS_dEta",       20, -5.0, 5.0, SignedDEta(lep_vec, muSS_vec),     cat_evt_wgt );
	  BookAndFill(tupF, "lep_muSS_dR",         20,    0, 6.0, lep_vec.DeltaR(muSS_vec),          cat_evt_wgt );
	  BookAndFill(tupF, "lep_muOS_dEta",       20, -5.0, 5.0, SignedDEta( lep_vec, muOS_vec),    cat_evt_wgt );
	  BookAndFill(tupF, "muSS_muOS_dEta",      20, -5.0, 5.0, SignedDEta(muSS_vec, muOS_vec),    cat_evt_wgt );  // "Wrong combo" of lep_muOS_dEta
	  BookAndFill(tupF, "lep_muOS_dR",         20,    0, 6.0, lep_vec .DeltaR(muOS_vec),         cat_evt_wgt );
	  BookAndFill(tupF, "muSS_muOS_dR",        30,    0, 6.0, muSS_vec.DeltaR(muOS_vec),         cat_evt_wgt );  // "Wrong combo" of lep_muOS_dR
	  BookAndFill(tupF, "lep_MET_dPhi_abs",    16,    0, 3.2, abs( lep_vec.DeltaPhi(MET_vec)),   cat_evt_wgt );
	  BookAndFill(tupF, "muSS_MET_dPhi_abs",   16,    0, 3.2, abs(muSS_vec.DeltaPhi(MET_vec)),   cat_evt_wgt );  // "Wrong combo" of lep_MET_dPhi_abs
	  BookAndFill(tupF, "lep_MET_MT",          20,    0, 200, lep_MET_vec.M(),                   cat_evt_wgt );
	  BookAndFill(tupF, "muSS_MET_MT",         30,    0, 300, (muSS_vecT + MET_vec).M(),         cat_evt_wgt );  // "Wrong combo" of lep_MET_MT
	  BookAndFill(tupF, "lep_MHT_MT",          20,    0, 200, lep_MHT_vec.M(),                   cat_evt_wgt );
	  BookAndFill(tupF, "muSS_MHT_MT",         30,    0, 300, (muSS_vecT + MHT_vec).M(),         cat_evt_wgt );  // "Wrong combo" of lep_MHT_MT
	  BookAndFill(tupF, "lep_MHT_dPhi_abs",    16,    0, 3.2, abs( lep_vec.DeltaPhi(MHT_vec)),   cat_evt_wgt );
	  BookAndFill(tupF, "muSS_MHT_dPhi_abs",   16,    0, 3.2, abs(muSS_vec.DeltaPhi(MHT_vec)),   cat_evt_wgt );  // "Wrong combo" of lep_MHT_dPhi_abs
	  BookAndFill(tupF, "MHT_MET_dPhi_abs",    16,    0, 3.2, abs(MHT_vec.DeltaPhi(MET_vec)),    cat_evt_wgt );


	  // Extra variables especially for Xunwu's BDTs
	  BookAndFill( b_map_flt, out_tree, h_pre, "dimu_mass",           b_map_flt["H_pair_mass"]        );  // In Xunwu's code: dimu_vec.M()
	  BookAndFill( b_map_flt, out_tree, h_pre, "dimu_pt",             b_map_flt["H_pair_pt"]          );  // In Xunwu's code: dimu_vec.Pt()
	  BookAndFill( b_map_flt, out_tree, h_pre, "lep_pt",              b_map_flt["lep_pt"]             );  // In Xunwu's code: lep_vec.Pt()
	  BookAndFill( b_map_flt, out_tree, h_pre, "ldimu_abs_dEta",  abs(b_map_flt["lep_H_pair_dEta"])   );  // In Xunwu's code: abs(lep_vec.Eta() - dimu_vec.Eta())
	  BookAndFill( b_map_flt, out_tree, h_pre, "ldimu_dR",            b_map_flt["lep_H_pair_dR"]      );  // In Xunwu's code: lep_vec.DeltaR(dimu_vec)
	  BookAndFill( b_map_flt, out_tree, h_pre, "cts_lmuSS",           b_map_flt["lep_muSS_cosThStar"] );  // In Xunwu's code: CosThetaStar(lep_vec, mu1_vec)
	  BookAndFill( b_map_flt, out_tree, h_pre, "cts_lmuOS",           b_map_flt["lep_muOS_cosThStar"] );  // In Xunwu's code: CosThetaStar(lep_vec, mu2_vec)
	  BookAndFill( b_map_flt, out_tree, h_pre, "lmuSS_abs_dEta",  abs(b_map_flt["lep_muSS_dEta"])     );  // In Xunwu's code: abs(lep_vec.Eta() - mu2_vec.Eta())
	  BookAndFill( b_map_flt, out_tree, h_pre, "lmuSS_dR",            b_map_flt["lep_muSS_dR"]        );  // In Xunwu's code: lep_vec.DeltaR(mu2_vec)
	  BookAndFill( b_map_flt, out_tree, h_pre, "lmuOS_abs_dEta",  abs(b_map_flt["lep_muOS_dEta"])     );  // In Xunwu's code: abs(lep_vec.Eta() - mu1_vec.Eta())
	  BookAndFill( b_map_flt, out_tree, h_pre, "lmuOS_dR",            b_map_flt["lep_muOS_dR"]        );  // In Xunwu's code: lep_vec.DeltaR(mu1_vec)
	  BookAndFill( b_map_flt, out_tree, h_pre, "met_pt",              b_map_flt["MET"]                );  // In Xunwu's code: met_vec.Pt()
	  BookAndFill( b_map_flt, out_tree, h_pre, "mt_lmet",             b_map_flt["lep_MET_MT"]         );  // In Xunwu's code: lMET_vec.M()
	  BookAndFill( b_map_flt, out_tree, h_pre, "abs_dPhi_lmet",       b_map_flt["lep_MET_dPhi_abs"]   );  // In Xunwu's code: abs(lep_vec.DeltaPhi(met_vec))  );
	  BookAndFill( b_map_flt, out_tree, h_pre, "mht_pt",            (FourVec(*br.mht)).Pt()           );  // In Xunwu's code: mht_vec.Pt()
	  BookAndFill( b_map_flt, out_tree, h_pre, "mt_lmht",           (FourVec(*br.mht) + lep_vecT).M() );  // In Xunwu's code: lMHT_vec.M() - MAY BE BUGGY IN XUNWU'S TRAINING!!! (AWB 15.05.2019)

	  // Extra variables especially for Andrew's BDTs
	  BookAndFill( b_map_flt, out_tree, h_pre, "muH1_eta_abs", abs(b_map_flt["muH1_eta"]) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "muH2_eta_abs", abs(b_map_flt["muH2_eta"]) );

	  // Evaluate MVA output values
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_XWZ_noMass",      BDT_XWZ_noMass  .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_XWZ_withMass",    BDT_XWZ_withMass.Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_AWB_v2_noMass",   BDT_AWB_2_noMass  .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_AWB_v2_withMass", BDT_AWB_2_withMass.Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDTG_UF_v1",          b_map_flt["BDT_AWB_v2_noMass"] );  // Needed for re-trained noMassBDT vs. mass
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_AWB_v2_retrain",  BDT_AWB_2_retrain .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_AWB_v3_noMass",   BDT_AWB_3_noMass  .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_AWB_v3_withMass", BDT_AWB_3_withMass.Evaluate(b_map_flt, b_map_int) );
	  b_map_flt["BDTG_UF_v1"] = b_map_flt["BDT_AWB_v3_noMass"];  // Reset input value for v3 re-trained noMassBDT vs. mass
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_AWB_v3_retrain",  BDT_AWB_3_retrain .Evaluate(b_map_flt, b_map_int) );
	  b_map_flt["BDTG_UF_v1"] = b_map_flt["BDT_AWB_v2_noMass"];  // ... and re-reset to avoid run-time errors

	  // Combine non-mass and with-mass BDT distributions into one plot
	  // Compress and shift non-mass BDT distribution from [-0.8, 0.9] to [-1.0, 0.4]
	  double noMass = (b_map_flt["BDT_AWB_v3_noMass"] * 0.8) - 0.3;
	  noMass = std::max(-0.999, std::min(0.399, noMass));
	  double withMass = b_map_flt["BDT_AWB_v3_withMass"];
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_AWB_v3_combo", withMass > 0.4 ? withMass : noMass );

	  // Use only even-numbered MC events to fill histograms
	  if ( (b_map_int["event"] % 2) == 0 || isData ) {
	    float two = (isData ? 1.0 : 2.0);  // Scale MC by a factor of 2
	    BookAndFill( h_map_1D, h_pre+"BDT_XWZ_noMass",      20, -1, 1, b_map_flt["BDT_XWZ_noMass"],      cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_XWZ_withMass",    20, -1, 1, b_map_flt["BDT_XWZ_withMass"],    cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v2_noMass",   20, -1, 1, b_map_flt["BDT_AWB_v2_noMass"],   cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v2_withMass", 20, -1, 1, b_map_flt["BDT_AWB_v2_withMass"], cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v2_retrain",  20, -1, 1, b_map_flt["BDT_AWB_v2_retrain"],  cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_noMass",   20, -1, 1, b_map_flt["BDT_AWB_v3_noMass"],   cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass", 20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"], cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_retrain",  20, -1, 1, b_map_flt["BDT_AWB_v3_retrain"],  cat_evt_wgt*two );

	    BookAndFill( h_map_1D, h_pre+"BDT_XWZ_noMass_zoom",      100, -1, 1, b_map_flt["BDT_XWZ_noMass"],      cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_XWZ_withMass_zoom",    100, -1, 1, b_map_flt["BDT_XWZ_withMass"],    cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v2_noMass_zoom",   100, -1, 1, b_map_flt["BDT_AWB_v2_noMass"],   cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v2_withMass_zoom", 100, -1, 1, b_map_flt["BDT_AWB_v2_withMass"], cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v2_retrain_zoom",  100, -1, 1, b_map_flt["BDT_AWB_v2_retrain"],  cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_noMass_zoom",   100, -1, 1, b_map_flt["BDT_AWB_v3_noMass"],   cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_zoom", 100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"], cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_retrain_zoom",  100, -1, 1, b_map_flt["BDT_AWB_v3_retrain"],  cat_evt_wgt*two );

	    std::vector<float> combo_bins;
	    // Fill with 0.2 wide bins from -1.0 to 1.4, width 0.01 after
	    for (int i = 0; i < 201; i++) if (i > 140 || (i % 20) == 0) combo_bins.push_back(-1 + i*0.01);
	    int nBins = combo_bins.size() - 1;
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo",              20, -1, 1, b_map_flt["BDT_AWB_v3_combo"], cat_evt_wgt*two );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_zoom", nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"], cat_evt_wgt*two );

	    // Write histograms with 20% systematic shift on no-mass BDT shape
	    float sysUp   = (isData ? 1.0 : 1 + b_map_flt["BDT_AWB_v3_noMass"] * 0.2);
	    float sysDown = (isData ? 1.0 : 1 - b_map_flt["BDT_AWB_v3_noMass"] * 0.2);
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_UP",           20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"], cat_evt_wgt*two*sysUp );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_zoom_UP",     100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"], cat_evt_wgt*two*sysUp );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_UP",              20, -1, 1, b_map_flt["BDT_AWB_v3_combo"],    cat_evt_wgt*two*sysUp );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_zoom_UP", nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],    cat_evt_wgt*two*sysUp );

	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_DOWN",           20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"], cat_evt_wgt*two*sysDown );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_zoom_DOWN",     100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"], cat_evt_wgt*two*sysDown );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_DOWN",              20, -1, 1, b_map_flt["BDT_AWB_v3_combo"],    cat_evt_wgt*two*sysDown );
	    BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_zoom_DOWN", nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],    cat_evt_wgt*two*sysDown );

	    // Write histograms also with tighter selection cuts
	    if ( lep_vec.Pt() > 20 ) {
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_lepW_pt20",                20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_lepW_pt20_zoom",          100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_lepW_pt20_zoom",       nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],     cat_evt_wgt*two );
	    } if ( muH1_vec.Pt() > 20 && muH2_vec.Pt() > 20 && lep_vec.Pt() > 20 ) {
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_allLep_pt20",              20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_allLep_pt20_zoom",        100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_allLep_pt20_zoom",     nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],     cat_evt_wgt*two );
	    } if ( LepMVA(muSS1, YEAR, "T") && (MU ? LepMVA(muSS2, YEAR, "T") : LepMVA(ele, YEAR, "T")) ) {
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_lepMVA_tight_SS",          20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_lepMVA_tight_SS_zoom",    100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_lepMVA_tight_SS_zoom", nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],     cat_evt_wgt*two );
	    } if ( LepMVA(muOS, YEAR, "T") && LepMVA(muSS1, YEAR, "T") && (MU ? LepMVA(muSS2, YEAR, "T") : LepMVA(ele, YEAR, "T")) ) {
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_lepMVA_tight",             20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_lepMVA_tight_zoom",       100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_lepMVA_tight_zoom",    nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],     cat_evt_wgt*two );
	    } if ( (!MU) || (H_pair_vec.Pt() > OS_pair_vec.Pt()) ) {
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_hiPt_Higgs",               20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_hiPt_Higgs_zoom",         100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_hiPt_Higgs_zoom",      nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],     cat_evt_wgt*two );
	    } if ( (!MU) || (H_pair_vec.Pt() > OS_pair_vec.Pt()) || OS_pair_vec.Pt() < 110 || OS_pair_vec.Pt() > 160 ) {
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_hiPt2_Higgs",              20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_hiPt2_Higgs_zoom",        100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_hiPt2_Higgs_zoom",     nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],     cat_evt_wgt*two );
	    } if ( (!MU) || abs(OS_pair_vec.M() - 91) > 10 ) {
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_noZ10",                    20, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_withMass_noZ10_zoom",              100, -1, 1, b_map_flt["BDT_AWB_v3_withMass"],  cat_evt_wgt*two );
	      BookAndFill( h_map_1D, h_pre+"BDT_AWB_v3_combo_noZ10_zoom",           nBins, combo_bins, b_map_flt["BDT_AWB_v3_combo"],     cat_evt_wgt*two );
	    }

	    // Write dimuon mass histogram in different no-mass BDT bins
	         if ( b_map_flt["BDT_AWB_v3_noMass"] > 0.6 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_p06_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_AWB_v3_noMass"] > 0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_p02_p06_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_AWB_v3_noMass"] > -0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_n02_p02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_n10_n02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );

	         if ( b_map_flt["BDT_AWB_v3_noMass"] > 0.76 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_p076_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_AWB_v3_noMass"] > 0.68 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_p068_p076_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_AWB_v3_noMass"] > 0.60 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_p060_p068_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_n100_p060_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );

	  } // End conditional: if ( (b_map_int["event"] % 2) == 0 || isData )


	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)

      // Fill branches in output tuple tree, if event passed at least one category
      if (pass_opt_cuts && pass_cat_cuts) {

	if (n_map_int != -99 && n_map_int != b_map_int.size()) {
	  std::cout << "\n\nBizzare error!!! b_map_int = " << n_map_int << ", b_map_int.size() = " << b_map_int.size() << std::endl;
	  for (std::map<TString, int>::iterator it = b_map_int.begin(); it != b_map_int.end(); ++it) {
	    std::cout << it->first << " = " << it->second << std::endl;
	  }
	  ASSERT(false, "false");
	}
	if (n_map_flt != -99 && n_map_flt != b_map_flt.size()) {
	  std::cout << "\n\nBizzare error!!! b_map_flt = " << n_map_flt << ", b_map_flt.size() = " << b_map_flt.size() << std::endl;
	  for (std::map<TString, float>::iterator it = b_map_flt.begin(); it != b_map_flt.end(); ++it) {
	    std::cout << it->first << " = " << it->second << std::endl;
	  }
	  ASSERT(false, "false");
	}
	if (n_map_str != -99 && n_map_str != b_map_str.size()) {
	  std::cout << "\n\nBizzare error!!! b_map_str = " << n_map_str << ", b_map_str.size() = " << b_map_str.size() << std::endl;
	  for (std::map<TString, std::string>::iterator it = b_map_str.begin(); it != b_map_str.end(); ++it) {
	    std::cout << it->first << " = " << it->second << std::endl;
	  }
	  ASSERT(false, "false");
	}

	n_map_int = b_map_int.size();
	n_map_flt = b_map_flt.size();
	n_map_str = b_map_str.size();

	out_tree->Fill();
      }

    } // End loop: for (int iPair = 0; iPair < 2; iPair++)

  } // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;


  std::cout << "\n******* normalizing histos, weight =  " << samp_weight << " *******" << std::endl;
  if (h_map_1D.empty()) std::cout << "h_map_1D is empty" << std::endl;
  else {
    for (std::map<TString, TH1*>::iterator it_term = h_map_1D.begin() ; it_term != h_map_1D.end() ; it_term++) {
        it_term->second->Scale(samp_weight);
    }
  }
  if (h_map_2D.empty()) std::cout << "h_map_2D is empty" << std::endl;
  else {
    for (std::map<TString, TH2*>::iterator it_term = h_map_2D.begin() ; it_term != h_map_2D.end() ; it_term++) {
        it_term->second->Scale(samp_weight);
    }
  }

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
	if ( StartsWith(h_name, (sample+"_"+optCatStr+"_").Data()) ) { // Histogram name starts with SAMP+OPT+CAT
	  // Remove optional selection and category cuts from histogram names
	  h_name.erase( h_name.find(optCatStr+"_"), optCatStr.length() + 1 );
	  it->second->SetName(h_name.c_str());
	  it->second->SetTitle(h_name.c_str());
	  // Re-scale shape-systematic shifted histograms to nominal yield
	  if ( EndsWith(h_name, "_UP") || EndsWith(h_name, "_DOWN") ) { // Histogram name ends with "UP" or "DOWN"
	    for (std::map<TString, TH1*>::iterator it2 = h_map_1D.begin(); it2 != h_map_1D.end(); ++it2) {
	      if (it->first == it2->first+"_UP" || it->first == it2->first+"_DOWN")
		it->second->Scale( it2->second->Integral() / it->second->Integral() );
	    }
	  }
	  // std::cout << "  * Writing 1D histogram " << it->second->GetName() << std::endl;
	  it->second->Write();
	}
      }
      // Write out 2D histograms
      for (std::map<TString, TH2*>::iterator it = h_map_2D.begin(); it != h_map_2D.end(); ++it) {
	std::string h_name = it->second->GetName();
	if ( StartsWith(h_name, (sample+"_"+optCatStr+"_").Data()) ) { // Histogram name starts with SAMP+OPT+CAT
	  // Remove optional selection and category cuts from histogram names
	  h_name.erase( h_name.find(optCatStr+"_"), optCatStr.length() + 1 );
	  it->second->SetName(h_name.c_str());
	  it->second->SetTitle(h_name.c_str());
	  // std::cout << "  * Writing 2D histogram " << it->second->GetName() << std::endl;
	  it->second->Write();
	}
      }
  
      out_file->Write();
      out_file->Close();
      std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
      
    } // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
  } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)


  // Only write output tree file if there are some events
  if (out_tree->GetEntries() > 0) {
    std::cout << "\n******* Writing output tuple file " << out_tuple_name.Data()
	      << " with " << out_tree->GetEntries() << " events *******" << std::endl;
    out_tuple->cd();
    out_tree->Write();
    out_tuple->Write();
    out_tuple->Close();
  } else {
    std::cout << "\n******* NOT writing output tuple file " << out_tuple_name.Data()
	      << " - " << out_tree->GetEntries() << " events! *******" << std::endl;
  }


  std::cout << "\nExiting WH_lep()\n";
  
} // End void WH_lep()
