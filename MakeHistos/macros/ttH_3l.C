
///////////////////////////////////////////////////
///   Macro to study ttH-->lvjj+mumu category   ///
///                                             ///
///        Andrew Brinkerhoff  19.01.2019       ///
///////////////////////////////////////////////////

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

const std::string YEAR = "2017";
// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// const std::string YEAR = "2018";
// // Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
// R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc700/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)


// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;    // Minimum index of input files to process
const int MAX_FILE = 1;    // Maximum index of input files to process
const int MAX_EVT  = 1000; // Maximum number of events to process
const int PRT_EVT  = 100;  // Print every N events
const float SAMP_WGT = 1.0;
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information

const TString IN_DIR = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8/H2Mu_ttH_125";
const TString SAMPLE = "H2Mu_ttH_125";
// const TString IN_DIR = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/tt_ll_MG";
// const TString SAMPLE = "tt_ll_MG";
// const TString IN_DIR = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/ttW";
// const TString SAMPLE = "ttW";
// const TString IN_DIR = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
// const TString SAMPLE = "SingleMu";

// const TString IN_DIR = "/eos/cms/store/user/bortigno/h2mm/ntuples/2018/102X/prod-v18.1.6.skim3l/ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8/H2Mu_ttH_125/190528_111755/0000";
// const TString SAMPLE = "H2Mu_ttH_125";
// const TString IN_DIR = "/eos/cms/store/user/bortigno/h2mm/ntuples/2018/102X/prod-v18.1.6.skim3l/SingleMuon/SingleMu_2018D/190528_111415/0000";
// const TString SAMPLE = "SingleMu";


const std::string SLIM  = (YEAR == "2018" ? "notSlim" : "Slim");  // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR   = "plots";
const TString HIST_TREE = "HistTree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both

// Cuts which every event must pass, applied in sequence
const std::vector<std::string> SEL_CUTS = {"PreselRun2"};
// Multiple selection cuts, applied independently in parallel
// const std::vector<std::string> OPT_CUTS = {"3lep", "3lep_allMass"};
const std::vector<std::string> OPT_CUTS = {"3lep"};
// Category selection cuts, also applied in parallel
// *** IMPORTANT!!! No category name may start with a sub-string which is identical to another entire category name! ***
const std::vector<std::string> CAT_CUTS = { "looseLepMVA_ge2j_btag",
					    "medLepMVA_noZ10_ge2j_btag",
                                            "hiPt_lepW20_medLepMVA_noZ10_ge2j_btag" };


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void ttH_3l( TString sample = "", TString in_dir = "", TString out_dir = "",
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
    if (YEAR == "2018") {
      for (int i = MIN_FILE; i <= MAX_FILE; i++) {
	in_file_name.Form("%s/tuple_%d.root", in_dir.Data(), i);
	std::cout << "Adding file " << in_file_name.Data() << std::endl;
	in_file_names.push_back(in_file_name.Data());
      }
    } else {
      in_file_name.Form("%s/NTuple_0.root", in_dir.Data());
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

  evt_sel.muPair_mass_min = 12;  // Allow masses down to 12 GeV (instead of 60 GeV) for background studies

  obj_sel.Print();
  evt_sel.Print();
  evt_wgt.Print();

  ObjectSelectionConfig obj_sel_orig = obj_sel;  // Store original object selection (may change for different CAT_CUTS)

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later

  std::cout << "\n******* About to load 2D LepMVA efficiency scale factor histograms *******" << std::endl;
  std::map<std::string, TH2F *> lepSF;
  lepSF["mu_T"]  = LoadSFsLepMVA(YEAR,  "mu", "T");
  lepSF["mu_M"]  = LoadSFsLepMVA(YEAR,  "mu", "M");
  lepSF["mu_L"]  = LoadSFsLepMVA(YEAR,  "mu", "L");
  lepSF["ele_T"] = LoadSFsLepMVA(YEAR, "ele", "T");
  lepSF["ele_M"] = LoadSFsLepMVA(YEAR, "ele", "M");
  lepSF["ele_L"] = LoadSFsLepMVA(YEAR, "ele", "L");

  std::cout << "\n******* About to load XML files for ttH reconstruction BDT *******" << std::endl;
  MVA BDT_ttH_reco( "data/XMLs/ttH_3l/Charlie/2017_01_25/ttH_HtoWW_SS2l_reconstruction/",
		    "weights/TMVAClassification_bloose_BDTG.weights.xml",
		    "BDTG" );
  
  MVA BDT_v1_all_withMass( "data/XMLs/ttH_3l/Andrew/2019_07_04/f_Opt_2019_07_04_varsAll_withMass_all_sig_all_bkg_withMass/",
			   "weights/f_Opt_2019_07_04_varsAll_withMass_all_sig_all_bkg_withMass_BDTG_AWB_lite.weights.xml",
			   "BDTG_AWB_lite" );
  
  MVA BDT_v1_all( "data/XMLs/ttH_3l/Andrew/2019_07_04/f_Opt_2019_07_04_varsAll_resWgt_all_sig_all_bkg_resWgt/",
	       "weights/f_Opt_2019_07_04_varsAll_resWgt_all_sig_all_bkg_resWgt_BDTG_AWB_lite.weights.xml",
	       "BDTG_AWB_lite" );

  MVA BDT_v1_med( "data/XMLs/ttH_3l/Andrew/2019_07_04/f_Opt_2019_07_04_varsM_resWgt_all_sig_all_bkg_resWgt/",
	       "weights/f_Opt_2019_07_04_varsM_resWgt_all_sig_all_bkg_resWgt_BDTG_AWB_lite.weights.xml",
	       "BDTG_AWB_lite" );
  
  MVA BDT_v1_med_noBDT( "data/XMLs/ttH_3l/Andrew/2019_07_04/f_Opt_2019_07_04_varsM_noBDT_resWgt_all_sig_all_bkg_resWgt/",
		     "weights/f_Opt_2019_07_04_varsM_noBDT_resWgt_all_sig_all_bkg_resWgt_BDTG_AWB_lite.weights.xml",
		     "BDTG_AWB_lite" );
  
  MVA BDT_v1_tight( "data/XMLs/ttH_3l/Andrew/2019_07_04/f_Opt_2019_07_04_varsT_resWgt_all_sig_all_bkg_resWgt/",
		 "weights/f_Opt_2019_07_04_varsT_resWgt_all_sig_all_bkg_resWgt_BDTG_AWB_lite.weights.xml",
		 "BDTG_AWB_lite" );
  
  MVA BDT_v1_tight_noBDT( "data/XMLs/ttH_3l/Andrew/2019_07_04/f_Opt_2019_07_04_varsT_noBDT_resWgt_all_sig_all_bkg_resWgt/",
		       "weights/f_Opt_2019_07_04_varsT_noBDT_resWgt_all_sig_all_bkg_resWgt_BDTG_AWB_lite.weights.xml",
		       "BDTG_AWB_lite" );
  


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

    // For original 2016 and some 2017 and 2018 NTuples, convert "SlimJets" collection into regular jets
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

	// Reset object selection to original state
	obj_sel = obj_sel_orig;

	// Choose Higgs mass (110 - 160 GeV) or all mass (> 12 GeV) requirements
	bool allMass = ( OPT_CUT.find("allMass") != std::string::npos );
	// If two muon pairs are in the signal mass window, choose the higher-pT pair
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
	// All opposite-charge lepton pairs must have mass > 12 GeV
	if (true) {
	  if ( SelectedEles(obj_sel, br).size() == 0 ) {
	    if ( SelectedMuPairs(obj_sel, br).size() != 2 ) continue;
	    if ( FourVec( SelectedMuPairs(obj_sel, br).at(0), PTC).M() < 12 ) continue;
	    if ( FourVec( SelectedMuPairs(obj_sel, br).at(1), PTC).M() < 12 ) continue;
	    MU = true;
	  } else if ( SelectedEles(obj_sel, br).size() == 1 ) {
	    if ( SelectedMuPairs(obj_sel, br).size() != 1 ) continue;
	    if ( iPair != 0 ) continue;  // Only one Higgs pair possible, hence iPair must be 0
	    int iMuOS = (SelectedMuons(obj_sel, br).at(0).charge != SelectedEles(obj_sel, br).at(0).charge ? 0 : 1);
	    if ( ( FourVec(SelectedMuons(obj_sel, br).at(iMuOS), PTC) + FourVec(SelectedEles(obj_sel, br).at(0)) ).M() < 12 ) continue;
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
	  if ( MU && hiPt && !allMass ) {
	    // Choose the Higgs candidate pair with the higher vector pT value
	    if ( FourVec(SelectedMuPairs(obj_sel, br).at((iPair+1) % 2), PTC).M()  > 110 &&
		 FourVec(SelectedMuPairs(obj_sel, br).at((iPair+1) % 2), PTC).M()  < 160 &&
		 FourVec(SelectedMuPairs(obj_sel, br).at((iPair+1) % 2), PTC).Pt() > H_pair_vec.Pt() ) continue;
	  }
	}

	// Special cuts for 3 muon category
	if ( StartsWith(OPT_CUT, "3mu") ) {
	  if (!MU) continue;
	}
	else if ( StartsWith(OPT_CUT, "e2mu") ) {
	  if (MU) continue;
	} else ASSERT( StartsWith(OPT_CUT, "3lep"), "StartsWith(OPT_CUT, '3lep')" );
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

	  // Reset object selection to original state
	  obj_sel = obj_sel_orig;

	  // Get selected muons and electrons in event
	  MuonInfos muons = SelectedMuons(obj_sel, br);
	  EleInfos  eles  = SelectedEles(obj_sel, br);
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
	      if ( IsGenMatched( muPair, *br.muons, *br.genMuons, *br.genParents, "Z") ||
		   IsGenMatched( muPair, *br.muons, *br.genMuons, *br.genParents, "gamma") ||
		   IsGenMatched( muPair, *br.muons, *br.genMuons, *br.genParents, "light_quark") ) {
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
	    ASSERT( MuonPt(muW, PTC) >= obj_sel.mu_pt_min, "MuonPt(muW, PTC) >= obj_sel.mu_pt_min" ); // We should always find a W candidate
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
	  
	  ///////////////////////////////////////////
	  ///  Apply the category selection cuts  ///
	  ///////////////////////////////////////////

	  std::string CAT_CUT   = CAT_CUTS.at(iCat);
	  std::string CAT_UNCUT = CAT_CUT; // Track what sub-strings don't match any known cuts
	  bool pass_cat_cut = true;

	  if ( CAT_CUT.find("noZ10") != std::string::npos ) {
	    if ( (abs(OS_pair_vec.M() - 91) < 10 && MU) ||
		  abs(H_pair_vec.M()  - 91) < 10        )         { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("noZ10"), std::string("noZ10").length() );
	  }
	  if ( CAT_CUT.find("noZ5") != std::string::npos ) {
	    if ( (abs(OS_pair_vec.M() - 91) < 5 && MU) ||
		  abs(H_pair_vec.M()  - 91) < 5        )          { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("noZ5"), std::string("noZ5").length() );
	  }
	  if ( CAT_CUT.find("onZ10") != std::string::npos ) {
	    if ( (abs(OS_pair_vec.M() - 91) > 10 || !MU) &&
		  abs(H_pair_vec.M()  - 91) > 10          )       { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("onZ10"), std::string("onZ10").length() );
	  }
	  if ( CAT_CUT.find("hiPt") != std::string::npos ) {
	    if ( MU && H_pair_vec.Pt() < OS_pair_vec.Pt() &&
		 OS_pair_vec.M() > 110 && OS_pair_vec.M() < 160 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("hiPt"), std::string("hiPt").length() );
	  }
	  if ( CAT_CUT.find("onlyHiPt") != std::string::npos ) {
	    if ( MU && H_pair_vec.Pt() < OS_pair_vec.Pt() )       { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("onlyHiPt"), std::string("onlyHiPt").length() );
	  }
	  if ( CAT_CUT.find("ele20") != std::string::npos ) {
	    obj_sel.ele_pt_min = 20;
	    if ( !MU && SelectedEles(obj_sel, br).size() != 1 )   { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("ele20"), std::string("ele20").length() );
	  }
	  if ( CAT_CUT.find("lepW20") != std::string::npos ) {
	    obj_sel.ele_pt_min = 20;
	    if ( !MU && SelectedEles(obj_sel, br).size() != 1 )   { pass_cat_cut = false; continue; }
	    if (  MU && lep_vec.Pt() < 20 )                       { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("lepW20"), std::string("lepW20").length() );
	  }
	  if ( CAT_CUT.find("lep20") != std::string::npos ) {
	    obj_sel.mu_pt_min  = 20;
	    obj_sel.ele_pt_min = 20;
	    if ( SelectedMuons(obj_sel, br).size() + SelectedEles(obj_sel, br).size() != 3 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("lep20"), std::string("lep20").length() );
	  }
	  if ( CAT_CUT.find("looseLepMVA") != std::string::npos ) {
	    obj_sel.mu_MVA_min  = -0.4;
	    obj_sel.ele_MVA_min = -0.4;
	    if ( SelectedMuons(obj_sel, br).size() + SelectedEles(obj_sel, br).size() != 3 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("looseLepMVA"), std::string("looseLepMVA").length() );
	  }
	  if ( CAT_CUT.find("medLepMVA") != std::string::npos ) {
	    obj_sel.mu_MVA_min  = 0.4;
	    obj_sel.ele_MVA_min = 0.4;
	    if ( SelectedMuons(obj_sel, br).size() + SelectedEles(obj_sel, br).size() != 3 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("medLepMVA"), std::string("medLepMVA").length() );
	  }
	  if ( CAT_CUT.find("muMeleT_lepMVA") != std::string::npos ) {
	    obj_sel.mu_MVA_min  = 0.4;
	    obj_sel.ele_MVA_min = 0.8;
	    if ( SelectedMuons(obj_sel, br).size() + SelectedEles(obj_sel, br).size() != 3 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("muMeleT_LepMVA"), std::string("muMeleT_LepMVA").length() );
	  }
	  if ( CAT_CUT.find("tightLepMVA") != std::string::npos ) {
	    obj_sel.mu_MVA_min  = 0.8;
	    obj_sel.ele_MVA_min = 0.8;
	    if ( SelectedMuons(obj_sel, br).size() + SelectedEles(obj_sel, br).size() != 3 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("tightLepMVA"), std::string("tightLepMVA").length() );
	  }
	  if ( CAT_CUT.find("tightLepCut") != std::string::npos ) {
	    if (        muH1_vec.Pt() < 20 || !muH1.isTightID || muH1.relIso > 0.12 ) { pass_cat_cut = false; continue; };
	    if (        muH2_vec.Pt() < 20 || !muH2.isTightID || muH2.relIso > 0.12 ) { pass_cat_cut = false; continue; };
	    if (  MU && (lep_vec.Pt() < 20 || !muW.isTightID  || muW.relIso > 0.12) ) { pass_cat_cut = false; continue; };
	    if ( !MU && (lep_vec.Pt() < 20 || !ele.isTightID  || ele.relIso > 0.12) ) { pass_cat_cut = false; continue; };
	    CAT_UNCUT.erase( CAT_UNCUT.find("tightLepCut"), std::string("tightLepCut").length() );
	  }
	  if ( CAT_CUT.find("ge2j") != std::string::npos ) {
	    if ( SelectedJets(obj_sel, br, "Central").size() < 2 )    { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("ge2j"), std::string("ge2j").length() );
	  }
	  if ( CAT_CUT.find("ge3j") != std::string::npos ) {
	    if ( SelectedJets(obj_sel, br, "Central").size() < 3 )    { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("ge3j"), std::string("ge3j").length() );
	  }
	  if ( CAT_CUT.find("btag") != std::string::npos ) {
	    if ( SelectedJets(obj_sel, br, "BTagMedium").size() < 1 &&
		 SelectedJets(obj_sel, br, "BTagLoose").size()  < 2 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("btag"), std::string("btag").length() );
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

	  // Compute lepMVA muon and electron efficiency scale factors
	  if ( CAT_CUT.find("LepMVA") != std::string::npos ) {
	    std::string SFm, SFe;
	    if ( CAT_CUT.find("looseLepMVA")    != std::string::npos ) { SFm = "mu_L"; SFe = "ele_L"; }
	    if ( CAT_CUT.find("medLepMVA")      != std::string::npos ) { SFm = "mu_M"; SFe = "ele_M"; }
	    if ( CAT_CUT.find("muMeleT_lepMVA") != std::string::npos ) { SFm = "mu_M"; SFe = "ele_T"; }
	    if ( CAT_CUT.find("tightLepMVA")    != std::string::npos ) { SFm = "mu_T"; SFe = "ele_T"; }

	    lepMVA_wgts[OPT_CUT+"_"+CAT_CUT]  = LepMVASF(lepSF[SFm], muH1.pt, muH1.eta);
	    lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] *= LepMVASF(lepSF[SFm], muH2.pt, muH2.eta);
	    if (MU) lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] *= LepMVASF(lepSF[SFm],  muW.pt, muW.eta);
	    else    lepMVA_wgts[OPT_CUT+"_"+CAT_CUT] *= LepMVASF(lepSF[SFe], ele.pt, ele.eta);
	  }

	  // Store the lepMVA weight separately (varies depending on category selection)
	  float lepMVA_wgt = (isData ? 1.0 : lepMVA_wgts[OPT_CUT+"_"+CAT_CUT]);
	  b_map_flt["lepMVA_wgt_OPT_"+OPT_CUT+"_CAT_"+CAT_CUT] = lepMVA_wgt;


	  //////////////////////////////////////
	  ///  Compute a few more variables  ///
	  //////////////////////////////////////

	  // Full event variables
	  JetInfos jets     = SelectedJets(obj_sel, br);
	  JetInfos jetsCent = SelectedJets(obj_sel, br, "Central");

	  TLorentzVector MET_vec     = FourVec(*br.met);
	  TLorentzVector lep_MET_vec = lep_vecT + MET_vec;
	  TLorentzVector evt_vec     = FourVec(muons, PTC, eles, jets);
	  TLorentzVector MHT_vec     = FourVec(muons, PTC, eles, jets, "T");
	  MHT_vec.RotateZ(TMath::Pi());
	  TLorentzVector lep_MHT_vec = lep_vecT + MHT_vec;

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

	  // Rank jets by highest b-tag score
	  int   bjet1_idx = -99;
	  int   bjet2_idx = -99;
	  float bjet1_CSV = -99;
	  float bjet2_CSV = -99;
	  TLorentzVector bjet1_vec;
	  TLorentzVector bjet2_vec;

	  for (int i = 0; i < jets.size(); i++) {
	    if (JetCSV(jets.at(i)) > bjet1_CSV) {
	      bjet2_idx = bjet1_idx;
	      bjet2_CSV = bjet1_CSV;
	      bjet2_vec = bjet1_vec;
	      bjet1_idx = i;
	      bjet1_CSV = max(JetCSV(jets.at(i)), float(-0.1));
	      bjet1_vec = FourVec(jets.at(i));
	    } else if (JetCSV(jets.at(i)) > bjet2_CSV) {
	      bjet2_idx = i;
	      bjet2_CSV = max(JetCSV(jets.at(i)), float(-0.1));
	      bjet2_vec = FourVec(jets.at(i));
	    }
	  }

	  ////////////////////////////////////////////////////////////////////
	  ///  Manual reconstruction of ttbar system in ttH signal events  ///
	  ////////////////////////////////////////////////////////////////////

	  // Guess a plausible ttbar reconstruction in ttH signal events
	  int MAN_top1_b_idx = -99;
	  int MAN_top2_b_idx = -99;
	  int MAN_topH_b_idx = -99;
	  int MAN_W_jet1_idx = -99;
	  int MAN_W_jet2_idx = -99;

	  float MAN_top1_b_CSV = -99;
	  float MAN_top2_b_CSV = -99;
	  float MAN_topH_b_CSV = -99;

	  TLorentzVector MAN_top1_b_vec;
	  TLorentzVector MAN_top2_b_vec;
	  TLorentzVector MAN_topH_b_vec;
	  TLorentzVector MAN_W_jet1_vec;
	  TLorentzVector MAN_W_jet2_vec;
	  TLorentzVector MAN_W_jj_vec;
	  TLorentzVector MAN_topH_vec;

	  // Associate loose b-tagged jets with top quark decays
	  if (bjet1_CSV > obj_sel.jet_btag_cuts.at(0)) {
	    MAN_top1_b_idx = bjet1_idx;
	    MAN_top1_b_CSV = bjet1_CSV;
	    MAN_top1_b_vec = FourVec(jets.at(MAN_top1_b_idx));
	  }
	  if (bjet2_CSV > obj_sel.jet_btag_cuts.at(0)) {
	    MAN_top2_b_idx = bjet2_idx;
	    MAN_top2_b_CSV = bjet2_CSV;
	    MAN_top2_b_vec = FourVec(jets.at(MAN_top2_b_idx));
	  }

	  // Pick the closest remaining pair in dR as the W candidate
	  float MAN_W_jj_dR = 99;
	  for (int i = 0; i < jetsCent.size(); i++) {
	    if (i == MAN_top1_b_idx || i == MAN_top2_b_idx)
	      continue;
	    TLorentzVector jv1 = FourVec(jetsCent.at(i));
	    for (int j = i+1; j < jetsCent.size(); j++) {
	      if (j == MAN_top1_b_idx || j == MAN_top2_b_idx)
		continue;
	      TLorentzVector jv2 = FourVec(jetsCent.at(j));
	      if (jv1.DeltaR(jv2) < MAN_W_jj_dR) {
		MAN_W_jj_dR    = jv1.DeltaR(jv2);
		MAN_W_jet1_idx = i;
		MAN_W_jet2_idx = j;
		MAN_W_jet1_vec = jv1;
		MAN_W_jet2_vec = jv2;
		MAN_W_jj_vec   = jv1+jv2;
	      }
	    }
	  }

	  // Pick the best bj(j) combination for an hadronic top
	  float top_mass = -9999;
	  // First consider the case with no selected di-jet W pair
	  if (MAN_W_jet1_idx < 0) {
	    for (int i = 0; i < jetsCent.size(); i++) {
	      if (i == MAN_top1_b_idx || i == MAN_top2_b_idx)
		continue;
	      TLorentzVector jv1 = FourVec(jetsCent.at(i));
	      if ( MAN_top1_b_idx > 0 && abs((MAN_top1_b_vec+jv1).M() - 173) < abs(top_mass - 173) ) {
		MAN_W_jet1_idx = i;
		MAN_W_jet1_vec = jv1;
		MAN_topH_b_idx = MAN_top1_b_idx;
		MAN_topH_b_CSV = MAN_top1_b_CSV;
		MAN_topH_b_vec = MAN_top1_b_vec;
		MAN_topH_vec   = MAN_top1_b_vec+MAN_W_jet1_vec;
	      }
	      if ( MAN_top2_b_idx > 0 && abs((MAN_top2_b_vec+jv1).M() - 173) < abs(top_mass - 173) ) {
		MAN_W_jet1_idx = i;
		MAN_W_jet1_vec = jv1;
		MAN_topH_b_idx = MAN_top2_b_idx;
		MAN_topH_b_CSV = MAN_top2_b_CSV;
		MAN_topH_b_vec = MAN_top2_b_vec;
		MAN_topH_vec   = MAN_top2_b_vec+MAN_W_jet1_vec;
	      }
	    }
	  } else { // Also consider the fully-reconstructed t --> bjj system
	    if ( MAN_top1_b_idx > 0 && abs((MAN_top1_b_vec+MAN_W_jj_vec).M() - 173) < abs(top_mass - 173) ) {
	      MAN_topH_b_idx = MAN_top1_b_idx;
	      MAN_topH_b_CSV = MAN_top1_b_CSV;
	      MAN_topH_b_vec = MAN_top1_b_vec;
	      MAN_topH_vec   = MAN_top1_b_vec+MAN_W_jj_vec;
	    }
	    if ( MAN_top2_b_idx > 0 && abs((MAN_top2_b_vec+MAN_W_jj_vec).M() - 173) < abs(top_mass - 173) ) {
	      MAN_topH_b_idx = MAN_top2_b_idx;
	      MAN_topH_b_CSV = MAN_top2_b_CSV;
	      MAN_topH_b_vec = MAN_top2_b_vec;
	      MAN_topH_vec   = MAN_top2_b_vec+MAN_W_jj_vec;
	    }
	  }


	  /////////////////////////////////////////////////////////////////
	  ///  BDT reconstruction of ttbar system in ttH signal events  ///
	  /////////////////////////////////////////////////////////////////

	  float BDT_score_jj_blv  = -1.1;  // Missing b-jet from hadronic top
	  float BDT_score_bjj_lv  = -1.1;  // Missing b-jet from leptonic top
	  float BDT_score_bj_blv  = -1.1;  // Missing one jet from hadronic W
	  float BDT_score_bjj_blv = -1.1;  // Full reconstruction
	  float BDT_score_max     = -1.1;  // Best BDT score out of all reconstructions
	  int   BDT_best_match    =    0;  // Best reconstruction of ttH system

	  // Find the most likely ttbar reconstruction in ttH signal events
	  int BDT_topL_b_idx = -1;
	  int BDT_topH_b_idx = -1;
	  int BDT_W_jet1_idx = -1;
	  int BDT_W_jet2_idx = -1;
	  int BDT_H_lep_idx  = -1; // 1 for muH_1, 2 for muH_2

	  // Input variables to reconstruction BDT
	  b_map_flt["topL_b_CSV"] = -1.0;
	  b_map_flt["topH_b_CSV"] = -1.0;
	  b_map_flt["topH_pt"]    =  0.0;
	  b_map_flt["W_jj_mass"]  =  0.0;
	  b_map_flt["topH_mass"]  =  0.0;
	  b_map_flt["top_Higgs_lep_pt_ratio"] = 0.0;
	  b_map_flt["lep_topL_b_dR"]  = -1.0;
	  b_map_flt["lep_topH_b_dR"]  = -1.0;
	  b_map_flt["muH_topL_b_dR"] = -1.0;

	  int max_iter = min(int(jetsCent.size()), 5);  // Only consider the 5 highest-pT central jets
	  for (int ibL = -1; ibL < max_iter; ibL++) {
	    if (ibL >= 0 && JetCSV(jetsCent.at(ibL)) < obj_sel.jet_btag_cuts.at(0)) continue;  // Only consider loose b-tagged jets
	    for (int ibH = -1; ibH < max_iter; ibH++) {
	      if (ibH >= 0 && ibH == ibL) continue;
	      if (ibH >= 0 && JetCSV(jetsCent.at(ibH)) < obj_sel.jet_btag_cuts.at(0)) continue;  // Only consider loose b-tagged jets
	      for (int ij1 = -1; ij1 < max_iter; ij1++) {
		if (ij1 >= 0 && (ij1 == ibL || ij1 == ibH)) continue;
		if (ij1 >= 0 && JetCSV(jetsCent.at(ij1)) > obj_sel.jet_btag_cuts.at(1)) continue;  // Only consider non-medium b-tagged jets
		for (int ij2 = -1; ij2 < max_iter; ij2++) {
		  if (ij2 >= 0 && (ij2 == ibL || ij2 == ibH || ij2 == ij1)) continue;
		  if (ij2 >= 0 && JetCSV(jetsCent.at(ij2)) > obj_sel.jet_btag_cuts.at(1)) continue;  // Only consider non-medium b-tagged jets
		  if (ij2 >= 0 && ij1 < 0) continue;                                                 // Second jet from W cannot be present if first is absent
		  if ((ibL < 0) + (ibH < 0) + (ij1 < 0) + (ij2 < 0) > 1) continue;                   // Require at least 3 matched jets

		  TLorentzVector W_jj_vec;
		  if (ij1 >= 0) W_jj_vec += FourVec(jetsCent.at(ij1));
		  if (ij2 >= 0) W_jj_vec += FourVec(jetsCent.at(ij2));
		  TLorentzVector topH_vec;
		  if (ibH >= 0) topH_vec += FourVec(jetsCent.at(ibH));
		  if (ij1 >= 0) topH_vec += W_jj_vec;

		  b_map_flt["topL_b_CSV"] = (ibL < 0 ? -1.0 : max(JetCSV(jetsCent.at(ibL), "CSVv2"), float(-1)));  // Use CSVv2 algorithm
		  b_map_flt["topH_b_CSV"] = (ibH < 0 ? -1.0 : max(JetCSV(jetsCent.at(ibH), "CSVv2"), float(-1)));  // Use CSVv2 algorithm
		  b_map_flt["topH_pt"]    = topH_vec.Pt();
		  b_map_flt["W_jj_mass"]  = (ij1 < 0 ?  0.0 : W_jj_vec.M());
		  b_map_flt["topH_mass"]  = topH_vec.M();
		  b_map_flt["lep_topL_b_dR"] = (ibL < 0 ? -1.0 : lep_vec.DeltaR( FourVec(jetsCent.at(ibL)) ) );
		  b_map_flt["lep_topH_b_dR"] = (ibH < 0 ? -1.0 : lep_vec.DeltaR( FourVec(jetsCent.at(ibH)) ) );
		  b_map_flt["muH_topL_b_dR"] = (ibL < 0 ? -1.0 : lep_vec.DeltaR( muH1_vec ) );
		  b_map_flt["top_Higgs_lep_pt_ratio"] = (lep_vec.Pt() / muH1_vec.Pt());

		  float BDT_score_1 = BDT_ttH_reco.Evaluate(b_map_flt, b_map_int);
		  b_map_flt["muH_topL_b_dR"] = (ibL < 0 ? -1.0 : lep_vec.DeltaR( muH2_vec ) );
		  b_map_flt["top_Higgs_lep_pt_ratio"] = (lep_vec.Pt() / muH2_vec.Pt());
		  float BDT_score_2 = BDT_ttH_reco.Evaluate(b_map_flt, b_map_int);
		  float BDT_score   = max(BDT_score_1, BDT_score_2);
		  int   BDT_match   = -99;

		  if (ibH >= 0 && ij1 >= 0 && ij2  < 0 && ibL >= 0 && BDT_score > BDT_score_bj_blv ) { BDT_score_bj_blv  = BDT_score; BDT_match = 2; }
		  if (ibH  < 0 && ij1 >= 0 && ij2 >= 0 && ibL >= 0 && BDT_score > BDT_score_jj_blv ) { BDT_score_jj_blv  = BDT_score; BDT_match = 4; }
		  if (ibH >= 0 && ij1 >= 0 && ij2 >= 0 && ibL  < 0 && BDT_score > BDT_score_bjj_lv ) { BDT_score_bjj_lv  = BDT_score; BDT_match = 6; }
		  if (ibH >= 0 && ij1 >= 0 && ij2 >= 0 && ibL >= 0 && BDT_score > BDT_score_bjj_blv) { BDT_score_bjj_blv = BDT_score; BDT_match = 8; }

		  if (BDT_score > BDT_score_max) {
		    ASSERT(BDT_match == 2 || BDT_match == 4 || BDT_match == 6 || BDT_match == 8, "BDT_match = 2, 4, 6, 8");
		    BDT_score_max  = BDT_score;
		    BDT_best_match = BDT_match + (BDT_score_2 > BDT_score_1);
		    BDT_topL_b_idx = ibL;
		    BDT_topH_b_idx = ibH;
		    BDT_W_jet1_idx = ij1;
		    BDT_W_jet2_idx = ij2;
		    BDT_H_lep_idx  = (BDT_score_1 > BDT_score_2 ? 1 : 2); // Tracks which lepton from the Higgs is better matched
		  } // End conditional: if (BDT_score > BDT_score_max)

		} // End loop: for (int ij2 = -1; ij2 < max_iter; ij2++)
	      } // End loop: for (int ij1 = -1; ij1 < max_iter; ij1++)
	    } // End loop: for (int ibH = -1; ibH < max_iter; ibH++)
	  } // End loop: for (int ibL = -1; ibL < max_iter; ibL++)

	  float BDT_topL_b_CSV = (BDT_topL_b_idx < 0 ? -99 : JetCSV(jetsCent.at(BDT_topL_b_idx)));
	  float BDT_topH_b_CSV = (BDT_topH_b_idx < 0 ? -99 : JetCSV(jetsCent.at(BDT_topH_b_idx)));

	  TLorentzVector zero_vec;
	  TLorentzVector BDT_topL_b_vec = (BDT_topL_b_idx < 0 ? zero_vec : FourVec(jetsCent.at(BDT_topL_b_idx)));
	  TLorentzVector BDT_topH_b_vec = (BDT_topH_b_idx < 0 ? zero_vec : FourVec(jetsCent.at(BDT_topH_b_idx)));
	  TLorentzVector BDT_W_jet1_vec = (BDT_W_jet1_idx < 0 ? zero_vec : FourVec(jetsCent.at(BDT_W_jet1_idx)));
	  TLorentzVector BDT_W_jet2_vec = (BDT_W_jet2_idx < 0 ? zero_vec : FourVec(jetsCent.at(BDT_W_jet2_idx)));
	  TLorentzVector BDT_W_jj_vec   = BDT_W_jet1_vec + BDT_W_jet2_vec;
	  TLorentzVector BDT_topH_vec   = BDT_topH_b_vec + BDT_W_jj_vec;


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

	  BookAndFill(tupI, "nJets",      13, -0.5, 12.5, SelectedJets(obj_sel, br).size(),               cat_evt_wgt );
	  BookAndFill(tupI, "nBJetsLoose", 6, -0.5,  5.5, SelectedJets(obj_sel, br, "BTagLoose").size(),  cat_evt_wgt );
	  BookAndFill(tupI, "nBJetsMed",   5, -0.5,  4.5, SelectedJets(obj_sel, br, "BTagMedium").size(), cat_evt_wgt );
	  BookAndFill(tupI, "nBJetsTight", 4, -0.5,  3.5, SelectedJets(obj_sel, br, "BTagTight").size(),  cat_evt_wgt );
	  BookAndFill(tupI, "nJetsCent",  11, -0.5, 10.5, SelectedJets(obj_sel, br, "Central").size(),    cat_evt_wgt );
	  BookAndFill(tupI, "nJetsFwd",    6, -0.5,  5.5, SelectedJets(obj_sel, br, "Forward").size(),    cat_evt_wgt );

	  obj_sel.jet_pt_min = 30;
	  BookAndFill(tupI, "nJetsCent30",  11, -0.5, 10.5, SelectedJets(obj_sel, br, "Central").size(),    cat_evt_wgt );
	  obj_sel.jet_pt_min = 40;
	  BookAndFill(tupI, "nJetsCent40",  11, -0.5, 10.5, SelectedJets(obj_sel, br, "Central").size(),    cat_evt_wgt );
	  obj_sel.jet_pt_min = 50;
	  BookAndFill(tupI, "nJetsCent50",  11, -0.5, 10.5, SelectedJets(obj_sel, br, "Central").size(),    cat_evt_wgt );
	  obj_sel.jet_pt_min = 20;

	  BookAndFill(tupF, "jet1_pt",  40,    0, 400, (jets.size() > 0 ? jets.at(0).pt                        : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet2_pt",  30,    0, 300, (jets.size() > 1 ? jets.at(1).pt                        : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet3_pt",  40,    0, 200, (jets.size() > 2 ? jets.at(2).pt                        : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet4_pt",  30,    0, 150, (jets.size() > 3 ? jets.at(3).pt                        : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet1_eta", 20, -5.0, 5.0, (jets.size() > 0 ? jets.at(0).eta                       : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet2_eta", 20, -5.0, 5.0, (jets.size() > 1 ? jets.at(1).eta                       : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet3_eta", 20, -5.0, 5.0, (jets.size() > 2 ? jets.at(2).eta                       : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet4_eta", 20, -5.0, 5.0, (jets.size() > 3 ? jets.at(3).eta                       : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet1_CSV", 22, -0.1, 1.0, (jets.size() > 0 ? max(JetCSV(jets.at(0)), float(-0.1)) : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet2_CSV", 22, -0.1, 1.0, (jets.size() > 1 ? max(JetCSV(jets.at(1)), float(-0.1)) : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet3_CSV", 22, -0.1, 1.0, (jets.size() > 2 ? max(JetCSV(jets.at(2)), float(-0.1)) : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "jet4_CSV", 22, -0.1, 1.0, (jets.size() > 3 ? max(JetCSV(jets.at(3)), float(-0.1)) : -9), cat_evt_wgt, false );

	  BookAndFill(tupF, "bjet1_pt",  40,    0, 200, (bjet1_idx > 0 ? bjet1_vec.Pt() : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "bjet2_pt",  40,    0, 200, (bjet2_idx > 0 ? bjet2_vec.Pt() : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "bjet1_CSV", 22, -0.1, 1.0, (bjet1_idx > 0 ? bjet1_CSV      : -9), cat_evt_wgt, false );
	  BookAndFill(tupF, "bjet2_CSV", 22, -0.1, 1.0, (bjet2_idx > 0 ? bjet2_CSV      : -9), cat_evt_wgt, false );

	  BookAndFill(tupF, "MET",  20, 0,  200, MET_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "MHT",  20, 0,  200, MHT_vec.Pt(), cat_evt_wgt );
	  BookAndFill(tupF, "MASS", 20, 0, 2000, evt_vec.M(),  cat_evt_wgt );

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

	  BookAndFill(tupF, "lep1_eta", 24, -2.4, 2.4, lep1_vec.Eta(), cat_evt_wgt );
	  BookAndFill(tupF, "lep2_eta", 24, -2.4, 2.4, lep2_vec.Eta(), cat_evt_wgt );
	  BookAndFill(tupF, "lep3_eta", 24, -2.4, 2.4, lep3_vec.Eta(), cat_evt_wgt );

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

	  // More useful variables for ttH vs. Z+jets and ttH vs. ttZ
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

	  // Manual top quark reconstruction variables in ttH signal events
	  BookAndFill(tupF, "MAN_top1_b_CSV", 22, -0.1,  1.0, (MAN_top1_b_idx > 0   ? MAN_top1_b_CSV      : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_top2_b_CSV", 22, -0.1,  1.0, (MAN_top2_b_idx > 0   ? MAN_top2_b_CSV      : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_topH_b_CSV", 22, -0.1,  1.0, (MAN_topH_b_idx > 0   ? MAN_topH_b_CSV      : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_top1_b_pt",  40,    0,  200, (MAN_top1_b_idx > 0   ? MAN_top1_b_vec.Pt() : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_top2_b_pt",  40,    0,  200, (MAN_top2_b_idx > 0   ? MAN_top1_b_vec.Pt() : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_topH_b_pt",  40,    0,  200, (MAN_topH_b_idx > 0   ? MAN_top1_b_vec.Pt() : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_W_jet1_pt",  30,    0,  300, (MAN_W_jet1_idx > 0   ? MAN_W_jet1_vec.Pt() : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_W_jet2_pt",  30,    0,  300, (MAN_W_jet2_idx > 0   ? MAN_W_jet2_vec.Pt() : -9), cat_evt_wgt, false);

	  BookAndFill(tupF, "MAN_topH_mass", 100,    0, 1000, (MAN_topH_vec.M() > 0 ? MAN_topH_vec.M()    : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_topH_pt",    50,    0,  500, (MAN_topH_vec.M() > 0 ? MAN_topH_vec.Pt()   : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_W_jj_mass", 100,    0, 1000, (MAN_W_jj_vec.M() > 0 ? MAN_W_jj_vec.M()    : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_W_jj_pt",    50,    0,  500, (MAN_W_jj_vec.M() > 0 ? MAN_W_jj_vec.Pt()   : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "MAN_W_jj_dR",    20,    0,  6.0, (MAN_W_jj_vec.M() > 0 ? MAN_W_jj_dR         : -9), cat_evt_wgt, false);

	  // BDT-based top quark reconstruction variables in ttH signal events
	  BookAndFill(tupF, "BDT_topL_b_CSV", 22, -0.1,  1.0, (BDT_topL_b_idx > 0   ? BDT_topL_b_CSV      : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_topH_b_CSV", 22, -0.1,  1.0, (BDT_topH_b_idx > 0   ? BDT_topH_b_CSV      : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_topL_b_pt",  40,    0,  200, (BDT_topL_b_idx > 0   ? BDT_topL_b_vec.Pt() : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_topH_b_pt",  40,    0,  200, (BDT_topH_b_idx > 0   ? BDT_topH_b_vec.Pt() : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_W_jet1_pt",  30,    0,  300, (BDT_W_jet1_idx > 0   ? BDT_W_jet1_vec.Pt() : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_W_jet2_pt",  30,    0,  300, (BDT_W_jet2_idx > 0   ? BDT_W_jet2_vec.Pt() : -9), cat_evt_wgt, false);

	  BookAndFill(tupF, "BDT_topH_mass", 100,    0, 1000, (BDT_topH_vec.M() > 0 ? BDT_topH_vec.M()    : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_topH_pt",    50,    0,  500, (BDT_topH_vec.M() > 0 ? BDT_topH_vec.Pt()   : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_W_jj_mass", 100,    0, 1000, (BDT_W_jj_vec.M() > 0 ? BDT_W_jj_vec.M()    : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_W_jj_pt",    50,    0,  500, (BDT_W_jj_vec.M() > 0 ? BDT_W_jj_vec.Pt()   : -9), cat_evt_wgt, false);

	  BookAndFill(tupF, "BDT_W_jj_dR",       20, 0, 6.0, (BDT_W_jet2_idx > 0 ? BDT_W_jet1_vec.DeltaR(BDT_W_jet2_vec) : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_lep_topL_b_dR", 20, 0, 6.0, (BDT_topL_b_idx > 0 ? lep_vec       .DeltaR(BDT_topL_b_vec) : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_lep_topH_b_dR", 20, 0, 6.0, (BDT_topH_b_idx > 0 ? lep_vec       .DeltaR(BDT_topH_b_vec) : -9), cat_evt_wgt, false);
	  if (BDT_H_lep_idx == 1)
	    BookAndFill(tupF, "BDT_muH_topL_b_dR", 20, 0, 6.0, (BDT_topL_b_idx > 0 ? muH1_vec    .DeltaR(BDT_topL_b_vec) : -9), cat_evt_wgt, false);
	  else if (BDT_H_lep_idx == 2)
	    BookAndFill(tupF, "BDT_muH_topL_b_dR", 20, 0, 6.0, (BDT_topL_b_idx > 0 ? muH2_vec    .DeltaR(BDT_topL_b_vec) : -9), cat_evt_wgt, false);
	  else
	    BookAndFill(tupF, "BDT_muH_topL_b_dR", 20, 0, 6.0,                                                             -99, cat_evt_wgt, false);

	  BookAndFill(tupF, "BDT_score_jj_blv",  22, -1.1, 1.0, BDT_score_jj_blv,  cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_score_bjj_lv",  22, -1.1, 1.0, BDT_score_bjj_lv,  cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_score_bj_blv",  22, -1.1, 1.0, BDT_score_bj_blv,  cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_score_bjj_blv", 22, -1.1, 1.0, BDT_score_bjj_blv, cat_evt_wgt, false);
	  BookAndFill(tupF, "BDT_score_max",     22, -1.1, 1.0, BDT_score_max,     cat_evt_wgt);
	  BookAndFill(tupI, "BDT_best_match",    11, -1.5, 9.5, BDT_best_match,    cat_evt_wgt);

	  // Top quark reconstruction variables in ttbar background events
	  BookAndFill(tupF, "bjet1_muSS_dR",   20, 0, 6.0, (bjet1_idx > 0 ? muSS_vec.DeltaR(bjet1_vec) : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "bjet1_lep_dR",    20, 0, 6.0, (bjet1_idx > 0 ?  lep_vec.DeltaR(bjet1_vec) : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "bjet1_muSS_pt",   30, 0, 300, (bjet1_idx > 0 ? (muSS_vec+bjet1_vec).Pt()  : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "bjet1_lep_pt",    30, 0, 300, (bjet1_idx > 0 ?  (lep_vec+bjet1_vec).Pt()  : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "bjet1_muSS_mass", 30, 0, 300, (bjet1_idx > 0 ? (muSS_vec+bjet1_vec).M()   : -9), cat_evt_wgt, false);
	  BookAndFill(tupF, "bjet1_lep_mass",  30, 0, 300, (bjet1_idx > 0 ?  (lep_vec+bjet1_vec).M()   : -9), cat_evt_wgt, false);

          // Extra variables especially for Andrew's BDTs
          BookAndFill( b_map_flt, out_tree, h_pre, "muH1_eta_abs", abs(b_map_flt["muH1_eta"]) );
          BookAndFill( b_map_flt, out_tree, h_pre, "muH2_eta_abs", abs(b_map_flt["muH2_eta"]) );

	  // Evaluate MVA output values
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_v1_all_withMass", BDT_v1_all_withMass.Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_v1_all",          BDT_v1_all         .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_v1_med",          BDT_v1_med         .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_v1_med_noBDT",    BDT_v1_med_noBDT   .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_v1_tight",        BDT_v1_tight       .Evaluate(b_map_flt, b_map_int) );
	  BookAndFill( b_map_flt, out_tree, h_pre, "BDT_v1_tight_noBDT",  BDT_v1_tight_noBDT .Evaluate(b_map_flt, b_map_int) );

	  // Use only even-numbered MC events to fill histograms
	  if ( (b_map_int["event"] % 2) == 0 || isData ) {
            float two = (isData ? 1.0 : 2.0);  // Scale MC by a factor of 2

	    BookAndFill(tupF, "BDT_v1_all_withMass", 40, -1, 1, b_map_flt["BDT_v1_all_withMass"], cat_evt_wgt*two);
	    BookAndFill(tupF, "BDT_v1_all",          20, -1, 1, b_map_flt["BDT_v1_all"],          cat_evt_wgt*two);
	    BookAndFill(tupF, "BDT_v1_med",          20, -1, 1, b_map_flt["BDT_v1_med"],          cat_evt_wgt*two);
	    BookAndFill(tupF, "BDT_v1_med_noBDT",    20, -1, 1, b_map_flt["BDT_v1_med_noBDT"],    cat_evt_wgt*two);
	    BookAndFill(tupF, "BDT_v1_tight",        20, -1, 1, b_map_flt["BDT_v1_tight"],        cat_evt_wgt*two);
	    BookAndFill(tupF, "BDT_v1_tight_noBDT",  20, -1, 1, b_map_flt["BDT_v1_tight_noBDT"],  cat_evt_wgt*two);

	    // Write dimuon mass histogram in different no-mass BDT bins
	    if      ( b_map_flt["BDT_v1_all"] >  0.5 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_all_p05_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_all"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_all_p02_p05_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    if      ( b_map_flt["BDT_v1_all"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_all_p02_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_all"] > -0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_all_n02_p02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_all_n10_n02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );

	    if      ( b_map_flt["BDT_v1_med"] >  0.5 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_p05_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_med"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_p02_p05_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    if      ( b_map_flt["BDT_v1_med"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_p02_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_med"] > -0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_n02_p02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_n10_n02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );

	    if      ( b_map_flt["BDT_v1_med_noBDT"] >  0.5 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_noBDT_p05_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_med_noBDT"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_noBDT_p02_p05_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    if      ( b_map_flt["BDT_v1_med_noBDT"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_noBDT_p02_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_med_noBDT"] > -0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_noBDT_n02_p02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_med_noBDT_n10_n02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );

	    if      ( b_map_flt["BDT_v1_tight"] >  0.5 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_p05_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_tight"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_p02_p05_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    if      ( b_map_flt["BDT_v1_tight"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_p02_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_tight"] > -0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_n02_p02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_n10_n02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );

	    if      ( b_map_flt["BDT_v1_tight_noBDT"] >  0.5 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_noBDT_p05_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_tight_noBDT"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_noBDT_p02_p05_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    if      ( b_map_flt["BDT_v1_tight_noBDT"] >  0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_noBDT_p02_p10_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else if ( b_map_flt["BDT_v1_tight_noBDT"] > -0.2 )
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_noBDT_n02_p02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );
	    else
	      BookAndFill( h_map_1D, h_pre+"H_pair_mass_BDT_v1_tight_noBDT_n10_n02_zoomH", 100, 110, 160, b_map_flt["H_pair_mass"], cat_evt_wgt*two );

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


  std::cout << "\nExiting ttH_3l()\n";
  
} // End void ttH_3l()
