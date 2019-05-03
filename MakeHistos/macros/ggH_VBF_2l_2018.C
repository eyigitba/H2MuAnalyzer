
/////////////////////////////////////////////////
///   Macro to study ggH and VBF categories   ///
///                                           ///
///      Andrew Brinkerhoff  04.16.2019       /// 
/////////////////////////////////////////////////

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

// #include "H2MuAnalyzer/MakeHistos/interface/SampleDatabase2016.h" // Input data and MC samples

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc700/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;      // Minimum index of input files to process
const int MAX_FILE = 1;      // Maximum index of input files to process
const int MAX_EVT  = 1000; // Maximum number of events to process
const int PRT_EVT  = 10;   // Print every N events
const float SAMP_WGT = 1.0;
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information

const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/SingleMuon/SingleMu_2018C/190426_125030/0000";
const TString SAMPLE   = "SingleMu";

const std::string YEAR  = "2018";
const std::string SLIM  = "notSlim";  // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 and 2018 NTuples are "Slim"
const TString OUT_DIR   = "plots";
const TString HIST_TREE = "HistTree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both

const std::vector<std::string> SEL_CUTS = {"Presel2018"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"2mu"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = { "NONE"};


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void ggH_VBF_2l_2018( TString sample = "", TString in_dir = "", TString out_dir = "",
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
    // in_file_name.Form("%s/NTuple_0.root", in_dir.Data());
    // std::cout << "Adding file " << in_file_name.Data() << std::endl;
    // in_file_names.push_back(in_file_name.Data());
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

  gROOT->cd(); // Navigate to "local" memory, so all histograms are not saved in out_tuple


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

    // For original 2016 and some 2017 or 2018 NTuples, convert "SlimJets" collection into regular jets
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

    // For each event, we set all the branch values to their default value (-99)
    for (std::map<TString, float>::iterator       it = b_map_flt.begin(); it != b_map_flt.end(); ++it) b_map_flt[it->first] = -99.0;
    for (std::map<TString, int>::iterator         it = b_map_int.begin(); it != b_map_int.end(); ++it) b_map_int[it->first] = -99;
    for (std::map<TString, std::string>::iterator it = b_map_str.begin(); it != b_map_str.end(); ++it) b_map_str[it->first] = "-99";


    // Get event weight for MC, defined in src/EventWeight.cc
    bool isData = sample.Contains("SingleMu");
    float event_wgt = ( isData ? 1.0 : EventWeight(br, evt_wgt, verbose) );

    // Initialize the selected Higgs candidate dimuon pair
    MuPairInfo        H_pair;
    TLorentzVector    H_pair_vec;


    /////////////////////////////////////////////////////////////////////////////////////////
    ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
    /////////////////////////////////////////////////////////////////////////////////////////
    
    bool pass_opt_cuts = false; // Check if event passes at least one optional cut
    bool pass_cat_cuts = false; // Check if event passes at least one category cut
    
    for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
      std::string OPT_CUT = OPT_CUTS.at(iOpt);
      assert(OPT_CUT == "2mu");

      // Assume event fails the optional cut
      BookAndFill(b_map_int, out_tree, sample, "OPT_"+OPT_CUT, 0 );

      if (OPT_CUT == "2mu") {

	// Require exactly 2 selected opposite-charge muons
	if ( SelectedMuPairs(obj_sel, br).size() != 1 ) continue;

	// Also require no extra muons or electrons used for the 3- or 4-lepton categories
	ConfigureObjectSelection(obj_sel, YEAR, "lepMVA");
	if ( SelectedMuPairs(obj_sel, br).size() > 1 ||
	     SelectedEles(obj_sel, br).size()    > 0 ) {
	  ConfigureObjectSelection(obj_sel, YEAR);
	  continue;
	} else {
	  ConfigureObjectSelection(obj_sel, YEAR);
	}
	
	// Require the Higgs candidate dimuon pair to fall inside the mass window
	H_pair     = SelectedMuPairs(obj_sel, br).at(0);
	H_pair_vec = FourVec(H_pair, PTC);
	if ( H_pair_vec.M() < 110 ||
	     H_pair_vec.M() > 160 ) continue;
      }
      pass_opt_cuts = true;
      
      // Now we know the event passes the optional cut
      b_map_int["OPT_"+OPT_CUT] = 1;
      
      
      //////////////////////////////////////////////////////////////////
      /// Loop through category cuts defined in src/CategoryCuts.cc  ///
      //////////////////////////////////////////////////////////////////
      
      for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {

	// Assume event fails the category cut
	BookAndFill(b_map_int, out_tree, sample+"_"+OPT_CUT, "CAT_"+CAT_CUTS.at(iCat), 0 );
	
	//////////////////////////////////////////////////////
	///  Compute variables relevent for category cuts  ///
	//////////////////////////////////////////////////////
	
	// Get selected muons and electrons in event
	MuonInfos muons = SelectedMuons(obj_sel, br);
	EleInfos  eles  = SelectedEles(obj_sel, br);

	assert( muons.size() == 2 && eles.size() == 0 );
	assert( muons.at(0).charge + muons.at(1).charge == 0 );

	MuonInfo mu1 = muons.at(0);
	MuonInfo mu2 = muons.at(1);

	TLorentzVector mu1_vec = FourVec(mu1, PTC);
	TLorentzVector mu2_vec = FourVec(mu2, PTC);
	TLorentzVector MET_vec = FourVec(*br.met);

	JetInfo        jet1;
	JetInfo        jet2;
	TLorentzVector jet1_vec;
	TLorentzVector jet2_vec;
	JetPairInfo    dijet;
	TLorentzVector dijet_vec;

	// Fake "jet pairs" made up of Higgs candidate pair and leading jet or dijet
	JetPairInfo H_jet1;
	JetPairInfo H_dijet;

	int nJets     = SelectedJets(obj_sel, br).size();
	if (nJets > 0) {
	  jet1     = SelectedJets(obj_sel, br).at(0);
	  jet1_vec = FourVec(jet1);
	  H_jet1   = MakeJetPair(H_pair_vec, jet1_vec);
	}
	if (nJets > 1) {
	  jet2      = SelectedJets(obj_sel, br).at(1);
	  jet2_vec  = FourVec(jet2);
	  dijet     = MakeJetPair(jet1_vec, jet2_vec);
	  dijet_vec = jet1_vec + jet2_vec;
	  H_dijet   = MakeJetPair(H_pair_vec, dijet_vec);
	}

	///////////////////////////////////////////
	///  Apply the category selection cuts  ///
	///////////////////////////////////////////

	std::string CAT_CUT   = CAT_CUTS.at(iCat);
	std::string CAT_UNCUT = CAT_CUT; // Track what sub-strings don't match any known cuts
	bool pass_cat_cut = true;

	// if ( CAT_CUT.find("MyCut") != std::string::npos ) {
	//   if ( not myCut ) { pass_cat_cut = false; continue; }
	//   CAT_UNCUT.erase( CAT_UNCUT.find("myCut"), std::string("myCut").length() );
	// }
	
	// Remove "_" characters left over after all known category sub-strings have been removed
	while ( CAT_UNCUT.find("_") != std::string::npos ) {
	  CAT_UNCUT.erase( CAT_UNCUT.find("_"), std::string("_").length() );
	}
	// If CAT_CUT contains unknown elements, look for it in the central repository
	if ( CAT_UNCUT.length() > 0 ) pass_cat_cut = InCategory(obj_sel, br, CAT_CUT, verbose);
	
	if (not pass_cat_cut) continue;
	
	if (verbose) std::cout << "\nPassed cut " << OPT_CUT << ", in category " << CAT_CUT << std::endl;
	std::string h_pre = (std::string) sample + "_"+OPT_CUT+"_"+CAT_CUT+"_";

	// Now we know the event passes the category cut
	b_map_int["CAT_"+CAT_CUT] = 1;
	pass_cat_cuts = true;


	///////////////////////////////////////////////////
	///  Generate and fill branches and histograms  ///
	///////////////////////////////////////////////////

	// Tuple containing maps and common options for "BookAndFill"
	std::tuple< const TString, std::map<TString, float> &, TTree *, std::map<TString, TH1*> &, const TString >
	  tupF{ hist_tree, b_map_flt, out_tree, h_map_1D, h_pre };
	std::tuple< const TString, std::map<TString, int> &,   TTree *, std::map<TString, TH1*> &, const TString >
	  tupI{ hist_tree, b_map_int, out_tree, h_map_1D, h_pre };
	
	// Store sample and event information
	BookAndFill(b_map_str, out_tree, h_pre, "sample", sample );
	BookAndFill(b_map_int, out_tree, h_pre, "event",  br.event->event );
	
	// Store event weights
	BookAndFill(tupF, "PU_wgt",    40, -2, 2, isData ? 1.0 : br.PU_wgt  );
	BookAndFill(tupF, "muon_wgt",  40, -2, 2, isData ? 1.0 : MuonWeight(br, evt_wgt, verbose) );
	BookAndFill(tupF, "GEN_wgt",   40, -2, 2, isData ? 1.0 : br.GEN_wgt );
	BookAndFill(tupF, "event_wgt", 40, -2, 2, isData ? 1.0 : event_wgt  );

	
	// Plot jet and full event kinematics
	BookAndFill(tupI, "nJets",       8, -0.5, 7.5, SelectedJets(obj_sel, br).size(),               event_wgt );
	BookAndFill(tupI, "nBJetsLoose", 6, -0.5, 5.5, SelectedJets(obj_sel, br, "BTagLoose").size(),  event_wgt );
	BookAndFill(tupI, "nBJetsMed",   4, -0.5, 3.5, SelectedJets(obj_sel, br, "BTagMedium").size(), event_wgt );
	BookAndFill(tupI, "nBJetsTight", 4, -0.5, 3.5, SelectedJets(obj_sel, br, "BTagTight").size(),  event_wgt );
	BookAndFill(tupI, "nJetsCent",   5, -0.5, 4.5, SelectedJets(obj_sel, br, "Central").size(),    event_wgt );
	BookAndFill(tupI, "nJetsFwd",    5, -0.5, 4.5, SelectedJets(obj_sel, br, "Forward").size(),    event_wgt );
	
	BookAndFill(tupF, "MET", 40, 0, 200, br.met->pt, event_wgt );
	BookAndFill(tupF, "MHT", 40, 0, 200, br.mht->pt, event_wgt );

	BookAndFill(tupF, "jet1_pt", 400, 0, 200, (nJets < 1 ? -1 : jet1_vec.Pt()), event_wgt, false );
	BookAndFill(tupF, "jet2_pt", 200, 0, 100, (nJets < 2 ? -1 : jet2_vec.Pt()), event_wgt, false );

	BookAndFill(tupF, "jet1_eta", 50, -5, 5, (nJets < 1 ? -6 : jet1_vec.Eta()), event_wgt, false );
	BookAndFill(tupF, "jet2_eta", 50, -5, 5, (nJets < 2 ? -6 : jet2_vec.Eta()), event_wgt, false );

	BookAndFill(tupF, "dijet_mass",  40,    0, 2000, (nJets < 2 ?  -1 : dijet.mass), event_wgt, false );
	BookAndFill(tupF, "dijet_pt",    50,    0,  500, (nJets < 2 ?  -1 : dijet.pt  ), event_wgt, false );
	BookAndFill(tupF, "dijet_eta",   50,   -5,    5, (nJets < 2 ?  -6 : dijet.eta ), event_wgt, false );
	BookAndFill(tupF, "dijet_phi",   32, -3.2,  3.2, (nJets < 2 ?  -4 : dijet.phi ), event_wgt, false );
	BookAndFill(tupF, "dijet_dR",   100,  0.0, 10.0, (nJets < 2 ?  -1 : dijet.dR  ), event_wgt, false );
	BookAndFill(tupF, "dijet_dEta", 100,  -10, 10.0, (nJets < 2 ? -11 : dijet.dEta), event_wgt, false );
	BookAndFill(tupF, "dijet_dPhi",  64, -3.2,  3.2, (nJets < 2 ?  -7 : dijet.dPhi), event_wgt, false );
	
	// Plot muon and di-muon kinematics
	BookAndFill(tupF, "mu1_pt", 60, 0, 300, mu1_vec.Pt(), event_wgt );
	BookAndFill(tupF, "mu2_pt", 40, 0, 200, mu2_vec.Pt(), event_wgt );
	
	BookAndFill(tupF, "mu1_eta", 48, -2.4, 2.4, mu1_vec.Eta(), event_wgt );
	BookAndFill(tupF, "mu2_eta", 48, -2.4, 2.4, mu2_vec.Eta(), event_wgt );
	
	BookAndFill(tupI, "mu1_tightID", 2, -0.5, 1.5, mu1.isTightID, event_wgt );
	BookAndFill(tupI, "mu2_tightID", 2, -0.5, 1.5, mu2.isTightID, event_wgt );

	BookAndFill(tupF, "mu1_iso", 30, 0, 0.3, mu1.relIso, event_wgt );
	BookAndFill(tupF, "mu2_iso", 30, 0, 0.3, mu2.relIso, event_wgt );
	
	BookAndFill(tupF, "mu1_d0", 100, -0.05, 0.05, mu1.d0_PV, event_wgt );
	BookAndFill(tupF, "mu2_d0", 100, -0.05, 0.05, mu2.d0_PV, event_wgt );

	BookAndFill(tupF, "mu1_dz", 100, -0.1, 0.1, mu1.dz_PV, event_wgt );
	BookAndFill(tupF, "mu2_dz", 100, -0.1, 0.1, mu2.dz_PV, event_wgt );

	BookAndFill(tupF, "mu1_miniIso", 30, 0, 0.3, mu1.miniIso, event_wgt );
	BookAndFill(tupF, "mu2_miniIso", 30, 0, 0.3, mu2.miniIso, event_wgt );
	
	BookAndFill(tupF, "mu1_SIP", 100, 0, 10, mu1.SIP_3D, event_wgt );
	BookAndFill(tupF, "mu2_SIP", 100, 0, 10, mu2.SIP_3D, event_wgt );
	
	BookAndFill(tupF, "mu1_segCompat", 100, 0, 1, mu1.segCompat, event_wgt );
	BookAndFill(tupF, "mu2_segCompat", 100, 0, 1, mu2.segCompat, event_wgt );

	BookAndFill(tupF, "mu1_lepMVA", 40, -1, 1, mu1.lepMVA, event_wgt );
	BookAndFill(tupF, "mu2_lepMVA", 40, -1, 1, mu2.lepMVA, event_wgt );
	
	BookAndFill(tupF, "H_mass_zoom", 50, 110, 160, H_pair_vec.M(),  event_wgt );
	BookAndFill(tupF, "H_mass_on",   20, 110, 160, H_pair_vec.M(),  event_wgt );
	BookAndFill(tupF, "H_mass_err",  50,   0,  10, H_pair.massErr,  event_wgt );
	BookAndFill(tupF, "H_pt",        60,   0, 300, H_pair_vec.Pt(), event_wgt );

	BookAndFill(tupF, "dimu_eta",   40, -8.0, 8.0, H_pair_vec.Eta(),      event_wgt );
	BookAndFill(tupF, "dimu_phi",   32, -3.2, 3.2, H_pair_vec.Phi(),      event_wgt );
	BookAndFill(tupF, "dimu_rapid", 30, -3.0, 3.0, H_pair_vec.Rapidity(), event_wgt );

	BookAndFill(tupF, "dimu_dR",   60,  0.0, 6.0, H_pair.dR,   event_wgt );
	BookAndFill(tupF, "dimu_dEta", 48, -4.8, 4.8, H_pair.dEta, event_wgt );
	BookAndFill(tupF, "dimu_dPhi", 32, -3.2, 3.2, H_pair.dPhi, event_wgt );

	BookAndFill(tupF, "dimu_cos_theta_star", 40, -1, 1, CosThetaStar(mu1_vec, mu2_vec), event_wgt );
	BookAndFill(tupF, "dimu_cos_CS_theta",   40, -1, 1, Cos_CS_Theta(mu1_vec, mu2_vec), event_wgt );
	BookAndFill(tupF, "dimu_cos_CS_phi",     40, -1, 1, Cos_CS_Phi  (mu1_vec, mu2_vec), event_wgt );
	BookAndFill(tupF, "dimu_sin_CS_phi",     40,  0, 1, Sin_CS_Phi  (mu1_vec, mu2_vec), event_wgt );

	// Plot dimuon-jet kinematics
	BookAndFill(tupF, "H_jet1_dR",   40,  0.0, 8.0, (nJets < 1 ? -1 : H_jet1.dR  ), event_wgt, false );
	BookAndFill(tupF, "H_jet1_dEta", 40, -8.0, 8.0, (nJets < 1 ? -9 : H_jet1.dEta), event_wgt, false );
	BookAndFill(tupF, "H_jet1_dPhi", 32, -3.2, 3.2, (nJets < 1 ? -4 : H_jet1.dPhi), event_wgt, false );

	BookAndFill(tupF, "H_jet1_cos_theta_star", 40, -1, 1, (nJets < 1 ? -2 : CosThetaStar(H_pair_vec, jet1_vec)), event_wgt, false );
	BookAndFill(tupF, "H_jet1_cos_CS_theta",   40, -1, 1, (nJets < 1 ? -2 : Cos_CS_Theta(H_pair_vec, jet1_vec)), event_wgt, false );
	BookAndFill(tupF, "H_jet1_cos_CS_phi",     40, -1, 1, (nJets < 1 ? -2 : Cos_CS_Phi  (H_pair_vec, jet1_vec)), event_wgt, false );
	BookAndFill(tupF, "H_jet1_sin_CS_phi",     40,  0, 1, (nJets < 1 ? -1 : Sin_CS_Phi  (H_pair_vec, jet1_vec)), event_wgt, false );

	BookAndFill(tupF, "H_dijet_dR",   40,  0.0, 8.0, (nJets < 2 ? -1 : H_dijet.dR  ), event_wgt, false );
	BookAndFill(tupF, "H_dijet_dEta", 40, -8.0, 8.0, (nJets < 2 ? -9 : H_dijet.dEta), event_wgt, false );
	BookAndFill(tupF, "H_dijet_dPhi", 32, -3.2, 3.2, (nJets < 2 ? -4 : H_dijet.dPhi), event_wgt, false );

	BookAndFill(tupF, "H_dijet_cos_theta_star", 40, -1, 1, (nJets < 2 ? -2 : CosThetaStar(H_pair_vec, dijet_vec)), event_wgt, false );
	BookAndFill(tupF, "H_dijet_cos_CS_theta",   40, -1, 1, (nJets < 2 ? -2 : Cos_CS_Theta(H_pair_vec, dijet_vec)), event_wgt, false );
	BookAndFill(tupF, "H_dijet_cos_CS_phi",     40, -1, 1, (nJets < 2 ? -2 : Cos_CS_Phi  (H_pair_vec, dijet_vec)), event_wgt, false );
	BookAndFill(tupF, "H_dijet_sin_CS_phi",     40,  0, 1, (nJets < 2 ? -1 : Sin_CS_Phi  (H_pair_vec, dijet_vec)), event_wgt, false );
	
	
      } // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
    } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)

    // Fill branches in output tuple tree, if event passed at least one category
    if (pass_opt_cuts && pass_cat_cuts) {
      out_tree->Fill();
    }

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
	if (h_name.find(optCatStr+"_") != std::string::npos) {
	  // Remove optional selection and category cuts from histogram names
	  h_name.erase( h_name.find(optCatStr+"_"), optCatStr.length() + 1 );
	  it->second->SetName(h_name.c_str());
	  it->second->SetTitle(h_name.c_str());
	  // std::cout << "  * Writing 1D histogram " << it->second->GetName() << std::endl;
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
  
  
  std::cout << "\nExiting ggH_VBF_2l_2018()\n";
  
} // End void ggH_VBF_2l_2018()
