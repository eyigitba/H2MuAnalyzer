/////////////////////////////////////////////////////////////
///   Macro to make minintuple for faster BDT training    ///
///         make minintuple in selected categories        ///
///                Xunwu Zuo  08.11.2018                  ///
/////////////////////////////////////////////////////////////


//make map of trees, each contain different branches targeted for training in that category
//each tree writen in corresponding out root file
//everything defined before entry loop, calculated in entry loop, cat loop


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
#include "H2MuAnalyzer/MakeHistos/interface/SampleID.h"           // assign number ID to samples
#include "H2MuAnalyzer/MakeHistos/interface/ReadMVA.h"            // Read and evaluate XMLs for MVA

// #include "H2MuAnalyzer/MakeHistos/interface/SampleDatabase2016.h" // Input data and MC samples

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 1;     // Maximum index of input files to process
const int MAX_EVT  = 100000;    // Maximum number of events to process
const int PRT_EVT  = 1000;  // Print every N events
const float SAMP_WGT = 1.0;
const float LUMI = 41000; // pb-1   36814 for 2016, 41000 for 2017
const bool verbose = false; // Print extra information


const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/WminusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8/H2Mu_WH_neg_125/190115_144038/0000";
const TString SAMPLE   = "H2Mu_WH_neg_125";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos/170315_105045/0000";
//const TString SAMPLE   = "H2Mu_WH_pos";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WZ_3l_AMC";
//const TString SAMPLE   = "WZ_3l_AMC";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
//const TString SAMPLE   = "SingleMu";
const std::string YEAR = "2017";
const std::string SLIM = "Slim"; // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR  = "plots";
const TString HIST_TREE = "Tree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both. Not in use in this macro

const std::vector<std::string> SEL_CUTS = {"Presel2017"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"WH_3l_mu"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void MiniNTupliser( TString sample = "", TString in_dir = "", TString out_dir = "",
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
      //in_file_name.Form("%s/NTuple_0.root", in_dir.Data());
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

  // creating output file

  TString out_file_name;
  if (out_file_str.Length() > 0) out_file_name.Form( "%s/histos_%s_%s.root",    out_dir.Data(), sample.Data(), out_file_str.Data() );
  else                           out_file_name.Form( "%s/histos_%s_%d_%d.root", out_dir.Data(), sample.Data(), MIN_FILE, MAX_FILE );
  std::cout << "\nCreating output file " << out_file_name.Data() << std::endl;
  TFile * Out_File = TFile::Open( out_file_name, "RECREATE" );

  // plant TTree and add branches
  TTree * Out_Tree = PlantTree("tree", "flat_tree_minintuple");    

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

  evt_sel.muPair_mass_min = 105; // Require at least one Higgs candidate pair, default 60
// use default 60 GeV for WZ validation
  obj_sel.mu_pt_min       =  10; // Lower muon pT threshold for muons not from Higgs, default 20
  obj_sel.mu_iso_max	  = 0.4;
  obj_sel.muPair_Higgs = "sort_OS_dimuon_pt"; // alternate selection, used by Hamburg group
  //obj_sel.muPair_Higgs = "sort_WH_3_mu_v1"; //for 3mu category only

  //obj_sel.ele_pt_min = 20;
  obj_sel.ele_ID_cut = "loose";
  obj_sel.ele_iso_max = 0.4;

  if (verbose) obj_sel.Print();
  if (verbose) evt_sel.Print();
  if (verbose) evt_wgt.Print();

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later

  std::cout << "\n******* About to load XML files for signal-background BDTs *******" << std::endl;
  MVA::MVA BDT_noMass_old( "data/XMLs/WH_3l/Xunwu/2019_04_30/2017_WH_ele_without_mass_all_sig_all_bkg_ge0j/",
                           "weights/2017_WH_ele_without_mass_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml",
                           "BDTG_UF_v2" );
  MVA::MVA BDT_mass_old  ( "data/XMLs/WH_3l/Xunwu/2019_04_30/2017_WH_ele_with_mass_all_sig_all_bkg_ge0j/",
                           "weights/2017_WH_ele_with_mass_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml",
                           "BDTG_UF_v2" );

  MVA::MVA BDT_noMass_new( "data/XMLs/WH_3l/Xunwu/2019_05_12/2017_WH_ele_against_inclu_trimvar_all_sig_all_bkg_ge0j/",
                           "weights/2017_WH_ele_against_inclu_trimvar_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml",
                           "BDTG_UF_v2" );
  MVA::MVA BDT_mass_new  ( "data/XMLs/WH_3l/Xunwu/2019_05_12/2017_WH_ele_against_inclu_trimvar_with_mass_all_sig_all_bkg_ge0j/",
                           "weights/2017_WH_ele_against_inclu_trimvar_with_mass_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml",
                           "BDTG_UF_v2" );


  std::cout << "\n******* About to enter the event loop " << in_chain->GetEntries() << " events *******" << std::endl;
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

    if (pass_sel_cuts) {
      
      // Get event weight for MC, defined in src/EventWeight.cc
      float event_wgt = ( sample.Contains("SingleMu") ? 1.0 : EventWeight(br, evt_wgt, verbose) );
      float xsec_norm = samp_weight;
      int   Sample_ID = getSampleID(sample);

      /////////////////////////////////////////////////////////////////////////////////////////
      ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
      /////////////////////////////////////////////////////////////////////////////////////////
      ///  in this macro opt selection act as additional selection, can only have one   - XWZ 12.11.2018

      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	std::string OPT_CUT = OPT_CUTS.at(iOpt);
	std::string h_pre = (std::string) sample + "_"+OPT_CUT;
	///////////////////////////////////////////////////////////////
	///  select objects and get variables to be filled in tree  ///
	///////////////////////////////////////////////////////////////

	for (std::map<TString, float>::iterator       it = b_map_flt.begin(); it != b_map_flt.end(); ++it) b_map_flt[it->first] = -99.0;
        for (std::map<TString, int>::iterator         it = b_map_int.begin(); it != b_map_int.end(); ++it) b_map_int[it->first] = -99;
        for (std::map<TString, std::string>::iterator it = b_map_str.begin(); it != b_map_str.end(); ++it) b_map_str[it->first] = "-99";
	// To make sure that all branches in the tree are filled in the first event, 
	// it is safer to make all the BookAndFill method have the same indentation.
	// If some of them are in a condition block, make sure that condition(s) covers all possibilities
	// All variables to be filled needs to be properly initiated  - XWZ 14.05.2019


	MuPairInfo    dimu;	  // object initial values are all -999
  	MuonInfo      mu_1;
  	MuonInfo      mu_2;
  	MuonInfo      extra_mu;
  	EleInfo       ele;
  	JetPairInfo   dijet;
  	JetInfo       jet1;
  	JetInfo       jet2;
  	JetInfo       jet0;
  	MhtInfo       mht_info;

  	TLorentzVector dimu_vec;  // vector default initial is (0,0,0,0)
  	TLorentzVector mu1_vec;
  	TLorentzVector mu2_vec;
  	TLorentzVector extra_mu_vec;
  	TLorentzVector lep_vec;
	TLorentzVector lep_trans_vec;  
  	TLorentzVector met_vec;
  	TLorentzVector mht_vec;
  	TLorentzVector mlt_vec;
  	TLorentzVector lMET_vec;
  	TLorentzVector lMHT_vec;
  	TLorentzVector lMLT_vec;
	
	// select dimuon candidate
	bool have_Z_mass = false;
	for (const auto & muPair : SelectedMuPairs(obj_sel, br)) {
	    if (muPair.charge == 0 and muPair.mass > 81 and muPair.mass < 101) have_Z_mass = true;
	} 
	if (have_Z_mass) continue;

	dimu = SelectedCandPair(obj_sel, br);
        dimu_vec = FourVec( dimu, PTC);
        if ( dimu_vec.M() < 105 ||
             dimu_vec.M() > 160 ) continue;  // 70-110 for Z mass validation, 105-160 for analysis window

	MuonInfos muons = SelectedMuons(obj_sel, br);
        EleInfos  eles  = SelectedEles(obj_sel, br);
        JetInfos  jets  = SelectedJets(obj_sel, br);

	// ************ category cuts and fill lepton variables ***************
	int lep_ID = 0;  // default value is 0
	int lep_charge = 0;
	/////////////////////////
	///   WH_3l_ele cat   ///
	/////////////////////////
	if (OPT_CUT == "WH_3l_ele") {
	  if (muons.size() != 2 or eles.size() != 1) continue;
	  if ( not (ElePass(obj_sel, eles.at(0), br) && SelectedJets(obj_sel, br, "BTagMedium").size() == 0) ) continue;
          //if ( not (ElePass(obj_sel, eles.at(0), br) && nBJets_Med > 0) )  continue;   // for ttH 3l 

          ele  = eles.at(0);
	  lep_ID = 11;
	  lep_charge = ele.charge;
          lep_vec = FourVec(ele);
	  lep_trans_vec = FourVec(ele,"T");

	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep_pt", lep_vec.Pt() );
          BookAndFill( b_map_flt, Out_Tree, h_pre, "lep_abs_eta", abs(lep_vec.Eta()) );
          BookAndFill( b_map_flt, Out_Tree, h_pre, "lep_lepMVA", ele.lepMVA );
          BookAndFill( b_map_int, Out_Tree, h_pre, "lep_charge", ele.charge );	
	  BookAndFill( b_map_int, Out_Tree, h_pre, "lep_ID", lep_ID);
	} // if (OPT_CUT == "WH_3l_ele") 	

        /////////////////////////
        ///   WH_3l_mu cat    ///
        /////////////////////////
        if (OPT_CUT == "WH_3l_mu") {
          if (eles.size() != 0 or muons.size() != 3) continue;
	  bool have_muon_cand = false;
	  for (int imu = 0; imu < int((br.muons)->size()) ; imu++) {
	    if (have_muon_cand) continue;
	    if ( imu != dimu.iMu1 and imu != dimu.iMu2 and MuonPass(obj_sel, br.muons->at(imu), br) ) {
	      extra_mu = br.muons->at(imu);
	      have_muon_cand = true;
	    }
	  } 
	  if (not have_muon_cand) continue;
	  if ( not (MuonPass(obj_sel, extra_mu, br) && SelectedJets(obj_sel, br, "BTagMedium").size() == 0) ) continue;
	  //if ( not (MuonPass(obj_sel, extra_mu, br) && nBJets_Med > 0) )  continue; // for ttH 3l

	  lep_ID = 13;
	  lep_charge = extra_mu.charge;
	  lep_vec = FourVec(extra_mu, PTC);
	  lep_trans_vec = FourVec(extra_mu,PTC,"T");

	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep_pt", lep_vec.Pt() );
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep_abs_eta", abs(lep_vec.Eta()) );
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep_lepMVA", extra_mu.lepMVA );
	  BookAndFill( b_map_int, Out_Tree, h_pre, "lep_charge", extra_mu.charge );
	  BookAndFill( b_map_int, Out_Tree, h_pre, "lep_ID", lep_ID);;
	} // if (OPT_CUT == "WH_3l_mu")

	if (lep_charge == 0) {
	  std::cout << "weird case: lepton charge = 0" << std::endl;
	  continue;
	}
	
	// ********* categories done, lep variables filled *********
	// now must already be in one of the categories

	// ***************** auxiliary variables *****************
	BookAndFill( b_map_flt, Out_Tree, h_pre, "event_wgt",	event_wgt);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "xsec_norm",   xsec_norm);
	BookAndFill( b_map_int, Out_Tree, h_pre, "Sample_ID",   Sample_ID);
	BookAndFill( b_map_str, Out_Tree, h_pre, "Sample_name",    sample);

	// ****************** muon variables *******************
	int dimu_gen_ID = 0; // default value is 0
	if (not sample.Contains("SingleMu")) {
	    if ( IsGenMatched(dimu, *br.muons, *br.genMuons, "H") ) dimu_gen_ID = 25;
	    else if ( IsGenMatched(dimu, *br.muons, *br.genMuons, "Z") ) dimu_gen_ID = 23;
	}
	mu_1 = br.muons->at(dimu.iMu1);
        mu_2 = br.muons->at(dimu.iMu2);
        mu1_vec = FourVec(br.muons->at(dimu.iMu1), PTC);
        mu2_vec = FourVec(br.muons->at(dimu.iMu2), PTC);
	
	BookAndFill( b_map_int, Out_Tree, h_pre, "dimu_gen_ID",         dimu_gen_ID                     	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_mass", 		dimu_vec.M()				);
 	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_mass_err", 	dimu.massErr				);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_pt", 		dimu_vec.Pt()				);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_abs_eta", 	abs(dimu_vec.Eta())			);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_abs_dEta", 	abs(mu1_vec.Eta() - mu2_vec.Eta())	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_abs_dPhi", 	abs(mu1_vec.DeltaPhi(mu2_vec))		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_dR", 		mu1_vec.DeltaR(mu2_vec)			);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "mu1_pt", 		mu1_vec.Pt()		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mu1_abs_eta", 	abs(mu1_vec.Eta())	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mu1_lepMVA", 		mu_1.lepMVA		);
	BookAndFill( b_map_int, Out_Tree, h_pre, "mu1_charge", 		mu_1.charge		);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "mu2_pt", 		mu2_vec.Pt()		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mu2_abs_eta", 	abs(mu2_vec.Eta())	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mu2_lepMVA", 		mu_2.lepMVA		);
	BookAndFill( b_map_int, Out_Tree, h_pre, "mu2_charge", 		mu_2.charge		);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_mu1", 		CosThetaStar(mu1_vec, mu2_vec)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_mu_pos", 		mu_1.charge == 1 ? CosThetaStar(mu1_vec, mu2_vec) : CosThetaStar(mu2_vec, mu1_vec)  );

	// ***************** event variables ******************
	BookAndFill( b_map_int, Out_Tree, h_pre, "nMuons", 		muons.size()		);
	BookAndFill( b_map_int, Out_Tree, h_pre, "nEles",		eles.size()		);
	BookAndFill( b_map_int, Out_Tree, h_pre, "nJets", 		jets.size()		);
	BookAndFill( b_map_int, Out_Tree, h_pre, "nFwdJets", 		SelectedJets(obj_sel, br, "Forward").size()	);
	BookAndFill( b_map_int, Out_Tree, h_pre, "nCentJets", 		SelectedJets(obj_sel, br, "Central").size()	);
	BookAndFill( b_map_int, Out_Tree, h_pre, "nBJets_Loose", 	SelectedJets(obj_sel, br, "BTagLoose").size()	);
	BookAndFill( b_map_int, Out_Tree, h_pre, "nBJets_Med", 		SelectedJets(obj_sel, br, "BTagMedium").size()	);
	BookAndFill( b_map_int, Out_Tree, h_pre, "nBJets_Tight", 	SelectedJets(obj_sel, br, "BTagTight").size()	);

	// *************** jet variables *******************
	if (SelectedJetPairs(obj_sel, br).size() > 0) {
	  dijet = SelectedJetPairs(obj_sel, br).at(0);      // no need for jet vector since no PTC
	  jet1  = br.jets->at(dijet.iJet1);
	  jet2  = br.jets->at(dijet.iJet2);
	}
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_mass", 		dijet.mass	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_pt", 		dijet.pt	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_abs_eta", 	abs(dijet.eta)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_abs_dEta", 	abs(dijet.dEta)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_abs_dPhi", 	abs(dijet.dPhi)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_dR", 		dijet.dR);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "jet1_pt", 		jet1.pt		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "jet1_abs_eta", 	abs(jet1.eta)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "jet2_pt", 		jet2.pt		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "jet2_abs_eta", 	abs(jet2.eta)	);
	
        if (SelectedJets(obj_sel, br).size() > 0) {
	  jet0 = SelectedJets(obj_sel, br).at(0);
	}
	BookAndFill( b_map_flt, Out_Tree, h_pre, "jet0_pt", 		jet0.pt		);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "jet0_abs_eta",	abs(jet0.eta)	);

	// *************** ldimu variables ******************
	TLorentzVector ldimu_vec = lep_vec + dimu_vec;
	BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_ldimu", 		CosThetaStar(lep_vec, dimu_vec)		);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "ldimu_mass", 		ldimu_vec.M()				);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "ldimu_pt", 		ldimu_vec.Pt()				);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "ldimu_abs_eta", 	abs(ldimu_vec.Eta())			);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "ldimu_abs_dEta", 	abs(lep_vec.Eta() - dimu_vec.Eta())	);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "ldimu_abs_dPhi", 	abs(lep_vec.DeltaPhi(dimu_vec))		);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "ldimu_dR", 		lep_vec.DeltaR(dimu_vec)		);

	TLorentzVector lmu1_vec = lep_vec + mu1_vec;
        TLorentzVector lmu2_vec = lep_vec + mu2_vec;
        if (lep_charge == mu_1.charge) {
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_lmuSS", 	CosThetaStar(lep_vec, mu1_vec)		);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_lmuOS", 	CosThetaStar(lep_vec, mu2_vec)		);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_pt", 	lmu1_vec.Pt()				);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_abs_eta", 	abs(lmu1_vec.Eta())			);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_abs_dEta", 	abs(lep_vec.Eta() - mu1_vec.Eta())	);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_abs_dPhi", 	abs(lep_vec.DeltaPhi(mu1_vec))		);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_dR", 	lep_vec.DeltaR(mu1_vec)			);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_mass", 	lmu1_vec.M()				);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_pt", 	lmu2_vec.Pt()				);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_abs_eta", 	abs(lmu2_vec.Eta())			);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_abs_dEta", 	abs(lep_vec.Eta() - mu2_vec.Eta())	);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_abs_dPhi", 	abs(lep_vec.DeltaPhi(mu2_vec))		);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_dR", 	lep_vec.DeltaR(mu2_vec)			);
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_mass", 	lmu2_vec.M()				);
        }
        else {
	    BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_lmuSS",       CosThetaStar(lep_vec, mu2_vec)  	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_lmuOS",       CosThetaStar(lep_vec, mu1_vec)  	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_pt",        lmu2_vec.Pt()                   	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_abs_eta",   abs(lmu2_vec.Eta())                  	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_abs_dEta",  abs(lep_vec.Eta() - mu2_vec.Eta())   	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_abs_dPhi",  abs(lep_vec.DeltaPhi(mu2_vec))       	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_dR",        lep_vec.DeltaR(mu2_vec)         	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuSS_mass",      lmu2_vec.M()                    	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_pt",        lmu1_vec.Pt()                   	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_abs_eta",   abs(lmu1_vec.Eta())                  	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_abs_dEta",  abs(lep_vec.Eta() - mu1_vec.Eta())   	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_abs_dPhi",  abs(lep_vec.DeltaPhi(mu1_vec))       	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_dR",        lep_vec.DeltaR(mu1_vec)         	);
            BookAndFill( b_map_flt, Out_Tree, h_pre, "lmuOS_mass",      lmu1_vec.M()                    	);
        }

	// ***************** met variables ***************
    	met_vec  = FourVec(*br.met);
	mht_info = *br.mht;
	mht_vec  = FourVec(mht_info);

	mlt_vec  = -lep_trans_vec - FourVec(dimu,PTC,"T",*br.muons);
        lMET_vec = lep_trans_vec + met_vec;
        lMHT_vec = lep_trans_vec + mht_vec;
        lMLT_vec = lep_trans_vec + mlt_vec;

	BookAndFill( b_map_flt, Out_Tree, h_pre, "met_pt", 		met_vec.Pt()			);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mt_lmet", 		lMET_vec.M()			);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "abs_dPhi_lmet", 	abs(lep_vec.DeltaPhi(met_vec))	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mht_pt", 		mht_vec.Pt()			);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mht_mass", 		mht_info.mass			);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mt_lmht", 		lMHT_vec.M()			);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "abs_dPhi_lmht", 	abs(lep_vec.DeltaPhi(mht_vec))	);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mlt_pt", 		mlt_vec.Pt()			);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mt_lmlt", 		lMLT_vec.M()			);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "abs_dPhi_lmlt", 	abs(lep_vec.DeltaPhi(mlt_vec))	);

	// **************** BDT score ********************
	BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_noMass_old",	BDT_noMass_old.Evaluate(b_map_flt, b_map_int) );
	BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_mass_old",  	BDT_mass_old.Evaluate(b_map_flt, b_map_int)   );		
	BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_noMass_new",      BDT_noMass_new.Evaluate(b_map_flt, b_map_int) );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_mass_new",        BDT_mass_new.Evaluate(b_map_flt, b_map_int)   );

	//////////////////////////////////////////////////////////////////
	/// Loop through category cuts defined in src/CategoryCuts.cc  ///
	//////////////////////////////////////////////////////////////////
	
	  ///////////////////////////////////////////
	  ///  Apply the category selection cuts  ///
	  ///////////////////////////////////////////
        for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {  
	  std::string CAT_CUT = CAT_CUTS.at(iCat);
	  bool pass_cat_cut = true;

	  if (CAT_CUT == "something") {
	  }
	  else { pass_cat_cut = InCategory(obj_sel, br, CAT_CUT, verbose); }

	  if (not pass_cat_cut) continue;
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUT << ", in category " << CAT_CUT << std::endl;

	  
	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
	Out_Tree->Fill();  // again, only one opt sel
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
    } // End conditional: if (pass_sel_cuts)

  } // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;
     
  // Write output file
  Out_File->cd();
  Out_Tree->Write();
  Out_File->Write();
      
  
  std::cout << "\nExiting ()\n";
  
} // End void MiniNTupliser()
