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
#include "H2MuAnalyzer/MakeHistos/interface/MiniNTupleHelper.h"
#include "H2MuAnalyzer/MakeHistos/interface/SampleID.h"
#include "H2MuAnalyzer/MakeHistos/interface/KinematicAngles.h"
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
const bool verbose = false; // Print extra information

const bool CR_validation = false; // whether we are plotting the validation control region or not

//const TString IN_DIR    = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/prod-v18.1.6.skim3l/ZH_HToMuMu_ZToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8/H2Mu_ZH_125/190528_111720/0000";
//const TString SAMPLE    = "H2Mu_ZH_125";

const TString IN_DIR   = "/eos/cms/store/user/bortigno/h2mm/ntuples/2016/94X_v3/STR/ZH_HToMuMu_ZToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8/H2Mu_ZH_125/190625_204256/0000";
const TString SAMPLE   = "H2Mu_ZH_125";

//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/ZH_HToMuMu_ZToAll_M125_13TeV_powheg_pythia8/H2Mu_ZH_125";
//const TString SAMPLE   = "H2Mu_ZH_125";

//const TString IN_DIR    = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/prod-v18.1.6.skim3l/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/ZZ_4l/190528_113414/0000";
//const TString SAMPLE    = "ZZ_4l";

//const TString IN_DIR    = "/eos/cms/store/user/bortigno/h2mm/ntuples/2016/94X_v3/STR/ZZTo4L_13TeV-amcatnloFXFX-pythia8/ZZ_4l_amc/190625_205403/0000";
//const TString SAMPLE    = "ZZ_4l_amc";

//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/ZZTo4L_13TeV_powheg_pythia8/ZZ_4l";
//const TString SAMPLE   = "ZZ_4l";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
//const TString SAMPLE   = "SingleMu";
const std::string YEAR = "2016";
const std::string SLIM  = (YEAR == "2017" ? "Slim" : "notSlim");  // "Slim" or "notSlim" - Some 2017 NTuples are "Slim"
const TString OUT_DIR  = "plots";
const TString HIST_TREE = "Tree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both. Not in use in this macro

const std::vector<std::string> SEL_CUTS = {"PreselRun2"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"ZH_4l_ele", "ZH_4l_mu"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void MiniNTupliser_4l_cat( TString sample = "", TString in_dir = "", TString out_dir = "",
	     std::vector<TString> in_files = {}, TString out_file_str = "",
	     int max_evt = 0, int prt_evt = 0, float samp_weight = 1.0,
	     TString hist_tree = "", std::string SYS = "" ) {
  
  // Set variables to hard-coded values if they are not initialized
  if (sample.Length()      == 0) sample      = SAMPLE;
  if (in_dir.Length()      == 0) in_dir      = IN_DIR;
  if (out_dir.Length()     == 0) out_dir     = OUT_DIR;
  if (max_evt              == 0) max_evt     = MAX_EVT;
  if (prt_evt              == 0) prt_evt     = PRT_EVT;
  if (samp_weight          == 0) samp_weight = SAMP_WGT;
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
      if (YEAR == "2017") in_file_name.Form("%s/NTuple_0.root", in_dir.Data());
      else   		  in_file_name.Form("%s/tuple_%d.root", in_dir.Data(), i);
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
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM}); // Options in {} include "JES", "Flags", and "SFs"
    else
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM, "GEN", "Wgts"}, (SYS.find("JES_") != std::string::npos ? SYS : "noSys"), false); // Options in {} include "JES", "Flags", and "SFs"
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
  ConfigureObjectSelection(obj_sel, YEAR, "lepMVA");
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

  evt_wgt.muon_ID  = false;
  evt_wgt.muon_Iso = false;

  if (!CR_validation)  evt_sel.muPair_mass_min = 105; // Require at least one Higgs candidate pair, default 60
  // use default for ZZ validation
  obj_sel.mu_pt_min       =  20; // Lower muon pT threshold for muons not from Higgs, default 20
  obj_sel.mu_mIso_max      = 0.4;

  obj_sel.ele_pt_min = 20;
  obj_sel.ele_ID_cut = "loose";
  obj_sel.ele_mIso_max = 0.4;

  if (verbose) obj_sel.Print();
  if (verbose) evt_sel.Print();
  if (verbose) evt_wgt.Print();

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later


  std::cout << "\n******* About to load 2D LepMVA efficiency scale factor histograms *******" << std::endl;
  std::map<std::string, TH2F *> lepSF;
  lepSF["mu_T"]  = LoadSFsLepMVA(YEAR,  "mu", "T");
  lepSF["mu_M"]  = LoadSFsLepMVA(YEAR,  "mu", "M");
  lepSF["mu_L"]  = LoadSFsLepMVA(YEAR,  "mu", "L");
  lepSF["ele_T"] = LoadSFsLepMVA(YEAR, "ele", "T");
  lepSF["ele_M"] = LoadSFsLepMVA(YEAR, "ele", "M");
  lepSF["ele_L"] = LoadSFsLepMVA(YEAR, "ele", "L");

  std::cout << "\n******* About to load XML files for signal-background BDTs *******" << std::endl;
  MVA::MVA BDT_noMass_v3( "data/XMLs/ZH_4l/Xunwu/2019_05_22/ZH_2017_lep_against_inclu_sel_angle_without_mass_all_sig_all_bkg_ge0j/",
                       "weights/ZH_2017_lep_against_inclu_sel_angle_without_mass_all_sig_all_bkg_ge0j_BDTG_UF_v3.weights.xml",
                       "BDTG_UF_v3" );
  MVA::MVA BDT_mass_more( "data/XMLs/ZH_4l/Xunwu/2019_05_22/ZH_2017_lep_against_inclu_sel_angle_with_mass_all_sig_all_bkg_ge0j/",
                           "weights/ZH_2017_lep_against_inclu_sel_angle_with_mass_all_sig_all_bkg_ge0j_BDTG_UF_v3.weights.xml",
                           "BDTG_UF_v3" );
  MVA::MVA BDT_mass_min( "data/XMLs/ZH_4l/Xunwu/2019_05_22/ZH_2017_lep_against_inclu_min_with_mass_all_sig_all_bkg_ge0j/",
                           "weights/ZH_2017_lep_against_inclu_min_with_mass_all_sig_all_bkg_ge0j_BDTG_UF_v3.weights.xml",
                           "BDTG_UF_v3" );
  MVA::MVA BDT_and_mass( "data/XMLs/ZH_4l/Xunwu/2019_05_22/ZH_2017_lep_against_inclu_BDTv3_and_mass_all_sig_all_bkg_ge0j/",
                           "weights/ZH_2017_lep_against_inclu_BDTv3_and_mass_all_sig_all_bkg_ge0j_BDTG_UF_v2.weights.xml",
                           "BDTG_UF_v2" );



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
      float all_lepMVA_SF = 1.0;
      float xsec_norm = samp_weight;
      int Sample_ID = getSampleID(sample);

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

	MuPairInfo    dimu;         // object initial values are all -999
  	MuonInfo      mu_1;
  	MuonInfo      mu_2;
  	MuonInfo      mu_3;
  	MuonInfo      mu_4;
  	EleInfo       ele1;
  	EleInfo       ele2;
  	JetPairInfo   dijet;
  	JetInfo       jet1;
  	JetInfo       jet2;
  	JetInfo       jet0;
  	MhtInfo       mht_info;

  	TLorentzVector dimu_vec;    // vector default initial is (0,0,0,0)
  	TLorentzVector mu1_vec;
  	TLorentzVector mu2_vec;
  	TLorentzVector mu3_vec;
  	TLorentzVector mu4_vec;
  	TLorentzVector lep1_vec;
  	TLorentzVector lep2_vec;
  	TLorentzVector dilep_vec;
  	TLorentzVector quadlep_vec;
  	TLorentzVector dijet_vec;
  	TLorentzVector jet1_vec;
  	TLorentzVector jet2_vec;
  	TLorentzVector met_vec;
  	TLorentzVector mht_vec;


	// ************ category cuts and fill lepton variables ***************
	MuonInfos muons = SelectedMuons(obj_sel, br);
        EleInfos  eles  = SelectedEles(obj_sel, br);
        JetInfos  jets  = SelectedJets(obj_sel, br);

	int lep_ID = 0;  // default value is 0
	int lep1_charge = 0;
	int lep2_charge = 0;
	int dilep_gen_ID = 0;
	/////////////////////////
        ///    ZH_4l_ele cat    ///
        /////////////////////////
	if (OPT_CUT == "ZH_4l_ele") {
          if (muons.size() != 2 or SelectedMuPairs(obj_sel, br).size() != 1 or SelectedEles(obj_sel, br).size() != 2) continue;
	  if ( SelectedJets(obj_sel, br, "BTagMedium").size() != 0 or SelectedJets(obj_sel, br, "BTagLoose").size() > 1 ) continue;

	  for (const auto & electron : SelectedEles(obj_sel, br)) {  // this is for if more than 2 electrons, select highest pt ones
	     if (electron.lepMVA > -1.0) { //-1.0 as a place holder
		ele1 = electron;
		break;
	     }
	  }
	  //ele1 = eles.at(0);
	  for (const auto & electron : SelectedEles(obj_sel, br)) {
	     if (electron.charge + ele1.charge == 0 and electron.lepMVA > -1.0) { // this is for if more than 2 electrons, select highest pt ones
		ele2 = electron;
		break;
	     }
	  }
	  //ele2 = eles.at(1);
	  if (ele1.charge + ele2.charge != 0) continue; //make sure ele1 and ele2 are filled.

	  lep_ID = 11;
	  lep1_charge = ele1.charge;
	  lep2_charge = ele2.charge;
	  lep1_vec = FourVec(ele1);
	  lep2_vec = FourVec(ele2);
	  dilep_vec = lep1_vec + lep2_vec;
	  dilep_gen_ID = 0;
	  dimu = SelectedCandPair(obj_sel, br);

	  // veto events with electron pair from low mass resonances
	  if (dilep_vec.M() < 12) continue;

	  if ( not sample.Contains("SingleMu") ) {
	    all_lepMVA_SF *= LepMVASF(lepSF["ele_M"], lep1_vec.Pt(), lep1_vec.Eta(), (SYS.find("LepMVA_SF_") != std::string::npos ? SYS : "noSys") );
	    all_lepMVA_SF *= LepMVASF(lepSF["ele_M"], lep2_vec.Pt(), lep2_vec.Eta(), (SYS.find("LepMVA_SF_") != std::string::npos ? SYS : "noSys") );	  
	  }

	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep1_lepMVA", ele1.lepMVA);
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep1_charge", ele1.charge);
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep2_lepMVA", ele2.lepMVA);
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep2_charge", ele2.charge);
	  BookAndFill( b_map_int, Out_Tree, h_pre, "lep_ID", 	  lep_ID);
	} // end if (OPT_CUT == "ZH_4l_ele")


        /////////////////////////
        ///   ZH_4l_mu cat    ///
        /////////////////////////
        if (OPT_CUT == "ZH_4l_mu") {
	  if ( muons.size() != 4 or eles.size() != 0 ) continue;
	  if ( SelectedMuPairs(obj_sel, br).size() != 4 ) continue; //two + and two -, comment out if more than 4 muons
	  if ( SelectedJets(obj_sel, br, "BTagMedium").size() != 0 or SelectedJets(obj_sel, br, "BTagLoose").size() > 1 ) continue;

	  MuPairInfo Z_cand, H_cand;
	  Z_cand.init();
	  H_cand.init();
	  float Z_off = 10.0;  // allow Z candidate to be off shell by maximum 10 GeV

	  for (const auto & muPair1 : SelectedMuPairs(obj_sel, br)) { // loops through 4 pairs (2 combos X which is H)
	    MuPairInfo muPair2;  
	    for (const auto & muPair : SelectedMuPairs(obj_sel, br)) {
	      if (muPair.iMu1 != muPair1.iMu1 and muPair.iMu1 != muPair1.iMu2 and muPair.iMu2 != muPair1.iMu1 and muPair.iMu2 != muPair1.iMu2) {
		muPair2 = muPair;
		break;
	      }
	    } // found the other pair

	    // veto events with dimuon from low mass resonances
	    if ( muPair1.mass < 12 or muPair2.mass <12 ) {
	      Z_cand.init();
	      break;
	    }
	    // for analysis, veto events with two Z-mass pairs; for validation, use them
	    if ( abs(muPair1.mass-91)<10 and abs(muPair2.mass-91)<10 ) {
	      if (CR_validation) { // if making the validation plots
		H_cand = muPair1;  // pair 1 is higher pt
                Z_cand = muPair2;  // pair 2 is lower pt
	      }
	      else Z_cand.init();
	      break;
	    }
	    // for analysis, take the combination with the pair closest to Z mass; for validation, skip this step
	    if ( !CR_validation and muPair1.mass > 110 and muPair1.mass < 160 and abs(muPair2.mass - 91)<Z_off) { // this part is only for signal region
		H_cand = muPair1;   // muPair1 as H_cand, 70-110 for validation, 105-160 for signal window
                Z_cand = muPair2;   // muPair2 as Z_cand, closest to 91 GeV 
		Z_off = abs(Z_cand.mass - 91);
            } // if muPair1 fits for H, and muPair2 fits for Z
	  } // end for (const auto & muPair : SelectedMuPairs(obj_sel, br))

	  if (Z_cand.mass == -999) continue; // no Z pair found

	  mu_3 = br.muons->at(Z_cand.iMu1);
	  mu_4 = br.muons->at(Z_cand.iMu2);
	  lep_ID = 13;
	  lep1_charge = mu_3.charge;
	  lep2_charge = mu_4.charge;
	  lep1_vec = FourVec(mu_3, PTC);
	  lep2_vec = FourVec(mu_4, PTC);
	  dilep_vec = lep1_vec + lep2_vec;
	  dimu = H_cand;

	  if (not sample.Contains("SingleMu")) {
            if ( IsGenMatched(Z_cand, *br.muons, *br.genMuons, *br.genParents, "H") ) dilep_gen_ID = 25;
            else if ( IsGenMatched(Z_cand, *br.muons, *br.genMuons, *br.genParents, "Z") ) dilep_gen_ID = 23;
	    else if ( IsGenMatched(Z_cand, *br.muons, *br.genMuons, *br.genParents, "gamma") ) dilep_gen_ID = 22;
	    else if ( IsGenMatched(Z_cand, *br.muons, *br.genMuons, *br.genParents, "tau") ) dilep_gen_ID = 15;
	    else if ( IsGenMatched(Z_cand, *br.muons, *br.genMuons, *br.genParents, "light_quark") ) dilep_gen_ID = 1;
          }

	  if ( not sample.Contains("SingleMu") ) {
            all_lepMVA_SF *= LepMVASF(lepSF["mu_M"], lep1_vec.Pt(), abs(lep1_vec.Eta()), (SYS.find("LepMVA_SF_") != std::string::npos ? SYS : "noSys") );
            all_lepMVA_SF *= LepMVASF(lepSF["mu_M"], lep2_vec.Pt(), abs(lep2_vec.Eta()), (SYS.find("LepMVA_SF_") != std::string::npos ? SYS : "noSys") );
          }

	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep1_lepMVA", mu_3.lepMVA);
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep1_charge", mu_3.charge);
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep2_lepMVA", mu_4.lepMVA);
	  BookAndFill( b_map_flt, Out_Tree, h_pre, "lep2_charge", mu_4.charge);
	  BookAndFill( b_map_int, Out_Tree, h_pre, "lep_ID", 	  lep_ID);
	} // end if (OPT_CUT == "ZH_4l_mu")

	if ( lep1_charge == 0 or lep1_charge + lep2_charge != 0) {
	    std::cout << "weird case: no proper lepton found" << std::endl;
	    continue;
	}

	// ********* categories done, lep variables filled *********
	// now must already be in one of the categories

	// **************** remaining lepton variables ******************
	BookAndFill( b_map_flt, Out_Tree, h_pre, "lep1_pt", 		lep1_vec.Pt()		);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "lep1_abs_eta", 	abs(lep1_vec.Eta())	);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "lep2_pt", 		lep2_vec.Pt() 		);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "lep2_abs_eta", 	abs(lep2_vec.Eta())	);
 
	// ***************** muon variables ********************
        dimu_vec = FourVec( dimu, PTC);
	if (CR_validation) {if ( dimu_vec.M() < 81  || dimu_vec.M() > 101 ) continue;}  // 70-110 for ZZ validation,
	else 		   {if ( dimu_vec.M() < 110 || dimu_vec.M() > 160 ) continue;}  // 110-160 for signal window  

	int dimu_gen_ID = 0;
        if (not sample.Contains("SingleMu")) {
          if ( IsGenMatched(dimu, *br.muons, *br.genMuons, *br.genParents, "H") ) dimu_gen_ID = 25;
          else if ( IsGenMatched(dimu, *br.muons, *br.genMuons, *br.genParents, "Z") ) dimu_gen_ID = 23;
	  else if ( IsGenMatched(dimu, *br.muons, *br.genMuons, *br.genParents, "gamma") ) dimu_gen_ID = 22;
	  else if ( IsGenMatched(dimu, *br.muons, *br.genMuons, *br.genParents, "tau") ) dimu_gen_ID = 15;
	  else if ( IsGenMatched(dimu, *br.muons, *br.genMuons, *br.genParents, "light_quark") ) dimu_gen_ID = 1;
        }

	mu_1 = br.muons->at(dimu.iMu1);
        mu_2 = br.muons->at(dimu.iMu2);
        mu1_vec = FourVec(br.muons->at(dimu.iMu1), PTC);
        mu2_vec = FourVec(br.muons->at(dimu.iMu2), PTC);

	if ( not sample.Contains("SingleMu") ) {
          all_lepMVA_SF *= LepMVASF(lepSF["mu_M"], mu1_vec.Pt(), abs(mu1_vec.Eta()), (SYS.find("LepMVA_SF_") != std::string::npos ? SYS : "noSys") );
          all_lepMVA_SF *= LepMVASF(lepSF["mu_M"], mu2_vec.Pt(), abs(mu2_vec.Eta()), (SYS.find("LepMVA_SF_") != std::string::npos ? SYS : "noSys") );
        }

	BookAndFill( b_map_int, Out_Tree, h_pre, "dimu_gen_ID", 	dimu_gen_ID				);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_mass", 		dimu_vec.M()				);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_mass_err", 	dimu.massErr				);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_pt", 		dimu_vec.Pt()				);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_abs_eta",        abs(dimu_vec.Eta())                     );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_abs_dEta",       abs(mu1_vec.Eta() - mu2_vec.Eta())      );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_abs_dPhi",       abs(mu1_vec.DeltaPhi(mu2_vec))          );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dimu_dR",             mu1_vec.DeltaR(mu2_vec)                 );

	BookAndFill( b_map_flt, Out_Tree, h_pre, "mu1_pt",              mu1_vec.Pt()            );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mu1_abs_eta",         abs(mu1_vec.Eta())      );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mu1_lepMVA",          mu_1.lepMVA             );
        BookAndFill( b_map_int, Out_Tree, h_pre, "mu1_charge",          mu_1.charge             );

        BookAndFill( b_map_flt, Out_Tree, h_pre, "mu2_pt",              mu2_vec.Pt()            );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mu2_abs_eta",         abs(mu2_vec.Eta())      );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mu2_lepMVA",          mu_2.lepMVA             );
        BookAndFill( b_map_int, Out_Tree, h_pre, "mu2_charge",          mu_2.charge             );

	BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_mu1",             CosThetaStar(mu1_vec, mu2_vec)  );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_mu_pos",          mu_1.charge == 1 ? CosThetaStar(mu1_vec, mu2_vec) : CosThetaStar(mu2_vec, mu1_vec)  );

	// ******************* dilep variables ************************
	BookAndFill( b_map_int, Out_Tree, h_pre, "dilep_gen_ID",	dilep_gen_ID				);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dilep_mass", 		dilep_vec.M()				);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dilep_pt", 		dilep_vec.Pt()				);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dilep_abs_eta", 	abs(dilep_vec.Eta())			);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dilep_abs_dEta", 	abs(lep1_vec.Eta() - lep2_vec.Eta())	);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dilep_abs_dPhi", 	abs(lep1_vec.DeltaPhi(lep2_vec))	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dilep_dR", 		lep1_vec.DeltaR(lep2_vec)		);

        BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_lep1", 		CosThetaStar(lep1_vec, lep2_vec)	);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_lep_pos", 	lep1_charge == 1 ? CosThetaStar(lep1_vec, lep2_vec) : CosThetaStar(lep2_vec, lep1_vec));

	// ******************** quadlep variables *********************
	quadlep_vec = dimu_vec + dilep_vec;
	KinAngles  kinematic_4l(mu1_vec, mu2_vec, lep1_vec, lep2_vec);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "quadlep_mass", 	quadlep_vec.M()		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "quadlep_pt", 		quadlep_vec.Pt()	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "quadlep_abs_eta", 	abs(quadlep_vec.Eta()) 	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dipair_dEta_H", 	dimu_vec.Eta() - dilep_vec.Eta()	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dipair_dPhi_H", 	dimu_vec.DeltaPhi(dilep_vec)		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "dipair_dR_H", 	dimu_vec.DeltaR(dilep_vec)		);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "cts_dipair_H", 	CosThetaStar(dimu_vec, dilep_vec)	);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "cs_costheta", 	Cos_CS_Theta(dimu_vec, dilep_vec)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "cs_cosphi", 		Cos_CS_Phi(dimu_vec, dilep_vec)		);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "cs_sinphi", 		Sin_CS_Phi(dimu_vec, dilep_vec)		);

	BookAndFill( b_map_flt, Out_Tree, h_pre, "cos_theta1", 		kinematic_4l.MELA_Cos_Theta1()		);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "cos_phiH", 		kinematic_4l.MELA_Cos_PhiH()		);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "cos_phi1", 		kinematic_4l.MELA_Cos_Phi1()		);

	// ********************* event variables *********************
	BookAndFill( b_map_int, Out_Tree, h_pre, "nMuons",              muons.size()            );
        BookAndFill( b_map_int, Out_Tree, h_pre, "nEles",               eles.size()             );
        BookAndFill( b_map_int, Out_Tree, h_pre, "nJets",               jets.size()             );
        BookAndFill( b_map_int, Out_Tree, h_pre, "nFwdJets",            SelectedJets(obj_sel, br, "Forward").size()     );
        BookAndFill( b_map_int, Out_Tree, h_pre, "nCentJets",           SelectedJets(obj_sel, br, "Central").size()     );
        BookAndFill( b_map_int, Out_Tree, h_pre, "nBJets_Loose",        SelectedJets(obj_sel, br, "BTagLoose").size()   );
        BookAndFill( b_map_int, Out_Tree, h_pre, "nBJets_Med",          SelectedJets(obj_sel, br, "BTagMedium").size()  );
        BookAndFill( b_map_int, Out_Tree, h_pre, "nBJets_Tight",        SelectedJets(obj_sel, br, "BTagTight").size()   );

	// ******************** met variables ************************
        met_vec  = FourVec(*br.met);
        mht_info = *br.mht;
        mht_vec  = FourVec(mht_info);
	float R_met_4l = met_vec.Pt() / quadlep_vec.Pt();
	float R_mht_4l = mht_vec.Pt() / quadlep_vec.Pt();

	BookAndFill( b_map_flt, Out_Tree, h_pre, "met_pt",              met_vec.Pt()                    );
	BookAndFill( b_map_flt, Out_Tree, h_pre, "abs_dPhi_4lPt_met", 	abs(quadlep_vec.DeltaPhi(met_vec))	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "met_Rtrans_4lPt", 	max(min(R_met_4l * cos(quadlep_vec.DeltaPhi(met_vec)),1.5),-1.5) 	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "met_Rlongi_4lPt",    	max(min(R_met_4l * sin(quadlep_vec.DeltaPhi(met_vec)),1.5),-1.5)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mht_pt",              mht_vec.Pt()                    );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mht_mass",            mht_info.mass                   );
	BookAndFill( b_map_flt, Out_Tree, h_pre, "abs_dPhi_4lPt_mht",   abs(quadlep_vec.DeltaPhi(mht_vec))	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "mht_Rtrans_4lPt",    	max(min(R_mht_4l * cos(quadlep_vec.DeltaPhi(mht_vec)),2.5),-1.5)	);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "mht_Rlongi_4lPt",    	max(min(R_mht_4l * sin(quadlep_vec.DeltaPhi(mht_vec)),1.5),-1.5)	);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "abs_dPhi_met_mht",	abs(met_vec.DeltaPhi(mht_vec))		);
	// ratio of transverse met to quadlep_pt?

	// ******************* jet variables ************************
	if (SelectedJetPairs(obj_sel, br).size() > 0) {
          dijet = SelectedJetPairs(obj_sel, br).at(0);      // no need for jet vector since no PTC
          jet1  = br.jets->at(dijet.iJet1);
          jet2  = br.jets->at(dijet.iJet2);
        }
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_mass",          dijet.mass      );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_pt",            dijet.pt        );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_abs_eta",       abs(dijet.eta)  );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_abs_dEta",      abs(dijet.dEta) );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_abs_dPhi",      abs(dijet.dPhi) );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "dijet_dR",            dijet.dR);

        BookAndFill( b_map_flt, Out_Tree, h_pre, "jet1_pt",             jet1.pt         );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "jet1_abs_eta",        abs(jet1.eta)   );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "jet2_pt",             jet2.pt         );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "jet2_abs_eta",        abs(jet2.eta)   );

        if (SelectedJets(obj_sel, br).size() > 0) {
          jet0 = SelectedJets(obj_sel, br).at(0);
        }
        BookAndFill( b_map_flt, Out_Tree, h_pre, "jet0_pt",             jet0.pt         );
        BookAndFill( b_map_flt, Out_Tree, h_pre, "jet0_abs_eta",        abs(jet0.eta)   );

	// ****************** BDT score **************************
	BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_noMass_v3",      BDT_noMass_v3.Evaluate(b_map_flt, b_map_int) );
	BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_mass_more",      BDT_mass_more.Evaluate(b_map_flt, b_map_int) );
	BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_mass_min",       BDT_mass_min.Evaluate(b_map_flt, b_map_int) );
	BookAndFill( b_map_flt, Out_Tree, h_pre, "BDT_and_mass",       BDT_and_mass.Evaluate(b_map_flt, b_map_int) );

	// ***************** auxiliary variables *****************
	BookAndFill( b_map_flt, Out_Tree, h_pre, "event_wgt",   	event_wgt);
        BookAndFill( b_map_flt, Out_Tree, h_pre, "xsec_norm",   	xsec_norm);
	BookAndFill( b_map_flt, Out_Tree, h_pre, "all_lepMVA_SF", 	all_lepMVA_SF);
        BookAndFill( b_map_int, Out_Tree, h_pre, "Sample_ID",   	Sample_ID);
        BookAndFill( b_map_str, Out_Tree, h_pre, "Sample_name",    	sample);

	//} // end for (const auto & muPair : SelectedMuPairs(obj_sel, br))


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
  
} // End void MiniNTupliser_4l_cat()
