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

// #include "H2MuAnalyzer/MakeHistos/interface/SampleDatabase2016.h" // Input data and MC samples

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 1;     // Maximum index of input files to process
const int MAX_EVT  = 10000000;    // Maximum number of events to process
const int PRT_EVT  = 1000;  // Print every N events
const float SAMP_WGT = 1.0;
const float LUMI = 41000; // pb-1   36814 for 2016, 41000 for 2017
const bool verbose = false; // Print extra information

//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/ZH_HToMuMu_ZToAll_M125_13TeV_powheg_pythia8/H2Mu_ZH_125";
//const TString SAMPLE   = "H2Mu_ZH_125";

const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/ZZTo4L_13TeV_powheg_pythia8/ZZ_4l";
const TString SAMPLE   = "ZZ_4l";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
//const TString SAMPLE   = "SingleMu";
const std::string YEAR = "2017";
const std::string SLIM = "Slim"; // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR  = "plots";
const TString HIST_TREE = "Tree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both. Not in use in this macro

const std::vector<std::string> SEL_CUTS = {"Presel2017"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"ZH_4l_ele"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void MiniNTupliser_4l_cat( TString sample = "", TString in_dir = "", TString out_dir = "",
	     std::vector<TString> in_files = {}, TString out_file_str = "",
	     int max_evt = 0, int prt_evt = 0, float samp_weight = 1.0,
	     TString hist_tree = "" ) {
  
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
      //in_file_name.Form("%s/tuple_%d.root", in_dir.Data(), i);
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

  // Initialize set of pointers to all branches in tree
  NTupleBranches br;

  // Initialize empty map of histogram names to histograms
  TFile * Out_File;
  TTree * Out_Tree;
  // declare variables that could be filled in tree. Not all of them are used.
  // branches will be added in cat_sel

  int		nMuons;
  int		nEles;
  int		nJets;
  int		nBJets_Loose;
  int           nBJets_Med;
  int           nBJets_Tight;
  int		nFwdJets;
  int 		nCentJets;

  float		dimu_mass;
  float		dimu_mass_err;
  float		dimu_pt;
  float		dimu_eta;    // in BDT training, will use abs value.
  float 	dimu_dEta;   // Keep this for convenience of producing data/MC from miniNTuple  -- XWZ 19.11.2019
  float		dimu_dPhi;
  float		dimu_dR;
  int 		dimu_gen_ID;
  float 	mu1_pt;
  float		mu1_eta;
  float		mu1_lepMVA;
  int		mu1_charge;
  float		mu2_pt;
  float		mu2_eta;
  float		mu2_lepMVA;
  int		mu2_charge;

  float		lep1_pt;
  float		lep1_eta;
  float		lep1_lepMVA;
  int 		lep1_charge;
  float         lep2_pt;
  float         lep2_eta;
  float         lep2_lepMVA;
  int		lep2_charge;

  float 	cts_mu1;
  float		cts_mu_pos;

  float		dilep_mass;
  float		dilep_mass_err;
  float		dilep_pt;
  float		dilep_eta;
  float 	dilep_dEta;
  float		dilep_dPhi;
  float		dilep_dR;
  int		dilep_gen_ID;

  float		cts_lep1;
  float		cts_lep_pos;

  float		quadlep_mass;
  float		quadlep_pt;
  float		quadlep_eta;
  float		dipair_dEta_H;
  float		dipair_dPhi_H;
  float		dipair_dR_H;
  float         dipair_dEta_pt;
  float         dipair_dPhi_pt;
  float         dipair_dR_pt;

  float		cts_dipair_H;
  float		cts_dipair_pt;

  float 	cs_costheta;
  float 	cs_cosphi;
  float		cs_sinphi;

  float 	cos_theta1;
  float		cos_phi1;
  float		cos_phiH;

  float		dijet_mass;
  float		dijet_pt;
  float		dijet_eta;
  float		dijet_dEta;
  float		dijet_dPhi;
  float		dijet_dR;

  float		jet1_pt;
  float		jet1_eta;
  float		jet2_pt;
  float 	jet2_eta;
  float		jet0_pt;
  float		jet0_eta;

  float		met_pt;
  float		mht_pt;
  float		mht_mass;

  float		event_wgt;
  float		xsec_norm;
  bool 		lep_is_ele;
  bool		lep_is_mu;
  int 		Sample_ID;
  TString 	Sample_name = "";// for identifying sample in minintuple

  MuPairInfo 	dimu;
  MuonInfo 	mu_1;
  MuonInfo	mu_2;
  MuonInfo	mu_3;
  MuonInfo      mu_4;
  EleInfo 	ele1;
  EleInfo       ele2;
  JetPairInfo	dijet;
  JetInfo	jet1;
  JetInfo	jet2;
  JetInfo	jet0;
  MhtInfo	mht_info;

  TLorentzVector dimu_vec;  
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

  // declaration end.
  // creating output file

  TString out_file_name;
  if (out_file_str.Length() > 0) out_file_name.Form( "%s/histos_%s_%s.root",    out_dir.Data(), sample.Data(), out_file_str.Data() );
  else                           out_file_name.Form( "%s/histos_%s_%d_%d.root", out_dir.Data(), sample.Data(), MIN_FILE, MAX_FILE );
  std::cout << "\nCreating output file " << out_file_name.Data() << std::endl;
  Out_File = TFile::Open( out_file_name, "RECREATE" );

  // plant TTree and add branches
  Out_Tree = PlantTree("tree", "flat_tree_minintuple");    

  //muon var
  Out_Tree->Branch("mu1_pt", 		& mu1_pt, 		"mu1_pt/F");
  Out_Tree->Branch("mu1_eta", 		& mu1_eta,		"mu1_eta/F");	
  Out_Tree->Branch("mu1_lepMVA",	& mu1_lepMVA,		"mu1_lepMVA/F");
  Out_Tree->Branch("mu1_charge",	& mu1_charge,		"mu1_charge/I");
  Out_Tree->Branch("mu2_pt",    	& mu2_pt,       	"mu2_pt/F");
  Out_Tree->Branch("mu2_eta",   	& mu2_eta,      	"mu2_eta/F");
  Out_Tree->Branch("mu2_lepMVA",	& mu2_lepMVA,		"mu2_lepMVA/F");
  Out_Tree->Branch("mu2_charge",	& mu2_charge,		"mu2_charge/I");

  Out_Tree->Branch("dimu_mass",		& dimu_mass,		"dimu_mass/F");
  Out_Tree->Branch("dimu_mass_err",	& dimu_mass_err,	"dimu_mass_err/F");
  Out_Tree->Branch("dimu_pt", 		& dimu_pt, 		"dimu_pt/F");
  Out_Tree->Branch("dimu_eta", 		& dimu_eta, 		"dimu_eta/F");
  Out_Tree->Branch("dimu_dEta",		& dimu_dEta,		"dimu_dEta/F");
  Out_Tree->Branch("dimu_dPhi",		& dimu_dPhi,		"dimu_dPhi/F");
  Out_Tree->Branch("dimu_dR",		& dimu_dR,		"dimu_dR/F");
  Out_Tree->Branch("dimu_gen_ID",       & dimu_gen_ID,          "dimu_gen_ID/I");  

  Out_Tree->Branch("cts_mu1",           & cts_mu1,              "cts_mu1/F"  );
  Out_Tree->Branch("cts_mu_pos",        & cts_mu_pos,           "cts_mu_pos/F"  );

  //lep var
  Out_Tree->Branch("lep1_pt", 		& lep1_pt,		"lep1_pt/F");
  Out_Tree->Branch("lep1_eta", 		& lep1_eta, 		"lep1_eta/F");
  Out_Tree->Branch("lep1_lepMVA",	& lep1_lepMVA,		"lep1_lepMVA/F");
  Out_Tree->Branch("lep1_charge",	& lep1_charge,		"lep1_charge/I");
  Out_Tree->Branch("lep2_pt",           & lep2_pt,              "lep2_pt/F");
  Out_Tree->Branch("lep2_eta",          & lep2_eta,             "lep2_eta/F");
  Out_Tree->Branch("lep2_lepMVA",       & lep2_lepMVA,          "lep2_lepMVA/F");
  Out_Tree->Branch("lep2_charge",	& lep2_charge,		"lep2_charge/I");

  Out_Tree->Branch("dilep_mass",        & dilep_mass,           "dilep_mass/F");
  Out_Tree->Branch("dilep_mass_err",    & dilep_mass_err,       "dilep_mass_err/F");
  Out_Tree->Branch("dilep_pt",          & dilep_pt,             "dilep_pt/F");
  Out_Tree->Branch("dilep_eta",         & dilep_eta,            "dilep_eta/F");
  Out_Tree->Branch("dilep_dEta",        & dilep_dEta,           "dilep_dEta/F");
  Out_Tree->Branch("dilep_dPhi",        & dilep_dPhi,           "dilep_dPhi/F");
  Out_Tree->Branch("dilep_dR",          & dilep_dR,             "dilep_dR/F");
  Out_Tree->Branch("dilep_gen_ID",      & dilep_gen_ID,         "dilep_gen_ID/I"); 

  Out_Tree->Branch("cts_lep1",          & cts_lep1,             "cts_lep1/F"  );
  Out_Tree->Branch("cts_lep_pos",       & cts_lep_pos,          "cts_lep_pos/F"  );

  //quadlep var
  Out_Tree->Branch("quadlep_mass",	& quadlep_mass,		"quadlep_mass/F");
  Out_Tree->Branch("quadlep_pt",	& quadlep_pt,   	"quadlep_pt/F");
  Out_Tree->Branch("quadlep_eta",	& quadlep_eta,  	"quadlep_eta/F");
  Out_Tree->Branch("dipair_dEta_H",	& dipair_dEta_H,  	"dipair_dEta_H/F");
  Out_Tree->Branch("dipair_dPhi_H",	& dipair_dPhi_H,  	"dipair_dPhi_H/F");
  Out_Tree->Branch("dipair_dR_H",	& dipair_dR_H,    	"dipair_dR_H/F");
  Out_Tree->Branch("dipair_dEta_pt",    & dipair_dEta_pt,       "dipair_dEta_pt/F");
  Out_Tree->Branch("dipair_dPhi_pt",    & dipair_dPhi_pt,       "dipair_dPhi_pt/F");
  Out_Tree->Branch("dipair_dR_pt",      & dipair_dR_pt,         "dipair_dR_pt/F");

  Out_Tree->Branch("cts_dipair_H",	& cts_dipair_H,   	"cts_dipair_H/F");
  Out_Tree->Branch("cts_dipair_pt",     & cts_dipair_pt,        "cts_dipair_pt/F");

  Out_Tree->Branch("cs_costheta",	& cs_costheta,		"cs_costheta/F");
  Out_Tree->Branch("cs_cosphi",		& cs_cosphi,		"cs_cosphi/F");
  Out_Tree->Branch("cs_sinphi",         & cs_sinphi,            "cs_sinphi/F");

  Out_Tree->Branch("cos_theta1",	& cos_theta1,		"cos_theta1/F");
  Out_Tree->Branch("cos_phi1",		& cos_phi1,		"cos_phi1/F");
  Out_Tree->Branch("cos_phiH",		& cos_phiH,		"cos_phiH/F");

  //jet var
  Out_Tree->Branch("dijet_mass",	& dijet_mass,		"dijet_mass/F");
  Out_Tree->Branch("dijet_pt",		& dijet_pt,		"dijet_pt/F");
  Out_Tree->Branch("dijet_eta",		& dijet_eta,		"dijet_eta/F");
  Out_Tree->Branch("dijet_dEta",	& dijet_dEta,		"dijet_dEta/F");
  Out_Tree->Branch("dijet_dPhi",	& dijet_dPhi,		"dijet_dPhi/F");
  Out_Tree->Branch("dijet_dR",		& dijet_dR,		"dijet_dR/F");

  Out_Tree->Branch("jet1_pt",		& jet1_pt,		"jet1_pt/F");
  Out_Tree->Branch("jet1_eta",		& jet1_eta,		"jet1_eta/F");
  Out_Tree->Branch("jet2_pt",		& jet2_pt,		"jet2_pt/F");
  Out_Tree->Branch("jet2_eta",		& jet2_eta,		"jet2_eta/F");
  Out_Tree->Branch("jet0_pt",           & jet0_pt,              "jet0_pt/F");
  Out_Tree->Branch("jet0_eta",      	& jet0_eta,         	"jet0_eta/F");

  //met var
  Out_Tree->Branch("met_pt",		& met_pt,       	"met_pt/F");
  Out_Tree->Branch("mht_pt",            & mht_pt,               "mht_pt/F");
  Out_Tree->Branch("mht_mass",          & mht_mass,             "mht_mass/F");

  //event var
  Out_Tree->Branch("nJets",		& nJets,		"nJets/I");
  Out_Tree->Branch("nBJets_Loose",	& nBJets_Loose,		"nBJets_Loose/I");
  Out_Tree->Branch("nBJets_Med",        & nBJets_Med,           "nBJets_Med/I");
  Out_Tree->Branch("nBJets_Tight",      & nBJets_Tight,         "nBJets_Tight/I");
  Out_Tree->Branch("nFwdJets",		& nFwdJets,		"nFwdJets/I");
  Out_Tree->Branch("nCentJets",		& nCentJets,		"nCentJets/I");
  Out_Tree->Branch("nMuons",		& nMuons,		"nMuons/I");
  Out_Tree->Branch("nEles",		& nEles,		"nEles/I");  
  Out_Tree->Branch("lep_is_ele",	& lep_is_ele,		"lep_is_ele/O");
  Out_Tree->Branch("lep_is_mu",         & lep_is_mu,            "lep_is_mu/O");

  //weights and spectator 
  Out_Tree->Branch("event_wgt", 	& event_wgt, 		"event_wgt/F");
  Out_Tree->Branch("xsec_norm",		& xsec_norm,		"xsec_norm/F");
  Out_Tree->Branch("Sample_ID",		& Sample_ID,		"Sample_ID/I");
  Out_Tree->Branch("Sample_name",	& Sample_name,		"Sample_name/C");
//      Out_Tree->Branch();

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

  // Configuration for object selection, event selection, and object weighting
  ObjectSelectionConfig obj_sel;
  EventSelectionConfig  evt_sel;
  EventWeightConfig     evt_wgt;
  ConfigureObjectSelection(obj_sel, YEAR);
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

//  evt_sel.muPair_mass_min = 105; // Require at least one Higgs candidate pair, default 60
// use default for ZZ validation
  obj_sel.mu_pt_min       =  10; // Lower muon pT threshold for muons not from Higgs, default 20
  obj_sel.mu_mIso_max      = 0.4;

  //obj_sel.ele_pt_min = 20;
  obj_sel.ele_ID_cut = "loose";
  obj_sel.ele_mIso_max = 0.4;

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

    // For 2016 NTuples, convert "SlimJets" collection into regular jets
    JetInfos jets_tmp;
    if (SLIM == "Slim") {
      jets_tmp = ConvertSlimJets(*(br.slimJets));
      br.jets  = &jets_tmp;
    }


    /////////////////////////////////
    ///  initialize all branches  ///
    /////////////////////////////////
    event_wgt		= 0;
    xsec_norm		= 0;
    Sample_ID		= -999;
    Sample_name 	= "none";
    lep_is_ele		= false;
    lep_is_mu		= false;

    nMuons		= -999;
    nEles		= -999;
    nJets		= -999;
    nBJets_Loose	= -999;
    nBJets_Med  	= -999;
    nBJets_Tight	= -999;
    nFwdJets		= -999;
    nCentJets		= -999;

    dimu_mass		= -999;
    dimu_mass_err	= -999;
    dimu_pt		= -999;
    dimu_eta		= -999;
    dimu_dEta		= -999;
    dimu_dPhi		= -999;
    dimu_dR		= -999;
    dimu_gen_ID		= 0;
    cts_mu1             = -999;
    cts_mu_pos          = -999;

    mu1_pt		= -999;
    mu1_eta		= -999;
    mu1_lepMVA		= -999;
    mu1_charge		= -999;
    mu2_pt		= -999;
    mu2_eta		= -999;
    mu2_lepMVA		= -999;
    mu2_charge		= -999;

    lep1_pt             = -999;
    lep1_eta         	= -999;
    lep1_lepMVA		= -999;
    lep1_charge		= -999;
    lep2_pt             = -999;
    lep2_eta            = -999;
    lep2_lepMVA         = -999;
    lep2_charge		= -999;

    dilep_mass		= -999;
    dilep_mass_err	= -999; 
    dilep_pt		= -999; 
    dilep_eta		= -999; 
    dilep_dEta		= -999; 
    dilep_dPhi		= -999; 
    dilep_dR		= -999; 
    dilep_gen_ID	= 0;
                 
    cts_lep1	 	= -999; 
    cts_lep_pos		= -999; 
                 
    quadlep_mass	= -999; 
    quadlep_pt		= -999; 
    quadlep_eta		= -999; 
    dipair_dEta_H	= -999; 
    dipair_dPhi_H	= -999; 
    dipair_dR_H		= -999; 
    dipair_dEta_pt      = -999;
    dipair_dPhi_pt      = -999;
    dipair_dR_pt        = -999;
                 
    cts_dipair_H	= -999;
    cts_dipair_pt	= -999; 
 
    cs_costheta		= -999;
    cs_cosphi		= -999;
    cs_sinphi		= -999;

    cos_theta1		= -999;
    cos_phi1		= -999;
    cos_phiH		= -999;

    dijet_mass		= -999;
    dijet_pt		= -999;
    dijet_eta		= -999;
    dijet_dEta		= -999;
    dijet_dPhi		= -999;
    dijet_dR		= -999;

    jet1_pt		= -999;
    jet1_eta		= -999;
    jet2_pt		= -999;
    jet2_eta		= -999;
    jet0_pt		= -999;
    jet0_eta		= -999;

    met_pt		= -999;
    mht_pt		= -999;
    mht_mass		= -999;

    dimu	.init();
    mu_1	.init();
    mu_2	.init();
    mu_3	.init();
    mu_4	.init();
    ele1	.init();
    ele2	.init();
    dijet	.init();
    jet1	.init();
    jet2	.init();
    jet0	.init();
    mht_info	.init();

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
      event_wgt = ( sample.Contains("SingleMu") ? 1.0 : EventWeight(br, evt_wgt, verbose) );
      Sample_name = sample;
      xsec_norm = samp_weight;
      Sample_ID = getSampleID(sample);

      /////////////////////////////////////////////////////////////////////////////////////////
      ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
      /////////////////////////////////////////////////////////////////////////////////////////
      ///  in this macro opt selection act as additional selection, can only have one   - XWZ 12.11.2018

      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	std::string OPT_CUT = OPT_CUTS.at(iOpt);
	
	///////////////////////////////////////////////////////////////
	///  select objects and get variables to be filled in tree  ///
	///////////////////////////////////////////////////////////////

	MuonInfos muons = SelectedMuons(obj_sel, br);
        EleInfos  eles  = SelectedEles(obj_sel, br);
        JetInfos  jets  = SelectedJets(obj_sel, br);

	/////////////////////////
        ///    ZH_4l_ele cat    ///
        /////////////////////////
	if (OPT_CUT == "ZH_4l_ele") {

	  //obj_sel.muPair_Higgs == "sort_OS_dimuon_mass";
	  EleInfos  eles  = SelectedEles(obj_sel, br);
          if (muons.size() != 2 or SelectedMuPairs(obj_sel, br).size() != 1 or SelectedEles(obj_sel, br).size() != 2) continue;
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
	  if (ele1.charge + ele2.charge != 0) continue;
	  lep1_vec = FourVec(ele1);
	  lep1_lepMVA = ele1.lepMVA;
	  lep1_charge = ele1.charge;
	  lep2_vec = FourVec(ele2);
	  lep2_lepMVA = ele2.lepMVA;
	  lep2_charge = ele2.charge;
	  dilep_vec = lep1_vec + lep2_vec;
	  dimu = SelectedCandPair(obj_sel, br);
	  
	  if (not sample.Contains("SingleMu")) {
	    if ( IsGenMatched(dimu, *br.muons, *br.genMuons, "H") ) dimu_gen_ID = 25;
            else if ( IsGenMatched(dimu, *br.muons, *br.genMuons, "Z") ) dimu_gen_ID = 23;
	  }

	  if ( SelectedJets(obj_sel, br, "BTagMedium").size() == 0 ) lep_is_ele = true;
	  else continue;
	} // end if (OPT_CUT == "ZH_4l_ele")


        /////////////////////////
        ///   ZH_4l_mu cat    ///
        /////////////////////////
        if (OPT_CUT == "ZH_4l_mu") {
	  if ( muons.size() != 4 or eles.size() != 0 ) continue;
	  int mu1_id = -999, mu2_id = -999, mu3_id = -999, mu4_id = -999;
	  for (int imu = 0; imu < int((br.muons)->size()) ; imu++) {
	    if ( MuonPass(obj_sel, br.muons->at(imu)) ) {
		mu1_id = imu;
		break;
	    }
	  } 
	  for (int imu = mu1_id+1 ; imu < int((br.muons)->size()) ; imu++) {
	    if ( MuonPass(obj_sel, br.muons->at(imu)) ) {
                mu2_id = imu;
                break;
	    }
          }
	  for (int imu = mu2_id+1 ; imu < int((br.muons)->size()) ; imu++) {
	    if ( MuonPass(obj_sel, br.muons->at(imu)) ) {
                mu3_id = imu;
                break;
            }
          }
	  for (int imu = mu3_id+1 ; imu < int((br.muons)->size()) ; imu++) {
	    if ( MuonPass(obj_sel, br.muons->at(imu)) ) {
                mu4_id = imu;
                break;
            }
          } // finish getting muons
	  if ( SelectedMuPairs(obj_sel, br).size() != 4 ) continue; //so we have two + and two - good muons
	  //comment out if more than 4 muons
	  MuPairInfo Z_cand, H_cand;
	  Z_cand.init();
	  H_cand.init();
	  for (const auto & muPair1 : SelectedMuPairs(obj_sel, br)) {
	    MuPairInfo muPair2;
	    for (const auto & muPair : SelectedMuPairs(obj_sel, br)) {
	      if (muPair.iMu1 != muPair1.iMu1 and muPair.iMu1 != muPair1.iMu2 and muPair.iMu2 != muPair1.iMu1 and muPair.iMu2 != muPair1.iMu2) {
		muPair2 = muPair;
		break;
	      }
	    } // end for (const auto & muPair : SelectedMuPairs(obj_sel, br))
	    //if ( muPair1.mass > 81 and muPair1.mass < 101 and muPair2.mass > 105 and muPair2.mass < 160) {
	    //	Z_cand = muPair1;
	    //	H_cand = muPair2;   // with break, select highest dimu pt Z_cand; without break, select lowest dimu pt Z_cand
	    //	break;
	    //}
	    //else continue;
	    if ( muPair2.mass > 81 and muPair2.mass < 101 and muPair1.mass > 70 and muPair1.mass < 110) {  // 70-110 for validation, 105-160 for signal window
                Z_cand = muPair2;
                H_cand = muPair1;  // with break, select highest dimu pt H_cand; without break, select lowest dimu pt H_cand
                break;		   // combination is selecting highest passing pair(Z and H)
            }
	  } // end for (const auto & muPair : SelectedMuPairs(obj_sel, br))
	  if (Z_cand.mass == -999) continue; // no Z pair found
	  //if (Z_cand.mass < 86 or Z_cand.mass > 96) continue; //cut on Z mass
	  mu_3 = br.muons->at(Z_cand.iMu1);
	  mu_4 = br.muons->at(Z_cand.iMu2);

	  if (not sample.Contains("SingleMu")) {
            if ( IsGenMatched(H_cand, *br.muons, *br.genMuons, "H") ) dimu_gen_ID = 25;
            else if ( IsGenMatched(H_cand, *br.muons, *br.genMuons, "Z") ) dimu_gen_ID = 23;

	    if ( IsGenMatched(Z_cand, *br.muons, *br.genMuons, "H") ) dilep_gen_ID = 25;
            else if ( IsGenMatched(Z_cand, *br.muons, *br.genMuons, "Z") ) dilep_gen_ID = 23;
          }
	
	  dimu = H_cand;
	  lep1_vec = FourVec(mu_3, PTC);
	  lep1_lepMVA = mu_3.lepMVA;
	  lep1_charge = mu_3.charge;
	  lep2_vec = FourVec(mu_4, PTC);
          lep2_lepMVA = mu_4.lepMVA;
          lep2_charge = mu_4.charge;
	  dilep_vec = lep1_vec + lep2_vec;

	  if ( SelectedJets(obj_sel, br, "BTagMedium").size() == 0 ) lep_is_mu = true;
          else continue;
	} // end if (OPT_CUT == "ZH_4l_mu")


        dimu_vec = FourVec( dimu, PTC);
        if ( dimu_vec.M() < 70 ||        // 70-110 for ZZ validation, 105-160 for signal window
             dimu_vec.M() > 110 ) continue;    

	mu_1 = br.muons->at(dimu.iMu1);
        mu_2 = br.muons->at(dimu.iMu2);
        mu1_vec = FourVec(br.muons->at(dimu.iMu1), PTC);
        mu2_vec = FourVec(br.muons->at(dimu.iMu2), PTC);

	//muon vars
	mu1_pt = mu1_vec.Pt();
        mu1_eta = mu1_vec.Eta();
	mu1_lepMVA = mu_1.lepMVA;
	mu1_charge = mu_1.charge;
        mu2_pt = mu2_vec.Pt();
        mu2_eta = mu2_vec.Eta();
  	mu2_lepMVA = mu_2.lepMVA;
	mu2_charge = mu_2.charge;

	dimu_mass 	= dimu_vec.M();                	// all variables filled with vec in order to  
	dimu_mass_err 	= dimu.massErr;
        dimu_pt 	= dimu_vec.Pt();		// work with different PTC
        dimu_eta 	= dimu_vec.Eta();
	dimu_dEta	= mu1_vec.Eta() - mu2_vec.Eta();   
	dimu_dPhi	= mu1_vec.DeltaPhi(mu2_vec);  
	dimu_dR		= mu1_vec.DeltaR(mu2_vec);

	cts_mu1 = CosThetaStar(mu1_vec, mu2_vec);
        cts_mu_pos = mu_1.charge == 1 ? CosThetaStar(mu1_vec, mu2_vec) : CosThetaStar(mu2_vec, mu1_vec);

	//lep vars
	lep1_pt = lep1_vec.Pt();
	lep1_eta = lep1_vec.Eta();
	lep2_pt = lep2_vec.Pt();
        lep2_eta = lep2_vec.Eta();

	dilep_mass	= dilep_vec.M();
	dilep_pt	= dilep_vec.Pt();
	dilep_eta       = dilep_vec.Eta();
        dilep_dEta      = lep1_vec.Eta() - lep2_vec.Eta();
        dilep_dPhi      = lep1_vec.DeltaPhi(lep2_vec);
        dilep_dR        = lep1_vec.DeltaR(lep2_vec);
	
	cts_lep1 = CosThetaStar(lep1_vec, lep2_vec);
        cts_lep_pos = lep1_charge == 1 ? CosThetaStar(lep1_vec, lep2_vec) : CosThetaStar(lep2_vec, lep1_vec);

	//quadlep vars
	quadlep_vec = dimu_vec + dilep_vec;
	quadlep_mass    = quadlep_vec.M();
        quadlep_pt      = quadlep_vec.Pt();
        quadlep_eta     = quadlep_vec.Eta();

        dipair_dEta_H   = dimu_vec.Eta() - dilep_vec.Eta();
        dipair_dPhi_H   = dimu_vec.DeltaPhi(dilep_vec);
        dipair_dR_H     = dimu_vec.DeltaR(dilep_vec);
	dipair_dEta_pt  = (dimu_vec.Pt() > dilep_vec.Pt()) ? (dimu_vec.Eta() - dilep_vec.Eta()) : (dilep_vec.Eta() - dimu_vec.Eta());
        dipair_dPhi_pt  = (dimu_vec.Pt() > dilep_vec.Pt()) ? dimu_vec.DeltaPhi(dilep_vec) : dilep_vec.DeltaPhi(dimu_vec);
        dipair_dR_pt    = (dimu_vec.Pt() > dilep_vec.Pt()) ? dimu_vec.DeltaR(dilep_vec) : dilep_vec.DeltaR(dimu_vec);
	
	cts_dipair_H	= CosThetaStar(dimu_vec, dilep_vec);
	cts_dipair_pt	= (dimu_vec.Pt() > dilep_vec.Pt()) ? CosThetaStar(dimu_vec, dilep_vec) : CosThetaStar(dilep_vec, dimu_vec);

	cs_costheta     = Cos_CS_Theta(dimu_vec, dilep_vec);
  	cs_cosphi	= Cos_CS_Phi(dimu_vec, dilep_vec);
	cs_sinphi	= Sin_CS_Phi(dimu_vec, dilep_vec);

  	KinAngles  kinematic_4l(mu1_vec, mu2_vec, lep1_vec, lep2_vec);

	cos_theta1	= kinematic_4l.MELA_Cos_Theta1();
	cos_phiH        = kinematic_4l.MELA_Cos_PhiH();
	cos_phi1        = kinematic_4l.MELA_Cos_Phi1();
//	cos_theta1	= MELA_Cos_Theta1( mu1_vec, mu2_vec, lep1_vec, lep2_vec );
//	cos_phiH	= MELA_Cos_PhiH( mu1_vec, mu2_vec, lep1_vec, lep2_vec );
// 	cos_phi1	= MELA_Cos_Phi1( mu1_vec, mu2_vec, lep1_vec, lep2_vec );

	//event vars and jet vars
	nMuons = muons.size();
        nEles = eles.size();
        nJets = jets.size();
        nFwdJets     = SelectedJets(obj_sel, br, "Forward").size();
        nCentJets    = SelectedJets(obj_sel, br, "Central").size();
        nBJets_Loose = SelectedJets(obj_sel, br, "BTagLoose").size();
        nBJets_Med   = SelectedJets(obj_sel, br, "BTagMedium").size();
        nBJets_Tight = SelectedJets(obj_sel, br, "BTagTight").size();
        met_vec  = FourVec(*br.met);
        mht_info = *br.mht;
        mht_vec  = FourVec(mht_info);
	met_pt	 = met_vec.Pt();
	mht_pt	 = mht_vec.Pt();
	mht_mass = mht_vec.M();

	if (SelectedJetPairs(obj_sel, br).size() > 0) {
	  dijet = SelectedJetPairs(obj_sel, br).at(0);      // no need for jet vector since no PTC
	  jet1  = br.jets->at(dijet.iJet1);
	  jet2  = br.jets->at(dijet.iJet2);

	  dijet_mass = dijet.mass;
	  dijet_pt   = dijet.pt;
	  dijet_eta  = dijet.eta;
	  dijet_dEta = dijet.dEta;
	  dijet_dPhi = dijet.dPhi;
	  dijet_dR   = dijet.dR;

	  jet1_pt  = jet1.pt;
	  jet1_eta = jet1.eta;
	  jet2_pt  = jet2.pt;
	  jet2_eta = jet2.eta;
	}
        if (SelectedJets(obj_sel, br).size() > 0) {
	  jet0 = SelectedJets(obj_sel, br).at(0);
	  jet0_pt  = jet0.pt;
	  jet0_eta = jet0.eta;
	}
	Out_Tree->Fill(); 

	//} // end for (const auto & muPair : SelectedMuPairs(obj_sel, br))
	//} // end if (OPT_CUT == "ZH_4l_mu")


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
	//Out_Tree->Fill();  // again, only one opt sel
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
