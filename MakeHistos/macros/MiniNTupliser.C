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

const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos/170315_105045/0000";
const TString SAMPLE   = "H2Mu_WH_pos";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WZ_3l_AMC";
//const TString SAMPLE   = "WZ_3l_AMC";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
//const TString SAMPLE   = "SingleMu";
const std::string YEAR = "2017";
const std::string SLIM = "Slim"; // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {"Presel2017"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void MiniNTupliser( TString sample = "", TString in_dir = "", TString out_dir = "",
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
  float		dimu_pt;
  float		dimu_eta;    // in BDT training, will use abs value.
  float 	dimu_dEta;   // Keep this for convenience of producing data/MC from miniNTuple  -- XWZ 19.11.2019
  float		dimu_dPhi;
  float		dimu_dR;
  float 	mu1_pt;
  float		mu1_eta;
  float		mu2_pt;
  float		mu2_eta;

  float		ele_pt;
  float		ele_eta;
  float		extra_mu_pt;
  float		extra_mu_eta;

  float 	cts_mu1;
  float		cts_mu_pos;
  float		cts_edimu;
  float		cts_emuSS;
  float		cts_emuOS;
  
  float         edimu_mass;
  float         edimu_pt;
  float  	edimu_eta;
  float 	edimu_dEta;
  float		edimu_dPhi;
  float		edimu_dR;

  float 	emuSS_pt;	
  float 	emuSS_eta;
  float 	emuSS_dEta;
  float 	emuSS_dPhi;	
  float 	emuSS_dR;		
  float 	emuOS_pt;		
  float 	emuOS_eta;	
  float 	emuOS_dEta;
  float 	emuOS_dPhi;	
  float 	emuOS_dR;		

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
  float		mt_emet;
  float		dPhi_emet;
  float		mht_pt;
  float		mht_mass;
  float		mt_emht;
  float		dPhi_emht;
  float 	mlt_pt;
  float		mt_emlt;
  float		dPhi_emlt;

  float		event_wgt;
  float		xsec_norm;
  int 		In_cat_WH;  // need to change by hand based on what category is under study
  int 		Sample_ID;
  TString 	Sample_name = "";// for identifying sample in minintuple

  MuPairInfo 	dimu;
  MuonInfo 	mu_1;
  MuonInfo	mu_2;
  MuonInfo	extra_mu;
  EleInfo 	ele;
  JetPairInfo	dijet;
  JetInfo	jet1;
  JetInfo	jet2;
  JetInfo	jet0;
  MhtInfo	mht_info;

  TLorentzVector dimu_vec;  
  TLorentzVector mu1_vec;
  TLorentzVector mu2_vec;
  TLorentzVector extra_mu_vec;
  TLorentzVector ele_vec;
  TLorentzVector dijet_vec;
  TLorentzVector jet1_vec;
  TLorentzVector jet2_vec;
  TLorentzVector met_vec;
  TLorentzVector mht_vec;
  TLorentzVector mlt_vec;

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
  Out_Tree->Branch("mu2_pt",    	& mu2_pt,       	"mu2_pt/F");
  Out_Tree->Branch("mu2_eta",   	& mu2_eta,      	"mu2_eta/F");

  Out_Tree->Branch("dimu_mass",		& dimu_mass,		"dimu_mass/F");
  Out_Tree->Branch("dimu_pt", 		& dimu_pt, 		"dimu_pt/F");
  Out_Tree->Branch("dimu_eta", 		& dimu_eta, 		"dimu_eta/F");
  Out_Tree->Branch("dimu_dEta",		& dimu_dEta,		"dimu_dEta/F");
  Out_Tree->Branch("dimu_dPhi",		& dimu_dPhi,		"dimu_dPhi/F");
  Out_Tree->Branch("dimu_dR",		& dimu_dR,		"dimu_dR/F");
  
  Out_Tree->Branch("cts_mu1",           & cts_mu1,              "cts_mu1/F"  );
  Out_Tree->Branch("cts_mu_pos",        & cts_mu_pos,           "cts_mu_pos/F"  );

  Out_Tree->Branch("extra_mu_pt",       & extra_mu_pt,          "extra_mu_pt/F");
  Out_Tree->Branch("extra_mu_eta",  	& extra_mu_eta,     	"extra_mu_abs_eta/F");

  //ele var
  Out_Tree->Branch("ele_pt", 		& ele_pt,		"ele_pt/F");
  Out_Tree->Branch("ele_eta", 		& ele_eta, 		"ele_eta/F");

  Out_Tree->Branch("cts_edimu", 	& cts_edimu,		"cts_edimu/F" );
  Out_Tree->Branch("edimu_mass",	& edimu_mass,		"edimu_mass/F");  // how much is it correlated with dimu_mass, need to check
  Out_Tree->Branch("edimu_pt",		& edimu_pt,		"edimu_pt/F");
  Out_Tree->Branch("edimu_eta",		& edimu_eta,		"edimu_eta/F");
  Out_Tree->Branch("edimu_dEta",	& edimu_dEta,		"edimu_dEta/F");
  Out_Tree->Branch("edimu_dPhi",	& edimu_dPhi,		"edimu_dPhi/F");
  Out_Tree->Branch("edimu_dR",		& edimu_dR,		"edimu_dR/F");

  //ele_muon var
  Out_Tree->Branch("cts_emuSS",         & cts_emuSS,            "cts_emuSS/F" );
  Out_Tree->Branch("cts_emuOS",         & cts_emuOS,            "cts_emuOS/F" );
  Out_Tree->Branch("emuSS_pt",		& emuSS_pt,		"emuSS_pt/F");
  Out_Tree->Branch("emuSS_eta",		& emuSS_eta,		"emuSS_eta/F");
  Out_Tree->Branch("emuSS_dEta",	& emuSS_dEta,		"emuSS_dEta/F");
  Out_Tree->Branch("emuSS_dPhi",	& emuSS_dPhi,		"emuSS_dPhi/F");
  Out_Tree->Branch("emuSS_dR",		& emuSS_dR,		"emuSS_dR");
  Out_Tree->Branch("emuOS_pt",		& emuOS_pt,		"emuOS_pt/F");
  Out_Tree->Branch("emuOS_eta",		& emuOS_eta,		"emuOS_eta/F");
  Out_Tree->Branch("emuOS_dEta",	& emuOS_dEta,		"emuOS_dEta/F");
  Out_Tree->Branch("emuOS_dPhi",	& emuOS_dPhi,		"emuOS_dPhi/F");
  Out_Tree->Branch("emuOS_dR",		& emuOS_dR,		"emuOS_dR");

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
  Out_Tree->Branch("mt_emet",		& mt_emet,   		"mt_emet/F");
  Out_Tree->Branch("dPhi_emet",		& dPhi_emet, 		"dPhi_emet/F");
  Out_Tree->Branch("mht_pt",		& mht_pt,       	"mht_pt/F");
  Out_Tree->Branch("mht_mass",		& mht_mass,  		"mht_mass/F");
  Out_Tree->Branch("mt_emht",		& mt_emht,   		"mt_emht/F");
  Out_Tree->Branch("dPhi_emht",		& dPhi_emht, 		"dPhi_emht/F");
  Out_Tree->Branch("mlt_pt",		& mlt_pt,       	"mlt_pt/F");
  Out_Tree->Branch("mt_emlt",		& mt_emlt,   		"mt_emlt/F");
  Out_Tree->Branch("dPhi_emlt",		& dPhi_emlt, 		"dPhi_emlt/F");

  //event var
  Out_Tree->Branch("nJets",		& nJets,		"nJets/I");
  Out_Tree->Branch("nBJets_Loose",	& nBJets_Loose,		"nBJets_Loose/I");
  Out_Tree->Branch("nBJets_Med",        & nBJets_Med,           "nBJets_Med/I");
  Out_Tree->Branch("nBJets_Tight",      & nBJets_Tight,         "nBJets_Tight/I");
  Out_Tree->Branch("nFwdJets",		& nFwdJets,		"nFwdJets/I");
  Out_Tree->Branch("nCentJets",		& nCentJets,		"nCentJets/I");
  Out_Tree->Branch("nMuons",		& nMuons,		"nMuons/I");
  Out_Tree->Branch("nEles",		& nEles,		"nEles/I");  

  //weights and spectator 
  Out_Tree->Branch("event_wgt", 	& event_wgt, 		"event_wgt/F");
  Out_Tree->Branch("xsec_norm",		& xsec_norm,		"xsec_norm/F");
  Out_Tree->Branch("In_cat_WH",		& In_cat_WH,		"In_cat_WH/I");
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

  evt_sel.muPair_mass_min = 105; // Require at least one Higgs candidate pair
  obj_sel.mu_pt_min       =  10; // Lower muon pT threshold for muons not from Higgs

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
    In_cat_WH		= 0;
    Sample_ID		= -999;
    Sample_name 	= "none";

    nMuons		= -999;
    nEles		= -999;
    nJets		= -999;
    nBJets_Loose	= -999;
    nBJets_Med  	= -999;
    nBJets_Tight	= -999;
    nFwdJets		= -999;
    nCentJets		= -999;

    dimu_mass		= -999;
    dimu_pt		= -999;
    dimu_eta		= -999;
    dimu_dEta		= -999;
    dimu_dPhi		= -999;
    dimu_dR		= -999;
    cts_mu1             = -999;
    cts_mu_pos          = -999;

    mu1_pt		= -999;
    mu1_eta		= -999;
    mu2_pt		= -999;
    mu2_eta		= -999;
    extra_mu_pt		= -999;
    extra_mu_eta	= -999;

    ele_pt              = -999;
    ele_eta         	= -999;
    cts_edimu		= -999;
    edimu_mass		= -999;
    edimu_pt		= -999;		
    edimu_eta		= -999;	
    edimu_dEta		= -999;	
    edimu_dPhi		= -999;
    edimu_dR		= -999;		

    cts_emuSS 		= -999;
    cts_emuOS		= -999;
    emuSS_pt		= -999;		
    emuSS_eta		= -999;
    emuSS_dEta		= -999;	
    emuSS_dPhi		= -999;	
    emuSS_dR		= -999;		
    emuOS_pt		= -999;		
    emuOS_eta		= -999;	
    emuOS_dEta		= -999;	
    emuOS_dPhi		= -999;	
    emuOS_dR		= -999;		

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
    mt_emet		= -999;
    dPhi_emet		= -999;
    mht_pt		= -999;
    mht_mass		= -999;
    mt_emht		= -999;
    dPhi_emht		= -999;
    mlt_pt		= -999;
    mt_emlt		= -999;
    dPhi_emlt		= -999;

    dimu	.init();
    mu_1	.init();
    mu_2	.init();
    extra_mu	.init();
    ele		.init();
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

	dimu = SelectedCandPair(obj_sel, br);     //  must have this object for the next line to work, somehow
        dimu_vec = FourVec( SelectedCandPair(obj_sel, br), *br.muons, PTC);
        if ( dimu_vec.M() < 105 ||
             dimu_vec.M() > 160 ) continue;

	MuonInfos muons = SelectedMuons(obj_sel, br);
        EleInfos  eles  = SelectedEles(obj_sel, br);
        JetInfos  jets  = SelectedJets(obj_sel, br);

	mu_1 = br.muons->at(dimu.iMu1);
        mu_2 = br.muons->at(dimu.iMu2);
        mu1_vec = FourVec(br.muons->at(dimu.iMu1), PTC);
        mu2_vec = FourVec(br.muons->at(dimu.iMu2), PTC);

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

	dimu_mass 	= dimu_vec.M();                	// all variables filled with vec in order to  
        dimu_pt 	= dimu_vec.Pt();		// work with different PTC
        dimu_eta 	= dimu_vec.Eta();
	dimu_dEta	= mu1_vec.Eta() - mu2_vec.Eta();   
	dimu_dPhi	= mu1_vec.DeltaPhi(mu2_vec);  
	dimu_dR		= mu1_vec.DeltaR(mu2_vec);

	mu1_pt = mu1_vec.Pt();
	mu1_eta = mu1_vec.Eta();
   	mu2_pt = mu2_vec.Pt();
	mu2_eta = mu2_vec.Eta();

	cts_mu1 = CosThetaStar(mu1_vec, mu2_vec);
        cts_mu_pos = mu_1.charge == 1 ? CosThetaStar(mu1_vec, mu2_vec) : CosThetaStar(mu2_vec, mu1_vec);


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

	/////////////////////////
	///  toy WH_3l_e cat  ///
	/////////////////////////
	if (true || OPT_CUT == "WH_3l_e") {
	  obj_sel.ele_pt_min = 20;
	  obj_sel.ele_ID_cut = "tight";
	  obj_sel.ele_iso_max = 0.12;
	  EleInfos  eles  = SelectedEles(obj_sel, br);
	  if (SelectedEles(obj_sel, br).size() > 0) {
            ele  = eles.at(0);
            ele_vec = FourVec(ele);
            ele_pt = ele_vec.Pt();
            ele_eta = ele_vec.Eta();

            TLorentzVector edimu_temp = ele_vec + dimu_vec;
            cts_edimu  	= CosThetaStar(ele_vec, dimu_vec);
            edimu_mass 	= edimu_temp.M();
	    edimu_pt	= edimu_temp.Pt(); 
	    edimu_eta	= edimu_temp.Eta();
	    edimu_dEta	= ele_vec.Eta() - dimu_vec.Eta();
	    edimu_dPhi	= ele_vec.DeltaPhi(dimu_vec);
	    edimu_dR	= ele_vec.DeltaR(dimu_vec);

	    TLorentzVector emu1_vec = ele_vec + mu1_vec;
            TLorentzVector emu2_vec = ele_vec + mu2_vec;
            if (ele.charge == mu_1.charge) {
                cts_emuSS   = CosThetaStar(ele_vec, mu1_vec);
                cts_emuOS   = CosThetaStar(ele_vec, mu2_vec);
		emuSS_pt    = emu1_vec.Pt();	
                emuSS_eta   = emu1_vec.Eta();   
                emuSS_dEta  = ele_vec.Eta() - mu1_vec.Eta();   
                emuSS_dPhi  = ele_vec.DeltaPhi(mu1_vec);   
                emuSS_dR    = ele_vec.DeltaR(mu1_vec);        
                emuOS_pt    = emu2_vec.Pt();   
                emuOS_eta   = emu2_vec.Eta();   
                emuOS_dEta  = ele_vec.Eta() - mu2_vec.Eta();   
                emuOS_dPhi  = ele_vec.DeltaPhi(mu2_vec);   
                emuOS_dR    = ele_vec.DeltaR(mu2_vec);         
            }
            else {
                cts_emuSS   = CosThetaStar(ele_vec, mu2_vec);
                cts_emuOS   = CosThetaStar(ele_vec, mu1_vec);
		emuSS_pt    = emu2_vec.Pt();                   	
                emuSS_eta   = emu2_vec.Eta();                      
                emuSS_dEta  = ele_vec.Eta() - mu2_vec.Eta();        
                emuSS_dPhi  = ele_vec.DeltaPhi(mu2_vec);            
                emuSS_dR    = ele_vec.DeltaR(mu2_vec);              
                emuOS_pt    = emu1_vec.Pt();	          	
                emuOS_eta   = emu1_vec.Eta();                     
                emuOS_dEta  = ele_vec.Eta() - mu1_vec.Eta();       
                emuOS_dPhi  = ele_vec.DeltaPhi(mu1_vec);           
                emuOS_dR    = ele_vec.DeltaR(mu1_vec);             
            }
          }
	  
	  TLorentzVector mlt_vec  = -FourVec(ele,"T") - FourVec(dimu,*br.muons,PTC,"T");
	  TLorentzVector eMET_temp = FourVec(ele,"T") + met_vec;
	  TLorentzVector eMHT_temp = FourVec(ele,"T") + mht_vec;
	  TLorentzVector eMLT_temp = FourVec(ele,"T") + mlt_vec;

	  met_pt	= met_vec.Pt();	  
          mt_emet	= eMET_temp.M();
          dPhi_emet	= ele_vec.DeltaPhi(met_vec);
          mht_pt	= mht_vec.Pt();
          mht_mass	= mht_info.mass;
          mt_emht	= eMHT_temp.M();
          dPhi_emht	= ele_vec.DeltaPhi(mht_vec);
          mlt_pt	= mlt_vec.Pt();
          mt_emlt	= eMLT_temp.M();
          dPhi_emlt	= ele_vec.DeltaPhi(mlt_vec);

	  if ( ElePass(obj_sel, ele) && nBJets_Med == 0 && mt_emet < 150 )  In_cat_WH = true;
	  else continue;
	}	

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
