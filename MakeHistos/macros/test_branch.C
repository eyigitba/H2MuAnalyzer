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
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information

// const TString IN_DIR   = "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/H2Mu_WH_pos/170315_105045/0000";
// const TString SAMPLE   = "H2Mu_WH_pos";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WZ_3l_AMC";
// const TString SAMPLE   = "WZ_3l_AMC";
const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
const TString SAMPLE   = "SingleMu";
const std::string YEAR = "2016";
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {"Presel2016"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void test_branch( TString sample = "", TString in_dir = "", TString out_dir = "",
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
    in_file_name.Form("root://eoscms.cern.ch/%s/%s", in_dir.Data(), in_files.at(i).Data());
    std::cout << "Adding file " << in_file_name.Data() << std::endl;
    in_file_names.push_back(in_file_name);
  }
  if (in_files.size() == 0) {
    for (int i = MIN_FILE; i <= MAX_FILE; i++) {
      // in_file_name.Form("root://eoscms.cern.ch/%s/tuple_%d.root", in_dir.Data(), i);
      in_file_name.Form("root://eoscms.cern.ch/%s/NTuple_0.root", in_dir.Data());
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
  float		event_wgt;
  int		nMuons;
  int		nEles;
  int		nJets;
  int		nBJets_Loose;
  int           nBJets_Med;
  int           nBJets_Tight;
  int		nFwdJets;
  int 		nCentJets;

  float		dimu_pt;
  float		dimu_eta;
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

  float		dijet_mass;
  float		dijet_pt;
  float		dijet_eta;
  float		jet1_pt;
  float		jet1_eta;
  float		jet2_pt;
  float 	jet2_eta;
  float		jet0_pt;
  float		jet0_eta;

  int 		In_cat_WH;  // need to change by hand based on what category is under study
  int 		Sample_type;
  TString 	Sample_name = sample;// for identifying sample in minintuple

  MuPairInfo 	dimu;
  MuonInfo 	mu_1;
  MuonInfo	mu_2;
  MuonInfo	extra_mu;
  EleInfo 	ele;
  JetPairInfo	dijet;
  JetInfo	jet1;
  JetInfo	jet2;
  JetInfo	jet0;

  TLorentzVector dimu_vec;  
  TLorentzVector mu1_vec;
  TLorentzVector mu2_vec;
  TLorentzVector extra_mu_vec;
  TLorentzVector ele_vec;
  TLorentzVector dijet_vec;
  TLorentzVector jet1_vec;
  TLorentzVector jet2_vec;
  TLorentzVector met_vec;

  // declaration end.

  TString out_file_name;
  if (out_file_str.Length() > 0) out_file_name.Form( "%s/histos_%s_%s.root",    out_dir.Data(), sample.Data(), out_file_str.Data() );
  else                           out_file_name.Form( "%s/histos_%s_%d_%d.root", out_dir.Data(), sample.Data(), MIN_FILE, MAX_FILE );
  std::cout << "\nCreating output file " << out_file_name.Data() << std::endl;
  Out_File = TFile::Open( out_file_name, "RECREATE" );


  Out_Tree = PlantTree("tree", "flat_tree_minintuple");    // if trees in different cats share the same branch addresses. will they be filled with the same events?  XWZ - 09.11.2018

  Out_Tree->Branch("mu1_pt", 		& mu1_pt, 		"mu1_pt/F");
  Out_Tree->Branch("mu1_eta", 		& mu1_eta,		"mu1_eta/F");	
  Out_Tree->Branch("mu2_pt",    	& mu2_pt,       	"mu2_pt/F");
  Out_Tree->Branch("mu2_eta",   	& mu2_eta,      	"mu2_eta/F");
  Out_Tree->Branch("dimu_pt", 		& dimu_pt, 		"dimu_pt/F");
  Out_Tree->Branch("dimu_eta", 		& dimu_eta, 		"dimu_eta/F");
  Out_Tree->Branch("ele_pt", 		& ele_pt,		"ele_pt/F");
  Out_Tree->Branch("ele_eta", 		& ele_eta, 		"ele_eta/F");
  Out_Tree->Branch("extra_mu_pt",    	& extra_mu_pt,       	"mu2_pt/F");
  Out_Tree->Branch("extra_mu_eta",   	& extra_mu_eta,      	"mu2_eta/F");

  Out_Tree->Branch("cts_mu1",		& cts_mu1, 		"cts_mu1/F"  );
  Out_Tree->Branch("cts_mu_pos",        & cts_mu_pos,           "cts_mu_pos/F"  );
  Out_Tree->Branch("cts_edimu", 	& cts_edimu,		"cts_edimu/F" );
  Out_Tree->Branch("cts_emuSS", 	& cts_emuSS,		"cts_emuSS/F" );
  Out_Tree->Branch("cts_emuOS", 	& cts_emuOS,		"cts_emuOS/F" );

  Out_Tree->Branch("dijet_mass",	& dijet_mass,		"dijet_mass/F");
  Out_Tree->Branch("dijet_pt",		& dijet_pt,		"dijet_pt/F");
  Out_Tree->Branch("dijet_eta",		& dijet_eta,		"dijet_eta/F");
  Out_Tree->Branch("jet1_pt",		& jet1_pt,		"jet1_pt/F");
  Out_Tree->Branch("jet1_eta",		& jet1_eta,		"jet1_eta/F");
  Out_Tree->Branch("jet2_pt",		& jet2_pt,		"jet2_pt/F");
  Out_Tree->Branch("jet2_eta",		& jet2_eta,		"jet2_eta/F");
  Out_Tree->Branch("jet0_pt",           & jet0_pt,              "jet0_pt/F");
  Out_Tree->Branch("jet0_eta",          & jet0_eta,             "jet0_eta/F");


  Out_Tree->Branch("nJets",		& nJets,		"nJets/I");
  Out_Tree->Branch("nBJets_Loose",	& nBJets_Loose,		"nBJets_Loose/I");
  Out_Tree->Branch("nBJets_Med",        & nBJets_Med,           "nBJets_Med/I");
  Out_Tree->Branch("nBJets_Tight",      & nBJets_Tight,         "nBJets_Tight/I");
  Out_Tree->Branch("nFwdJets",		& nFwdJets,		"nFwdJets/I");
  Out_Tree->Branch("nCentJets",		& nCentJets,		"nCentJets/I");
  Out_Tree->Branch("nMuons",		& nMuons,		"nMuons/I");
  Out_Tree->Branch("nEles",		& nEles,		"nEles/I");  


  Out_Tree->Branch("event_wgt", 	& event_wgt, 		"event_wgt/F");
  Out_Tree->Branch("In_cat_WH",		& In_cat_WH,		"In_cat_WH/I");
  Out_Tree->Branch("Sample_type",	& Sample_type,		"Sample_type/I");
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
      SetBranchAddresses(*in_chain, br, {YEAR}, false); // Options in {} include "JES", "Flags", and "SFs"
    else
      SetBranchAddresses(*in_chain, br, {YEAR, "GEN", "Wgts"}, false); // Options in {} include "JES", "Flags", and "SFs"
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
    if (YEAR == "2016") {
      jets_tmp = ConvertSlimJets(*(br.slimJets));
      br.jets  = &jets_tmp;
    }


    /////////////////////////////////
    ///  initialize all branches  ///
    /////////////////////////////////
    In_cat_WH	= 0;
    Sample_type = -999;

    event_wgt	= -999;
    nMuons	= -999;
    nEles	= -999;
    nJets	= -999;
    nBJets_Loose= -999;
    nBJets_Med  = -999;
    nBJets_Tight= -999;
    nFwdJets	= -999;
    nCentJets	= -999;

    dimu_pt	= -999;
    dimu_eta	= -999;
    mu1_pt	= -999;
    mu1_eta	= -999;
    mu2_pt	= -999;
    mu2_eta	= -999;
    ele_pt	= -999;
    ele_eta	= -999;
    extra_mu_pt	= -999;
    extra_mu_eta= -999;

    cts_mu1	= -999;
    cts_mu_pos  = -999;
    cts_edimu	= -999;
    cts_emuSS 	= -999;
    cts_emuOS	= -999;

    dijet_mass	= -999;
    dijet_pt	= -999;
    dijet_eta	= -999;
    jet1_pt	= -999;
    jet1_eta	= -999;
    jet2_pt	= -999;
    jet2_eta	= -999;
    jet0_pt	= -999;
    jet0_eta	= -999;

    dimu	.init();
    mu_1	.init();
    mu_2	.init();
    extra_mu	.init();
    ele		.init();
    dijet	.init();
    jet1	.init();
    jet2	.init();
    jet0	.init();

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
      if ( sample.Contains("SingleMu") ) Sample_type = -1;
      else if (sample.Contains("H2Mu") ) Sample_type = 1;
      else Sample_type = 0;

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
        dimu_vec = FourVec( SelectedCandPair(obj_sel, br), PTC);
        dimu_pt = dimu_vec.Pt();
        dimu_eta = dimu_vec.Eta();
        if ( dimu_vec.M() < 105 ||
             dimu_vec.M() > 160 ) continue;

	MuonInfos muons = SelectedMuons(obj_sel, br);
	EleInfos  eles  = SelectedEles(obj_sel, br);
	JetInfos  jets  = SelectedJets(obj_sel, br);

	nMuons = muons.size();
	nEles = eles.size();
	nJets = jets.size();
 	nFwdJets     = SelectedJets(obj_sel, br, "Forward").size();
	nCentJets    = SelectedJets(obj_sel, br, "Central").size();
	nBJets_Loose = SelectedJets(obj_sel, br, "BTagLoose").size();
        nBJets_Med   = SelectedJets(obj_sel, br, "BTagMedium").size();
	nBJets_Tight = SelectedJets(obj_sel, br, "BTagTight").size();

	mu_1 = br.muons->at(dimu.iMu1);
        mu_2 = br.muons->at(dimu.iMu2);
	mu1_vec = FourVec(br.muons->at(dimu.iMu1), PTC);
        mu2_vec = FourVec(br.muons->at(dimu.iMu2), PTC);
	mu1_pt = mu1_vec.Pt();
	mu1_eta = mu1_vec.Eta();
   	mu2_pt = mu2_vec.Pt();
	mu2_eta = mu2_vec.Eta();

	met_vec = FourVec(*br.met);

    	cts_mu1 = CosThetaStar(mu1_vec, mu2_vec);
	cts_mu_pos = mu_1.charge == 1 ? CosThetaStar(mu1_vec, mu2_vec) : CosThetaStar(mu2_vec, mu1_vec);

	if (SelectedJetPairs(obj_sel, br).size() > 0) {
	  dijet = SelectedJetPairs(obj_sel, br).at(0);
	  jet1  = br.jets->at(dijet.iJet1);
	  jet2  = br.jets->at(dijet.iJet2);

	  dijet_mass = dijet.mass;
	  dijet_pt   = dijet.pt;
	  dijet_eta  = dijet.eta;
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

	if (SelectedEles(obj_sel, br).size() > 0) {
	  ele  = eles.at(0);
          ele_vec = FourVec(ele);
	  ele_pt = ele_vec.Pt();
	  ele_eta = ele_vec.Eta();

	  cts_edimu = CosThetaStar(ele_vec, dimu_vec);
	  if (ele.charge == mu_1.charge) {
	      cts_emuSS = CosThetaStar(ele_vec, mu1_vec);
	      cts_emuOS = CosThetaStar(ele_vec, mu2_vec);
	  }
	  else {
	      cts_emuSS = CosThetaStar(ele_vec, mu2_vec);
	      cts_emuOS = CosThetaStar(ele_vec, mu1_vec);
	  }
	}


	/////////////////////////
	///  toy WH_3l_e cat  ///
	/////////////////////////
	if (true || OPT_CUT == "WH_3l_e") {
	  obj_sel.ele_pt_min = 20;
	  obj_sel.ele_ID_cut = "tight";
	  obj_sel.ele_iso_max = 0.12;
	  TLorentzVector MT_temp = FourVec(ele,"T") + met_vec;
	  if ( ElePass(obj_sel, ele) && nBJets_Med == 0 && MT_temp.M() < 150 )  {
		std::cout << "in cat" << std::endl;
		In_cat_WH = true;
	  }
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


  // Place histograms made with separate selection and cateogry cuts into separate files
     
  // Write output file
  Out_File->cd();
  Out_Tree->Write();
  Out_File->Write();
      
  
  std::cout << "\nExiting ()\n";
  
} // End void WH_lep()
