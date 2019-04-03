
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

// #include "H2MuAnalyzer/MakeHistos/interface/SampleDatabase2016.h" // Input data and MC samples

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 1;     // Maximum index of input files to process
const int MAX_EVT  = -1;    // Maximum number of events to process
const int PRT_EVT  = 10000;  // Print every N events
const float SAMP_WGT = 1.0;
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information

const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1/ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8/H2Mu_ttH_125";
const TString SAMPLE   = "H2Mu_ttH_125";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/WZ_3l_AMC";
// const TString SAMPLE   = "WZ_3l_AMC";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM/SingleMuon";
// const TString SAMPLE   = "SingleMu";

const std::string YEAR = "2017";
const std::string SLIM = "Slim";  // "Slim" or "notSlim" - original 2016 NTuples were in "Slim" format, some 2017 NTuples are "Slim"
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {"Presel2017"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"3mu", "e2mu"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = { "NONE", "ge2j_btag_mass12", "noZ_ge3j_btag_mass12", "tightLepMVA_ge2j_btag_mass12",
					    "looseLepMVA_noZ_ge3j_btag_mass12",
					    "tightLepMVA_noZ_ge3j_btag_mass12",
					    "tightLepCut_noZ_ge3j_btag_mass12" }; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void ttH_3l( TString sample = "", TString in_dir = "", TString out_dir = "",
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

  // Initialize set of pointers to all branches in tree
  NTupleBranches br;

  // Initialize empty map of histogram names to histograms
  std::map<TString, TH1*> h_map_1D;
  std::map<TString, TH2*> h_map_2D;

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


  // Configuration for object selection, event selection, and object weighting
  ObjectSelectionConfig obj_sel;
  EventSelectionConfig  evt_sel;
  EventWeightConfig     evt_wgt;
  ConfigureObjectSelection(obj_sel, YEAR, "lepMVA");
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

  evt_sel.muPair_mass_min = 105; // Require at least one Higgs candidate pair
  // obj_sel.muPair_Higgs = "sort_WH_3_mu_v1"; // Choose Higgs candidate based on MT(W muon, MET)

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

    if (pass_sel_cuts) {
      
      // Get event weight for MC, defined in src/EventWeight.cc
      float event_wgt = ( sample.Contains("SingleMu") ? 1.0 : EventWeight(br, evt_wgt, verbose) );

      // Initialize the selected Higgs candidate dimuon pair
      std::vector<int> iCandPairs;  // Indices of Higgs candidate pairs (can have up to 2)
      MuPairInfo        H_pair;     // When filling histograms, have only one candidate pair at a time
      TLorentzVector    H_pair_vec;
      
      bool MU = false;  // Says whether event has 3 muons or 1 electron + 2 muons

      /////////////////////////////////////////////////////////////////////////////////////////
      ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
      /////////////////////////////////////////////////////////////////////////////////////////

      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	std::string OPT_CUT = OPT_CUTS.at(iOpt);
	if (OPT_CUT == "3mu") {
	  // Require exactly 0 selected electrons
	  if ( SelectedEles(obj_sel, br).size() != 0 ) continue;
	  // Require exactly 3 selected muons whose charge sums to +/-1
	  if ( SelectedMuPairs(obj_sel, br).size() != 2 ) continue;
	  // Require at least one Higgs candidate dimuon pair to fall inside mass window
	  for (int iPair = 0; iPair < 2; iPair++) {
	    H_pair     = SelectedMuPairs(obj_sel, br).at(iPair);
	    H_pair_vec = FourVec(H_pair, PTC);
	    if ( H_pair_vec.M() < 105 ||
		 H_pair_vec.M() > 160 ) continue;
	    iCandPairs.push_back(iPair);
	    MU = true;
	  }
	}
	else if (OPT_CUT == "e2mu") {
	  // Require exactly 1 selected electron
	  if ( SelectedEles(obj_sel, br).size() != 1 ) continue;
	  // Require exactly 2 selected muons whose charge sums to +/-1
	  if ( SelectedMuPairs(obj_sel, br).size() != 1 ) continue;
	  // Require selected Higgs candidate dimuon pair to fall inside mass window
	  H_pair     = SelectedMuPairs(obj_sel, br).at(0);
	  H_pair_vec = FourVec(H_pair, PTC);
	  if ( H_pair_vec.M() < 105 ||
	       H_pair_vec.M() > 160 ) continue;
	  iCandPairs.push_back(0);
	  MU = false;
	}
	else assert(OPT_CUT == "3mu" || OPT_CUT == "e2mu");
	// else if (not PassSelection(br, evt_sel, obj_sel, OPT_CUT, verbose)) continue;

	
	//////////////////////////////////////////////////////////////////
	/// Loop through category cuts defined in src/CategoryCuts.cc  ///
	//////////////////////////////////////////////////////////////////
	
	for (int iCat = 0; iCat < CAT_CUTS.size()*iCandPairs.size(); iCat++) {

	  // In case of two valid Higgs candidate pairs, choose based on even/odd iCat
	  if (iCandPairs.size() == 2 ) {
	    int iCandPair = iCat % 2;
	    H_pair     = SelectedMuPairs(obj_sel, br).at(iCandPair);
	    H_pair_vec = FourVec(H_pair, PTC);
	  } else {  // Otherwise pick the single pair that fell in the mass [105, 160] range
	    H_pair     = SelectedMuPairs(obj_sel, br).at(iCandPairs.at(0));
	    H_pair_vec = FourVec(H_pair, PTC);
	  }

	  //////////////////////////////////////////////////////
	  ///  Compute variables relevent for category cuts  ///
	  //////////////////////////////////////////////////////

	  // Get selected muons and electrons in event
	  MuonInfos muons = SelectedMuons(obj_sel, br);
	  EleInfos  eles  = SelectedEles(obj_sel, br);
	  int sum_lep_charge = 0;
	  if (MU) {
	    assert(muons.size() == 3 && eles.size() == 0);
	    sum_lep_charge = muons.at(0).charge + muons.at(1).charge + muons.at(2).charge;
	  } else {
	    assert(muon.size() == 2 && eles.size() == 1);
	    sum_lep_charge = muons.at(0).charge + muons.at(1).charge + eles.at(0).charge;
	  }
	  assert(abs(sum_lep_charge) == 1);
	  
	  MuPairInfo H_true;     // Dimuon pair matched to GEN-level Higgs decay
	  MuPairInfo Z_true;     // Dimuon pair matched to GEN-level Z boson decay
	  MuonInfo   H_mu1;      // Higher pT muon from chosen Higgs candidate pair
	  MuonInfo   H_mu2;      // Lower pT muon from chosen Higgs candidate pair
	  MuonInfo   nonH_mu;    // Muon not from chosen Higgs candidate pair
	  EleInfo    ele;        // Electron
	  MuPairInfo nonH_pair;  // OS dimuon pair not chosen as Higgs candidate
	  MuonInfo   SS_mu1;     // Same-sign muon with higher pT
	  MuonInfo   SS_mu2;     // Same-sign muon with lower pT
	  MuonInfo   OS_mu;      // Opposite-sign muon

	  TLorentzVector H_true_vec;
	  TLorentzVector Z_true_vec;
	  TLorentzVector H_mu1_vec;
	  TLorentzVector H_mu2_vec;
	  TLorentzVector nonH_mu_vec;
	  TLorentzVector ele_vec;
	  TLorentzVector nonH_lep_vec;   // Lepton not chosen to be in Higgs candidate pair
	  TLorentzVector nonH_lep_vecT;  // Transverse component only (eta = 0, pZ = 0, mass = 0)
	  TLorentzVector nonH_pair_vec;  // OS lepton pair not chosen as Higgs candidate (mu-mu or ele-mu)
	  TLorentzVector SS_mu1_vec;
	  TLorentzVector SS_mu2_vec;
	  TLorentzVector SS_lep1_vec;  // Same-sign lepton with higher pT
	  TLorentzVector SS_lep2_vec;  // Same-sign lepton with lower pT
	  TLorentzVector OS_mu_vec;

	  H_mu1 = br.muons->at(H_pair.iMu1);
	  H_mu2 = br.muons->at(H_pair.iMu2);
	  H_mu1_vec = FourVec(H_mu1, PTC);
	  H_mu2_vec = FourVec(H_mu2, PTC);

	  // Loop over selected dimuon pairs
	  for (const auto & muPair : SelectedMuPairs(obj_sel, br)) {
	    // Check if the pair matches a GEN muon pair from H
	    if (not sample.Contains("SingleMu")) {
	      if ( IsGenMatched( muPair, *br.muons, *br.genMuons, "H") ) {
		H_true     = muPair;
		H_true_vec = FourVec(muPair, PTC);
	      }
	      if ( IsGenMatched( muPair, *br.muons, *br.genMuons, "Z") ) {
		Z_true     = muPair;
		Z_true_vec = FourVec(muPair, PTC);
	      }
	    }
	    // Check if the pair does not match the selected Higgs candidate pair
	    if ( muPair.iMu1 != H_pair.iMu1 || muPair.iMu2 != H_pair.iMu2 ) {
	      nonH_pair     = muPair;
	      nonH_pair_vec = FourVec(muPair, PTC);
	    }
	  }

	  // For 3 muon events
	  if (MU) {
	    // Find the muon not in the H pair, and the same-sign and opposite-sign muons
	    for (const auto & mu : muons) {
	      if ( mu.pt != H_mu1.pt && mu.eta != H_mu1.eta &&
		   mu.pt != H_mu2.pt && mu.eta != H_mu2.eta ) {
		assert(nonH_mu.pt <= 0); // We should not have found a W candidate before
		nonH_mu = mu;
		nonH_mu_vec   = FourVec(mu, PTC);
		nonH_lep_vec  = nonH_mu_vec;
		nonH_lep_vecT = FourVec(mu, PTC, "T", *br.muons);
	      }
	      if (mu.charge != sum_lep_charge) {
		OS_mu     = mu;
		OS_mu_vec = FourVec(mu, PTC);
	      }
	      // Loop through all other muons
	      for (const auto & mu2 : muons) {
		// Skip if mu2 is the same as mu
		if (mu.charge == mu2.charge && (mu.pt != mu2.pt || mu.eta != mu2.eta)) {
		  SS_mu1 = ( MuonPt(mu, PTC) > MuonPt(mu2, PTC) ? mu : mu2 );
		  SS_mu2 = ( MuonPt(mu, PTC) > MuonPt(mu2, PTC) ? mu2 : mu );
		  SS_mu1_vec  = FourVec(SS_mu1, PTC);
		  SS_mu2_vec  = FourVec(SS_mu2, PTC);
		  SS_lep1_vec = SS_mu1_vec;
		  SS_lep2_vec = SS_mu2_vec;
		}
	      }
	    }
	    assert(nonH_mu.pt >= obj_sel.mu_pt_min); // We should always find a W candidate
	    assert(SS_mu1.charge == SS_mu2.charge && OS_mu.charge == -1*sum_lep_charge); // We should always find two SS and one OS muon
	  }
	  // For ele + 2 muon events
	  else {
	    ele           = eles.at(0);
	    ele_vec       = FourVec(ele);
	    nonH_lep_vec  = ele_vec;
	    nonH_lep_vecT = FourVec(ele, "T");
	    for (const auto & mu : muons) {
	      if (mu.charge == ele.charge) { SS_mu1 = mu;  SS_mu1_vec = FourVec(mu, PTC); }
	      else                         { OS_mu  = mu;  OS_mu_vec  = FourVec(mu, PTC); }
	    }
	    nonH_pair_vec = ele_vec + OS_mu_vec;
	    SS_lep1_vec   = (SS_mu1_vec.Pt() > ele_vec.Pt() ? SS_mu1_vec : ele_vec);
	    SS_lep2_vec   = (SS_mu1_vec.Pt() > ele_vec.Pt() ? ele_vec : SS_mu1_vec);
	    assert(SS_mu1.charge == ele.charge && OS_mu.charge == -1*sum_lep_charge); // We should always find one SS and one OS muon
	  }
	  
	  TLorentzVector MET_vec          = FourVec(*br.met);
	  TLorentzVector nonH_lep_MET_vec = nonH_lep_vecT + MET_vec;
	  
	  ///////////////////////////////////////////
	  ///  Apply the category selection cuts  ///
	  ///////////////////////////////////////////

	  std::string CAT_CUT   = CAT_CUTS.at(iCat / iCandPairs.size());
	  std::string CAT_UNCUT = CAT_CUT; // Track what sub-strings don't match any known cuts
	  bool pass_cat_cut = true;
	  if ( CAT_CUT.find("mass12") != std::string::npos ) {
	    if ( nonH_pair_vec.M() < 12 )                             { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("mass12"), std::string("mass12").length() );
	  }
	  if ( CAT_CUT.find("noZ") != std::string::npos ) {
	    if ( MU && abs(nonH_pair_vec.M() - 91) < 10 )             { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("noZ"), std::string("noZ").length() );
	  }
	  if ( CAT_CUT.find("ge2j") != std::string::npos ) {
	    if ( SelectedJets(obj_sel, br).size() < 2 )               { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("ge2j"), std::string("ge2j").length() );
	  }
	  if ( CAT_CUT.find("ge3j") != std::string::npos ) {
	    if ( SelectedJets(obj_sel, br).size() < 3 )               { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("ge3j"), std::string("ge3j").length() );
	  }
	  if ( CAT_CUT.find("btag") != std::string::npos ) {
	    if ( SelectedJets(obj_sel, br, "BTagMedium").size() < 1 &&
		 SelectedJets(obj_sel, br, "BTagLoose").size()  < 2 ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("btag"), std::string("btag").length() );
	  }
	  if ( CAT_CUT.find("looseLepMVA") != std::string::npos ) {
	    if ( OS_mu.lepMVA < 0.4 )                                 { pass_cat_cut = false; continue; }
	    if ( MU && (SS_mu1.lepMVA < 0.4 || SS_mu2.lepMVA < 0.4) ) { pass_cat_cut = false; continue; }
	    if (!MU && (SS_mu1.lepMVA < 0.4 ||    ele.lepMVA < 0.4) ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("looseLepMVA"), std::string("looseLepMVA").length() );
	  }
	  if ( CAT_CUT.find("tightLepMVA") != std::string::npos ) {
	    if ( OS_mu.lepMVA < 0.4 )                                 { pass_cat_cut = false; continue; }
	    if ( MU && (SS_mu1.lepMVA < 0.8 || SS_mu2.lepMVA < 0.8) ) { pass_cat_cut = false; continue; }
	    if (!MU && (SS_mu1.lepMVA < 0.8 ||    ele.lepMVA < 0.8) ) { pass_cat_cut = false; continue; }
	    CAT_UNCUT.erase( CAT_UNCUT.find("tightLepMVA"), std::string("tightLepMVA").length() );
	  }
	  if ( CAT_CUT.find("tightLepCut") != std::string::npos ) {
	    if (         MuonPt(OS_mu,  PTC) < 10 || !OS_mu.isMediumID || OS_mu.relIso  > 0.25  ) { pass_cat_cut = false; continue; };
	    if (  MU && (MuonPt(SS_mu1, PTC) < 20 || !SS_mu1.isTightID || SS_mu1.relIso > 0.12) ) { pass_cat_cut = false; continue; };
	    if (  MU && (MuonPt(SS_mu2, PTC) < 20 || !SS_mu2.isTightID || SS_mu2.relIso > 0.12) ) { pass_cat_cut = false; continue; };
	    if ( !MU && (             ele.pt < 20 ||    !ele.isTightID ||    ele.relIso > 0.12) ) { pass_cat_cut = false; continue; };
	    CAT_UNCUT.erase( CAT_UNCUT.find("tightLepCut"), std::string("tightLepCut").length() );
	  }
	  // Remove "_" characters left over after all known category sub-strings have been removed
	  while ( CAT_UNCUT.find("_") != std::string::npos ) {
	    CAT_UNCUT.erase( CAT_UNCUT.find("_"), std::string("_").length() );
	  }
	  // If CAT_CUT contains unknown elements, look for it in the central repository
	  if ( CAT_UNCUT.length() > 0 ) pass_cat_cut = InCategory(obj_sel, br, CAT_CUT, verbose);

	  ///  *** SPECIAL CUTS ***  ///
	  // Only keep signal events if the chosen Higgs candidate pair really comes from the Higgs
	  if (  sample.Contains("H2Mu") && H_true_vec.M() != H_pair_vec.M() ) pass_cat_cut = false;
	  // Throw away background events if the chosen Higgs candidate pair actually comes from the Higgs (e.g. ttH, tHq, etc.)
	  if ( !sample.Contains("H2Mu") && H_true_vec.M() == H_pair_vec.M() ) pass_cat_cut = false;


	  if (not pass_cat_cut) continue;
	  if (verbose) std::cout << "\nPassed cut " << OPT_CUT << ", in category " << CAT_CUT << std::endl;
	  std::string h_pre = (std::string)sample + "_"+OPT_CUT+"_"+CAT_CUT+"_";


	  //////////////////////////////////////
	  ///  Compute a few more variables  ///
	  //////////////////////////////////////

	  // Leptons ordered by pT
	  TLorentzVector lep1_vec;
	  TLorentzVector lep2_vec;
	  TLorentzVector lep3_vec;

	  if ( SS_lep1_vec.Pt() > OS_mu_vec.Pt() ) {
	    lep1_vec = SS_lep1_vec;
	    lep2_vec = (OS_mu_vec.Pt() > SS_lep2_vec.Pt() ? OS_mu_vec : SS_lep2_vec);
	    lep3_vec = (OS_mu_vec.Pt() > SS_lep2_vec.Pt() ? SS_lep2_vec : OS_mu_vec);
	  } else {
	    lep1_vec = OS_mu_vec;
	    lep2_vec = SS_lep1_vec;
	    lep3_vec = SS_lep2_vec;
	  }

	  TLorentzVector SS_pair_vec = SS_lep1_vec + SS_lep2_vec;
	  TLorentzVector trilep_vec  = OS_mu_vec + SS_pair_vec;

	  float SS_mu_MVA  = 0;
	  float SS_ele_MVA = 0;
	  float SS_lep1_iso;
	  float SS_lep2_iso;
	  float SS_lep1_MVA;
	  float SS_lep2_MVA;
	  float SS_lep1_SIP;
	  float SS_lep2_SIP;
	  float SS_lep1_seg;
	  float SS_lep2_seg;
	  bool  SS_lep1_tID;
	  bool  SS_lep2_tID;
	  if (MU) {
	    SS_lep1_iso = SS_mu1.relIso;
	    SS_lep2_iso = SS_mu2.relIso;
	    SS_lep1_MVA = SS_mu1.lepMVA;
	    SS_lep2_MVA = SS_mu2.lepMVA;
	    SS_lep1_SIP = SS_mu1.SIP_3D;
	    SS_lep2_SIP = SS_mu2.SIP_3D;
	    SS_lep1_seg = SS_mu1.segCompat;
	    SS_lep2_seg = SS_mu2.segCompat;
	    SS_lep1_tID = SS_mu1.isTightID;
	    SS_lep2_tID = SS_mu2.isTightID;
	  } else {
	    SS_mu_MVA   = SS_mu1.lepMVA;
	    SS_ele_MVA  = ele.lepMVA;
	    SS_lep1_iso = (SS_mu1_vec.Pt() > ele.pt ? SS_mu1.relIso : ele.relIso);
	    SS_lep2_iso = (SS_mu1_vec.Pt() > ele.pt ? ele.relIso : SS_mu1.relIso);
	    SS_lep1_MVA = (SS_mu1_vec.Pt() > ele.pt ? SS_mu1.lepMVA : ele.lepMVA);
	    SS_lep2_MVA = (SS_mu1_vec.Pt() > ele.pt ? ele.lepMVA : SS_mu1.lepMVA);
	    SS_lep1_SIP = (SS_mu1_vec.Pt() > ele.pt ? SS_mu1.SIP_3D : ele.SIP_3D);
	    SS_lep2_SIP = (SS_mu1_vec.Pt() > ele.pt ? ele.SIP_3D : SS_mu1.SIP_3D);
	    SS_lep1_seg = (SS_mu1_vec.Pt() > ele.pt ? SS_mu1.segCompat : -0.01);
	    SS_lep2_seg = (SS_mu1_vec.Pt() > ele.pt ? -0.01 : SS_mu1.segCompat);
	    SS_lep1_tID = (SS_mu1_vec.Pt() > ele.pt ? SS_mu1.isTightID : ele.isTightID);
	    SS_lep2_tID = (SS_mu1_vec.Pt() > ele.pt ? ele.isTightID : SS_mu1.isTightID);
	  }


	  /////////////////////////////////
	  ///  Generate and fill plots  ///
	  /////////////////////////////////

	  // Store event weights
	  if (not sample.Contains("SingleMu")) {
	    BookAndFill(h_map_1D, h_pre+"PU_wgt",    40, -2, 2, br.PU_wgt  );
	    BookAndFill(h_map_1D, h_pre+"muon_wgt",  40, -2, 2, event_wgt / (br.PU_wgt * br.GEN_wgt) );
	    BookAndFill(h_map_1D, h_pre+"GEN_wgt",   40, -2, 2, br.GEN_wgt );
	    BookAndFill(h_map_1D, h_pre+"event_wgt", 40, -2, 2, event_wgt  );
	  }


	  // Plot kinematic histograms
	  BookAndFill(h_map_1D, h_pre+"nJets",      13, -0.5, 12.5, SelectedJets(obj_sel, br).size(),               event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nBJetsLoose", 6, -0.5,  5.5, SelectedJets(obj_sel, br, "BTagLoose").size(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nBJetsMed",   5, -0.5,  4.5, SelectedJets(obj_sel, br, "BTagMedium").size(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nBJetsTight", 4, -0.5,  3.5, SelectedJets(obj_sel, br, "BTagTight").size(),  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nJetsCent",  11, -0.5, 10.5, SelectedJets(obj_sel, br, "Central").size(),    event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nJetsFwd",    6, -0.5,  5.5, SelectedJets(obj_sel, br, "Forward").size(),    event_wgt );

	  BookAndFill(h_map_1D, h_pre+"MET", 20, 0, 200, br.met->pt, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"MHT", 20, 0, 200, br.mht->pt, event_wgt );

	  BookAndFill(h_map_1D, h_pre+"lep1_pt", 30, 0, 300, lep1_vec.Pt(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"lep2_pt", 20, 0, 200, lep2_vec.Pt(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"lep3_pt", 20, 0, 100, lep3_vec.Pt(), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"H_mu1_pt",    30, 0, 300, H_mu1_vec.Pt(),    event_wgt );
	  BookAndFill(h_map_1D, h_pre+"H_mu2_pt",    20, 0, 200, H_mu2_vec.Pt(),    event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonH_lep_pt", 20, 0, 200, nonH_lep_vec.Pt(), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_pt",   20, 0, 200, OS_mu_vec.Pt(),   event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep1_pt", 20, 0, 200, SS_lep1_vec.Pt(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep2_pt", 20, 0, 100, SS_lep2_vec.Pt(), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"lep1_eta", 24, -2.4, 2.4, lep1_vec.Eta(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"lep2_eta", 24, -2.4, 2.4, lep2_vec.Eta(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"lep3_eta", 24, -2.4, 2.4, lep3_vec.Eta(), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"H_mu1_eta",    24, -2.4, 2.4, H_mu1_vec.Eta(),    event_wgt );
	  BookAndFill(h_map_1D, h_pre+"H_mu2_eta",    24, -2.4, 2.4, H_mu2_vec.Eta(),    event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonH_lep_eta", 24, -2.4, 2.4, nonH_lep_vec.Eta(), event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_iso",   15, 0, 0.3, OS_mu.relIso, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep1_iso", 15, 0, 0.3, SS_lep1_iso,  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep2_iso", 15, 0, 0.3, SS_lep2_iso,  event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_tightID",   2, -0.5, 1.5, OS_mu.isTightID, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep1_tightID", 2, -0.5, 1.5, SS_lep1_tID,     event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep2_tightID", 2, -0.5, 1.5, SS_lep2_tID,     event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_lepMVA",   40, -1, 1, OS_mu.lepMVA, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep1_lepMVA", 40, -1, 1, SS_lep1_MVA,  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep2_lepMVA", 40, -1, 1, SS_lep2_MVA,  event_wgt );
	  BookAndFill(h_map_2D, h_pre+"SS_lep2_vs_lep2_lepMVA", 40, -1, 1, 40, -1, 1, SS_lep1_MVA, SS_lep2_MVA,  event_wgt );

	  BookAndFill(h_map_1D, h_pre+"SS_mu_lepMVA",  40, -1, 1, SS_mu_MVA,  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_ele_lepMVA", 40, -1, 1, SS_ele_MVA, event_wgt );
	  BookAndFill(h_map_2D, h_pre+"SS_ele_vs_mu_lepMVA", 40, -1, 1, 40, -1, 1, SS_mu_MVA, SS_ele_MVA, event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_SIP",   50, 0, 10, OS_mu.SIP_3D, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep1_SIP", 50, 0, 10, SS_lep1_SIP,  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep2_SIP", 50, 0, 10, SS_lep2_SIP,  event_wgt );

	  BookAndFill(h_map_1D, h_pre+"OS_mu_segCompat",   50, 0, 1, OS_mu.segCompat, event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep1_segCompat", 50, 0, 1, SS_lep1_seg,  event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_lep2_segCompat", 50, 0, 1, SS_lep2_seg,  event_wgt );

	  if (not sample.Contains("SingleMu")) {
	    BookAndFill(h_map_1D, h_pre+"H_mass_true", 55, 105, 160, MuPairMass(H_true, PTC), event_wgt, false );  // Don't include overflow
	    BookAndFill(h_map_1D, h_pre+"Z_mass_true", 40,   1, 201, MuPairMass(Z_true, PTC), event_wgt, false );  // Don't include overflow
	  }
	  BookAndFill(h_map_1D, h_pre+"H_mass_zoom",        55, 105, 160, H_pair_vec.M(),       event_wgt );
	  BookAndFill(h_map_1D, h_pre+"H_mass_on",          11, 105, 160, H_pair_vec.M(),       event_wgt );
	  BookAndFill(h_map_1D, h_pre+"H_mass_off",         40,   0, 400, H_pair_vec.M(),       event_wgt );
	  BookAndFill(h_map_1D, h_pre+"H_pt",               30,   0, 300, H_pair_vec.Pt(),      event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonH_OS_dilep_mass", 40,   0, 400, nonH_pair_vec.M(),    event_wgt );
	  BookAndFill(h_map_1D, h_pre+"nonH_lep_MET_MT",    40,   0, 200, nonH_lep_MET_vec.M(), event_wgt );
	  BookAndFill(h_map_1D, h_pre+"SS_dilep_mass",      40,   0, 400, SS_pair_vec.M(),      event_wgt );
	  BookAndFill(h_map_1D, h_pre+"trilep_mass",        30,   0, 600, trilep_vec.M(),       event_wgt );

	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size()*iCandCuts.size(); iCat++)
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
    } // End conditional: if (pass_sel_cuts)


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

  // Concatenate basic selection cuts into string for output file name
  std::string selStr = "";
  for (int i = 0; i < SEL_CUTS.size(); i++) selStr = selStr+SEL_CUTS.at(i)+"_";
  selStr.pop_back();

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
	  std::cout << "  * Writing 1D histogram " << it->second->GetName() << std::endl;
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
	  std::cout << "  * Writing 2D histogram " << it->second->GetName() << std::endl;
	  it->second->Write();
	}
      }
  
      out_file->Write();
      std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
      
    } // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
  } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
  
  std::cout << "\nExiting ttH_3l()\n";
  
} // End void ttH_3l()
