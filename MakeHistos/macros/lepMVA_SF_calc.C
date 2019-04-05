
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
#include "TH3.h"

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
const int MAX_EVT  = 100000;    // Maximum number of events to process
const int PRT_EVT  = 1000;  // Print every N events
const float SAMP_WGT = 1.0;
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information

const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2018_12_13_LepMVA_2l_test_v1/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/181213_154326/0000";
const TString SAMPLE   = "ZJets_AMC";
//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2018_12_13_LepMVA_2l_test_v1/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/tt_ll_MG/181213_211517/0000";
//const TString SAMPLE   = "tt_ll_MG"; 

//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2018_12_13_LepMVA_2l_test_v1/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/tt_ll_POW/181213_154754/0000";
//const TString SAMPLE   = "tt_ll_POW";

//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_14_LepMVA_2l_test_v2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/tt_ll_MG/190114_180858/0000";
//const TString SAMPLE   = "tt_ll_MG";

//const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_14_LepMVA_2l_test_v2/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/tt_ljj_POW_1/190114_182014/0000";
//const TString SAMPLE = "tt_ljj_POW_1";


const std::string YEAR = "2017";
const std::string SLIM = "Slim";
const TString OUT_DIR  = "plots";

const std::vector<std::string> SEL_CUTS = {"NONE"}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel


// Command-line options for running in batch.  Running "root -b -l -q macros/ReadNTupleChain.C" will use hard-coded options above.
void lepMVA_SF_calc( TString sample = "", TString in_dir = "", TString out_dir = "",
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
//      in_file_name.Form("%s/NTuple_0.root", in_dir.Data());
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

  // Initialize histograms
  TH3D *mu_hist  = new TH3D( sample + "_muon", sample +"_muon", 50,0,500, 16,-2.4,2.4, 200,-1,1);
  TH3D *ele_hist = new TH3D( sample + "_ele", sample + "_ele", 50,0,500, 16,-2.4,2.4, 200,-1,1); 


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
  ConfigureObjectSelection(obj_sel, YEAR);
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

//  evt_sel.muPair_mass_min = 105; // Require at least one Higgs candidate pair
//  obj_sel.mu_pt_min       =  10; // Lower muon pT threshold for muons not from Higgs
//  obj_sel.muPair_Higgs = "sort_WH_3_mu_v1"; // Choose Higgs candidate based on MT(W muon, MET) 
  obj_sel.mu_pt_corr = "PF";

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


    ///////////////////////////////////////////
    ///  Apply selection and category cuts  ///
    ///////////////////////////////////////////

    // Check if event passes basic selection cuts defined in src/SelectionCuts.cc
    bool pass_sel_cuts = true;
    for (int i = 0; i < SEL_CUTS.size(); i++) {
      if (not PassSelection(br, evt_sel, obj_sel, SEL_CUTS.at(i), verbose)) pass_sel_cuts = false;
    }   // only selection in this macro is "NONE", pass_sel_cuts is always true

    if (pass_sel_cuts) {
      
      // Get event weight for MC, defined in src/EventWeight.cc
      float event_wgt = ( sample.Contains("SingleMu") ? 1.0 : EventWeight(br, evt_wgt, verbose) );


      /////////////////////////////////////////////////////////////////////////////////////////
      ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
      /////////////////////////////////////////////////////////////////////////////////////////

      for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
	std::string OPT_CUT = OPT_CUTS.at(iOpt);
	
	//////////////////////////////////////////////////////////////////
	/// Loop through category cuts defined in src/CategoryCuts.cc  ///
	//////////////////////////////////////////////////////////////////
	
	for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {

	  //////////////////////////////////////////////////////
	  ///  Compute variables relevent for category cuts  ///
	  //////////////////////////////////////////////////////

	  // Get selected muons and electrons in event
	  MuonInfos muons = *br.muons;
	  EleInfos  eles  = *br.eles;
	  bool baseline_pass = false;
	  // here is the real preselection for this study
 	  for (int imu=0 ; imu < muons.size(); imu++) {
	      MuonInfo muon = muons.at(imu);
	      if ( MuonTrig(muon, obj_sel.year) and muon.pt > 30 ) baseline_pass = true; 
 	  }
	  if (br.met->pt > 30) baseline_pass = false; // for non-prompt preselection, no need for this in prompt preselection
	  if ( !baseline_pass ) continue;   // only require trigger match and pt, no need for the trigger muon to be prompt

	  ///////////////////////////////////////////
	  ///  Apply the category selection cuts  ///
	  ///////////////////////////////////////////

	  std::string CAT_CUT = CAT_CUTS.at(iCat);
	  std::string h_pre = (std::string)sample + "_"+OPT_CUT+"_"+CAT_CUT+"_";

	  /////////////////////////////////
	  ///  Generate and fill plots  ///
	  /////////////////////////////////

	  // Plot kinematic histograms
	  if (muons.size() == 2) {
	    MuPairInfo & dimuon = br.muPairs->at(0);
	    if (dimuon.mass > 81 and dimuon.mass < 101 and dimuon.charge == 0) {    // prompt control region
	      for(int imu=0 ; imu < muons.size(); imu++) {
	        MuonInfo muon_tag   = muons.at(imu);
	        MuonInfo muon_probe = muons.at( (imu+1)%2 );
		if ( muon_probe.SIP_3D > 8 or muon_probe.segCompat < 0.3 or not MuonID(muon_probe, "medium")) continue;
		if (muon_tag.pt > 20 and muon_tag.relIso < 0.12 and MuonID(muon_tag, "tight") ) { // prompt control region
		    mu_hist->Fill( (muon_probe.pt<499)?muon_probe.pt:499, muon_probe.eta, muon_probe.lepMVA, event_wgt);
	            
		}
	      }
	    }
	  }
	
	  if ( muons.size() == 1 and eles.size() == 1 ) {
	      MuonInfo muon = muons.at(0);	      
              EleInfo ele = eles.at(0);
	      TLorentzVector mu_vec  = FourVec(muon, obj_sel.mu_pt_corr);
	      TLorentzVector ele_vec = FourVec(ele);
	      TLorentzVector mu_ele = mu_vec + ele_vec;
	      if ( ele.SIP_3D > 8 or muon.segCompat < 0.3 ) continue;
	      if ( muon.relIso < 0.12 and MuonID(muon, "tight") and muon.charge + ele.charge == 0 and br.nBMed == 2) { // prompt
		  ele_hist->Fill( (ele.pt<499)?ele.pt:499, ele.eta, ele.lepMVA, event_wgt);
                
	      }
	  }

	} // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
      } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
    } // End conditional: if (pass_sel_cuts)


  } // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;

  std::cout << "\n******* normalizing histos, weight =  " << samp_weight << " *******" << std::endl;
  

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

      mu_hist->Write();
      ele_hist->Write();
 
      out_file->Write();
      std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
      
    } // End loop: for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++)
  } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
  
  std::cout << "\nExiting lepMVA_SF_calc()\n";
  
} // End void WH_lep()
