
///////////////////////////////////////////////////
///     Re-train Linear Discriminant with       ///
///    non-mass BDT + dimuon mass as input      ///
///                                             ///
///       Andrew Brinkerhoff  15.05.2019        ///
///////////////////////////////////////////////////

// Basic ROOT includes to read and write files
#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

#include "H2MuAnalyzer/TrainMVA/interface/TMVA_helper.h" // Tools for TMVA


// Hard-coded options for running locally / manually
// Options passed in as arguments to TMVA_retrain_WH_lep when running in batch mode
const int MAX_EVT  = -1;    // Maximum number of events to process per sample
const int PRT_EVT  = 10000; // Print every N events

const bool verbose = false; // Print extra information

// const TString IN_DIR  = "/afs/cern.ch/work/x/xzuo/public/H2Mu/2017/Histograms/VH_selection_2019april/pt10_iso04/WH_ele_with_BDT/plots/";
const TString IN_DIR  = "/afs/cern.ch/user/a/abrinke1/HiggsToMuMu/2017/CMSSW_9_4_10/src/H2MuAnalyzer/TrainMVA/output/";
const TString OUT_DIR = "output";

// Prescales for data and MC: select 1/Xth of the events in each sample
const int SIG_PRESC  = 1;
const int BKG_PRESC  = 1;
const int DAT_PRESC  = 1;

const double PI = 3.14159265359;
const double BIT = 0.000001; // Tiny value or offset

const double LUMI = 41500; // pb-1

using namespace TMVA;


// Command-line options for running in batch.  Running "root -b -l -q macros/TMVA_retrain_WH_lep.C" will use hard-coded options above.
void TMVA_retrain_WH_lep( TString myMethodList = "", TString in_dir = "",
			   TString out_dir = "", TString out_file_str = "",
			   int max_evt = 0, int prt_evt = 0) {
  
  // Set variables to hard-coded values if they are not initialized
  if (in_dir.Length()  == 0) in_dir  = IN_DIR;
  if (out_dir.Length() == 0) out_dir = OUT_DIR;
  if (max_evt          == 0) max_evt = MAX_EVT;
  if (prt_evt          == 0) prt_evt = PRT_EVT;
  
  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // Linear Discriminant Analysis
  Use["LD"]             = 0;
  // Boosted Decision Trees                                                                                                                                                                                        
  Use["BDTG_UF_v1"]     = 1;
  Use["BDTG_AWB_lite"]  = 0;

  // ---------------------------------------------------------------
  std::cout << "\n==> Start TMVA_retrain_WH_lep" << std::endl;

  // Select methods (don't look at this code - not of interest)
  std::vector<TString> mlist;
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
    mlist = gTools().SplitString( myMethodList, ',' );
    for (int i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      
      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	std::cout << std::endl;
	return;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------
  
  // Here the preparation phase begins
  TString out_file_name;
  out_file_name.Form( "%s/TMVA_retrain_WH_lep_AWB_2019_05_20_v3.root", out_dir.Data() );
  TFile * out_file = TFile::Open( out_file_name, "RECREATE" );

  ///////////////////////////////////////////////////////
  ///  Input samples: MC signal, MC background, data  ///
  ///////////////////////////////////////////////////////

  // Tuple of sample name and sample ID
  std::vector< std::tuple<TString, int> > sig_samps;
  std::vector< std::tuple<TString, int> > bkg_samps;
  std::vector< std::tuple<TString, int> > dat_samps;
  std::vector< std::tuple<TString, int> > all_samps;

  sig_samps.push_back( std::make_tuple("Sig", 0) );
  bkg_samps.push_back( std::make_tuple("Bkg", 1) );
  
  all_samps.insert( all_samps.end(), sig_samps.begin(), sig_samps.end() );
  all_samps.insert( all_samps.end(), bkg_samps.begin(), bkg_samps.end() );
  all_samps.insert( all_samps.end(), dat_samps.begin(), dat_samps.end() );

  
  //////////////////////////////////////////////////////////////////
  ///  Factories: Use different sets of variables, weights, etc. ///
  //////////////////////////////////////////////////////////////////
  
  TString fact_set = "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
  //      fact_set = "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification";
          fact_set = "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification";

  std::vector<TString> var_names; // Holds names of variables for a given factory and permutation
  std::vector<double> var_vals; // Holds values of variables for a given factory and permutation
  TMVA::Factory* nullF = new TMVA::Factory("NULL", out_file, fact_set); // Placeholder factory
  TMVA::DataLoader* nullL = new TMVA::DataLoader("NULL");                 // Placeholder loader
  
  // Tuple is defined by the factory and dataloader,  followed by a name,
  // var name and value vectors, and hex bit masks for input variables.
  // Each hex bit represents four variables, e.g. 0x1 would select only the 1st variable,
  // 0xf the 1st 4, 0xff the 1st 8, 0xa the 2nd and 4th, 0xf1 the 1st and 5th-8th, etc.
  // The last three strings indicate which signal, background, and nJets to use.
  std::vector< std::tuple<TMVA::Factory*, TMVA::DataLoader*, TString, std::vector<TString>, std::vector<double>,
    int, int, int, TString, TString, TString> > factories;

  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_noMassBDT_mass_v3", var_names, var_vals,
  					0x3, 0x0, 0x0, "all", "all", "ge0j") );
  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_withMassBDT_nEles_v3", var_names, var_vals,
  					0x5, 0x0, 0x0, "all", "all", "ge0j") );


  // Initialize factories and dataloaders
  for (int iFact = 0; iFact < factories.size(); iFact++) {
    TString fact_name;
    fact_name.Form( "%s_%s_sig_%s_bkg_%s", std::get<2>(factories.at(iFact)).Data(), std::get<8>(factories.at(iFact)).Data(),
		    std::get<9>(factories.at(iFact)).Data(), std::get<10>(factories.at(iFact)).Data() );
    std::get<0>(factories.at(iFact)) = new TMVA::Factory( fact_name, out_file, fact_set );
    std::get<1>(factories.at(iFact)) = new TMVA::DataLoader( fact_name );
    std::get<2>(factories.at(iFact)) = fact_name;
  }

  std::cout << "Initialized factories" << std::endl;

  // Defined in interface/MVA_helper.h
  // TMVA_var(TString name, TString descr, TString unit, TString type, double def_val)
  std::vector<TMVA_var> in_vars;    // Input variables
  std::vector<TMVA_var> spec_vars;  // Spectator variables

  /////////////////////////////////////////////////////////
  ///  Input variables: used in BDT to estimate the pT  ///
  /////////////////////////////////////////////////////////

  // in_vars.push_back( TMVA_var( "BDTG_UF_v2", "BDT output",   "", 'F', -88 ) ); // 0x0001
  // in_vars.push_back( TMVA_var( "dimu_mass",  "mass(#mu#mu)", "", 'F', -88 ) ); // 0x0002
  // in_vars.push_back( TMVA_var( "mu1_lepMVA",  "LepMVA(#mu1)", "", 'F', -88 ) ); // 0x0004
  in_vars.push_back( TMVA_var( "BDTG_UF_v1",  "BDT output",     "", 'F', -88 ) ); // 0x0001
  in_vars.push_back( TMVA_var( "dimu_mass",   "mass(#mu#mu)",   "", 'F', -88 ) ); // 0x0002
  in_vars.push_back( TMVA_var( "nEles",       "# of electrons", "", 'F', -88 ) ); // 0x0004

  /////////////////////////////
  ///  Spectator variables  ///
  /////////////////////////////
  
  spec_vars.push_back( TMVA_var( "classID",    "0 for Sig, 1 for Bkg",    "", 'I', -77 ) );
  spec_vars.push_back( TMVA_var( "samp_ID",    "Sample ID",               "", 'I', -77 ) );
  spec_vars.push_back( TMVA_var( "event",      "Event number",            "", 'I', -77 ) );
  // spec_vars.push_back( TMVA_var( "event_wgt",  "Per-event scale factors", "", 'F', -77 ) );
  // spec_vars.push_back( TMVA_var( "samp_wgt",   "Xsec x lumi for sample",  "", 'F', -77 ) );
  // spec_vars.push_back( TMVA_var( "lepMVA_wgt", "lepMVA efficiency SF",    "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "total_wgt",  "Overall weight",          "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "res_wgt",    "Mass resolution weight",  "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "dimu_mass",  "mass(#mu#mu)",         "GeV", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "nEles",      "Number of electrons",     "", 'I', -77 ) );


  // Fill each factory with the correct set of variables
  for (int iFact = 0; iFact < factories.size(); iFact++) {
    std::cout << "\n*** Factory " << std::get<2>(factories.at(iFact)) << " variables ***" << std::endl;

    std::cout << "*** Input variables ***" << std::endl;
    for (int i = 0; i < in_vars.size(); i++) {
      if ( 0x1 & (std::get<5>(factories.at(iFact)) >> i) ) { // Hex bit mask for in_vars
	TMVA_var v = in_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Spectator variables ***" << std::endl;
    for (int i = 0; i < spec_vars.size(); i++) {
      TMVA_var v = spec_vars.at(i);
      std::cout << v.name << std::endl;
      std::get<1>(factories.at(iFact))->AddSpectator( v.name, v.descr, v.unit, v.type );
      std::get<3>(factories.at(iFact)).push_back( v.name );
      std::get<4>(factories.at(iFact)).push_back( v.def_val );
    }

  } // End loop: for (int iFact = 0; iFact < factories.size(); iFact++)


  std::cout << "\n******* About to loop over samples *******" << std::endl;
  int nEvt_tot = 0;
  int nEvt_sig = 0;
  int nEvt_bkg = 0;
  std::vector<int> nTrain_sig;
  std::vector<int> nTrain_bkg;
  std::vector<int> nTest_sig;
  std::vector<int> nTest_bkg;

  for (int iFact = 0; iFact < factories.size(); iFact++) {
    nTrain_sig.push_back(0);
    nTrain_bkg.push_back(0);
    nTest_sig.push_back(0);
    nTest_bkg.push_back(0);
  }

  // TChain * ch_train_noMass = new TChain("2017_WH_ele_against_inclu_trimvar_all_sig_all_bkg_ge0j/TrainTree");
  // TChain * ch_test_noMass  = new TChain("2017_WH_ele_against_inclu_trimvar_all_sig_all_bkg_ge0j/TestTree");
  // TChain * ch_train_mass   = new TChain("2017_WH_ele_against_inclu_trimvar_with_mass_all_sig_all_bkg_ge0j/TrainTree");
  // TChain * ch_test_mass    = new TChain("2017_WH_ele_against_inclu_trimvar_with_mass_all_sig_all_bkg_ge0j/TestTree");

  // TChain * ch_train_noMass = new TChain("f_Opt_AWB_noMass_v2_all_sig_all_bkg_ge0j/TrainTree");
  // TChain * ch_test_noMass  = new TChain("f_Opt_AWB_noMass_v2_all_sig_all_bkg_ge0j/TestTree");
  // TChain * ch_train_mass   = new TChain("f_Opt_AWB_withMass_v2_all_sig_all_bkg_ge0j/TrainTree");
  // TChain * ch_test_mass    = new TChain("f_Opt_AWB_withMass_v2_all_sig_all_bkg_ge0j/TestTree");

  TChain * ch_train_noMass = new TChain("f_Opt_AWB_noMass_v3_resWgt_all_sig_all_bkg_resWgt/TrainTree");
  TChain * ch_test_noMass  = new TChain("f_Opt_AWB_noMass_v3_resWgt_all_sig_all_bkg_resWgt/TestTree");
  TChain * ch_train_mass   = new TChain("f_Opt_AWB_withMass_v3_all_sig_all_bkg_/TrainTree");
  TChain * ch_test_mass    = new TChain("f_Opt_AWB_withMass_v3_all_sig_all_bkg_/TestTree");

  std::vector<TChain *> all_chains = {ch_train_noMass, ch_test_noMass, ch_train_mass, ch_test_mass};

  // TString in_file = in_dir+"2017_WH_ele_against_inclu_lepMVA04.root";
  TString in_file = in_dir+"TMVA_BDT_2017_WH_lep_all_vs_all_v2.root";
  TFile *file_tmp(0);
  std::cout << "\nTrying to open file " << in_file << std::endl;
  if ( !gSystem->AccessPathName( in_file ) )
    file_tmp = TFile::Open( in_file ); // Check if file exists
  if (!file_tmp) {
    std::cout << "ERROR: could not open input file " << in_file << std::endl;
    return;
  }
  ch_train_noMass->Add( in_file );
  ch_test_noMass ->Add( in_file );
  ch_train_mass  ->Add( in_file );
  ch_test_mass   ->Add( in_file );
  
  // Maps of branch names to branch addresses
  std::map<TString, std::string *> br_str;
  std::map<TString, float>         br_flt;
  std::map<TString, int>           br_int;

  // ************************* //
  //  Looping over all chains  //
  // ************************* //

  for (int iCh = 0; iCh < all_chains.size(); iCh++) {
    TChain * in_chain = all_chains.at(iCh);
    
    std::cout << "Loading all branches from tree:" << std::endl;
    TObjArray * br_list = in_chain->GetListOfBranches();
    
    // Map each branch to a pointer of the appropriate type, and set branch address
    for (int i = 0; i < br_list->GetEntries(); i++) {
      TString br_name  = br_list->At(i)->GetName();
      TString br_title = br_list->At(i)->GetTitle();
      std::cout << "\nLoading branch " << br_title << " with type: ";
      
      if ( br_title.EqualTo("sample") ) {
	// std::cout << "STRING" << std::endl;
	br_str.insert( std::pair<TString, std::string *>(br_name, new std::string) );
	in_chain->SetBranchAddress(br_name, &br_str.at(br_name));
      } else if ( br_title.EndsWith("/F") ) {
	// std::cout << "FLOAT" << std::endl;
	br_flt.insert( std::pair<TString, float>(br_name, -88.) );
	in_chain->SetBranchAddress(br_name, &br_flt.at(br_name));
      } else if ( br_title.EndsWith("/I") ) {
	// std::cout << "INT" << std::endl;
	br_int.insert( std::pair<TString, int>(br_name, -88) );
	in_chain->SetBranchAddress(br_name, &br_int.at(br_name));
      } else {
	std::cout << "\n\nMAJOR BUG!!! Branch " << br_title << "type not recognized!!!\n\n" << std::endl;
      }
    } // End loop: for (int i = 0; i < br_list->GetEntries(); i++)
  
    int nEvents = in_chain->GetEntries();
    int nEvt_pass = 0;
    std::cout << "\nLooping over the " << nEvents << " events in chain #" << iCh << std::endl;
    int every_N = 1;
    if (max_evt > 0 && max_evt < nEvents) {
      every_N = int(nEvents / max_evt);
    }
    std::cout << "  * In order to process only " << max_evt << ", will skip " << every_N-1 << "/" << every_N << " of events" << std::endl;

    for (int iEvt = 0; iEvt < nEvents; iEvt++) {

      // Skip events to only process max_evt total
      if (iEvt % every_N != 0)
	continue;
    
      if (iEvt % (PRT_EVT*every_N) == 0)
	std::cout << "Looking at event " << iEvt << " / " << nEvents << " (" << nEvt_pass << " passed so far, " << nEvt_tot << " in all samples)" << std::endl;
    
      in_chain->GetEntry(iEvt);

      // See if event is from a known sample, is
      // unprescaled, and passes selection cuts
      bool keep_event = false;

      TString samp_name;
      int     samp_ID = br_flt.at("samp_ID");
      int     classID;
      for (int iSamp = 0; iSamp < all_samps.size(); iSamp++) {
	TString iSamp_name = std::get<0>(all_samps.at(iSamp));
	int     iClassID   = std::get<1>(all_samps.at(iSamp));
	
	// Check sample ID
	if ( iClassID == br_int.at("classID") ) {
	  keep_event = true;
	  samp_name  = iSamp_name;
	  classID    = iClassID;
	} else continue;
      }
      if (not keep_event) continue;
      
      // // Select "odd" events after prescale
      // int presc = 1;
      // if (classID == 0) presc = SIG_PRESC * every_N;
      // if (classID == 1) presc = BKG_PRESC * every_N;
      // // if (samp_ID == 0) presc = DAT_PRESC * every_N;
      // if ( (iEvt % presc) != (presc - 1) ) {
      // 	continue;
      // }

      
      ///////////////////////////////////////////////
      ///  Weight signal events by cross section  ///
      ///////////////////////////////////////////////
      
      double event_wgt  = br_flt.at("event_wgt");
      double samp_wgt   = br_flt.at("samp_wgt");
      double lepMVA_wgt = br_flt.at("lepMVA_wgt");
      double total_wgt  = br_flt.at("total_wgt");
      double res_wgt    = (classID == 0 ? br_flt.at("res_wgt") : 1.0);
      
      if (verbose) std::cout << "Event " << iEvt << ", sample " << samp_name << " has sample weight " << samp_wgt
			     << ", event weight = " << event_wgt << ", and res_wgt = " << res_wgt
			     << " (total = " << samp_wgt*event_wgt << ")" << std::endl;


      /////////////////////////////////////////////////////
      ///  Loop over factories and set variable values  ///
      /////////////////////////////////////////////////////
      for (int iFact = 0; iFact < factories.size(); iFact++) {

	if ((iCh == 0 || iCh == 1) && iFact != 0) continue;  // First and second chain and first factory are no-mass
	if ((iCh == 2 || iCh == 3) && iFact != 1) continue;  // Third and fourth chain and second factory are with-mass
	
	TString sig_name = std::get<8>(factories.at(iFact));
	TString bkg_name = std::get<9>(factories.at(iFact));
	TString jet_cut = std::get<10>(factories.at(iFact));
      
	if (verbose) std::cout << "\n  * For factory " << iFact << ", sig_name = " << sig_name << ", bkg_name = " << bkg_name
			       << ", jet_cut = " << jet_cut << std::endl;
	
	if (verbose) std::cout << "  * Event passed the selection" << std::endl;
	
	// Set vars equal to default vector of variables for this factory
	var_names = std::get<3>(factories.at(iFact));
	var_vals = std::get<4>(factories.at(iFact));
	
	// Fill all variables
	for (int iVar = 0; iVar < var_names.size(); iVar++) {
	  TString vName = var_names.at(iVar);
	  
	  /////////////////////////////
	  ///  Spectator variables  ///
	  /////////////////////////////
	  
    	  if (verbose) std::cout << "Variable " << vName;
	  
	  if      ( vName == "samp_ID" )
	    var_vals.at(iVar) = samp_ID;
	  else if ( vName == "classID" )
	    var_vals.at(iVar) = classID;
	  else if ( vName == "samp_wgt" )
	    var_vals.at(iVar) = samp_wgt;
	  else if ( vName == "event_wgt" )
	    var_vals.at(iVar) = event_wgt;
	  else if ( vName == "lepMVA_wgt" )
	    var_vals.at(iVar) = lepMVA_wgt;
	  else if ( vName == "total_wgt" )
	    var_vals.at(iVar) = total_wgt;
	  
	  /////////////////////////
	  ///  Input variables  ///
	  /////////////////////////
	  
	  else if ( vName == "BDTG_UF_v1" )
	    var_vals.at(iVar) = br_flt.at("BDTG_UF_v1");
	  else if ( vName == "BDTG_UF_v2" )
	    var_vals.at(iVar) = br_flt.at("BDTG_UF_v2");
	  else if ( vName == "H_pair_mass" )
	    var_vals.at(iVar) = br_flt.at("H_pair_mass");
	  else if ( vName == "dimu_mass" )
	    var_vals.at(iVar) = br_flt.at("dimu_mass");
	  else if ( vName == "mu1_lepMVA" )
	    var_vals.at(iVar) = br_flt.at("mu1_lepMVA");
	  else if ( vName == "nEles" )
	    var_vals.at(iVar) = br_flt.at("nEles");

    	  if (verbose) std::cout << "filled with value " << var_vals.at(iVar) << std::endl;
	  
	} // End loop: for (int iVar = 0; iVar < var_names.size(); iVar++)
	

	// // Unweighted signal and background
	// double sig_evt_weight = 1.0;
	// double bkg_evt_weight = 1.0;
	
	// Weight by expected sample normalization
	double sig_evt_weight = total_wgt;
	double bkg_evt_weight = total_wgt;
	
	// Load values into event
	if (classID == 0) { // Signal MC
	  if ( (iCh % 2) == 0 ) { // Chains 0 and 2 are from the TrainTree
	    std::get<1>(factories.at(iFact))->AddSignalTrainingEvent( var_vals, sig_evt_weight );
	    nTrain_sig.at(iFact) += 1;
	  } else {                // Chains 1 and 3 are from the TestTree
	    std::get<1>(factories.at(iFact))->AddSignalTestEvent( var_vals, sig_evt_weight );
	    nTest_sig.at(iFact) += 1;
	  }
	}
	if (classID == 1) { // Background MC
	  if ( (iCh % 2) == 0 ) { // Chains 0 and 2 are from the TrainTree
	    std::get<1>(factories.at(iFact))->AddBackgroundTrainingEvent( var_vals, bkg_evt_weight );
	    nTrain_bkg.at(iFact) += 1;
	  } else {                // Chains 1 and 3 are from the TestTree
	    std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, bkg_evt_weight );
	    nTest_bkg.at(iFact) += 1;
	  }
	}
	// if (samp_ID == 0) { // Data
	// 	std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, bkg_evt_weight );
	// 	nTest_bkg.at(iFact) += 1;
	// }
	
      } // End loop: for (int iFact = 0; iFact < factories.size(); iFact++)
      
      nEvt_pass += 1;
      nEvt_tot  += 1;
      if (classID == 0) nEvt_sig += 1;
      if (classID == 1) nEvt_bkg += 1;
      
    } // End loop: for (int iEvt = 0; iEvt < nEvents; iEvt++)
  
    std::cout << "\n******* Made it out of the event loop *******\n\n" << std::endl;

  } // End loop: for (int iCh = 0; iCh < all_chains.size(); iCh++)
  
  std::cout << "\n******* Made it out of the loop over chains *******\n\n" << std::endl;
  
  // Run all the factories
  for (int iFact = 0; iFact < factories.size(); iFact++) {
    
    std::cout << "Factory " << iFact << " has " << nTrain_sig.at(iFact) << "sig / " << nTrain_bkg.at(iFact)
	      << " bkg training events, " << nTest_sig.at(iFact) << " sig / " << nTest_bkg.at(iFact) << " bkg testing" << std::endl;
    
    std::string NTrS;
    std::string NTrB;
    std::ostringstream convertTrS;
    convertTrS << nTrain_sig.at(iFact);
    NTrS = convertTrS.str();
    std::ostringstream convertTrB;
    convertTrB << nTrain_bkg.at(iFact);
    NTrB = convertTrB.str();
    
    std::string numTrainStr = "nTrain_Signal="+NTrS+":nTrain_Background="+NTrB+":";
    
    std::cout << "For factory " << iFact << ", loading " << numTrainStr << std::endl;
    
    // // global event weights per tree (see below for setting event-wise weights)
    // double regWeight  = 1.0;
    
    TMVA::Factory* factX = std::get<0>(factories.at(iFact));
    TMVA::DataLoader* loadX = std::get<1>(factories.at(iFact));
    
    // // You can add an arbitrary number of regression trees
    // loadX->AddRegressionTree( regTree, regWeight );
    
    // // This would set individual event weights (the variables defined in the
    // // expression need to exist in the original TTree)
    // loadX->SetWeightExpression( "var1", "Regression" );
    loadX->SetWeightExpression( 1.0 );
    
    // // Apply additional cuts on the signal and background samples (can be different)
    // TCut mycut = "( abs(muon.eta[0]) > 1.25 && abs(muon.eta[1]) < 2.4 )"; // && track.mode[0] == 15 )";
    
    // Tell the dataloader how to use the training and testing events
    loadX->PrepareTrainingAndTestTree( "", "", numTrainStr+"SplitMode=Random:NormMode=NumEvents:!V" );
    // loadX->PrepareTrainingAndTestTree( mycut, "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );
    
    // If no numbers of events are given, half of the events in the tree are used
    // for training, and the other half for testing:
    //
    //     loadX->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
    
    // Book MVA methods
    //
    // Please lookup the various method configuration options in the corresponding cxx files, eg:
    // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
    // it is possible to preset ranges in the option string in which the cut optimisation should be done:
    // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
    
    // Linear discriminant
    if (Use["LD"])
      factX->BookMethod( loadX,  TMVA::Types::kLD, "LD",
			 "!H:!V:VarTransform=None" );
    
    // Boosted Decision Trees

    if (Use["BDTG_UF_v1"]) // Optimized settings - AWB 04.04.2017
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v1", (std::string)
			 "!H:!V:NTrees=500::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5" );
    
    if (Use["BDTG_AWB_lite"]) // Fast, simple BDT
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_lite", (std::string)
			 "!H:!V:NTrees=50::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
    
    // --------------------------------------------------------------------------------------------------
    // Now you can tell the factory to train, test, and evaluate the MVAs
    
    // Train MVAs using the set of training events
    factX->TrainAllMethods();
    
    // Evaluate all MVAs using the set of test events
    factX->TestAllMethods();
    
    // Evaluate and compare performance of all configured MVAs
    factX->EvaluateAllMethods();
    // --------------------------------------------------------------
    
  } // End loop: for (int iFact = 0; iFact < factories.size(); iFact++)                                                                                                                                                                                          
  
  // Save the output
  out_file->Close();
  
  std::cout << "==> Wrote root file: " << out_file->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  
  // delete factory;
  // delete dataloader;
  
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( out_file_name );
}


int main( int argc, char** argv ) {
  // Select methods (don't look at this code - not of interest)
  TString methodList;
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(",");
    methodList += regMethod;
  }
  TMVA_retrain_WH_lep(methodList);
  return 0;
}
