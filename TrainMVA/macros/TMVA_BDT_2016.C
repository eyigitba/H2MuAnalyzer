
/////////////////////////////////////////////////////
///   General H2Mu vs. background BDT from 2016   ///
///        Copied from UFDimuAnalysis file        ///
///        bin/TMVAClassification_H2Mu.cxx        ///
///                                               ///
///        Andrew Brinkerhoff  01.10.2018         ///
/////////////////////////////////////////////////////

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
#include "H2MuAnalyzer/MakeHistos/interface/SampleDatabase2016.h" // Input data and MC samples

#include "H2MuAnalyzer/TrainMVA/interface/TMVA_helper.h" // Tools for TMVA


// Load the library of the local, compiled H2MuAnalyzer/MakeHistos and TrainMVA directories
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/TrainMVA/src/H2MuAnalyzerTrainMVA/libH2MuAnalyzerTrainMVA.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to TMVA_BDT_2016 when running in batch mode
const int MAX_EVT  =    -1;  // Maximum number of events to process per sample
const int PRT_EVT  = 10000;  // Print every N events
// const double SAMP_WGT = 1.0;
const bool verbose = false; // Print extra information

// const TString IN_DIR   = "";
// const TString SAMPLE   = "H2Mu_WH_pos";
const std::string YEAR = "2016";
const TString OUT_DIR  = "output";

const std::vector<std::string> SEL_CUTS = {"Presel2016"}; // Cuts which every event must pass
// const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
// const std::vector<std::string> CAT_CUTS = {"NONE"}; // Event selection categories, also applied in parallel


// Prescales for data and MC: select 1/Xth of the events in each sample
const int SIG_PRESC  = 1;
const int BKG_PRESC  = 1;
const int DAT_PRESC  = 1;

const bool MULTICLASS = false;

const double PI = 3.14159265359;
const double BIT = 0.000001; // Tiny value or offset

const double LUMI = 36814; // pb-1

using namespace TMVA;


// Command-line options for running in batch.  Running "root -b -l -q macros/TMVA_BDT_2016.C" will use hard-coded options above.
void TMVA_BDT_2016( TString myMethodList = "",
		    TString sample = "", TString in_dir = "", TString out_dir = "",
		    std::vector<TString> in_files = {}, TString out_file_str = "",
		    int max_evt = 0, int prt_evt = 0, double samp_weight = 1.0) {

  // Set variables to hard-coded values if they are not initialized
  // if (sample.Length()  == 0) sample  	 = SAMPLE;
  // if (in_dir.Length()  == 0) in_dir  	 = IN_DIR;
  if (out_dir.Length() == 0) out_dir 	 = OUT_DIR;
  if (max_evt          == 0) max_evt 	 = MAX_EVT;
  if (prt_evt          == 0) prt_evt 	 = PRT_EVT;
  // if (samp_weight      == 0) samp_weight = SAMP_WGT;
  
  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]          = 0;
  Use["PDEFoam"]        = 0;
  Use["KNN"]            = 0;
  // Linear Discriminant Analysis
  Use["LD"]             = 0;
  // Function Discriminant Analysis
  Use["FDA_GA"]         = 0;
  Use["FDA_MC"]         = 0;
  Use["FDA_MT"]         = 0;
  Use["FDA_GAMT"]       = 0;
  // Neural Network
  Use["MLP"]            = 0;
  Use["DNN"]            = 0;
  // Support Vector Machine
  Use["SVM"]            = 0;
  // Boosted Decision Trees                                                                                                                                                                                                                                         
  Use["BDT"]            = 0;
  Use["BDTG_default"]   = 0;
  Use["BDTG_UF_v1"]     = 1;
  Use["BDTG_AWB"]       = 0;
  Use["BDTG_AWB_lite"]  = 0;
  Use["BDTG_Carnes"]    = 0;

  // ---------------------------------------------------------------
  std::cout << "\n==> Start TMVA_BDT_2016" << std::endl;

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
  out_file_name.Form( "%s/TMVA_BDT_2016_AMC_2j.root", out_dir.Data() );
  TFile * out_file = TFile::Open( out_file_name, "RECREATE" );

  ///////////////////////////////////////////////////////
  ///  Input samples: MC signal, MC background, data  ///
  ///////////////////////////////////////////////////////

  std::map<TString, Sample*> samples;
  std::cout << "\nAbout to get the samples" << std::endl;

  // Load samples into our map. "UF" if you are on the UF servers
  // or "CERN" if you at CERN. "ALL" specifies that we want to load the Data
  // and all of the MC samples. Can loop through and remove the ones you don't want
  // to use if you desire or just grab the ones you care about from the map.
  GetSamples2016(samples, "CERN_hiM", "SIGNAL" );
  GetSamples2016(samples, "CERN_hiM", "ZJets_AMC");
  // GetSamples2016(samples, "CERN_hiM", "ZJets_MG");
  GetSamples2016(samples, "CERN_hiM", "tt_ll_MG");
  GetSamples2016(samples, "CERN_hiM", "singleTop");
  // GetSamples2016(samples, "CERN_hiM", "VV");
  // GetSamples2016(samples, "CERN_hiM", "ttX");
  // GetSamples2016(samples, "CERN_hiM", "DATA");
  
  std::cout << std::endl << "\nGot the samples" << std::endl;
  
  // Tuple of sample and sample ID
  std::vector< std::tuple<Sample*, int> > sig_samps;
  std::vector< std::tuple<Sample*, int> > bkg_samps;
  std::vector< std::tuple<Sample*, int> > dat_samps;
  std::vector< std::tuple<Sample*, int> > all_samps;

  sig_samps.push_back( std::make_tuple(samples["H2Mu_gg"],     -1) );
  sig_samps.push_back( std::make_tuple(samples["H2Mu_VBF"],    -2) );
  sig_samps.push_back( std::make_tuple(samples["H2Mu_ZH"],     -3) );
  sig_samps.push_back( std::make_tuple(samples["H2Mu_WH_pos"], -4) );
  sig_samps.push_back( std::make_tuple(samples["H2Mu_WH_neg"], -5) );
  // sig_samps.push_back( std::make_tuple(samples["H2Mu_ttH"],    -6) );
  
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_AMC"], + 1) );
  bkg_samps.push_back( std::make_tuple(samples["ZJets_AMC_0j"], + 1) );
  bkg_samps.push_back( std::make_tuple(samples["ZJets_AMC_1j"], + 2) );
  bkg_samps.push_back( std::make_tuple(samples["ZJets_AMC_2j"], + 3) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG"],              + 1) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_70_100"],    + 2) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_100_200"],   + 3) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_200_400"],   + 4) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_400_600"],   + 5) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_600_800"],   + 6) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_800_1200"],  + 7) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_1200_2500"], + 8) );
  // bkg_samps.push_back( std::make_tuple(samples["ZJets_MG_HT_2500_inf"],  + 9) );

  bkg_samps.push_back( std::make_tuple(samples["tt_ll_MG"], +10) );
  bkg_samps.push_back( std::make_tuple(samples["tW_pos"],   +11) );
  bkg_samps.push_back( std::make_tuple(samples["tW_neg"],   +12) );
  // bkg_samps.push_back( std::make_tuple(samples["tZq"],      +13) );
  
  // bkg_samps.push_back( std::make_tuple(samples["WW"],       +14) );
  // bkg_samps.push_back( std::make_tuple(samples["WZ_2l"],    +15) );
  // bkg_samps.push_back( std::make_tuple(samples["WZ_3l"],    +16) );
  // bkg_samps.push_back( std::make_tuple(samples["ZZ_2l_2q"], +17) );
  // bkg_samps.push_back( std::make_tuple(samples["ZZ_4l"],    +18) );
  
  // bkg_samps.push_back( std::make_tuple(samples["ttW"], +19) );
  // bkg_samps.push_back( std::make_tuple(samples["ttZ"], +20) );
  // // bkg_samps.push_back( std::make_tuple(samples["ttH"], +21) );
  
  // dat_samps.push_back( std::make_tuple(samples["RunB"],   0) );
  // dat_samps.push_back( std::make_tuple(samples["RunC"],   0) );
  // dat_samps.push_back( std::make_tuple(samples["RunD"],   0) );
  // dat_samps.push_back( std::make_tuple(samples["RunE"],   0) );
  // dat_samps.push_back( std::make_tuple(samples["RunF_1"], 0) );
  // dat_samps.push_back( std::make_tuple(samples["RunF_2"], 0) );
  // dat_samps.push_back( std::make_tuple(samples["RunG"],   0) );
  // dat_samps.push_back( std::make_tuple(samples["RunH"],   0) );
  // // dat_samps.push_back( std::make_tuple(samples["RunAll"], 0) );
  
  all_samps.insert( all_samps.end(), sig_samps.begin(), sig_samps.end() );
  all_samps.insert( all_samps.end(), bkg_samps.begin(), bkg_samps.end() );
  all_samps.insert( all_samps.end(), dat_samps.begin(), dat_samps.end() );

  // // Get branches, set addresses
  // // Tells the TTree that it should load the event information into samp->vars
  // for (int iSamp = 0; iSamp < all_samps.size(); iSamp++) {
  //   std::cout << "  * Setting branches for " << std::get<0>(all_samps.at(iSamp))->name << " (i = " << iSamp << ")" << std::endl;
  //   std::get<0>(all_samps.at(iSamp))->setBranchAddresses(2);
  //   std::get<0>(all_samps.at(iSamp))->calculateNoriginal();
  // }

  // std::cout << std::endl << "\nGot the branches, set addresses" << std::endl;
  
  //////////////////////////////////////////////////////////////////
  ///  Factories: Use different sets of variables, weights, etc. ///
  //////////////////////////////////////////////////////////////////
  
  TString         fact_set = "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
  fact_set                 = "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification";
  if (MULTICLASS) fact_set = "!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=multiclass";

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

  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_v1_AMC_2j", var_names, var_vals,
					0x0650, 0x033c, 0x001e, "all", "all", "ge0j") );

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
  std::vector<TMVA_var> mu_vars;    // Muon input variables
  std::vector<TMVA_var> jet_vars;   // Jet input variables
  std::vector<TMVA_var> evt_vars;   // Global event / combined object input variables
  std::vector<TMVA_var> spec_vars;  // All spectator variables

  /////////////////////////////////////////////////////////
  ///  Input variables: used in BDT to estimate the pT  ///
  /////////////////////////////////////////////////////////

  // Muon variables
  mu_vars.push_back( TMVA_var( "mu1_pt",        "p_{T}(#mu1)",       "GeV", 'F', -88 ) ); // 0x0001
  mu_vars.push_back( TMVA_var( "mu2_pt",        "p_{T}(#mu2)",       "GeV", 'F', -88 ) ); // 0x0002
  mu_vars.push_back( TMVA_var( "mu1_abs_eta",   "|#eta(#mu1)|",         "", 'F', -88 ) ); // 0x0004
  mu_vars.push_back( TMVA_var( "mu2_abs_eta",   "|#eta(#mu2)|",         "", 'F', -88 ) ); // 0x0008

  mu_vars.push_back( TMVA_var( "dimu_pt",       "p_{T}(#mu#mu)",     "GeV", 'F', -88 ) ); // 0x0010
  mu_vars.push_back( TMVA_var( "dimu_dMass",    "#sigma M(#mu#mu)",  "GeV", 'F', -88 ) ); // 0x0020
  mu_vars.push_back( TMVA_var( "dimu_eta",      "#eta(#mu#mu)",         "", 'F', -88 ) ); // 0x0040
  mu_vars.push_back( TMVA_var( "dimu_rapid",    "rapid(#mu#mu)",        "", 'F', -88 ) ); // 0x0080

  mu_vars.push_back( TMVA_var( "dimu_dR",       "dR(#mu#mu)",           "", 'F', -88 ) ); // 0x0100
  mu_vars.push_back( TMVA_var( "dimu_abs_dEta", "|d#eta(#mu#mu)|",      "", 'F', -88 ) ); // 0x0200
  mu_vars.push_back( TMVA_var( "dimu_abs_dPhi", "|d#phi(#mu#mu)|",      "", 'F', -88 ) ); // 0x0400
  mu_vars.push_back( TMVA_var( "dimu_dPhiStar", "d#phi*(#mu#mu)",       "", 'F', -88 ) ); // 0x0800

  // Jet variables
  jet_vars.push_back( TMVA_var( "jet1_pt",         "p_{T}(jet1)",       "GeV", 'F', -88 ) ); // 0x0001
  jet_vars.push_back( TMVA_var( "jet2_pt",         "p_{T}(jet2)",       "GeV", 'F', -88 ) ); // 0x0002
  jet_vars.push_back( TMVA_var( "jet1_eta",        "#eta(jet1)",           "", 'F', -88 ) ); // 0x0004
  jet_vars.push_back( TMVA_var( "jet2_eta",        "#eta(jet2)",           "", 'F', -88 ) ); // 0x0008

  jet_vars.push_back( TMVA_var( "dijet1_mass",     "1^{st} M(jj)",      "GeV", 'F', -88 ) ); // 0x0010
  jet_vars.push_back( TMVA_var( "dijet2_mass",     "2^{nd} M(jj)",      "GeV", 'F', -88 ) ); // 0x0020
  jet_vars.push_back( TMVA_var( "dijet3_mass",     "3^{rd} M(jj)",      "GeV", 'F', -88 ) ); // 0x0040
  jet_vars.push_back( TMVA_var( "dijet4_mass",     "4^{th} M(jj)",      "GeV", 'F', -88 ) ); // 0x0080

  jet_vars.push_back( TMVA_var( "dijet1_abs_dEta", "1^{st} |d#eta(jj)|",   "", 'F', -88 ) ); // 0x0100
  jet_vars.push_back( TMVA_var( "dijet2_abs_dEta", "2^{nd} |d#eta(jj)|",   "", 'F', -88 ) ); // 0x0200
  jet_vars.push_back( TMVA_var( "dijet3_abs_dEta", "3^{rd} |d#eta(jj)|",   "", 'F', -88 ) ); // 0x0400
  jet_vars.push_back( TMVA_var( "dijet4_abs_dEta", "4^{th} |d#eta(jj)|",   "", 'F', -88 ) ); // 0x0800

  // Global event variables
  evt_vars.push_back( TMVA_var( "nJets",       "# of jets",            "", 'I', -88 ) ); // 0x0001
  evt_vars.push_back( TMVA_var( "nJetsCent",   "# of central jets",    "", 'I', -88 ) ); // 0x0002
  evt_vars.push_back( TMVA_var( "nJetsFwd",    "# of forward jets",    "", 'I', -88 ) ); // 0x0004
  evt_vars.push_back( TMVA_var( "nBMed",       "# of medium b-tags",   "", 'I', -88 ) ); // 0x0008
  // evt_vars.push_back( TMVA_var( "nBLoose",     "# of loose b-tags",    "", 'I', -88 ) ); // 0x0008

  evt_vars.push_back( TMVA_var( "MET",         "MET",               "GeV", 'F', -88 ) ); // 0x0010
  evt_vars.push_back( TMVA_var( "MHT",         "MHT",               "GeV", 'F', -88 ) ); // 0x0020
  evt_vars.push_back( TMVA_var( "MT_had",      "M_{T} of jets",     "GeV", 'F', -88 ) ); // 0x0040
  evt_vars.push_back( TMVA_var( "mass_had",    "Mass of jets",      "GeV", 'F', -88 ) ); // 0x0080


  /////////////////////////////////////////////////////////////////////////////
  ///  Spectator variables: not used in training, but saved in output tree  ///
  /////////////////////////////////////////////////////////////////////////////

  spec_vars.push_back( TMVA_var( "samp_ID",   "Sample ID",            "", 'I', -77 ) );
  spec_vars.push_back( TMVA_var( "samp_wgt",  "Sample weight",        "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "res_wgt",   "Resolution weight",    "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "LHE_HT",    "Sample weight",     "GeV", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "dimu_mass", "mass(#mu#mu)",      "GeV", 'F', -77 ) );

  // Fill each factory with the correct set of variables
  for (int iFact = 0; iFact < factories.size(); iFact++) {
    std::cout << "\n*** Factory " << std::get<2>(factories.at(iFact)) << " variables ***" << std::endl;

    std::cout << "*** Input muon variables ***" << std::endl;
    for (int i = 0; i < mu_vars.size(); i++) {
      if ( 0x1 & (std::get<5>(factories.at(iFact)) >> i) ) { // Hex bit mask for mu_vars
	TMVA_var v = mu_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Input jet variables ***" << std::endl;
    for (int i = 0; i < jet_vars.size(); i++) {
      if ( 0x1 & (std::get<6>(factories.at(iFact)) >> i) ) { // Hex bit mask for jet_vars
	TMVA_var v = jet_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Input event variables ***" << std::endl;
    for (int i = 0; i < evt_vars.size(); i++) {
      if ( 0x1 & (std::get<7>(factories.at(iFact)) >> i) ) { // Hex bit mask for evt_vars
	TMVA_var v = evt_vars.at(i);
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


  // For inclusive ggH signal, Viktor gets the following triple gaussian (GluGlu__NoCats__125__ExpGaus__Separate__TripleGaus.png)
  // 0.73*Gaus(124.8, 1.52) + 0.23*Gaus(122.8, 4.24) + 0.04*Gaus(126, 2.1)
  TString tripGausExpr = "(0.73*TMath::Gaus(x, 124.8, 1.52, 1) + 0.23*TMath::Gaus(x, 122.8, 4.24, 1) + 0.04*TMath::Gaus(x, 126, 2.1, 1))";
  TF1* tripGaus   = new TF1("tripGaus", tripGausExpr, 113.8, 147.8);
  TF1* tripGausSq = new TF1("tripGaus", tripGausExpr+" * "+tripGausExpr, 113.8, 147.8);
  double tripGausNorm = tripGausSq->Integral(-1000, 1000);
  std::cout << "Triple gaussian has a normalization of " << tripGaus->Integral(-1000, 1000) << ", squared is " << tripGausNorm << std::endl;


  std::cout << "\n******* About to loop over samples *******" << std::endl;
  int nEvt_tot = 0;
  int nEvt_sig = 0;
  int nEvt_bkg = 0;
  std::vector<int> nTrain_sig;
  std::vector<int> nTrain_bkg;
  std::vector<int> nTest_sig;
  std::vector<int> nTest_bkg;

  std::vector<int> nTrain_ggH;
  std::vector<int> nTrain_VBF;
  std::vector<int> nTrain_VH;
  std::vector<int> nTrain_EWK;
  std::vector<int> nTrain_TOP;

  std::vector<int> nTest_ggH;
  std::vector<int> nTest_VBF;
  std::vector<int> nTest_VH;
  std::vector<int> nTest_EWK;
  std::vector<int> nTest_TOP;

  for (int iFact = 0; iFact < factories.size(); iFact++) {
    nTrain_sig.push_back(0);
    nTrain_bkg.push_back(0);
    nTest_sig.push_back(0);
    nTest_bkg.push_back(0);

    nTrain_ggH.push_back(0);
    nTrain_VBF.push_back(0);
    nTrain_VH.push_back(0);
    nTrain_EWK.push_back(0);
    nTrain_TOP.push_back(0);

    nTest_ggH.push_back(0);
    nTest_VBF.push_back(0);
    nTest_VH.push_back(0);
    nTest_EWK.push_back(0);
    nTest_TOP.push_back(0);
  }

  for (int iSamp = 0; iSamp < all_samps.size(); iSamp++) {
    Sample * samp = std::get<0>(all_samps.at(iSamp));
    int   samp_ID = std::get<1>(all_samps.at(iSamp));

    int prev_evt = -99;
    int this_evt = -88;

    int presc = 1;
    if (samp_ID  < 0) presc = SIG_PRESC;
    if (samp_ID  > 0) presc = BKG_PRESC;
    if (samp_ID == 0) presc = DAT_PRESC;

    /////////////////////////////////////////////////////////////
    ///  Block from macros in H2MuAnalyzer/MakeHistos/macros  ///
    /////////////////////////////////////////////////////////////

    // Initialize empty file to access each file in the list
    TFile *file_tmp(0);
    // Open all input files
    for (int i = 0; i < samp->filenames.size(); i++) {
      if ( !gSystem->AccessPathName(samp->filenames.at(i)) )
	file_tmp = TFile::Open( samp->filenames.at(i) ); // Check if file exists
      if (!file_tmp) {
	std::cout << "ERROR: could not open data file " << samp->filenames.at(i) << std::endl;
	return;
      }
    }

    // Initialize set of pointers to all branches in tree
    NTupleBranches br;
    // Add trees from the input files to the TChain
    TChain * in_chain = new TChain("dimuons/tree");
    for (UInt_t i = 0; i < samp->filenames.size(); i++) {
      in_chain->Add( samp->filenames.at(i) );
      // Set branch addresses, from H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h
      SetBranchAddresses(*in_chain, br, {YEAR, "Wgts"}, false); // Options in {} include "JES", "Flags", "GEN", "SFs", and "Wgts"
    }
    
    // Configuration for object selection, event selection, and object weighting
    ObjectSelectionConfig obj_sel;
    EventSelectionConfig  evt_sel;
    EventWeightConfig     evt_wgt;
    ConfigureObjectSelection(obj_sel, YEAR);
    ConfigureEventSelection (evt_sel, YEAR);
    ConfigureEventWeight    (evt_wgt, YEAR);
    obj_sel.Print();
    evt_sel.Print();
    evt_wgt.Print();

    /////////////////////////////////////////////////////////////////
    ///  End block from macros in H2MuAnalyzer/MakeHistos/macros  ///
    /////////////////////////////////////////////////////////////////


    int nEvents = in_chain->GetEntries();
    int nEvt_pass = 0;
    std::cout << "\nLooping over the " << nEvents << " events in sample " << samp->name << "\n" << std::endl;
    for (int iEvt = 0; iEvt < nEvents; iEvt++) {
      
      if (max_evt > 0 && iEvt > max_evt) break;
      
      if (iEvt % PRT_EVT == 0)
	std::cout << "Looking at event " << iEvt << " / " << nEvents << " (" << nEvt_pass << " passed so far, " << nEvt_tot << " in all samples)" << std::endl;
      
      in_chain->GetEntry(iEvt);

      
      // For 2016 NTuples, convert "SlimJets" collection into regular jets
      JetInfos jets_tmp;
      if (YEAR == "2016") {
	jets_tmp = ConvertSlimJets(*(br.slimJets));
	br.jets  = &jets_tmp;
      }


      ///////////////////////////////////////////
      ///  Apply selection and category cuts  ///
      ///////////////////////////////////////////
      
      // Throw away events that fail MC specific cuts
      if ( samp->name == "ZJets_MG" && br.LHE_HT > 70 )
	continue;
      if ( (samp->name == "ttH" || samp->name == "H2Mu_ttH") && br.genParents->size() == 0 ) {
	continue;
      }
      if ( samp->name == "ttH" && br.genParents->at(0).daughter_1_ID == 13 )
	continue;
      if ( samp->name == "H2Mu_ttH" && br.genParents->at(0).daughter_1_ID != 13 )
	continue;
      
      // Check if event passes basic selection cuts defined in src/SelectionCuts.cc
      bool pass_sel_cuts = true;
      for (int i = 0; i < SEL_CUTS.size(); i++) {
	if (not PassSelection(br, evt_sel, obj_sel, SEL_CUTS.at(i), verbose)) pass_sel_cuts = false;
      }
      if (not pass_sel_cuts) continue;
      
      // // Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc
      // for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
      //   if (not PassSelection(br, evt_sel, obj_sel, OPT_CUTS.at(iOpt), verbose)) continue;
      
      // // Loop through category cuts defined in src/CategoryCuts.cc
      // for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
      //   if (not InCategory(br, CAT_CUTS.at(iCat), verbose)) continue;
      
      // Get event weight for MC, defined in src/EventWeight.cc
      float event_wgt = ( samp_ID == 0 ? 1.0 : EventWeight(br, evt_wgt, verbose) );
      
      this_evt = br.event->event;
      if ( this_evt == prev_evt ) {
	std::cout << "\niEvt " << iEvt << " / " << nEvents << " (Run " << br.event->run << ", event " << this_evt
		  << ") is identical to the previous event (" << prev_evt << ")" << std::endl;
	std::cout << 100. - 100.0*iEvt/nEvents << " percent of sample " << samp->name << " lost" << std::endl;
	// continue;
	break;
	// return;
      }
      prev_evt = this_evt;
      
      if (samp->name == "RunF_1" && br.event->run > 278801)
	continue;
      if (samp->name == "RunF_2" && br.event->run < 278802)
	continue;
      
      // Discard half of signal events for use in limit-setting
      if (samp_ID < 0 && (this_evt % 2) == 0)
	continue;
      
      // Set the Higgs dimuon candidate pair
      MuPairInfo candPair = SelectedCandPair(obj_sel, br);
      float candMass = MuPairMass(candPair, obj_sel.mu_pt_corr);
      float candPt   = MuPairPt  (candPair, obj_sel.mu_pt_corr);
      
      // In data, ~the same number events in [113.8, 120], [120, 130], and [130, 147.8] GeV
      if ( candMass < 113.8 || candMass > 147.8 )
	continue;
      
      ///////////////////////////////////////////////
      ///  Set objects and variables of interest  ///
      ///////////////////////////////////////////////
      
      MuonInfo mu1        = br.muons->at(candPair.iMu1);
      MuonInfo mu2        = br.muons->at(candPair.iMu2);
      JetInfos jets       = SelectedJets(obj_sel, br);
      JetPairInfos dijets = SelectedJetPairs(obj_sel, br);
      int nJets           = NumJets(obj_sel, br);
      int nJetPairs       = dijets.size();
      int nBMed           = NumJets(obj_sel, br, "BTagMedium");
      int nBLoose         = NumJets(obj_sel, br, "BTagLoose");
      float MET           = br.met->pt;
      
      assert(nJets == jets.size());
      
      
      //////////////////////////////////////////////////////////////////////////////////
      ///  Weight signal events by cross section and inclusive H2Mu mass resolution  ///
      //////////////////////////////////////////////////////////////////////////////////
      
      double samp_wgt = 1.0;
      double lumi_SF  = samp->getLumiScaleFactor(LUMI);
      if (samp_ID != 0) // Half of signal / background MC events go into training, half into testing
	samp_wgt = 2.0 * event_wgt * lumi_SF;
      if (samp_ID != 0) // Reweight for prescale
	samp_wgt *= (1.0 * presc / DAT_PRESC);
      if (samp_ID < 0)
	samp_wgt *= 2.0; // Reserve even-numbered signal events for limit-setting
      
      double LHE_HT = br.LHE_HT;
      
      // // Scale for missing events in MC and data
      // // Related to duplicate events, i.e. "prev_evt" above? Or not?
      // // Seems that duplicate event problem does not appear in H2MuAnalyzer framework - AWB 03.10.2018
      // if (samp->name == "ZJets_MG_HT_100_200") samp_wgt *= (1.0 / (1.0 - 0.1217));
      // if (samp->name == "ZJets_MG_HT_200_400") samp_wgt *= (1.0 / (1.0 - 0.1008));
      // if (samp->name == "ZJets_MG_HT_400_600") samp_wgt *= (1.0 / (1.0 - 0.1651));
      // if (samp->name == "tt_ll_MG") samp_wgt *= (1.0 / (1.0 - 0.2164));
      // if (samp->name == "tW_pos")   samp_wgt *= (1.0 / (1.0 - 0.6218));
      // if (samp->name == "tW_neg")   samp_wgt *= (1.0 / (1.0 - 0.6254));
      // if (samp->name == "ttW")      samp_wgt *= (1.0 / (1.0 - 0.4098));
      // if (samp_ID == 0) samp_wgt *= (36814. / (36814. - 0.9751*9014.)); // 97.5% of RunH lost
      
      double res_wgt = tripGaus->Eval(candMass) / tripGausNorm;
      
      if (verbose) std::cout << "Event " << iEvt << ", sample " << samp->name << " has weight " << event_wgt
			     << ", lumi scale factor " << lumi_SF << " (total = " << samp_wgt << "), and LHE HT = " << LHE_HT << std::endl;
      
      
      //////////////////////////////
      ///  Fill event variables  ///
      //////////////////////////////
      
      if (verbose) std::cout << "Dimuon mass = " << candMass << ", pT = " << candPt << ", eta = " << candPair.eta << std::endl;
      if (verbose) std::cout << "Muon 1 pT = " << mu1.pt << ", eta = " << mu1.eta << std::endl;
      if (verbose) std::cout << "Muon 2 pT = " << mu2.pt << ", eta = " << mu2.eta << std::endl;
      
      /////////////////////////////////////////////////////
      ///  Loop over factories and set variable values  ///
      /////////////////////////////////////////////////////
      for (int iFact = 0; iFact < factories.size(); iFact++) {
	
	TString sig_name = std::get<8>(factories.at(iFact));
	TString bkg_name = std::get<9>(factories.at(iFact));
	TString jet_cut = std::get<10>(factories.at(iFact));
	
	if (verbose) std::cout << "\n  * For factory " << iFact << ", sig_name = " << sig_name << ", bkg_name = " << bkg_name
			       << ", jet_cut = " << jet_cut << ", nJets = " << nJets << std::endl;
	
	if (sig_name == "ggH" && samp_ID < 0 && samp->name != "H2Mu_gg") continue;
	if (sig_name == "VBF" && samp_ID < 0 && samp->name != "H2Mu_VBF") continue;
	if (sig_name ==  "VH" && samp_ID < 0 && samp->name != "H2Mu_ZH" &&
	    samp->name != "H2Mu_WH_neg" && samp->name != "H2Mu_WH_pos") continue;
	
	if (bkg_name == "ZJets" && samp_ID > 0 && (samp_ID <  1 || samp_ID >  9)) continue;
	if (bkg_name == "ttbar" && samp_ID > 0 && (samp_ID < 10 || samp_ID > 13)) continue;
	
	if (jet_cut.Contains("eq0j") && nJets != 0) continue;
	if (jet_cut.Contains("eq1j") && nJets != 1) continue;
	if (jet_cut.Contains("eq2j") && nJets != 2) continue;
	if (jet_cut.Contains("eq3j") && nJets != 3) continue;
	
	if (jet_cut.Contains("ge1j") && nJets < 1) continue;
	if (jet_cut.Contains("ge2j") && nJets < 2) continue;
	if (jet_cut.Contains("ge3j") && nJets < 3) continue;
	if (jet_cut.Contains("ge4j") && nJets < 4) continue;
	
	if (jet_cut.Contains("le1j") && nJets > 1) continue;
	
	if (jet_cut.Contains("eq0b") && nBMed != 0) continue;
	if (jet_cut.Contains("eq1b") && nBMed != 1) continue;
	if (jet_cut.Contains("ge1b") && nBMed  < 1) continue;
	if (jet_cut.Contains("ge2b") && nBMed  < 2) continue;
	
	if (jet_cut.Contains("met80") && MET*(1 - 0.05*nJets) > 80) continue;
	
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
	  
	  if      ( vName == "samp_ID" )
	    var_vals.at(iVar) = samp_ID;
	  else if ( vName == "samp_wgt" )
	    var_vals.at(iVar) = samp_wgt;
	  else if ( vName == "res_wgt" )
	    var_vals.at(iVar) = res_wgt;
	  else if ( vName == "LHE_HT" )
	    var_vals.at(iVar) = LHE_HT;
	  else if ( vName == "dimu_mass" )
	    var_vals.at(iVar) = candMass;
	  
	  /////////////////////////////////////////////////////
	  ///  Variables automatically set in lib/VarSet.h  ///
	  /////////////////////////////////////////////////////
	  
	  else { // In 2016, variables were automatically set in lib/VarSet.h
	    
	    // var_vals.at(iVar) = samp->vars.getValue(vName.Data());
	    
	    // Muon variables
	    if      (vName == "mu1_pt")      var_vals.at(iVar) = MuonPt(mu1, obj_sel.mu_pt_corr);
	    else if (vName == "mu2_pt")      var_vals.at(iVar) = MuonPt(mu2, obj_sel.mu_pt_corr);
	    else if (vName == "mu1_abs_eta") var_vals.at(iVar) = abs(mu1.eta);
	    else if (vName == "mu2_abs_eta") var_vals.at(iVar) = abs(mu2.eta);
	    
	    else if (vName == "dimu_pt")    var_vals.at(iVar) = MuPairPt(candPair, obj_sel.mu_pt_corr);
	    else if (vName == "dimu_dMass") var_vals.at(iVar) = MuPairMassErr(candPair, obj_sel.mu_pt_corr);
	    else if (vName == "dimu_eta")   var_vals.at(iVar) = candPair.eta;
	    else if (vName == "dimu_rapid") var_vals.at(iVar) = candPair.rapid;
	    
	    else if (vName == "dimu_dR")       var_vals.at(iVar) = candPair.dR;
	    else if (vName == "dimu_abs_dEta") var_vals.at(iVar) = abs(candPair.dEta);
	    else if (vName == "dimu_abs_dPhi") var_vals.at(iVar) = abs(candPair.dPhi);
	    else if (vName == "dimu_dPhiStar") var_vals.at(iVar) = abs(candPair.dPhiStar);  // TODO: check proper computation - AWB 2018.10.02
	    
	    // Jet variables
	    else if (vName == "jet1_pt")  var_vals.at(iVar) = (nJets > 0 ? jets.at(0).pt  : -10);
	    else if (vName == "jet2_pt")  var_vals.at(iVar) = (nJets > 1 ? jets.at(1).pt  : -10);
	    else if (vName == "jet1_eta") var_vals.at(iVar) = (nJets > 0 ? jets.at(0).eta : -10);
	    else if (vName == "jet2_eta") var_vals.at(iVar) = (nJets > 1 ? jets.at(1).eta : -10);
	    
	    else if (vName == "dijet1_mass") var_vals.at(iVar) = (nJetPairs > 0 ? dijets.at(0).mass : -10);
	    else if (vName == "dijet2_mass") var_vals.at(iVar) = (nJetPairs > 1 ? dijets.at(1).mass : -10);
	    else if (vName == "dijet3_mass") var_vals.at(iVar) = (nJetPairs > 2 ? dijets.at(2).mass : -10);
	    else if (vName == "dijet4_mass") var_vals.at(iVar) = (nJetPairs > 3 ? dijets.at(3).mass : -10);
	    
	    else if (vName == "dijet1_abs_dEta") var_vals.at(iVar) = (nJetPairs > 0 ? abs(dijets.at(0).dEta) : -10);
	    else if (vName == "dijet2_abs_dEta") var_vals.at(iVar) = (nJetPairs > 1 ? abs(dijets.at(1).dEta) : -10);
	    else if (vName == "dijet3_abs_dEta") var_vals.at(iVar) = (nJetPairs > 2 ? abs(dijets.at(2).dEta) : -10);
	    else if (vName == "dijet4_abs_dEta") var_vals.at(iVar) = (nJetPairs > 3 ? abs(dijets.at(3).dEta) : -10);
	    
	    // Global event variables
	    else if (vName == "nJets")     var_vals.at(iVar) = nJets;
	    else if (vName == "nJetsCent") var_vals.at(iVar) = NumJets(obj_sel, br, "Central");
	    else if (vName == "nJetsFwd")  var_vals.at(iVar) = NumJets(obj_sel, br, "Forward");
	    else if (vName == "nBMed")     var_vals.at(iVar) = nBMed;
	    // evt_else if (vName == "nBLoose") var_vals.at(iVar) = nBLoose;
	    
	    else if (vName == "MET")      var_vals.at(iVar) = MET;
	    else if (vName == "MHT")      var_vals.at(iVar) = br.mht->pt;
	    else if (vName == "MT_had")   var_vals.at(iVar) = br.mht->MT_had;
	    else if (vName == "mass_had") var_vals.at(iVar) = br.mht->mass_had;
	    
	    if (verbose) std::cout << "  * Filled variable " << vName << " with value " << var_vals.at(iVar) << std::endl;
	  }
	  
	} // End loop: for (int iVar = 0; iVar < var_names.size(); iVar++)
	
	TString multi_str;
	if      (samp_ID == -1)
	  multi_str = "ggH";
	else if (samp_ID == -2)
	  multi_str = "VBF";
	else if (samp_ID < -2 && samp_ID > -6)
	  multi_str = "VH";
	else if ( (samp_ID > 0 && samp_ID < 10) || (samp_ID > 13 && samp_ID < 19) )
	  multi_str = "EWK";
	else if ( (samp_ID > 9 && samp_ID < 14) || (samp_ID > 18) )
	  multi_str = "TOP";
	else if (samp_ID == 0)
	  multi_str = "DATA";
	else
	  multi_str = "OTHER";
	
	// // Unweighted signal and background
	// double sig_evt_weight = 1.0;
	// double bkg_evt_weight = 1.0;
	
	// Weight by expected sample normalization
	double sig_evt_weight = samp_wgt * 1000.;
	double bkg_evt_weight = samp_wgt;
	if (MULTICLASS) sig_evt_weight *= 0.001; // Don't weight signal for MultiClass
	
	// // Weight by expected sample normalization x signal resolution
	// double sig_evt_weight = samp_wgt * res_wgt * 1000.;
	// double bkg_evt_weight = samp_wgt;
	
	// Load values into event
	if (samp_ID < 0) { // Signal MC
	  if ( (iEvt % (2*presc)) == 0 ) {
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddSignalTrainingEvent( var_vals, sig_evt_weight );
	      nTrain_sig.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTrainingEvent( multi_str, var_vals, sig_evt_weight );
	      if (multi_str == "ggH") nTrain_ggH.at(iFact) += 1;
	      if (multi_str == "VBF") nTrain_VBF.at(iFact) += 1;
	      if (multi_str == "VH")  nTrain_VH.at(iFact)  += 1;
	    }
	  } else {
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddSignalTestEvent( var_vals, sig_evt_weight );
	      nTest_sig.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTestEvent( multi_str, var_vals, sig_evt_weight );
	      if (multi_str == "ggH") nTest_ggH.at(iFact) += 1;
	      if (multi_str == "VBF") nTest_VBF.at(iFact) += 1;
	      if (multi_str == "VH")  nTest_VH.at(iFact)  += 1;
	    }
	  }
	}
	if (samp_ID > 0) { // Background MC
	  if ( (iEvt % (2*presc)) == 0 ) {
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddBackgroundTrainingEvent( var_vals, bkg_evt_weight );
	      nTrain_bkg.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTrainingEvent( multi_str, var_vals, bkg_evt_weight );
	      if (multi_str == "EWK") nTrain_EWK.at(iFact) += 1;
	      if (multi_str == "TOP") nTrain_TOP.at(iFact) += 1;
	    }
	  } else {
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, bkg_evt_weight );
	      nTest_bkg.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTestEvent( multi_str, var_vals, bkg_evt_weight );
	      if (multi_str == "EWK") nTest_EWK.at(iFact) += 1;
	      if (multi_str == "TOP") nTest_TOP.at(iFact) += 1;
	    }
	  }
	}
	if (samp_ID == 0) { // Data
	  if (!MULTICLASS) {
	    std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, bkg_evt_weight );
	    nTest_bkg.at(iFact) += 1;
	  } else {
	    std::get<1>(factories.at(iFact))->AddTestEvent( multi_str, var_vals, bkg_evt_weight );
	  }
	}
	
      } // End loop: for (int iFact = 0; iFact < factories.size(); iFact++)
      
      nEvt_pass += 1;
      nEvt_tot  += 1;
      if (samp_ID < 0) nEvt_sig += 1;
      if (samp_ID > 0) nEvt_bkg += 1;
    } // End loop: for (int iEvt = 0; iEvt < samp->GetEntries(); iEvt++)
  } // End loop: for (int iSamp = 0; iSamp < all_samps.size(); iSamp++)

  
  std::cout << "\n******* Made it out of the event loop *******\n\n" << std::endl;
  
  
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
    
    // Settings for MultiClass
    std::string NTr_ggH;
    std::string NTr_VBF;
    std::string NTr_VH;
    std::string NTr_EWK;
    std::string NTr_TOP;
    
    std::ostringstream convertTr_ggH;
    convertTr_ggH << nTrain_ggH.at(iFact);
    NTr_ggH = convertTr_ggH.str();
    
    std::ostringstream convertTr_VBF;
    convertTr_VBF << nTrain_VBF.at(iFact);
    NTr_VBF = convertTr_VBF.str();
    
    std::ostringstream convertTr_VH;
    convertTr_VH << nTrain_VH.at(iFact);
    NTr_VH = convertTr_VH.str();
    
    std::ostringstream convertTr_EWK;
    convertTr_EWK << nTrain_EWK.at(iFact);
    NTr_EWK = convertTr_EWK.str();
    
    std::ostringstream convertTr_TOP;
    convertTr_TOP << nTrain_TOP.at(iFact);
    NTr_TOP = convertTr_TOP.str();
    
    std::string numTrainStr = "nTrain_Signal="+NTrS+":nTrain_Background="+NTrB+":";
    if (MULTICLASS) {
      numTrainStr = "";
      if (nTrain_ggH.at(iFact) > 0) numTrainStr += "nTrain_ggH="+NTr_ggH+":";
      if (nTrain_VBF.at(iFact) > 0) numTrainStr += "nTrain_VBF="+NTr_VBF+":";
      if (nTrain_VH.at(iFact)  > 0) numTrainStr += "nTrain_VH=" +NTr_VH +":";
      if (nTrain_EWK.at(iFact) > 0) numTrainStr += "nTrain_EWK="+NTr_EWK+":";
      if (nTrain_TOP.at(iFact) > 0) numTrainStr += "nTrain_TOP="+NTr_TOP+":";
    }
    
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
    
    // Neural network (MLP)
    if (Use["MLP"])
      factX->BookMethod( loadX,  TMVA::Types::kMLP, "MLP", (std::string)
			 "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:"+
			 "TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:"+
			 "ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );
    
    if (Use["DNN"])
      {
	// TString layoutString ("Layout=TANH|(N+100)*2,LINEAR");
	// TString layoutString ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|20,LINEAR");
	// TString layoutString ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");
	// TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
	// TString layoutString ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
	// TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
	// TString layoutString ("Layout=TANH|100,TANH|30,LINEAR");
	
	TString layoutString ("Layout=TANH|100,LINEAR");
	
	TString training0 ( (std::string) "LearningRate=1e-5,Momentum=0.5,Repetitions=1,"+
			    "ConvergenceSteps=500,BatchSize=50,TestRepetitions=7,WeightDecay=0.01,"+
			    "Regularization=NONE,DropConfig=0.5+0.5+0.5+0.5,DropRepetitions=2");
	TString training1 ( (std::string) "LearningRate=1e-5,Momentum=0.9,Repetitions=1,"+
			    "ConvergenceSteps=170,BatchSize=30,TestRepetitions=7,WeightDecay=0.01,"+
			    "Regularization=L2,DropConfig=0.1+0.1+0.1,DropRepetitions=1");
	TString training2 ( (std::string) "LearningRate=1e-5,Momentum=0.3,Repetitions=1,ConvergenceSteps=150,"+
			    "BatchSize=40,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE");
	TString training3 ( (std::string) "LearningRate=1e-6,Momentum=0.1,Repetitions=1,ConvergenceSteps=500,"+
			    "BatchSize=100,TestRepetitions=7,WeightDecay=0.0001,Regularization=NONE");
	
	TString trainingStrategyString ("TrainingStrategy=");
	trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;
	
	
	// TString trainingStrategyString ( (std::string) "TrainingStrategy=LearningRate=1e-1,Momentum=0.3,"+
	//                               "Repetitions=3,ConvergenceSteps=20,BatchSize=30,TestRepetitions=7,"+
	//                               "WeightDecay=0.0,L1=false,DropFraction=0.0,DropRepetitions=5");
	
	TString nnOptions ("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM");
	// TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
	nnOptions.Append (":"); nnOptions.Append (layoutString);
	nnOptions.Append (":"); nnOptions.Append (trainingStrategyString);
	
	factX->BookMethod(loadX, TMVA::Types::kDNN, "DNN", nnOptions ); // NN
      }
    
    // Support Vector Machine
    if (Use["SVM"])
      factX->BookMethod( loadX,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
    
    // Boosted Decision Trees
    if (Use["BDT"])
      factX->BookMethod( loadX,  TMVA::Types::kBDT, "BDT", (std::string)
			 "!H:!V:NTrees=100:MinNodeSize=1.0%:BoostType=AdaBoostR2:"+
			 "nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );
    
    // Default TMVA settings
    if (Use["BDTG_default"])
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_default", (std::string)
			 "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
    
    if (Use["BDTG_UF_v1"]) // Optimized settings - AWB 04.04.2017
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v1", (std::string)
			 "!H:!V:NTrees=500::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5" );
    
    if (Use["BDTG_AWB"]) // Optimized settings from EMTF pT assignment
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB", (std::string)
			 "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001" );
    
    if (Use["BDTG_AWB_lite"]) // Fast, simple BDT
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_lite", (std::string)
			 "!H:!V:NTrees=40::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.000001" );
    
    // Factory settings from Andrew Carnes ... what do they do? - AWB 04.01.2017
    if (Use["BDTG_Carnes"])
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_Carnes", (std::string)
			 "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:"+
			 "NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning" );
    
    
    // --------------------------------------------------------------------------------------------------
    
    // Now you can tell the factory to train, test, and evaluate the MVAs
    
    // Train MVAs using the set of training events
    
    factX->TrainAllMethods();
    
    // Evaluate all MVAs using the set of test events
    factX->TestAllMethods();
    
    // Evaluate and compare performance of all configured MVAs
    factX->EvaluateAllMethods();
    
    // // Instead of "EvaluateAllMethods()", just write out the training and testing trees
    // // Skip unnecessary evaluation histograms, which take time on large datasets
    // // Code gleaned from original "EvaluateAllMethods()" function in tmva/tmva/src/Factory.cxx - AWB 31.01.2017
    // if ( factX->fMethodsMap.empty() )
    //   std::cout << "factX->fMethodsMap is empty" << std::endl;
    
    // std::map<TString, std::vector<IMethod*>*>::iterator itrMap;
    // for (itrMap = factX->fMethodsMap.begin(); itrMap != factX->fMethodsMap.end(); itrMap++) {
    
    //   std::vector<IMethod*> *methods = itrMap->second;
    //   std::list<TString> datasets;
    //   int nmeth_used[2] = {int(mlist.size()), 1};
    
    //   for (int k = 0; k < 2; k++) {
    //          for (int i = 0; i < nmeth_used[k]; i++) {
    //            MethodBase* theMethod = dynamic_cast<MethodBase*>((*methods)[i]);
    //            if (theMethod == 0) {
    //              std::cout << "For k = " << k << ", i = " << i << ", no valid method" << std::endl;
    //              continue;
    //            }
    //            TDirectory* RootBaseDir = (TDirectory*) out_file;
    //            RootBaseDir->cd( std::get<2>(factories.at(iFact)) );
    //              theMethod->Data()->GetTree(Types::kTesting)->Write( "", TObject::kOverwrite );
    //          } // End loop: for (int i = 0; i < nmeth_used[k]; i++)
    //   } // End loop: for (int k = 0; k < 2; k++)
    // } // End loop: for (itrMap = factX->fMethodsMap.begin(); itrMap != factX->fMethodsMap.end(); itrMap++)
    
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
  TMVA_BDT_2016(methodList);
  return 0;
}
