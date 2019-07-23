
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

#include "H2MuAnalyzer/TrainMVA/interface/TMVA_helper.h" // Tools for TMVA


// Hard-coded options for running locally / manually
// Options passed in as arguments to TMVA_BDT_2017_ggH_hiPt when running in batch mode
const int MAX_EVT =  -1;     // Maximum number of events to process per sample
const int PRT_EVT =  10000;  // Print every N events

const bool verbose = false; // Print extra information

const std::vector<TString> IN_DIRS  = {"/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/WH_lep_AWB_2019_06_25_v1/files/",
				       "/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms/WH_lep_AWB_2019_06_25_v1/files/"};
// const TString OUT_DIR = "output";
const TString OUT_DIR = "/afs/cern.ch/work/a/abrinke1/public/H2Mu/Run2/TrainMVA/output";

// Prescales for data and MC: select 1/Xth of the events in each sample
const int SIG_PRESC  = 1;
const int BKG_PRESC  = 1;
const int DAT_PRESC  = 1;

const bool MULTICLASS = false;

const double PI = 3.14159265359;
const double BIT = 0.000001; // Tiny value or offset

// Events must fall into one of the allowed categories
const std::vector<TString> CATS = {"OPT_3lep_CAT_looseLepMVA_noZ5_noBtag"};


using namespace TMVA;


// Command-line options for running in batch.  Running "root -b -l -q macros/Run2_WH_lep.C" will use hard-coded options above.
void Run2_WH_lep( TString myMethodList = "", std::vector<TString> in_dirs = {},
		  TString out_dir = "", TString out_file_str = "",
		  int max_evt = 0, int prt_evt = 0) {
  
  // Set variables to hard-coded values if they are not initialized
  if (in_dirs.size()   == 0) in_dirs = IN_DIRS;
  if (out_dir.Length() == 0) out_dir = OUT_DIR;
  if (max_evt          == 0) max_evt = MAX_EVT;
  if (prt_evt          == 0) prt_evt = PRT_EVT;
  
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
  Use["BDTG_UF_v1a"]    = 0;
  Use["BDTG_UF_v1b"]    = 0;
  Use["BDTG_UF_v1c"]    = 0;
  Use["BDTG_UF_v2"]     = 1;
  Use["BDTG_UF_v3"]     = 1;
  Use["BDTG_AWB"]       = 0;
  Use["BDTG_AWB_lite"]  = 0;
  Use["BDTG_Carnes"]    = 0;

  // ---------------------------------------------------------------
  std::cout << "\n==> Start Run2_WH_lep" << std::endl;

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
  out_file_name.Form( "%s/Run2_WH_lep_all_vs_all_2019_06_26_v1.root", out_dir.Data() );
  TFile * out_file = TFile::Open( out_file_name, "RECREATE" );

  ///////////////////////////////////////////////////////
  ///  Input samples: MC signal, MC background, data  ///
  ///////////////////////////////////////////////////////

  // Tuple of sample name and sample ID
  std::vector< std::tuple<TString, float> > sig_samps;
  std::vector< std::tuple<TString, float> > bkg_samps;
  std::vector< std::tuple<TString, float> > dat_samps;
  std::vector< std::tuple<TString, float> > all_samps;

  // sig_samps.push_back( std::make_tuple("H2Mu_gg_125_NLO",  -1) );
  // sig_samps.push_back( std::make_tuple("H2Mu_VBF",         -2) );
  sig_samps.push_back( std::make_tuple("H2Mu_ZH_125",      -3) );
  sig_samps.push_back( std::make_tuple("H2Mu_WH_pos_125",  -4) );
  sig_samps.push_back( std::make_tuple("H2Mu_WH_neg_125",  -5) );
  sig_samps.push_back( std::make_tuple("H2Mu_ttH_125",     -6) );
  
  // sig_samps.push_back( std::make_tuple("H2Mu_gg_120_NLO", -11) );
  // sig_samps.push_back( std::make_tuple("H2Mu_VBF_120",    -12) );
  sig_samps.push_back( std::make_tuple("H2Mu_ZH_120",     -13) );
  sig_samps.push_back( std::make_tuple("H2Mu_WH_pos_120", -14) );
  sig_samps.push_back( std::make_tuple("H2Mu_WH_neg_120", -15) );
  sig_samps.push_back( std::make_tuple("H2Mu_ttH_120",    -16) );
  
  // sig_samps.push_back( std::make_tuple("H2Mu_gg_130_NLO", -21) );
  // sig_samps.push_back( std::make_tuple("H2Mu_VBF_130",    -22) );
  sig_samps.push_back( std::make_tuple("H2Mu_ZH_130",     -23) );
  sig_samps.push_back( std::make_tuple("H2Mu_WH_pos_130", -24) );
  sig_samps.push_back( std::make_tuple("H2Mu_WH_neg_130", -25) );
  sig_samps.push_back( std::make_tuple("H2Mu_ttH_130",    -26) );

  // bkg_samps.push_back( std::make_tuple("ZJets_MG",      +1) );
  // bkg_samps.push_back( std::make_tuple("ZJets_AMC",     +2) );
  bkg_samps.push_back( std::make_tuple("ZJets_hiM_MG",  +1) );
  bkg_samps.push_back( std::make_tuple("ZJets_hiM_AMC", +2) );

  bkg_samps.push_back( std::make_tuple("tt_ll_POW", +4) );
  bkg_samps.push_back( std::make_tuple("tt_ll_MG",  +5) );
  bkg_samps.push_back( std::make_tuple("tW_pos",    +6) );
  bkg_samps.push_back( std::make_tuple("tW_neg",    +7) );
  bkg_samps.push_back( std::make_tuple("tZq",       +8) );
  bkg_samps.push_back( std::make_tuple("tZW",       +9) );
  
  bkg_samps.push_back( std::make_tuple("ttW",      +10) );
  bkg_samps.push_back( std::make_tuple("ttZ",      +11) );
  bkg_samps.push_back( std::make_tuple("ttZ_lowM", +12) );
  bkg_samps.push_back( std::make_tuple("ttH",      +13) );
  bkg_samps.push_back( std::make_tuple("tHq",      +14) );
  bkg_samps.push_back( std::make_tuple("tHW",      +15) );
  bkg_samps.push_back( std::make_tuple("ttWW",     +16) );

  bkg_samps.push_back( std::make_tuple("WZ_3l",          +20) );
  bkg_samps.push_back( std::make_tuple("ZZ_4l",          +21) );
  bkg_samps.push_back( std::make_tuple("ZZ_4l_gg_2e2mu", +22) );
  bkg_samps.push_back( std::make_tuple("ggZZ_2e2mu",     +22.2) );
  bkg_samps.push_back( std::make_tuple("ggZZ_2mu2tau",   +22.4) );
  bkg_samps.push_back( std::make_tuple("ZZ_4l_gg_4mu",   +23) );
  bkg_samps.push_back( std::make_tuple("ggZZ_4mu",       +23.2) );
  bkg_samps.push_back( std::make_tuple("ggZZ_4tau",      +23.4) );
  
  bkg_samps.push_back( std::make_tuple("WWW", +25) );
  bkg_samps.push_back( std::make_tuple("WWZ", +26) );
  bkg_samps.push_back( std::make_tuple("WZZ", +27) );
  bkg_samps.push_back( std::make_tuple("ZZZ", +28) );

  bkg_samps.push_back( std::make_tuple("H2W_ZH_125",       +30) );
  bkg_samps.push_back( std::make_tuple("H2W_WH_neg_125",   +31) );
  bkg_samps.push_back( std::make_tuple("H2W_WH_pos_125",   +32) );
  bkg_samps.push_back( std::make_tuple("H2Tau_WH_neg_125", +33) );
  bkg_samps.push_back( std::make_tuple("H2Tau_WH_pos_125", +34) );
  bkg_samps.push_back( std::make_tuple("H2Tau_ZH_125",     +35) );
  bkg_samps.push_back( std::make_tuple("H2Z_WH_neg_125",   +36) );
  bkg_samps.push_back( std::make_tuple("H2Z_WH_pos_125",   +37) );
  bkg_samps.push_back( std::make_tuple("H2Z_ZH_125",       +38) );
  
  // dat_samps.push_back( std::make_tuple("RunB",   0) );
  // dat_samps.push_back( std::make_tuple("RunC",   0) );
  // dat_samps.push_back( std::make_tuple("RunD",   0) );
  // dat_samps.push_back( std::make_tuple("RunE",   0) );
  // dat_samps.push_back( std::make_tuple("RunF_1", 0) );
  // dat_samps.push_back( std::make_tuple("RunF_2", 0) );
  // dat_samps.push_back( std::make_tuple("RunG",   0) );
  // dat_samps.push_back( std::make_tuple("RunH",   0) );
  // dat_samps.push_back( std::make_tuple("RunAll", 0) );
  
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
  // Default set of Transformations gives the following error: <GetSeparation> signal and background histograms have different or invalid dimensions
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

  // // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_XWZ_noMass_v2", var_names, var_vals,
  // // 					0xfffe, 0x0000, 0x0000, "all", "all", "ge0j") );
  // // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_XWZ_withMass_v2", var_names, var_vals,
  // // 					0xffff, 0x0000, 0x0000, "all", "all", "ge0j") );
  // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v2", var_names, var_vals,
  // 					0xfffe, 0xffff, 0x1, "all", "all", "ge0j") );
  // // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_withMass_v2", var_names, var_vals,
  // // 					0xffff, 0xffff, 0x1, "all", "all", "ge0j") );

  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_withMass_v3", var_names, var_vals,
  					0xffff, 0xffff, 0x1, "all", "all", "") );
  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v3", var_names, var_vals,
  					0xfffe, 0xffff, 0x1, "all", "all", "") );
  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v3_resWgt", var_names, var_vals,
  					0xfffe, 0xffff, 0x7, "all", "all", "resWgt") );

  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v3_WH_vs_WZ_trueZ", var_names, var_vals,
  					0xfffe, 0xffff, 0x1, "WH", "WZ_trueZ", "") );
  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v3_WH_vs_WZ_nonZ", var_names, var_vals,
  					0xfffe, 0xffff, 0x0, "WH", "WZ_nonZ", "") );
  factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v3_WH_vs_ZJets", var_names, var_vals,
  					0xfffe, 0xffff, 0x1, "WH", "ZJets", "") );

  // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v2_WH_vs_WZ_trueZ", var_names, var_vals,
  // 					0xfffe, 0xffff, 0x0, "WH_e2mu", "WZ_e2mu_trueZ", "ge0j") );
  // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v2_WH_vs_WZ_nonZ", var_names, var_vals,
  // 					0xfffe, 0xffff, 0x0, "WH_e2mu", "WZ_e2mu_nonZ", "ge0j") );
  // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v2_WH_vs_ZJets", var_names, var_vals,
  // 					0xfffe, 0xffff, 0x0, "WH_e2mu", "ZJets_e2mu", "ge0j") );

  // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v2_WH_vs_WZ_trueZ", var_names, var_vals,
  // 					0xfffe, 0xffff, 0x0, "WH_3mu", "WZ_3mu_trueZ", "ge0j") );
  // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v2_WH_vs_WZ_nonZ", var_names, var_vals,
  // 					0xfffe, 0xffff, 0x0, "WH_3mu", "WZ_3mu_nonZ", "ge0j") );
  // factories.push_back( std::make_tuple( nullF, nullL, "f_Opt_AWB_noMass_v2_WH_vs_ZJets", var_names, var_vals,
  // 					0xfffe, 0xffff, 0x0, "WH_3mu", "ZJets_3mu", "ge0j") );


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
  std::vector<TMVA_var> xwz_vars;   // Variables originally from Xunwu
  std::vector<TMVA_var> awb_vars;   // Variables introduced by Andrew
  std::vector<TMVA_var> opt_vars;   // Optional extras
  std::vector<TMVA_var> spec_vars;  // Spectator variables

  /////////////////////////////////////////////////////////
  ///  Input variables: used in BDT to estimate the pT  ///
  /////////////////////////////////////////////////////////

  // Variables originally from Xunwu
  xwz_vars.push_back( TMVA_var( "H_pair_mass",     "H_pair_mass",     "", 'F', -88 ) ); //0x0001
  xwz_vars.push_back( TMVA_var( "H_pair_pt",       "H_pair_pt",       "", 'F', -88 ) ); //0x0002
  xwz_vars.push_back( TMVA_var( "lep_pt",          "lep_pt",          "", 'F', -88 ) ); //0x0004
  xwz_vars.push_back( TMVA_var( "lep_H_pair_dEta", "lep_H_pair_dEta", "", 'F', -88 ) ); //0x0008

  xwz_vars.push_back( TMVA_var( "lep_H_pair_dR",      "lep_H_pair_dR",      "", 'F', -88 ) ); //0x0010
  xwz_vars.push_back( TMVA_var( "lep_muSS_cosThStar", "lep_muSS_cosThStar", "", 'F', -88 ) ); //0x0020
  xwz_vars.push_back( TMVA_var( "lep_muOS_cosThStar", "lep_muOS_cosThStar", "", 'F', -88 ) ); //0x0040
  xwz_vars.push_back( TMVA_var( "lep_muSS_dEta",      "lep_muSS_dEta",      "", 'F', -88 ) ); //0x0080

  xwz_vars.push_back( TMVA_var( "lep_muSS_dR",   "lep_muSS_dR",   "", 'F', -88 ) ); //0x0100
  xwz_vars.push_back( TMVA_var( "lep_muOS_dEta", "lep_muOS_dEta", "", 'F', -88 ) ); //0x0200
  xwz_vars.push_back( TMVA_var( "lep_muOS_dR",   "lep_muOS_dR",   "", 'F', -88 ) ); //0x0400
  xwz_vars.push_back( TMVA_var( "MET",           "MET",           "", 'F', -88 ) ); //0x0800

  xwz_vars.push_back( TMVA_var( "lep_MET_MT",       "lep_MET_MT",       "", 'F', -88 ) ); //0x1000
  xwz_vars.push_back( TMVA_var( "lep_MET_dPhi_abs", "lep_MET_dPhi_abs", "", 'F', -88 ) ); //0x2000
  xwz_vars.push_back( TMVA_var( "MHT",              "MHT",              "", 'F', -88 ) ); //0x4000
  xwz_vars.push_back( TMVA_var( "lep_MHT_MT",       "lep_MHT_MT",       "", 'F', -88 ) ); //0x8000

  // Variables added by Andrew
  awb_vars.push_back( TMVA_var( "OS_pair_mass",      "OS_pair_mass",      "", 'F', -88 ) ); //0x0001
  awb_vars.push_back( TMVA_var( "OS_pair_pt",        "OS_pair_pt",        "", 'F', -88 ) ); //0x0002
  awb_vars.push_back( TMVA_var( "muSS_pt",           "muSS_pt",           "", 'F', -88 ) ); //0x0004
  awb_vars.push_back( TMVA_var( "muSS_OS_pair_dEta", "muSS_OS_pair_dEta", "", 'F', -88 ) ); //0x0008

  awb_vars.push_back( TMVA_var( "muSS_OS_pair_dR",     "muSS_OS_pair_dR",     "", 'F', -88 ) ); //0x0010
  awb_vars.push_back( TMVA_var( "muSS_lep_cosThStar",  "muSS_lep_cosThStar",  "", 'F', -88 ) ); //0x0020
  awb_vars.push_back( TMVA_var( "muSS_muOS_cosThStar", "muSS_muOS_cosThStar", "", 'F', -88 ) ); //0x0040
  awb_vars.push_back( TMVA_var( "muSS_muOS_dEta",      "muSS_muOS_dEta",      "", 'F', -88 ) ); //0x0080

  awb_vars.push_back( TMVA_var( "muSS_muOS_dR",      "muSS_muOS_dR",      "", 'F', -88 ) ); //0x0100
  awb_vars.push_back( TMVA_var( "muSS_MET_MT",       "muSS_MET_MT",       "", 'F', -88 ) ); //0x0200
  awb_vars.push_back( TMVA_var( "muSS_MET_dPhi_abs", "muSS_MET_dPhi_abs", "", 'F', -88 ) ); //0x0400
  awb_vars.push_back( TMVA_var( "muSS_MHT_MT",       "muSS_MHT_MT",       "", 'F', -88 ) ); //0x0800

  awb_vars.push_back( TMVA_var( "lep_MHT_dPhi_abs",  "lep_MHT_dPhi_abs",  "", 'F', -88 ) ); //0x1000
  awb_vars.push_back( TMVA_var( "muSS_MHT_dPhi_abs", "muSS_MHT_dPhi_abs", "", 'F', -88 ) ); //0x2000
  awb_vars.push_back( TMVA_var( "MHT_MET_dPhi_abs",  "MHT_MET_dPhi_abs",  "", 'F', -88 ) ); //0x4000
  awb_vars.push_back( TMVA_var( "lep_charge",        "lep_charge",        "", 'I', -88 ) ); //0x8000

  // Optional extras
  opt_vars.push_back( TMVA_var( "nEles",          "Number of electrons",    "", 'I', -88 ) ); //0x0001
  opt_vars.push_back( TMVA_var( "muH1_eta_abs",   "|#eta| #mu1 from Higgs", "", 'F', -88 ) ); //0x0002
  opt_vars.push_back( TMVA_var( "muH2_eta_abs",   "|#eta| #mu2 from Higgs", "", 'F', -88 ) ); //0x0004
  // opt_vars.push_back( TMVA_var( "", "", "", 'F', -88 ) ); //0x0008


  /////////////////////////////////////////////////////////////////////////////
  ///  Spectator variables: not used in training, but saved in output tree  ///
  /////////////////////////////////////////////////////////////////////////////

  spec_vars.push_back( TMVA_var( "samp_ID",          "Sample ID",               "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "event",            "Event number",            "", 'I', -77 ) );
  spec_vars.push_back( TMVA_var( "event_wgt",        "Per-event scale factors", "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "samp_wgt",         "Xsec x lumi for sample",  "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "lepMVA_wgt",       "lepMVA efficiency SF",    "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "total_wgt",        "Overall weight",          "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "res_wgt",          "Mass resolution weight",  "", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "dimu_mass",        "mass(#mu#mu)",         "GeV", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "nEles",            "Number of electrons",     "", 'I', -77 ) );
  // spec_vars.push_back( TMVA_var( "BDT_XWZ_noMass",   "Xwunwu's no-mass BDT",    "", 'F', -77 ) );
  // spec_vars.push_back( TMVA_var( "BDT_XWZ_withMass", "Xunwu's with-mass BDT",   "", 'F', -77 ) );
  // spec_vars.push_back( TMVA_var( "BDT_AWB_noMass",   "Andrew's no-mass BDT",    "", 'F', -77 ) );
  // spec_vars.push_back( TMVA_var( "BDT_AWB_withMass", "Andrew's with-mass BDT",  "", 'F', -77 ) );
  // spec_vars.push_back( TMVA_var( "BDT_AWB_retrain",  "AWB no-mass BDT + mass"   "", 'F', -77 ) );

  // Fill each factory with the correct set of variables
  for (int iFact = 0; iFact < factories.size(); iFact++) {
    std::cout << "\n*** Factory " << std::get<2>(factories.at(iFact)) << " variables ***" << std::endl;

    std::cout << "*** Variables originally from Xunwu ***" << std::endl;
    for (int i = 0; i < xwz_vars.size(); i++) {
      if ( 0x1 & (std::get<5>(factories.at(iFact)) >> i) ) { // Hex bit mask for xwz_vars
	TMVA_var v = xwz_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Variables added by Andrew ***" << std::endl;
    for (int i = 0; i < awb_vars.size(); i++) {
      if ( 0x1 & (std::get<6>(factories.at(iFact)) >> i) ) { // Hex bit mask for awb_vars
	TMVA_var v = awb_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Optional extras ***" << std::endl;
    for (int i = 0; i < opt_vars.size(); i++) {
      if ( 0x1 & (std::get<7>(factories.at(iFact)) >> i) ) { // Hex bit mask for opt_vars
	TMVA_var v = opt_vars.at(i);
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


  TChain * in_chain = new TChain("tree");

  std::vector<TString> in_files = {"HADD/tuples.root"};

  for (int i = 0; i < in_dirs.size(); i++) {
    for (int j = 0; j < in_files.size(); j++) {
      TString in_file = in_dirs.at(i)+in_files.at(j);
      // Open input file
      TFile *file_tmp(0);
      std::cout << "\nTrying to open file " << in_file << std::endl;
      if ( !gSystem->AccessPathName( in_file ) )
	file_tmp = TFile::Open( in_file ); // Check if file exists
      if (!file_tmp) {
	std::cout << "ERROR: could not open input file " << in_file << std::endl;
	return;
      }
      in_chain->Add( in_file );
    }
  }

  // Maps of branch names to branch addresses
  std::map<TString, std::string *> br_str;
  std::map<TString, float>         br_flt;
  std::map<TString, int>           br_int;

  std::cout << "Loading all branches from tree:" << std::endl;
  TObjArray * br_list = in_chain->GetListOfBranches();

  // Map each branch to a pointer of the appropriate type, and set branch address
  for (int i = 0; i < br_list->GetEntries(); i++) {
    TString br_name  = br_list->At(i)->GetName();
    TString br_title = br_list->At(i)->GetTitle();
    // std::cout << "\nLoading branch " << br_title << " with type: ";

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
  }


  std::vector<TString> used_samps;
  std::vector<TString> skipped_samps;

  int nEvents = in_chain->GetEntries();
  int nEvt_pass = 0;
  std::cout << "\nLooping over the " << nEvents << " events in the input file" << std::endl;
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
    float   samp_ID;
    for (int iSamp = 0; iSamp < all_samps.size(); iSamp++) {
      TString iSamp_name = std::get<0>(all_samps.at(iSamp));
      float   iSamp_ID   = std::get<1>(all_samps.at(iSamp));
      
      // Check sample name
      if ( iSamp_name.EqualTo(br_str.at("sample")->c_str()) ) {
	keep_event = true;
	samp_name  = iSamp_name;
	samp_ID    = iSamp_ID;
      } else continue;
    }

    // Skip event if it does not contain a known sample
    if (not keep_event) {
      bool already_skipped = false;
      for (const auto & skipped_samp : skipped_samps) {
	if (br_str.at("sample")->c_str() == skipped_samp) {
	  already_skipped = true;
	}
      }
      if (not already_skipped) {
	std::cout << "\nSkipping sample " << br_str.at("sample")->c_str() << " because it is not found in the list of samples\n" << std::endl;
	skipped_samps.push_back(br_str.at("sample")->c_str());
      }
      continue;
    }

    bool already_used = false;
    for (const auto & used_samp : used_samps) {
      if (br_str.at("sample")->c_str() == used_samp) {
	already_used = true;
      }
    }
    if (not already_used) {
      std::cout << "\nUsing sample " << br_str.at("sample")->c_str() << " which is found in the list of samples\n" << std::endl;
      used_samps.push_back(br_str.at("sample")->c_str());
    }


    // Select "odd" events after prescale
    int presc = 1;
    int event = br_int.at("event");
    if (samp_ID  < 0) presc = SIG_PRESC;
    if (samp_ID  > 0) presc = BKG_PRESC;
    if (samp_ID == 0) presc = DAT_PRESC;
    if ( (event % presc) != (presc - 1) ) continue;

    // Apply "optional" and "category" selection cuts
    bool pass_cat_cut = false;
    double lepMVA_wgt = 1.0;  // Lepton MVA efficiency SF weight depends on category
    for (int i = 0; i < CATS.size(); i++) {
      if ( br_int.at(CATS.at(i)) == 1 ) {
	pass_cat_cut = true;
	lepMVA_wgt = br_flt.at("lepMVA_wgt_"+CATS.at(i));
	// If using looseLepMVA to model medLepMVA, re-weight non-prompt background by 0.5
	if ( CATS.at(i).Contains("looseLepMVA") && ( samp_name.Contains("ZJets") ||
						     samp_name.Contains("tt_ll") ||
						     samp_name.Contains("tW_") ) ) lepMVA_wgt *= 0.5;
      }
    }
    if (!pass_cat_cut) continue;

    // For signal, only use the "true" di-muon pair from the Higgs
    if (samp_ID < 0 && br_flt.at("H_pair_mass") != br_flt.at("H_mass_true")) continue;
      

    //////////////////////////////////////////////////////////////////////////////////
    ///  Weight signal events by cross section and inclusive H2Mu mass resolution  ///
    //////////////////////////////////////////////////////////////////////////////////

    double event_wgt  = br_flt.at("event_wgt");
    double samp_wgt   = 1.0;
    double total_wgt  = 1.0;

    if (samp_ID != 0) // Half of signal / background MC events go into training, half into testing
      samp_wgt = 2.0 * br_flt.at("samp_wgt");
    if (samp_ID != 0) // Reweight for prescale
      samp_wgt *= (1.0 * presc / DAT_PRESC);
    
    // Special weights for processes with multiple MC samples
    if (samp_name == "ZJets_hiM_MG")  samp_wgt *= 0.5;
    if (samp_name == "ZJets_hiM_AMC") samp_wgt *= 0.5;
    if (samp_name == "tt_ll_POW")     samp_wgt *= 0.7;
    if (samp_name == "tt_ll_MG")      samp_wgt *= 0.3;    

    // Total weight combines all 3 weights
    total_wgt = samp_wgt*event_wgt*lepMVA_wgt;

    if (verbose) std::cout << "Event " << iEvt << ", sample " << samp_name << " has sample weight " << samp_wgt
    			   << ", event weight = " << event_wgt << ", and lepMVA SF = " << lepMVA_wgt
			   << " (total = " << total_wgt << ")" << std::endl;

    // For signal samples with mass 120 or 130, shift observed mass by 5 GeV
    double mass_shift = 0.0;
    if (samp_ID < 0 && samp_name.Contains("120")) mass_shift =  5.0;
    if (samp_ID < 0 && samp_name.Contains("130")) mass_shift = -5.0;


    /////////////////////////////////////////////////////
    ///  Loop over factories and set variable values  ///
    /////////////////////////////////////////////////////
    for (int iFact = 0; iFact < factories.size(); iFact++) {
      
      TString sig_name = std::get<8>(factories.at(iFact));
      TString bkg_name = std::get<9>(factories.at(iFact));
      TString options  = std::get<10>(factories.at(iFact));
      
      if (verbose) std::cout << "\n  * For factory " << iFact << ", sig_name = " << sig_name
			     << ", bkg_name = " << bkg_name << ", options = " << options << std::endl;

      // Select only events matching sample name
      if (sig_name.Contains("WH")    && samp_ID < 0 && !samp_name.Contains("WH"))    continue;
      if (bkg_name.Contains("WZ")    && samp_ID > 0 && !samp_name.Contains("WZ"))    continue;
      if (bkg_name.Contains("ZJets") && samp_ID > 0 && !samp_name.Contains("ZJets")) continue;
      // Select only events matching lepton flavor
      if (sig_name.Contains("e2mu") && samp_ID < 0 && br_int.at("nEles") != 1) continue;
      if (bkg_name.Contains("e2mu") && samp_ID > 0 && br_int.at("nEles") != 1) continue;
      if (sig_name.Contains("3mu")  && samp_ID < 0 && br_int.at("nEles") != 0) continue;
      if (bkg_name.Contains("3mu")  && samp_ID > 0 && br_int.at("nEles") != 0) continue;
      // Select only events with Higgs candidate matched to true Z (or not matched to true Z)
      if (bkg_name.Contains("trueZ") && samp_ID > 0 && br_flt.at("H_pair_mass") != br_flt.at("Z_mass_true")) continue;
      if (bkg_name.Contains("nonZ")  && samp_ID > 0 && br_flt.at("H_pair_mass") == br_flt.at("Z_mass_true")) continue;
      // If a 3-muon event, mask events within 5 GeV of the Z boson mass
      if (br_int.at("nEles") == 0 && abs(br_flt.at("OS_pair_mass") - 91) < 5) continue;
      
      
      // Apply additional selections or weightings
      double res_wgt = 1.0;
      if (options.Contains("resWgt")) {
	// Weight by avg. (mass_err/mass) / this pair's (mass_err/mass)
	// Computation of mass_err in the NTuples is too low by a factor of 2
	res_wgt = (1.1/125) / (br_flt.at("H_pair_mass_err") / br_flt.at("H_pair_mass"));
      }

      // Set vars equal to default vector of variables for this factory
      var_names = std::get<3>(factories.at(iFact));
      var_vals  = std::get<4>(factories.at(iFact));
      
      // Fill all variables
      for (int iVar = 0; iVar < var_names.size(); iVar++) {
    	TString vName = var_names.at(iVar);

	if (verbose) std::cout << "  * About to fill variable " << vName << " ... ";

    	/////////////////////////////
    	///  Spectator variables  ///
    	/////////////////////////////
	
    	if      ( vName == "samp_ID" )
    	  var_vals.at(iVar) = samp_ID;
    	else if ( vName == "event" )
    	  var_vals.at(iVar) = br_int.at("event");
    	else if ( vName == "samp_wgt" )
    	  var_vals.at(iVar) = samp_wgt;
    	else if ( vName == "event_wgt" )
    	  var_vals.at(iVar) = event_wgt;
    	else if ( vName == "lepMVA_wgt" )
    	  var_vals.at(iVar) = lepMVA_wgt;
    	else if ( vName == "total_wgt" )
    	  var_vals.at(iVar) = total_wgt;
    	else if ( vName == "res_wgt" )
    	  var_vals.at(iVar) = res_wgt;
    	else if ( vName == "dimu_mass" )
    	  var_vals.at(iVar) = br_flt.at("H_pair_mass") + mass_shift;
        else if ( vName == "nEles" )
	  var_vals.at(iVar) = br_int.at("nEles");
	
    	/////////////////////////////////////////////
    	///  Other variables with unusual values  ///
    	/////////////////////////////////////////////

        else if ( vName == "lep_charge" )
	  var_vals.at(iVar) = br_int.at("lep_charge");
    	else if ( vName == "H_pair_mass" )
    	  var_vals.at(iVar) = br_flt.at("H_pair_mass") + mass_shift;
    	else if ( vName == "muH1_eta_abs" )
    	  var_vals.at(iVar) = abs(br_flt.at("muH1_eta"));
    	else if ( vName == "muH2_eta_abs" )
    	  var_vals.at(iVar) = abs(br_flt.at("muH2_eta"));

    	////////////////////////////////////////////////////////////////////////////
    	///  Default behavior is to just get floating point value out of branch  ///
    	////////////////////////////////////////////////////////////////////////////

	else if ( vName.BeginsWith("n") )
	  var_vals.at(iVar) = br_int.at(vName);
	else
	  var_vals.at(iVar) = br_flt.at(vName);


	if (verbose) std::cout << "filled with value " << var_vals.at(iVar) << std::endl;
	  
      } // End loop: for (int iVar = 0; iVar < var_names.size(); iVar++)
	
      TString multi_str;
      if      (samp_ID <= -1 && samp_ID > -2)
    	multi_str = "ggH";
      else if (samp_ID <= -2 && samp_ID > -3)
    	multi_str = "VBF";
      else if (samp_ID < -2 && samp_ID > -6)
    	multi_str = "VH";
      else if ( (samp_ID > 0 && samp_ID < 10) || (samp_ID >= 14 && samp_ID < 19) )
    	multi_str = "EWK";
      else if ( (samp_ID >= 10 && samp_ID < 14) || (samp_ID >= 19) )
    	multi_str = "TOP";
      else if (samp_ID == 0)
    	multi_str = "DATA";
      else
    	multi_str = "OTHER";
      
      // Weight by expected sample normalization, plus di-muon mass resolution (if option "resWgt" used)
      double sig_evt_weight = total_wgt * res_wgt;
      double bkg_evt_weight = total_wgt;
	
      // Load values into event
      if (samp_ID < 0) { // Signal MC
	if ( (event % (2*presc)) == (2*presc - 1) ) {  // Use odd event numbers for training
	// if ( (iEvt % (2*every_N)) == 0 ) {
    	  if (!MULTICLASS) {
    	    std::get<1>(factories.at(iFact))->AddSignalTrainingEvent( var_vals, sig_evt_weight );
    	    nTrain_sig.at(iFact) += 1;
    	  } else {
    	    std::get<1>(factories.at(iFact))->AddTrainingEvent( multi_str, var_vals, sig_evt_weight );
    	    if (multi_str == "ggH") nTrain_ggH.at(iFact) += 1;
    	    if (multi_str == "VBF") nTrain_VBF.at(iFact) += 1;
    	    if (multi_str == "VH")  nTrain_VH.at(iFact)  += 1;
    	  }
    	} else {  // Use even event numbers for testing
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
	if ( (event % (2*presc)) == (2*presc - 1) ) {  // Use odd event numbers for training
        // if ( (iEvt % (2*every_N)) == 0 ) {
    	  if (!MULTICLASS) {
    	    std::get<1>(factories.at(iFact))->AddBackgroundTrainingEvent( var_vals, bkg_evt_weight );
    	    nTrain_bkg.at(iFact) += 1;
    	  } else {
    	    std::get<1>(factories.at(iFact))->AddTrainingEvent( multi_str, var_vals, bkg_evt_weight );
    	    if (multi_str == "EWK") nTrain_EWK.at(iFact) += 1;
    	    if (multi_str == "TOP") nTrain_TOP.at(iFact) += 1;
    	  }
    	} else {  // Use even event numbers for training
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

  } // End loop: for (int iEvt = 0; iEvt < nEvents; iEvt++)
  
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
    
    if (Use["BDTG_UF_v1a"]) // Exploratory settings - AWB 15.05.2018
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v1a", (std::string)
			 "!H:!V:NTrees=200::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=5" );
    
    if (Use["BDTG_UF_v1b"]) // Exploratory settings - AWB 15.05.2018
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v1b", (std::string)
			 "!H:!V:NTrees=500::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=10:MaxDepth=5" );
    
    if (Use["BDTG_UF_v1c"]) // Exploratory settings - AWB 15.05.2018
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v1c", (std::string)
			 "!H:!V:NTrees=500::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
    
    if (Use["BDTG_UF_v2"]) // Optimized settings - XWZ 15.05.2018
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v2", (std::string)
			 "!H:!V:NTrees=200::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=10:MaxDepth=3" );
    
    if (Use["BDTG_UF_v3"]) // Exploratory settings - AWB 15.05.2018
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v3", (std::string)
			 "!H:!V:NTrees=1000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=40:MaxDepth=6" );
    
    if (Use["BDTG_AWB"]) // Optimized settings from EMTF pT assignment
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB", (std::string)
			 "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001" );
    
    if (Use["BDTG_AWB_lite"]) // Fast, simple BDT
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_lite", (std::string)
			 "!H:!V:NTrees=50::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
    
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
  Run2_WH_lep(methodList);
  return 0;
}
