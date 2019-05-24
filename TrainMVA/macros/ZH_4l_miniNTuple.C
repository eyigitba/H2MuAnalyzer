
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
// R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/TrainMVA/src/H2MuAnalyzerTrainMVA/libH2MuAnalyzerTrainMVA.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to WH_e2mu_AWB_v1 when running in batch mode
const int MAX_EVT  =    -1;  // Maximum number of events to process per sample
const int PRT_EVT  = 1000;  // Print every N events
// const double SAMP_WGT = 1.0;
const bool verbose = false; // Print extra information

// const TString IN_DIR   = "";
// const TString SAMPLE   = "H2Mu_WH_pos";
const TString OUT_DIR  = "output";

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


// Command-line options for running in batch.  Running "root -b -l -q macros/ZH_4l_miniNTuple.C" will use hard-coded options above.
void ZH_4l_miniNTuple( TString myMethodList = "",
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
  Use["BDTG_UF_v2"]     = 1;
  Use["BDTG_UF_v3"]     = 1;
  Use["BDTG_AWB"]       = 0;
  Use["BDTG_AWB_lite"]  = 0;
  Use["BDTG_Carnes"]    = 0;

  // ---------------------------------------------------------------
  std::cout << "\n==> Start ZH_4l_miniNTuple" << std::endl;

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
  out_file_name.Form( "%s/2017_ZH_lep_against_inclu_with_mass.root", out_dir.Data() );
  TFile * out_file = TFile::Open( out_file_name, "RECREATE" );

  ///////////////////////////////////////////////////////
  ///  Input samples: MC signal, MC background, data  ///
  ///////////////////////////////////////////////////////
  std::map<TString, TString> MN_name; // <short_name, full path to miniNtuple>
  MN_name["signal"] 	= "/afs/cern.ch/work/x/xzuo/public/H2Mu/2017/Histograms/VH_selection_2019april/pt10_iso04/ZH_ele_BDT/sum_e_m_training.root";
  MN_name["bkg"]	= "/afs/cern.ch/work/x/xzuo/public/H2Mu/2017/Histograms/VH_selection_2019april/pt10_iso04/ZH_ele_BDT/sum_e_m_training.root";
 
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
    int, int, int, int, int, int, TString, TString, TString> > factories;


//  factories.push_back( std::make_tuple( nullF, nullL, "MNT_2017_inclusive", var_names, var_vals,
//                                        0x0fef, 0xffff, 0xffff, 0xffff, 0xffff, 0x0017, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_vs_ZZ_4l_training_with_mass_v1", var_names, var_vals,
//                                        0x0050, 0x003e, 0x000d, 0x0001, 0x0000, 0x0000, "all", "all", "ge0j") );


//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_ZZ_all_without_mass", var_names, var_vals,
//                                        0x0fef, 0xffff, 0xffff, 0x0fff, 0xffff, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_ZZ_test_without_mass", var_names, var_vals,
//                                        0x0fe0, 0xfff0, 0xfffe, 0x0fff, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_ZZ_trim_without_mass", var_names, var_vals,
//                                        0x0d60, 0xfd70, 0xffc8, 0x00d1, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_ZZ_sel_angle_without_mass", var_names, var_vals,
//                                        0x0360, 0xf370, 0x00d8, 0x00d1, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_ZZ_skim_without_mass", var_names, var_vals,
//                                        0x0360, 0xf370, 0x0058, 0x0010, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_ZZ_min_without_mass", var_names, var_vals,
//                                        0x0220, 0xf330, 0x0048, 0x0010, 0x0000, 0x0000, "all", "all", "ge0j") );

//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_all_without_mass", var_names, var_vals,
//                                        0x0fef, 0xffff, 0xffff, 0x0fff, 0xffff, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_test_without_mass", var_names, var_vals,
//                                        0x0fe0, 0xfff0, 0xfffe, 0x0fff, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_trim_without_mass", var_names, var_vals,
//                                        0x0d60, 0xfd70, 0xffc8, 0x00d1, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_sel_angle_without_mass", var_names, var_vals,
//                                        0x0360, 0xf370, 0x00d8, 0x00d1, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_skim_without_mass", var_names, var_vals,
//                                        0x0360, 0xf370, 0x0058, 0x0010, 0x0000, 0x0000, "all", "all", "ge0j") );
//  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_min_without_mass", var_names, var_vals,
//                                        0x0220, 0xf330, 0x0048, 0x0010, 0x0000, 0x0000, "all", "all", "ge0j") );


  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_sel_angle_with_mass", var_names, var_vals,
                                        0x0370, 0xf370, 0x00d8, 0x00d1, 0x0000, 0x0000, "all", "all", "ge0j") );
  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_skim_with_mass", var_names, var_vals,
                                        0x0370, 0xf370, 0x0058, 0x0010, 0x0000, 0x0000, "all", "all", "ge0j") );
  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_min_with_mass", var_names, var_vals,
                                        0x0230, 0xf330, 0x0048, 0x0010, 0x0000, 0x0000, "all", "all", "ge0j") );
  factories.push_back( std::make_tuple( nullF, nullL, "ZH_2017_lep_against_inclu_BDTv3_and_mass", var_names, var_vals,
                                        0x1010, 0x1000, 0x0000, 0x0000, 0x0000, 0x0000, "all", "all", "ge0j") );

  // Initialize factories and dataloaders
  for (int iFact = 0; iFact < factories.size(); iFact++) {
    TString fact_name;
    fact_name.Form( "%s_%s_sig_%s_bkg_%s", std::get<2>(factories.at(iFact)).Data(), std::get<11>(factories.at(iFact)).Data(),
		    std::get<12>(factories.at(iFact)).Data(), std::get<13>(factories.at(iFact)).Data() );
    std::get<0>(factories.at(iFact)) = new TMVA::Factory( fact_name, out_file, fact_set );
    std::get<1>(factories.at(iFact)) = new TMVA::DataLoader( fact_name );
    std::get<2>(factories.at(iFact)) = fact_name;
  }

  std::cout << "Initialized factories" << std::endl;

  // Defined in interface/MVA_helper.h
  // TMVA_var(TString name, TString descr, TString unit, TString type, double def_val)
  std::vector<TMVA_var> mu_vars;    // Muon input variables
  std::vector<TMVA_var> lep_vars;   // Electron input variables
  std::vector<TMVA_var> pair_vars;   // electron_muon kinematic variables
  std::vector<TMVA_var> met_vars;   // met related variables
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

  mu_vars.push_back( TMVA_var( "dimu_mass", 	"M(#mu#mu)",	     "GeV", 'F', -88 ) ); // 0x0010
  mu_vars.push_back( TMVA_var( "dimu_pt",       "p_{T}(#mu#mu)",     "GeV", 'F', -88 ) ); // 0x0020
  mu_vars.push_back( TMVA_var( "dimu_abs_eta",  "|#eta(#mu#mu)|",       "", 'F', -88 ) ); // 0x0040
  mu_vars.push_back( TMVA_var( "dimu_abs_dEta", "|d#eta(#mu#mu)|",      "", 'F', -88 ) ); // 0x0080 

  mu_vars.push_back( TMVA_var( "dimu_abs_dPhi", "|d#phi(#mu#mu)|",      "", 'F', -88 ) ); // 0x0100  
  mu_vars.push_back( TMVA_var( "dimu_dR",       "dR(#mu#mu)",           "", 'F', -88 ) ); // 0x0200 
  mu_vars.push_back( TMVA_var( "cts_mu1",       "cos(#theta*)(#mu1#mu)", "", 'F', -88 ) ); // 0x0400
  mu_vars.push_back( TMVA_var( "cts_mu_pos",    "cos(#theta*)(#mu+#mu)", "", 'F', -88 ) ); // 0x0800

  mu_vars.push_back( TMVA_var( "BDT_noMass_v3", "BDT v3 score",         "", 'F', -88 ) ); // 0x1000
  mu_vars.push_back( TMVA_var( "dimu_mass_err", "#sigma M(#mu#mu)",   "GeV", 'F', -88 ) ); // 0x2000  //not in use

  // lep variables
  lep_vars.push_back( TMVA_var( "lep1_pt", 	"p_{T}(lep1)", 	     "GeV",  'F', -88 ) );// 0x0001
  lep_vars.push_back( TMVA_var( "lep2_pt", 	"p_{T}(lep2)", 	     "GeV",  'F', -88 ) );// 0x0002
  lep_vars.push_back( TMVA_var( "lep1_abs_eta", "|#eta(lep1)|", 	"",  'F', -88 ) );// 0x0004
  lep_vars.push_back( TMVA_var( "lep2_abs_eta", "|#eta(lep2)|", 	"",  'F', -88 ) );// 0x0008
  
  lep_vars.push_back( TMVA_var( "dilep_mass", 	     "M(ll)",		"GeV", 'F', -88) );// 0x0010
  lep_vars.push_back( TMVA_var( "dilep_pt",	     "P_{T}(ll)",	"GeV", 'F', -88) );// 0x0020
  lep_vars.push_back( TMVA_var( "dilep_abs_eta",     "|#eta(ll)|",	"",    'F', -88) );// 0x0040
  lep_vars.push_back( TMVA_var( "dilep_abs_dEta",    "|d#eta(ll)|",   	"",    'F', -88) );// 0x0080

  lep_vars.push_back( TMVA_var( "dilep_abs_dPhi",    "|d#phi(ll)|",   	"",    'F', -88) );// 0x0100
  lep_vars.push_back( TMVA_var( "dilep_dR",	     "dR(ll)",      	"",    'F', -88) );// 0x0200
  lep_vars.push_back( TMVA_var( "cts_lep1",	     "cos(#theta*)(lep1)", "", 'F', -88) );// 0x0400
  lep_vars.push_back( TMVA_var( "cts_lep_pos",	     "cos(#theta*)(lep+)", "", 'F', -88) );// 0x0800

  lep_vars.push_back( TMVA_var( "lep_ID",       "lep_ID",       	"",    'I', -77 ) );// 0x1000

  // dipair variables
  pair_vars.push_back( TMVA_var( "quadlep_mass",     "M((#mu#mu),(ll))", 	 "GeV", 'F', -88) );// 0x0001
  pair_vars.push_back( TMVA_var( "quadlep_pt",	     "P_{T}((#mu#mu),(ll))",     "GeV", 'F', -88) );// 0x0002
  pair_vars.push_back( TMVA_var( "quadlep_abs_eta",  "|#eta((#mu#mu),(ll))|",  	 "", 	'F', -88) );// 0x0004
  pair_vars.push_back( TMVA_var( "dipair_dEta_H",    "d#eta((#mu#mu),(ll))",     "",    'F', -88) );// 0x0008

  pair_vars.push_back( TMVA_var( "dipair_dPhi_H",    "d#phi((#mu#mu),(ll))", 	 "", 	'F', -88) );// 0x0010
  pair_vars.push_back( TMVA_var( "dipair_dR_H",      "dR((#mu#mu),(ll))", 	 "", 	'F', -88) );// 0x0020
  pair_vars.push_back( TMVA_var( "cts_dipair_H",     "cos(#theta*)((#mu#mu),(ll))", "", 'F', -88) );// 0x0040
  pair_vars.push_back( TMVA_var( "cs_costheta",      "cos(#theta_{CS})", 	 "", 	'F', -88) );// 0x0080

  pair_vars.push_back( TMVA_var( "cs_cosphi", 	     "cos(#phi_{CS})", 	 	 "", 	'F', -88) );// 0x0100
  pair_vars.push_back( TMVA_var( "cs_sinphi", 	     "sin(#phi_{CS})", 		 "", 	'F', -88) );// 0x0200
  pair_vars.push_back( TMVA_var( "cos_theta1", 	     "cos(#theta_{1})", 	 "", 	'F', -88) );// 0x0400
  pair_vars.push_back( TMVA_var( "cos_phiH", 	     "cos(#phi_{H})", 		 "", 	'F', -88) );// 0x0800

  pair_vars.push_back( TMVA_var( "cos_phi1", 	     "cos(#phi_{1})", 		 "", 	'F', -88) );// 0x1000

   // met variables
  met_vars.push_back( TMVA_var( "met_pt", 		"p_{T}(met)",	      "GeV", 'F', -88 ) ); // 0x0001
  met_vars.push_back( TMVA_var( "abs_dPhi_4lPt_met",	"|d#phi(4l,met)|",	 "", 'F', -88 ) ); // 0x0002
  met_vars.push_back( TMVA_var( "met_Rtrans_4lPt", 	"(met/P_{T}(4l))_{T}",	 "", 'F', -88 ) ); // 0x0004
  met_vars.push_back( TMVA_var( "met_Rlongi_4lPt", 	"(met/P_{T}(4l))_{L}",   "", 'F', -88 ) ); // 0x0008

  met_vars.push_back( TMVA_var( "mht_pt",		"p_{T}(mht)",	      "GeV", 'F', -88 ) ); // 0x0010
  met_vars.push_back( TMVA_var( "abs_dPhi_4lPt_mht", 	"|d#phi(4l,mht)|", 	 "", 'F', -88 ) ); // 0x0020
  met_vars.push_back( TMVA_var( "mht_Rtrans_4lPt", 	"(mht/P_{T}(4l))_{T}",   "", 'F', -88 ) ); // 0x0040
  met_vars.push_back( TMVA_var( "mht_Rlongi_4lPt", 	"(mht/P_{T}(4l))_{L}", 	 "", 'F', -88 ) ); // 0x0080

  met_vars.push_back( TMVA_var( "mht_mass",		"M(mht)",	      "GeV", 'F', -88 ) ); // 0x0100
  met_vars.push_back( TMVA_var( "abs_dPhi_met_mht",	"|d#phi(met,mht)|",	 "", 'F', -88 ) ); // 0x0200

  // Jet variables
  jet_vars.push_back( TMVA_var( "dijet_mass",	 "M(jj)",	     "GeV", 'F', -88 ) ); // 0x0001
  jet_vars.push_back( TMVA_var( "dijet_pt",	 "p_{T}(jj)",	     "GeV", 'F', -88 ) ); // 0x0002
  jet_vars.push_back( TMVA_var( "dijet_abs_eta", "|#eta(jj)|",		"", 'F', -88 ) ); // 0x0004
  jet_vars.push_back( TMVA_var( "dijet_abs_dEta","|d#eta(jj)|",		"", 'F', -88 ) ); // 0x0008

  jet_vars.push_back( TMVA_var( "dijet_abs_dPhi","|d#phi(jj)|",		"", 'F', -88 ) ); // 0x0010
  jet_vars.push_back( TMVA_var( "dijet_dR",	 "dR(jj)",		"", 'F', -88 ) ); // 0x0020 
  jet_vars.push_back( TMVA_var( "jet1_pt",       "p_{T}(jet1)",      "GeV", 'F', -88 ) ); // 0x0040
  jet_vars.push_back( TMVA_var( "jet1_abs_eta",  "|#eta(jet1)|",        "", 'F', -88 ) ); // 0x0080

  jet_vars.push_back( TMVA_var( "jet2_pt",      "p_{T}(jet2)",       "GeV", 'F', -88 ) ); // 0x0100
  jet_vars.push_back( TMVA_var( "jet2_abs_eta", "|#eta(jet2)|",         "", 'F', -88 ) ); // 0x0200
  jet_vars.push_back( TMVA_var( "jet0_pt",	"p_{T}(jet0)",	     "GeV", 'F', -88 ) ); // 0x0400
  jet_vars.push_back( TMVA_var( "jet0_abs_eta", "|#eta(jet0)|",		"", 'F', -88 ) ); // 0x0800

  // Global event variables
  evt_vars.push_back( TMVA_var( "nJets",       "# of jets",            "", 'I', -88 ) ); // 0x0001
  evt_vars.push_back( TMVA_var( "nCentJets",   "# of central jets",    "", 'I', -88 ) ); // 0x0002
  evt_vars.push_back( TMVA_var( "nFwdJets",    "# of forward jets",    "", 'I', -88 ) ); // 0x0004
  evt_vars.push_back( TMVA_var( "nBJets_Med",  "# of medium b-tags",   "", 'I', -88 ) ); // 0x0008 // not in use

  evt_vars.push_back( TMVA_var( "nBJets_Loose",	"# of loose b-tags",   "", 'I', -88 ) ); // 0x0010 
  evt_vars.push_back( TMVA_var( "nBJets_Tight", "# of tight b-tags",   "", 'I', -88 ) ); // 0x0020 // not in use
  evt_vars.push_back( TMVA_var( "nMuons",	"# muons",	       "", 'I', -88 ) ); // 0x0040 // constant, not in use 
  evt_vars.push_back( TMVA_var( "nEles",	"# electrons",	       "", 'I', -88 ) ); // 0x0080 // constant, not in use 

  /////////////////////////////////////////////////////////////////////////////
  ///  Spectator variables: not used in training, but saved in output tree  ///
  /////////////////////////////////////////////////////////////////////////////

  spec_vars.push_back( TMVA_var( "event_wgt",   "event weight",        	"", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "xsec_norm",   "xsec normal",   	"", 'F', -77 ) );
  spec_vars.push_back( TMVA_var( "Sample_ID",   "Sample ID", 		"", 'I', -77 ) );
  spec_vars.push_back( TMVA_var( "lep_ID",	"the Z decays to",	"", 'I', -77 ) );
//  spec_vars.push_back( TMVA_var( "res_wgt",   "Resolution weight",    "", 'F', -77 ) );
//  spec_vars.push_back( TMVA_var( "LHE_HT",    "Sample weight",     "GeV", 'F', -77 ) );
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

    std::cout << "*** Input lep variables ***" << std::endl;
    for (int i = 0; i < lep_vars.size(); i++) {
      if ( 0x1 & (std::get<6>(factories.at(iFact)) >> i) ) { // Hex bit mask for ele_vars
	TMVA_var v = lep_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Input pair variables ***" << std::endl;
    for (int i = 0; i < pair_vars.size(); i++) {
      if ( 0x1 & (std::get<7>(factories.at(iFact)) >> i) ) { // Hex bit mask for emu_vars
	TMVA_var v = pair_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Input met variables ***" << std::endl;
    for (int i = 0; i < met_vars.size(); i++) {
      if ( 0x1 & (std::get<8>(factories.at(iFact)) >> i) ) { // Hex bit mask for met_vars
	TMVA_var v = met_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Input jet variables ***" << std::endl;
    for (int i = 0; i < jet_vars.size(); i++) {
      if ( 0x1 & (std::get<9>(factories.at(iFact)) >> i) ) { // Hex bit mask for jet_vars
	TMVA_var v = jet_vars.at(i);
	std::cout << v.name << std::endl;
	std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader
	std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
      }
    }

    std::cout << "*** Input event variables ***" << std::endl;
    for (int i = 0; i < evt_vars.size(); i++) {
      if ( 0x1 & (std::get<10>(factories.at(iFact)) >> i) ) { // Hex bit mask for evt_vars
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

  std::vector<int> nTrain_nonVH;
  std::vector<int> nTrain_ZH;
  std::vector<int> nTrain_WH;
  std::vector<int> nTrain_WZ;
  std::vector<int> nTrain_ZZ;
  std::vector<int> nTrain_ttX;

  std::vector<int> nTest_nonVH;
  std::vector<int> nTest_ZH;
  std::vector<int> nTest_WH;
  std::vector<int> nTest_WZ;
  std::vector<int> nTest_ZZ;
  std::vector<int> nTest_ttX;

  for (int iFact = 0; iFact < factories.size(); iFact++) {
    nTrain_sig.push_back(0);
    nTrain_bkg.push_back(0);
    nTest_sig.push_back(0);
    nTest_bkg.push_back(0);

    nTrain_nonVH.push_back(0);
    nTrain_ZH.push_back(0);
    nTrain_WH.push_back(0);
    nTrain_WZ.push_back(0);
    nTrain_ZZ.push_back(0);
    nTrain_ttX.push_back(0);

    nTest_nonVH.push_back(0);
    nTest_ZH.push_back(0);
    nTest_WH.push_back(0);
    nTest_WZ.push_back(0);
    nTest_ZZ.push_back(0);
    nTest_ttX.push_back(0);
  }

  int presc = 1;

  for (auto& MN: MN_name) {

    /////////////////////////////////////////////////////////////
    ///  Block from macros in H2MuAnalyzer/MakeHistos/macros  ///
    /////////////////////////////////////////////////////////////
    TChain * in_chain = new TChain("tree");
    in_chain->Add(MN.second);
    // Initialize set of pointers to all branches in tree
    
    int           nMuons;
    int           nEles;
    int           nJets;
    int           nBJets_Loose;
    int           nBJets_Med;
    int           nBJets_Tight;
    int           nFwdJets;
    int           nCentJets;
    float 	  BDT_noMass_v3;

    int 	  dimu_gen_ID;
    float         dimu_mass;
    float	  dimu_mass_err;
    float         dimu_pt;
    float         dimu_abs_eta;
    float	  dimu_abs_dEta;
    float	  dimu_abs_dPhi;
    float	  dimu_dR;
    float         mu1_pt;
    float         mu1_abs_eta;
    float	  mu1_lepMVA;
    float         mu2_pt;
    float         mu2_abs_eta;
    float	  mu2_lepMVA;

    float         cts_mu1;
    float         cts_mu_pos;

    float         lep1_pt;
    float         lep1_abs_eta;
    float         lep1_lepMVA;
    float         lep2_pt;
    float         lep2_abs_eta;
    float         lep2_lepMVA;

    int 	  dilep_gen_ID;
    float  	  dilep_mass; 
    float         dilep_pt; 
    float         dilep_abs_eta;
    float         dilep_abs_dEta;
    float         dilep_abs_dPhi;
    float         dilep_dR;
    float         cts_lep1;	 
    float         cts_lep_pos;	 
                          
    float         quadlep_mass;
    float         quadlep_pt; 
    float	  quadlep_abs_eta;
    float	  dipair_dEta_H;
    float	  dipair_dPhi_H;
    float	  dipair_dR_H;
    float         cts_dipair_H;

    float	  cs_costheta;
    float	  cs_cosphi;
    float	  cs_sinphi;
    float	  cos_theta1;
    float	  cos_phiH;
    float	  cos_phi1;

    float         dijet_mass;
    float         dijet_pt;
    float         dijet_abs_eta;
    float         dijet_abs_dEta;
    float         dijet_abs_dPhi;
    float         dijet_dR;

    float         jet1_pt;
    float         jet1_abs_eta;
    float         jet2_pt;
    float         jet2_abs_eta;
    float         jet0_pt;
    float         jet0_abs_eta;

    float         met_pt;
    float	  abs_dPhi_4lPt_met;
    float	  met_Rtrans_4lPt;
    float	  met_Rlongi_4lPt;
    float         mht_pt;
    float         mht_mass;
    float	  abs_dPhi_4lPt_mht;
    float	  mht_Rtrans_4lPt;
    float	  mht_Rlongi_4lPt;
    float	  abs_dPhi_met_mht;

    float	  event_wgt;
    float 	  xsec_norm;
    int           lep_ID; 
    int           Sample_ID;

    in_chain->SetBranchAddress("nMuons",       	& nMuons     	);
    in_chain->SetBranchAddress("nEles",        	& nEles      	);
    in_chain->SetBranchAddress("nJets",        	& nJets      	);
    in_chain->SetBranchAddress("nBJets_Loose", 	& nBJets_Loose	);
    in_chain->SetBranchAddress("nBJets_Med",   	& nBJets_Med 	);
    in_chain->SetBranchAddress("nBJets_Tight", 	& nBJets_Tight	);
    in_chain->SetBranchAddress("nFwdJets",     	& nFwdJets   	);
    in_chain->SetBranchAddress("nCentJets",    	& nCentJets  	);
    in_chain->SetBranchAddress("BDT_noMass_v3", & BDT_noMass_v3 );

    in_chain->SetBranchAddress("dimu_gen_ID",   & dimu_gen_ID   );
    in_chain->SetBranchAddress("dimu_mass",	& dimu_mass	);
    in_chain->SetBranchAddress("dimu_mass_err",	& dimu_mass_err );
    in_chain->SetBranchAddress("dimu_pt",    	& dimu_pt	);
    in_chain->SetBranchAddress("dimu_abs_eta",  & dimu_abs_eta	);
    in_chain->SetBranchAddress("dimu_abs_dEta", & dimu_abs_dEta	);
    in_chain->SetBranchAddress("dimu_abs_dPhi", & dimu_abs_dPhi	);
    in_chain->SetBranchAddress("dimu_dR",    	& dimu_dR	);
    in_chain->SetBranchAddress("mu1_pt",     	& mu1_pt	);
    in_chain->SetBranchAddress("mu1_abs_eta",   & mu1_abs_eta	);
    in_chain->SetBranchAddress("mu1_lepMVA",	& mu1_lepMVA	);
    in_chain->SetBranchAddress("mu2_pt",     	& mu2_pt	);
    in_chain->SetBranchAddress("mu2_abs_eta",   & mu2_abs_eta	);
    in_chain->SetBranchAddress("mu2_lepMVA",	& mu2_lepMVA	);
                
    in_chain->SetBranchAddress("cts_mu1",   	& cts_mu1   	);
    in_chain->SetBranchAddress("cts_mu_pos",	& cts_mu_pos	);

    in_chain->SetBranchAddress("dilep_gen_ID",		& dilep_gen_ID	);
    in_chain->SetBranchAddress("dilep_mass",		& dilep_mass	);
    in_chain->SetBranchAddress("dilep_pt",		& dilep_pt	);
    in_chain->SetBranchAddress("dilep_abs_eta",		& dilep_abs_eta );
    in_chain->SetBranchAddress("dilep_abs_dEta",	& dilep_abs_dEta);
    in_chain->SetBranchAddress("dilep_abs_dPhi",	& dilep_abs_dPhi);
    in_chain->SetBranchAddress("dilep_dR",		& dilep_dR	);

    in_chain->SetBranchAddress("lep1_pt",		& lep1_pt	);
    in_chain->SetBranchAddress("lep1_abs_eta",		& lep1_abs_eta	);
    in_chain->SetBranchAddress("lep1_lepMVA",		& lep1_lepMVA	);
    in_chain->SetBranchAddress("lep2_pt",		& lep2_pt	);
    in_chain->SetBranchAddress("lep2_abs_eta",		& lep2_abs_eta	);
    in_chain->SetBranchAddress("lep2_lepMVA",		& lep2_lepMVA	);
    in_chain->SetBranchAddress("cts_lep1",		& cts_lep1	);
    in_chain->SetBranchAddress("cts_lep_pos",		& cts_lep_pos	);

    in_chain->SetBranchAddress("quadlep_mass",		& quadlep_mass	 );
    in_chain->SetBranchAddress("quadlep_pt",		& quadlep_pt	 );
    in_chain->SetBranchAddress("quadlep_abs_eta",	& quadlep_abs_eta);
    in_chain->SetBranchAddress("dipair_dEta_H",		& dipair_dEta_H	 );
    in_chain->SetBranchAddress("dipair_dPhi_H",		& dipair_dPhi_H	 );
    in_chain->SetBranchAddress("dipair_dR_H",		& dipair_dR_H	 );
    in_chain->SetBranchAddress("cts_dipair_H",		& cts_dipair_H	 );
 
    in_chain->SetBranchAddress("cs_costheta",		& cs_costheta	);
    in_chain->SetBranchAddress("cs_cosphi",		& cs_cosphi	);
    in_chain->SetBranchAddress("cs_sinphi",		& cs_sinphi	);
    in_chain->SetBranchAddress("cos_theta1",		& cos_theta1	);
    in_chain->SetBranchAddress("cos_phiH",		& cos_phiH	);
    in_chain->SetBranchAddress("cos_phi1",		& cos_phi1	);
 
    in_chain->SetBranchAddress("dijet_mass",		& dijet_mass	);
    in_chain->SetBranchAddress("dijet_pt",  		& dijet_pt  	);
    in_chain->SetBranchAddress("dijet_abs_eta", 	& dijet_abs_eta );
    in_chain->SetBranchAddress("dijet_abs_dEta",	& dijet_abs_dEta);
    in_chain->SetBranchAddress("dijet_abs_dPhi",    	& dijet_abs_dPhi);
    in_chain->SetBranchAddress("dijet_dR",      	& dijet_dR      );

    in_chain->SetBranchAddress("jet1_pt",   	& jet1_pt   	);
    in_chain->SetBranchAddress("jet1_abs_eta",  & jet1_abs_eta  );
    in_chain->SetBranchAddress("jet2_pt",   	& jet2_pt   	);
    in_chain->SetBranchAddress("jet2_abs_eta",  & jet2_abs_eta  );
    in_chain->SetBranchAddress("jet0_pt",   	& jet0_pt   	);
    in_chain->SetBranchAddress("jet0_abs_eta",  & jet0_abs_eta  );

    in_chain->SetBranchAddress("met_pt",		& met_pt   		);
    in_chain->SetBranchAddress("abs_dPhi_4lPt_met", 	& abs_dPhi_4lPt_met	);
    in_chain->SetBranchAddress("met_Rtrans_4lPt",	& met_Rtrans_4lPt	);
    in_chain->SetBranchAddress("met_Rlongi_4lPt",	& met_Rlongi_4lPt	);
    in_chain->SetBranchAddress("mht_pt",     		& mht_pt		);
    in_chain->SetBranchAddress("mht_mass",   		& mht_mass		);
    in_chain->SetBranchAddress("abs_dPhi_4lPt_mht",	& abs_dPhi_4lPt_mht	);
    in_chain->SetBranchAddress("mht_Rtrans_4lPt",	& mht_Rtrans_4lPt	);
    in_chain->SetBranchAddress("mht_Rlongi_4lPt",	& mht_Rlongi_4lPt	);
    in_chain->SetBranchAddress("abs_dPhi_met_mht",	& abs_dPhi_met_mht	);

    in_chain->SetBranchAddress("event_wgt",     & event_wgt     );
    in_chain->SetBranchAddress("xsec_norm",     & xsec_norm     );
    in_chain->SetBranchAddress("lep_ID", 	& lep_ID 	);
    in_chain->SetBranchAddress("Sample_ID",	& Sample_ID	);


    /////////////////////////////////////////////////////////////////
    ///  End block from macros in H2MuAnalyzer/MakeHistos/macros  ///
    /////////////////////////////////////////////////////////////////


    int nEvents = in_chain->GetEntries();
    int nEvt_pass = 0;
    std::cout << "\nLooping over the " << nEvents << "\nSample_ID  " << Sample_ID << "\n" << std::endl;
    for (int iEvt = 0; iEvt < nEvents; iEvt++) {
      
      if (max_evt > 0 && iEvt > max_evt) break;
      
      if (iEvt % PRT_EVT == 0)
	std::cout << "Looking at event " << iEvt << " / " << nEvents << " (" << nEvt_pass << " passed so far, " << nEvt_tot << " in all samples)" << std::endl;
      
      in_chain->GetEntry(iEvt);
      if (iEvt % 1000 == 0)  std::cout << "Sample_ID   " << Sample_ID << "\nevent_wgt  " << event_wgt << "\nBDT_noMass_v3" << BDT_noMass_v3 << std::endl; 

 
      /////////////////////////////////////////
      ///  Get information about the event  ///
      /////////////////////////////////////////

      // Get event weight for MC, defined in src/EventWeight.cc
      
     
     
      // Discard half of signal events for use in limit-setting
//      if (Sample_ID > 0 && (iEvt % 2) == 0)  
//	continue;

      
      //////////////////////////////////////////////////////////////////////////////////
      ///  Weight signal events by cross section and inclusive H2Mu mass resolution  ///
      //////////////////////////////////////////////////////////////////////////////////
      
//      double samp_wgt = 1.0;
//      double lumi_SF  = samp->getLumiScaleFactor(LUMI);
//      if (samp_ID != 0) // Half of signal / background MC events go into training, half into testing
//	samp_wgt = 2.0 * event_wgt * lumi_SF;
//      if (samp_ID != 0) // Reweight for prescale
//	samp_wgt *= (1.0 * presc / DAT_PRESC);
//      if (samp_ID < 0)
//	samp_wgt *= 2.0; // Reserve even-numbered signal events for limit-setting
//      
     
//      double res_wgt = tripGaus->Eval(H_pair_vec.M()) / tripGausNorm;
      
      if (verbose) std::cout << "Event " << iEvt << ", sample " << Sample_ID << " has weight " << event_wgt << std::endl;
//			     << ", lumi scale factor " << lumi_SF << " (total = " << samp_wgt << "), and LHE HT = " << LHE_HT << std::endl;
      
      //////////////////////////////
      ///  Fill event variables  ///
      //////////////////////////////
      
//      if (verbose) std::cout << "Dimuon mass = " << H_pair_vec.M() << ", pT = " << H_pair_vec.Pt() << ", eta = " << H_pair.eta << std::endl;
      if (verbose) std::cout << "Muon 1 pT = " << mu1_pt << ", eta = " << mu1_abs_eta << std::endl;
      if (verbose) std::cout << "Muon 2 pT = " << mu2_pt << ", eta = " << mu2_abs_eta << std::endl;
  
      //  for 120 and 130 samples, shift mass by 5 GeV
      if (Sample_ID == 1202325) dimu_mass += 5.0;
      if (Sample_ID == 1302325) dimu_mass -= 5.0;
    
      /////////////////////////////////////////////////////
      ///  Loop over factories and set variable values  ///
      /////////////////////////////////////////////////////
      for (int iFact = 0; iFact < factories.size(); iFact++) {
	
	TString sig_name = std::get<11>(factories.at(iFact));
	TString bkg_name = std::get<12>(factories.at(iFact));
	TString jet_cut = std::get<13>(factories.at(iFact));
	
	if (verbose) std::cout << "\n  * For factory " << iFact << ", sig_name = " << sig_name << ", bkg_name = " << bkg_name
			       << ", jet_cut = " << jet_cut << ", nJets = " << nJets << std::endl;
	

	// Set vars equal to default vector of variables for this factory
	var_names = std::get<3>(factories.at(iFact));
	var_vals = std::get<4>(factories.at(iFact));
	
	// Fill all variables
	for (int iVar = 0; iVar < var_names.size(); iVar++) {
	  TString vName = var_names.at(iVar);

	  /////////////////////////////
	  ///  Spectator variables  ///
	  /////////////////////////////
	  
	  if      ( vName == "Sample_ID" )
	    var_vals.at(iVar) = Sample_ID;
	  else if ( vName == "event_wgt")
	    var_vals.at(iVar) = event_wgt;
	  else if ( vName == "xsec_norm")
	    var_vals.at(iVar) = xsec_norm;
	  else if ( vName == "lep_ID" )
	    var_vals.at(iVar) = lep_ID;
//	  else if ( vName == "dimu_mass" )
//	    var_vals.at(iVar) = dimu_mass;
	  
	  /////////////////////////////////////////////////////
	  ///  Variables automatically set in lib/VarSet.h  ///
	  /////////////////////////////////////////////////////
	  
	  else { // In 2016, variables were automatically set in lib/VarSet.h
	    
	    // var_vals.at(iVar) = samp->vars.getValue(vName.Data());
	    
	    // Muon variables
	    if      (vName == "mu1_pt")      	var_vals.at(iVar) = mu1_pt;
	    else if (vName == "mu2_pt")      	var_vals.at(iVar) = mu2_pt;
	    else if (vName == "mu1_abs_eta")    var_vals.at(iVar) = mu1_abs_eta;
	    else if (vName == "mu2_abs_eta")    var_vals.at(iVar) = mu2_abs_eta;
	    else if (vName == "mu1_lepMVA")     var_vals.at(iVar) = mu1_lepMVA;	    
	    else if (vName == "mu2_lepMVA")     var_vals.at(iVar) = mu2_lepMVA;

	    else if (vName == "dimu_mass")      var_vals.at(iVar) = dimu_mass;
	    else if (vName == "dimu_pt")    	var_vals.at(iVar) = dimu_pt;
	    else if (vName == "dimu_abs_eta")   var_vals.at(iVar) = dimu_abs_eta;
	    else if (vName == "dimu_abs_dEta") 	var_vals.at(iVar) = dimu_abs_dEta;
	    else if (vName == "dimu_abs_dPhi") 	var_vals.at(iVar) = dimu_abs_dPhi;
	    else if (vName == "dimu_dR")       	var_vals.at(iVar) = dimu_dR;
	    else if (vName == "dimu_mass_err")  var_vals.at(iVar) = dimu_mass_err;
	    else if (vName == "dimu_gen_ID")    var_vals.at(iVar) = dimu_gen_ID;

	    else if (vName == "cts_mu1")        var_vals.at(iVar) = cts_mu1;
            else if (vName == "cts_mu_pos")     var_vals.at(iVar) = cts_mu_pos;
	    
	    // lep variables
	    else if (vName == "lep1_pt")         var_vals.at(iVar) = lep1_pt;
            else if (vName == "lep2_pt")         var_vals.at(iVar) = lep2_pt;
            else if (vName == "lep1_abs_eta")    var_vals.at(iVar) = lep1_abs_eta;
            else if (vName == "lep2_abs_eta")    var_vals.at(iVar) = lep2_abs_eta;
            else if (vName == "lep1_lepMVA")     var_vals.at(iVar) = lep1_lepMVA;
            else if (vName == "lep2_lepMVA")     var_vals.at(iVar) = lep2_lepMVA;

	    else if (vName == "dilep_mass")	 var_vals.at(iVar) = dilep_mass;
	    else if (vName == "dilep_pt")	 var_vals.at(iVar) = dilep_pt;
	    else if (vName == "dilep_abs_eta")	 var_vals.at(iVar) = dilep_abs_eta;
	    else if (vName == "dilep_abs_dEta")	 var_vals.at(iVar) = dilep_abs_dEta;
	    else if (vName == "dilep_abs_dPhi")	 var_vals.at(iVar) = dilep_abs_dPhi;
	    else if (vName == "dilep_dR")	 var_vals.at(iVar) = dilep_dR;
	    else if (vName == "dilep_gen_ID")    var_vals.at(iVar) = dilep_gen_ID;

	    else if (vName == "cts_lep1")	 var_vals.at(iVar) = cts_lep1;
	    else if (vName == "cts_lep_pos")	 var_vals.at(iVar) = cts_lep_pos;

	    // pair variables
	    else if (vName == "quadlep_mass")	 var_vals.at(iVar) = quadlep_mass;
	    else if (vName == "quadlep_pt")	 var_vals.at(iVar) = quadlep_pt;
	    else if (vName == "quadlep_abs_eta") var_vals.at(iVar) = quadlep_abs_eta;
	    else if (vName == "dipair_dEta_H")	 var_vals.at(iVar) = dipair_dEta_H;
	    else if (vName == "dipair_dPhi_H")   var_vals.at(iVar) = dipair_dPhi_H;
	    else if (vName == "dipair_dR_H")     var_vals.at(iVar) = dipair_dR_H;
	    else if (vName == "cts_dipair_H")    var_vals.at(iVar) = cts_dipair_H;

            else if (vName == "cs_costheta")     var_vals.at(iVar) = cs_costheta;
	    else if (vName == "cs_cosphi")       var_vals.at(iVar) = cs_cosphi;
            else if (vName == "cs_sinphi")       var_vals.at(iVar) = cs_sinphi;
	    else if (vName == "cos_theta1")      var_vals.at(iVar) = cos_theta1;
            else if (vName == "cos_phiH")        var_vals.at(iVar) = cos_phiH;
	    else if (vName == "cos_phi1")        var_vals.at(iVar) = cos_phi1;

	    //met variables
            else if (vName == "met_pt")       		var_vals.at(iVar) = met_pt;
	    else if (vName == "abs_dPhi_4lPt_met")      var_vals.at(iVar) = abs_dPhi_4lPt_met;
            else if (vName == "met_Rtrans_4lPt")       	var_vals.at(iVar) = met_Rtrans_4lPt; 
	    else if (vName == "met_Rlongi_4lPt")       	var_vals.at(iVar) = met_Rlongi_4lPt;

            else if (vName == "mht_pt")       		var_vals.at(iVar) = mht_pt;
	    else if (vName == "mht_mass")       	var_vals.at(iVar) = mht_mass;
	    else if (vName == "abs_dPhi_4lPt_mht")      var_vals.at(iVar) = abs_dPhi_4lPt_mht;
            else if (vName == "mht_Rtrans_4lPt")       	var_vals.at(iVar) = mht_Rtrans_4lPt; 
	    else if (vName == "mht_Rlongi_4lPt")       	var_vals.at(iVar) = mht_Rlongi_4lPt;

            else if (vName == "abs_dPhi_met_mht")       var_vals.at(iVar) = abs_dPhi_met_mht; 

//            else if (vName == "")       var_vals.at(iVar) = ;
//            else if (vName == "")       var_vals.at(iVar) = ;
//            else if (vName == "")       var_vals.at(iVar) = ;
//            else if (vName == "")       var_vals.at(iVar) = ;

	    //dijet variables
	    else if (vName == "dijet_mass") 	var_vals.at(iVar) = dijet_mass;
	    else if (vName == "dijet_pt") 	var_vals.at(iVar) = dijet_pt;
	    else if (vName == "dijet_abs_eta") 	var_vals.at(iVar) = dijet_abs_eta;
	    else if (vName == "dijet_abs_dEta") var_vals.at(iVar) = dijet_abs_dEta;
            else if (vName == "dijet_abs_dPhi") var_vals.at(iVar) = dijet_abs_dPhi;
            else if (vName == "dijet_dR")       var_vals.at(iVar) = dijet_dR;

	    // Jet variables
	    else if (vName == "jet1_pt")  	var_vals.at(iVar) = jet1_pt;
	    else if (vName == "jet2_pt")  	var_vals.at(iVar) = jet2_pt;
	    else if (vName == "jet1_abs_eta") 	var_vals.at(iVar) = jet1_abs_eta;
	    else if (vName == "jet2_abs_eta") 	var_vals.at(iVar) = jet2_abs_eta;
	    else if (vName == "jet0_pt")  	var_vals.at(iVar) = jet0_pt;
	    else if (vName == "jet0_abs_eta") 	var_vals.at(iVar) = jet0_abs_eta;
 
	    // Global event variables
	    else if (vName == "nMuons")     	var_vals.at(iVar) = nMuons;     	
	    else if (vName == "nEles")     	var_vals.at(iVar) = nEles;      	
	    else if (vName == "nJets")     	var_vals.at(iVar) = nJets;      	
	    else if (vName == "nBJets_Loose")   var_vals.at(iVar) = nBJets_Loose;	
	    else if (vName == "nBJets_Med")     var_vals.at(iVar) = nBJets_Med; 	   
            else if (vName == "nBJets_Tight")   var_vals.at(iVar) = nBJets_Tight;	   
            else if (vName == "nFwdJets")     	var_vals.at(iVar) = nFwdJets;   	   
            else if (vName == "nCentJets")     	var_vals.at(iVar) = nCentJets;  	
	    else if (vName == "BDT_noMass_v3")	var_vals.at(iVar) = BDT_noMass_v3;	   
 
	    if (verbose) std::cout << "  * Filled variable " << vName << " with value " << var_vals.at(iVar) << std::endl;
	  }
	  
	} // End loop: for (int iVar = 0; iVar < var_names.size(); iVar++)
	
	TString multi_str = " " ;
//	if (samp_ID == -1 || samp_ID == -2)
//	  multi_str = "nonVH";
//	else if (samp_ID == -3)
//	  multi_str = "ZH";
//	else if (samp_ID == -4 || samp_ID == -5)
//	  multi_str = "WH";
//	else if (samp_ID == 16)
//	  multi_str = "WZ";
//	else if (samp_ID == 18)
//	  multi_str = "ZZ";
//	else if (samp_ID == 19 || samp_ID == 20)
//	  multi_str = "ttX";
//	else if (samp_ID == 0)
//	  multi_str = "DATA";
//	else
//	  multi_str = "OTHER";
	
	// // Unweighted signal and background
	// double sig_evt_weight = 1.0;
	// double bkg_evt_weight = 1.0;


	// Weight by expected sample normalization
//	double sig_evt_weight = samp_wgt * 1000.;
//	double bkg_evt_weight = samp_wgt;
//	if (MULTICLASS) sig_evt_weight *= 0.001; // Don't weight signal for MultiClass
	
	// // Weight by expected sample normalization x signal resolution
	// double sig_evt_weight = samp_wgt * res_wgt * 1000.;
	// double bkg_evt_weight = samp_wgt;

	double MVA_cut = 0.4;
//	if (mu1_lepMVA < MVA_cut or mu2_lepMVA < MVA_cut or lep1_lepMVA < MVA_cut or lep2_lepMVA < MVA_cut) continue;

	
	// Load values into event
	if ( MN.first == "signal" and (Sample_ID == 2325 or Sample_ID == 1202325 or Sample_ID == 1302325) ) { // Signal MC    Sample_ID > 0 , 2325 is ZH
	  if ( abs(dilep_mass-91)>10 or dimu_gen_ID!=25 ) continue;  // Z-cand mass range, dimu from H
	  if ( (iEvt % (2*presc)) != 1 ) {   // only iEvt % 2 == 1 is used for BDT (the other half for limits), now limited by stats.   -- XWZ 02.12.2018
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddSignalTrainingEvent( var_vals, event_wgt * xsec_norm ); // signal and bkg are automatically reweighted by TMVA
	      nTrain_sig.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTrainingEvent( multi_str, var_vals, event_wgt * xsec_norm );
	      if (multi_str == "nonVH") nTrain_nonVH.at(iFact) += 1;
	      if (multi_str == "ZH")    nTrain_ZH.at(iFact)    += 1;
	      if (multi_str == "WH")    nTrain_WH.at(iFact)    += 1;
	    }
	  } else {
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddSignalTestEvent( var_vals, event_wgt * xsec_norm );
	      nTest_sig.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTestEvent( multi_str, var_vals, event_wgt * xsec_norm );
	      if (multi_str == "nonVH") nTest_nonVH.at(iFact) += 1;
	      if (multi_str == "ZH")    nTest_ZH.at(iFact)    += 1;
	      if (multi_str == "WH")    nTest_WH.at(iFact)    += 1;
	    }
	  }
	}
	if ( MN.first == "bkg" and Sample_ID < 0 and Sample_ID != -999 ) { // Background MC       Sample_ID < 0
	  if ( abs(dilep_mass-91)>10 ) continue; // Z-cand mass range, dimu from Z
	  if ( (iEvt % (2*presc)) != 0 ) {
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddBackgroundTrainingEvent( var_vals, event_wgt * xsec_norm );
	      nTrain_bkg.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTrainingEvent( multi_str, var_vals, event_wgt * xsec_norm );
	      if (multi_str == "WZ")  nTrain_WZ.at(iFact)  += 1;
	      if (multi_str == "ZZ")  nTrain_ZZ.at(iFact)  += 1;
	      if (multi_str == "ttX") nTrain_ttX.at(iFact) += 1;
	    }
	  } else {
	    if (!MULTICLASS) {
	      std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, event_wgt * xsec_norm );
	      nTest_bkg.at(iFact) += 1;
	    } else {
	      std::get<1>(factories.at(iFact))->AddTestEvent( multi_str, var_vals, event_wgt * xsec_norm );
	      if (multi_str == "WZ")  nTest_WZ.at(iFact)  += 1;
	      if (multi_str == "ZZ")  nTest_ZZ.at(iFact)  += 1;
	      if (multi_str == "ttX") nTest_ttX.at(iFact) += 1;
	    }
	  }
	}
//	if (Sample_ID == 0) { // Data
//	  if (!MULTICLASS) {
//	    std::get<1>(factories.at(iFact))->AddBackgroundTestEvent( var_vals, event_wgt );
//	    nTest_bkg.at(iFact) += 1;
//	  } else {
//	    std::get<1>(factories.at(iFact))->AddTestEvent( multi_str, var_vals, event_wgt );
//	  }
//	}
	
      } // End loop: for (int iFact = 0; iFact < factories.size(); iFact++)
      
      nEvt_pass += 1;
      nEvt_tot  += 1;
      if (Sample_ID == 2325) nEvt_sig += 1;   //Sample_ID > 0
      if (Sample_ID == -2323) nEvt_bkg += 1;    //Sample_ID < 0
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
    std::string NTr_nonVH;
    std::string NTr_ZH;
    std::string NTr_WH;
    std::string NTr_WZ;
    std::string NTr_ZZ;
    std::string NTr_ttX;
    
    std::ostringstream convertTr_nonVH;
    convertTr_nonVH << nTrain_nonVH.at(iFact);
    NTr_nonVH = convertTr_nonVH.str();
    
    std::ostringstream convertTr_ZH;
    convertTr_ZH << nTrain_ZH.at(iFact);
    NTr_ZH = convertTr_ZH.str();
    
    std::ostringstream convertTr_WH;
    convertTr_WH << nTrain_WH.at(iFact);
    NTr_WH = convertTr_WH.str();
    
    std::ostringstream convertTr_WZ;
    convertTr_WZ << nTrain_WZ.at(iFact);
    NTr_WZ = convertTr_WZ.str();
    
    std::ostringstream convertTr_ZZ;
    convertTr_ZZ << nTrain_ZZ.at(iFact);
    NTr_ZZ = convertTr_ZZ.str();
    
    std::ostringstream convertTr_ttX;
    convertTr_ttX << nTrain_ttX.at(iFact);
    NTr_ttX = convertTr_ttX.str();
    
    std::string numTrainStr = "nTrain_Signal="+NTrS+":nTrain_Background="+NTrB+":";
    if (MULTICLASS) {
      numTrainStr = "";
      if (nTrain_nonVH.at(iFact) > 0) numTrainStr += "nTrain_nonVH="+NTr_nonVH+":";
      if (nTrain_ZH.at(iFact)    > 0) numTrainStr += "nTrain_ZH="+NTr_ZH+":";
      if (nTrain_WH.at(iFact)    > 0) numTrainStr += "nTrain_WH=" +NTr_WH +":";
      if (nTrain_WZ.at(iFact)    > 0) numTrainStr += "nTrain_WZ="+NTr_WZ+":";
      if (nTrain_ZZ.at(iFact)    > 0) numTrainStr += "nTrain_ZZ="+NTr_ZZ+":";
      if (nTrain_ttX.at(iFact)   > 0) numTrainStr += "nTrain_ttX="+NTr_ttX+":";
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
			 "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=10:MaxDepth=3" );
    
    if (Use["BDTG_UF_v2"]) // Test settings - AWB 29.10.2018
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v2", (std::string)
			 "!H:!V:NTrees=200::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			 "BaggedSampleFraction=0.5:nCuts=10:MaxDepth=3" );

    if (Use["BDTG_UF_v3"]) // Test settings - AWB 29.10.2018
      factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_UF_v3", (std::string)
                         "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
                         "BaggedSampleFraction=0.5:nCuts=10:MaxDepth=5" );   
 
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
  ZH_4l_miniNTuple(methodList);
  return 0;
}
