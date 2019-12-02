
///////////////////////////////////////////////////////////
///  Macro to study dependence of RECO - GEN pT vs. d0  ///
///                                                     ///
///            Andrew Brinkerhoff 24.08.2018            ///
///////////////////////////////////////////////////////////

// Basic ROOT includes to read and write files
#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h" // List of branches in the NTuple tree
#include "H2MuAnalyzer/MakeHistos/interface/HistoHelper.h"        // Use to book and fill output histograms
// #include "H2MuAnalyzer/MakeHistos/interface/SelectionCuts.h"      // Common event selections
// #include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"       // Common category definitions

// #include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h" // List of branches in the NTuple tree
// #include "H2MuAnalyzer/MakeHistos/interface/HistoHelper.h"        // Use to book and fill output histograms
// #include "H2MuAnalyzer/MakeHistos/interface/MassCalPlots.h"     // class and config for mass calibration study
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"    // Common object selections
#include "H2MuAnalyzer/MakeHistos/interface/EventSelection.h"     // Common event selections
#include "H2MuAnalyzer/MakeHistos/interface/EventWeight.h"        // Common event weights
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"       // Common category definitions
#include "H2MuAnalyzer/MakeHistos/interface/MiniNTupleHelper.h"   // "PlantTree" and "BookBranch" functions
#include "H2MuAnalyzer/MakeHistos/interface/ReadMVA.h"            // Read and evaluate XMLs for MVA

// Load the library of the local, compiled H2MuAnalyzer/MakeHistos directory
R__LOAD_LIBRARY(../../../tmp/slc6_amd64_gcc630/src/H2MuAnalyzer/MakeHistos/src/H2MuAnalyzerMakeHistos/libH2MuAnalyzerMakeHistos.so)

// Hard-coded options for running locally / manually
// Options passed in as arguments to ReadNTupleChain when running in batch mode
const int MIN_FILE = 1;     // Minimum index of input files to process
const int MAX_FILE = 5;     // Maximum index of input files to process
const int MAX_EVT  = 10000; // Maximum number of events to process
const int PRT_EVT  = 1000;  // Print every N events
const float SAMP_WGT = 1.0;
// const float LUMI = 36814; // pb-1
const bool verbose = false; // Print extra information


// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2016/94X_v3/prod-v16.0.7/SingleMuon/SingleMu_2016B/190714_182320/0000/";
const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2016/94X_v3/prod-v16.0.7/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/190714_182527/0000/";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/prod-v18-pre-tag/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/tt_ll_POW/190521_174438/0000/";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/prod-v18-pre-tag/ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8/H2Mu_ttH_125/190521_173336/0000/";
const TString SAMPLE   = "ZJets_AMC";
// const TString SAMPLE   = "H2Mu_ttH_125";
// const TString SAMPLE   = "tt_ll_POW";
// const TString SAMPLE   = "SingleMu_2016B";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/180802_165055/0000";
// const TString SAMPLE   = "ZJets_AMC";
// const TString OUT_DIR  = "plots";
const TString OUT_DIR = "/afs/cern.ch/work/e/eyigitba/public/H2Mu/2016/Histograms";
const std::string YEAR  = "2016";
const std::string SLIM  = "notSlim"; 
const TString HIST_TREE = "HistTree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both


const std::vector<std::string> SEL_CUTS = {}; // Cuts which every event must pass
const std::vector<std::string> OPT_CUTS = {}; // Cuts which every event must pass
// const std::vector<std::string> CAT_CUTS = {"nVtx_0_21", "nVtx_22_26", "nVtx_27_33", "nVtx_34_inf", "inclusive"};  // 5 nVertices bins: 20 - 35, 35 - 42.5, 42.5 - 50, > 50 GeV, and inclusive
const std::vector<std::string> CAT_CUTS = {"inclusive"};  // inclusive
const std::vector<std::string> ETA_CUTS = {"eta_0_0p9", "eta_0p9_1p7", "eta_1p7_inf", "eta_inc"};  // 5 muon pT bins: 20 - 35, 35 - 42.5, 42.5 - 50, > 50 GeV, and inclusive

// Command-line options for running in batch.  Running "root -b -l -q macros/GenRecoPtDiffVsD0VsPt.C" will use hard-coded options above.
void Compared0BSd0PV( TString sample = "", TString in_dir = "", TString out_dir = "",
			    std::vector<TString> in_files = {}, TString out_file_str = "",
			    int max_evt = 0, int prt_evt = 0, float samp_weight = 1.0,
          TString hist_tree = "", std::string SYS = "" ) {
  
  // Set variables to hard-coded values if they are not initialized
  if (sample.Length()  == 0) sample  = SAMPLE;
  if (in_dir.Length()  == 0) in_dir  = IN_DIR;
  if (out_dir.Length() == 0) out_dir = OUT_DIR;
  if (max_evt          == 0) max_evt = MAX_EVT;
  if (prt_evt          == 0) prt_evt = PRT_EVT;
  if (samp_weight      == 0) samp_weight = SAMP_WGT;
  // if (hist_tree.Length() == 0) hist_tree   = HIST_TREE;

  
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
      std::cout << "Adding file " << in_file_name.Data() << std::endl;
      in_file_names.push_back(in_file_name.Data());
    }
  }

  // Open all input files
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    if ( !gSystem->AccessPathName(in_file_names.at(i)) )
      file_tmp = TFile::Open( in_file_names.at(i) ); // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
      return;
    }
  }
  
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
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM}, "noSys", false); // Options in {} include "JES", "Flags", and "SFs"
    else
      SetBranchAddresses(*in_chain, br, {YEAR, SLIM, "GEN", "Wgts"}, "noSys", false); // Options in {} include "JES", "Flags", and "SFs"
  }

  gROOT->cd(); // Navigate to "local" memory, so all histograms are not saved in out_tuple

  if (verbose) std::cout << "Finished with SetBranchAddress calls" << std::endl;

  // Configuration for object selection, event selection, and object weighting
  ObjectSelectionConfig obj_sel;
  EventSelectionConfig  evt_sel;
  EventWeightConfig     evt_wgt;
  ConfigureObjectSelection(obj_sel, YEAR);
  ConfigureEventSelection (evt_sel, YEAR);
  ConfigureEventWeight    (evt_wgt, YEAR);

  evt_sel.muPair_mass_min = 12; // Allow masses down to 12 GeV (instead of 60 GeV) for background studies

  if (verbose) obj_sel.Print();
  if (verbose) evt_sel.Print();
  if (verbose) evt_wgt.Print();

  std::string PTC = obj_sel.mu_pt_corr; // Store muon pT correction in a shorter string; not changed later

  // MassCalPlot config add or not??
  
  std::cout << "\n******* About to enter the event loop *******" << std::endl;
  for (UInt_t iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++) {
    
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
    if (not pass_sel_cuts) continue;

    // Get event weight for MC, defined in src/EventWeight.cc
    bool isData = sample.Contains("SingleMu");
    float event_wgt = ( isData ? 1.0 : EventWeight(br, evt_wgt, verbose) );
    // float event_wgt = 1.0;

    int nVertices = br.nVertices;

 
    if ( SelectedMuPairs(obj_sel, br).size()!=1 or SelectedMuons(obj_sel, br).size()!=2 ) continue;
    MuPairInfo    dimu;
    dimu = SelectedCandPair(obj_sel, br);

    MuonInfo mu1, mu2, muP, muN;
    GenMuonInfo genMu1, genMu2, genMuDummy;
    mu1 = br.muons->at(dimu.iMu1);
    mu2 = br.muons->at(dimu.iMu2);



    TLorentzVector mu_vec1 = FourVec( mu1, "PF" );
    TLorentzVector mu_vec2 = FourVec( mu2, "PF" );
   
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Custom selection: require exactly one GEN Z --> mu-mu pair, matched to a pair of RECO muons //
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // for MC samples, find gen match
    GenParentInfo    gen_dimu;
    TLorentzVector gen_vec;
    gen_vec.SetPtEtaPhiM(0,0,0,-999);

    // Z parent matching
    if (not isData) {
      for (const auto & genPa : *br.genParents) {
        if (genPa.ID != 23 or genPa.daughter_1_ID + genPa.daughter_2_ID != 0) continue;  // legit Z 
        if (genPa.daughter_1_ID != 13 and genPa.daughter_1_ID != -13) continue;         // decays to dimuon
        if (genPa.daughter_1_idx < 0 or genPa.daughter_2_idx < 0) continue;             // present in genMuons collection
        TLorentzVector gen_vec1 = FourVec( br.genMuons->at(genPa.daughter_1_idx) );
        genMu1 = br.genMuons->at(genPa.daughter_1_idx);
        TLorentzVector gen_vec2 = FourVec( br.genMuons->at(genPa.daughter_2_idx) );
        genMu2 = br.genMuons->at(genPa.daughter_2_idx);
        if ( mu_vec1.DeltaR(gen_vec1)<0.05 and mu_vec2.DeltaR(gen_vec2)<0.05 ) gen_dimu = genPa;
        else if ( mu_vec1.DeltaR(gen_vec2)<0.05 and mu_vec2.DeltaR(gen_vec1)<0.05 ) {
          gen_dimu = genPa;
          genMuDummy = genMu1;
          genMu1 = genMu2;
          genMu2 = genMuDummy;
        }
      }
    }

    // H parent matching
    // if (not isData) {
    //   for (const auto & genPa : *br.genParents) {
    //     if (genPa.ID != 25 or genPa.daughter_1_ID + genPa.daughter_2_ID != 0) continue;  // legit H 
    //     if (genPa.daughter_1_ID != 13 and genPa.daughter_1_ID != -13) continue;         // decays to dimuon
    //     if (genPa.daughter_1_idx < 0 or genPa.daughter_2_idx < 0) continue;             // present in genMuons collection
    //     TLorentzVector gen_vec1 = FourVec( br.genMuons->at(genPa.daughter_1_idx) );
    //     genMu1 = br.genMuons->at(genPa.daughter_1_idx);
    //     TLorentzVector gen_vec2 = FourVec( br.genMuons->at(genPa.daughter_2_idx) );
    //     genMu2 = br.genMuons->at(genPa.daughter_2_idx);
    //     if ( mu_vec1.DeltaR(gen_vec1)<0.05 and mu_vec2.DeltaR(gen_vec2)<0.05 ) gen_dimu = genPa;
    //     else if ( mu_vec1.DeltaR(gen_vec2)<0.05 and mu_vec2.DeltaR(gen_vec1)<0.05 ) {
    //       gen_dimu = genPa;
    //       genMuDummy = genMu1;
    //       genMu1 = genMu2;
    //       genMu2 = genMuDummy;
    //     }
    //   }
    // }

    // for ttbar MC samples W parent matching
    // if (not isData) {
    //   for (const auto & genMuon : *br.genMuons) {
    //     if (abs(genMuon.mother_ID != 24)) continue; // comes from legit W
    //     gen_vec = FourVec( genMuon );
    //     if ( abs(mu_vec1.DeltaR(gen_vec)) < 0.05 and genMuon.charge == mu1.charge) genMu1 = genMuon;
    //     else if ( abs(mu_vec2.DeltaR(gen_vec)) < 0.05 and genMuon.charge == mu2.charge) genMu2 = genMuon;
    //   }
    // }

    // if ( gen_vec.M() == -999 and not isData) continue;
    if ( gen_dimu.mass == -999 and not isData) continue;

    // gen match found

    // Loop through category cuts defined in src/CategoryCuts.cc

    for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
      // TString optCatStr = OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat);
      
      double nVtx[2] = {0, 999};
      if        (CAT_CUTS.at(iCat).find("nVtx_0_21") != string::npos) {
        nVtx[0] = 0;
        nVtx[1] = 21;
      } 
      else if (CAT_CUTS.at(iCat).find("nVtx_22_26") != string::npos) {
        nVtx[0] = 22;
        nVtx[1] = 26;
      } 
      else if (CAT_CUTS.at(iCat).find("nVtx_27_33") != string::npos) {
        nVtx[0] = 27;
        nVtx[1] = 33;
      } 
      else if (CAT_CUTS.at(iCat).find("nVtx_34_inf") != string::npos) {
        nVtx[0] = 34;
        nVtx[1] = 999;
      }
      else if (CAT_CUTS.at(iCat).find("inclusive") != string::npos) {
        nVtx[0] = 0;
        nVtx[1] = 999;
      }

      for (int iEta = 0; iEta < ETA_CUTS.size(); iEta++) {
        double eta[2] = {0, 999};
        if        (ETA_CUTS.at(iEta).find("eta_0_0p9") != string::npos) {
          eta[0] = 0;
          eta[1] = 0.9;
        } 
        else if (ETA_CUTS.at(iEta).find("eta_0p9_1p7") != string::npos) {
          eta[0] = 0.9;
          eta[1] = 1.7;
        } 
        else if (ETA_CUTS.at(iEta).find("eta_1p7_inf") != string::npos) {
          eta[0] = 1.7;
          eta[1] = 3.0;
        } 
        else if (ETA_CUTS.at(iEta).find("eta_inc") != string::npos) {
          eta[0] = 0;
          eta[1] = 999;
        } 

        // std::string h_pre = "h_"+OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)+"_"+ETA_CUTS.at(iEta)+"_";
        std::string h_pre = "h_"+CAT_CUTS.at(iCat)+"_"+ETA_CUTS.at(iEta)+"_";

        /////////////////////////////////
        ///  Generate and fill plots  ///
        /////////////////////////////////
        for(int iMuon = 0; iMuon < 2; iMuon++){
          MuonInfo mu;
          GenMuonInfo genMu;
          if (iMuon == 0) {
            mu = mu1;
            genMu = genMu1;
          }
          else{
            mu = mu2;
            genMu = genMu2;
          }
          TLorentzVector mu_vec = FourVec( mu, "PF" );
          gen_vec = FourVec( genMu );
          if (abs(mu_vec.DeltaR(gen_vec)) > 0.05) continue;

          //See if muon passes the eta cuts
          if (eta[0] < fabs(mu.eta) && fabs(mu.eta) < eta[1]){
            float mu_pt_KinRoch  = mu.pt_kinfit * mu.pt_Roch / mu.pt;
            float mu_pt_KinKaMu  = mu.pt_kinfit * mu.pt_KaMu / mu.pt;
            float mu_ptCorr = 0.0;
            if ( fabs(mu.eta) < 0.9 ) {
              if ( YEAR == "2016" ) mu_ptCorr = mu.pt_Roch - 411.343/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
              if ( YEAR == "2017" ) mu_ptCorr = mu.pt_Roch - 650.839/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
              if ( YEAR == "2018" ) mu_ptCorr = mu.pt_Roch - 582.32/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
            }
            if ( 0.9 < fabs(mu.eta) && fabs(mu.eta) < 1.7){
              if ( YEAR == "2016" ) mu_ptCorr = mu.pt_Roch - 673.398/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
              if ( YEAR == "2017" ) mu_ptCorr = mu.pt_Roch - 988.369/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
              if ( YEAR == "2018" ) mu_ptCorr = mu.pt_Roch - 974.047/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
            }
            if ( 1.7 < fabs(mu.eta) ){
              if ( YEAR == "2016" ) mu_ptCorr = mu.pt_Roch - 1098.984/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
              if ( YEAR == "2017" ) mu_ptCorr = mu.pt_Roch - 1484.616/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
              if ( YEAR == "2018" ) mu_ptCorr = mu.pt_Roch - 1263.388/10000 * mu.pt_Roch * mu.pt_Roch * mu.d0_BS ;
            }
            // if (nVtx[0] <= nVertices && nVertices <= nVtx[1]){
              // BookAndFill( h_map_1D, h_pre+"d0_BS", 400, -0.01,  0.01, mu.d0_BS*mu.charge,event_wgt ); // d0_BS
              // BookAndFill( h_map_1D, h_pre+"d0_PV", 400, -0.01,  0.01, mu.d0_PV*mu.charge,event_wgt ); // d0_PV
              // BookAndFill( h_map_1D, h_pre+"d0_diff", 400, -0.01,  0.01, mu.d0_BS*mu.charge - mu.d0_PV*mu.charge,event_wgt ); // d0_BS - d0_PV
              // BookAndFill( h_map_2D, h_pre+"d0_BS_vs_d0_PV", 400, -0.01,  0.01, 400, -0.01,  0.01, mu.d0_PV*mu.charge, mu.d0_BS*mu.charge, event_wgt);

              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_vs_d0_PV", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_PV*mu.charge, 10000*(mu.pt      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_Roch_vs_d0_PV", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_PV*mu.charge, 10000*(mu.pt_Roch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KaMu_vs_d0_PV", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_PV*mu.charge, 10000*(mu.pt_KaMu      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_kinfit_vs_d0_PV", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_PV*mu.charge, 10000*(mu.pt_kinfit      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KinRoch_vs_d0_PV", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_PV*mu.charge,  10000*(mu_pt_KinRoch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KinKaMu_vs_d0_PV", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_PV*mu.charge, 10000*(mu_pt_KinKaMu      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );

              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_vs_d0_BS", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu.pt      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              BookAndFill( h_map_2D, h_pre+"dRelPt2p0_Roch_vs_d0_BS", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu.pt_Roch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              BookAndFill( h_map_2D, h_pre+"pt_Roch_vs_d0_BS", 400, -0.01,  0.01, 4000,    15, 300, mu.d0_BS*mu.charge, mu.pt_Roch,event_wgt );
              BookAndFill( h_map_2D, h_pre+"dRelPt1p0_Roch_vs_d0_BS", 400, -0.01,  0.01, 4000,    -100, 100, mu.d0_BS*mu.charge, 100*(mu.pt_Roch      - genMu.pt) /(genMu.pt),event_wgt );
              BookAndFill( h_map_2D, h_pre+"dPt_Roch_vs_d0_BS", 400, -0.01,  0.01, 4000,    -100, 100, mu.d0_BS*mu.charge, (mu.pt_Roch      - genMu.pt),event_wgt );
              
              BookAndFill( h_map_2D, h_pre+"dRelPt2p0_corr_vs_d0_BS", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu_ptCorr      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              BookAndFill( h_map_2D, h_pre+"pt_corr_vs_d0_BS", 400, -0.01,  0.01, 4000,    15, 300, mu.d0_BS*mu.charge, mu_ptCorr,event_wgt );
              BookAndFill( h_map_2D, h_pre+"dRelPt1p0_corr_vs_d0_BS", 400, -0.01,  0.01, 4000,    -100, 100, mu.d0_BS*mu.charge, 100*(mu_ptCorr      - genMu.pt) /(genMu.pt),event_wgt );
              BookAndFill( h_map_2D, h_pre+"dPt_corr_vs_d0_BS", 400, -0.01,  0.01, 4000,    -100, 100, mu.d0_BS*mu.charge, (mu_ptCorr      - genMu.pt),event_wgt );

              BookAndFill( h_map_2D, h_pre+"pTCorr_vs_d0_BS_vs_pt_Roch", 400, -0.01,  0.01, 4000,    15, 300, mu.d0_BS*mu.charge, mu.pt_Roch, fabs(mu_ptCorr - mu.pt_Roch) / mu.pt_Roch );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KaMu_vs_d0_BS", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu.pt_KaMu      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_kinfit_vs_d0_BS", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu.pt_kinfit      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KinRoch_vs_d0_BS", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge,  10000*(mu_pt_KinRoch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KinKaMu_vs_d0_BS", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu_pt_KinKaMu      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              // if (mu.charge > 0){
              //     BookAndFill( h_map_2D, h_pre+"dRelPt2p0_Roch_vs_d0_BS_muP", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu.pt_Roch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              //     BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KinRoch_vs_d0_BS_muP", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge,  10000*(mu_pt_KinRoch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );

              // }
              // else if (mu.charge < 0){
              //     BookAndFill( h_map_2D, h_pre+"dRelPt2p0_Roch_vs_d0_BS_muN", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge, 10000*(mu.pt_Roch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
              //     BookAndFill( h_map_2D, h_pre+"dRelPt2p0_KinRoch_vs_d0_BS_muN", 400, -0.01,  0.01, 4000,    -40, 40, mu.d0_BS*mu.charge,  10000*(mu_pt_KinRoch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );

              // }
            // }
          } // End eta cuts
        }// End iMuon loop
      } // End iEta loop
    } // End iCat loop
  } // // End loop: for (int iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // End custom selection: require exactly one GEN Z --> mu-mu pair, matched to a pair of RECO muons //
    /////////////////////////////////////////////////////////////////////////////////////////////////////
      
  
  std::cout << "\n******* Leaving the event loop *******" << std::endl;

  std::cout << "\n******* normalizing histos, weight =  " << samp_weight << " *******" << std::endl;
  if (h_map_1D.empty()) std::cout << "h_map_1D is empty" << std::endl;
  else {
    for (std::map<TString, TH1*>::iterator it_term = h_map_1D.begin() ; it_term != h_map_1D.end() ; it_term++) {
        it_term->second ->Scale(samp_weight);
    }
  }
  if (h_map_2D.empty()) std::cout << "h_map_2D is empty" << std::endl;
  else {
    for (std::map<TString, TH2*>::iterator it_term = h_map_2D.begin() ; it_term != h_map_2D.end() ; it_term++) {
        it_term->second ->Scale(samp_weight);
    }
  }


  // Concatenate basic selection cuts into string for output file name
  std::string selStr = "";
  for (int i = 0; i < SEL_CUTS.size(); i++) selStr = selStr+SEL_CUTS.at(i)+"_";
  selStr.pop_back();

  // Place histograms made with separate selection cuts into separate files
  // Create output file
  TString out_file_name;
  if (out_file_str.Length() > 0) out_file_name.Form( "%s/histos_%s_%s_%s_new.root",    out_dir.Data(), sample.Data(), selStr.c_str(), out_file_str.Data() );
  else                           out_file_name.Form( "%s/histos_%s_%s_%d_%d_new.root", out_dir.Data(), sample.Data(), selStr.c_str(), MIN_FILE, MAX_FILE );
  std::cout << "\nCreating output file " << out_file_name.Data() << std::endl;
  TFile* out_file = TFile::Open( out_file_name, "RECREATE" );
  
  // Write output file
  if (verbose) std::cout << "\nWriting output file " << out_file_name.Data() << std::endl;
  out_file->cd();
  
  // Write out 1D histograms
  for (std::map<TString, TH1*>::iterator it = h_map_1D.begin(); it != h_map_1D.end(); ++it) {
     std::cout << "  * Writing 1D histogram " << it->second->GetName() << std::endl;
     it->second->Write();
  }
  // Write out 2D histograms
  for (std::map<TString, TH2*>::iterator it = h_map_2D.begin(); it != h_map_2D.end(); ++it) {
  	std::cout << "  * Writing 2D histogram " << it->second->GetName() << std::endl;
  	it->second->Write();
  }
  
  out_file->Write();
  std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
        
  std::cout << "\nExiting GenRecoPtDiffVsD0VsPt()\n";
  
} // End void GenRecoPtDiffVsD0VsPt()
