
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

const TString SAMPLE   = "ZJets_AMC";
// const TString SAMPLE   = "SingleMu_2016B";
// const TString IN_DIR   = "/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZJets_AMC/180802_165055/0000";
// const TString SAMPLE   = "ZJets_AMC";
// const TString OUT_DIR  = "plots";
const TString OUT_DIR = "/afs/cern.ch/work/e/eyigitba/public/H2Mu/2016/Histograms";
const std::string YEAR  = "2016";
const std::string SLIM  = "notSlim"; 
const TString HIST_TREE = "HistTree"; // "Hist", "Tree", or "HistTree" to output histograms, trees, or both


const std::vector<std::string> SEL_CUTS = {}; // Cuts which every event must pass
// const std::vector<std::string> OPT_CUTS = {"NONE"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> OPT_CUTS = {"d0_m39","d0_m38","d0_m37","d0_m36","d0_m35","d0_m34","d0_m33","d0_m32","d0_m31","d0_m30","d0_m29","d0_m28","d0_m27",
                                          "d0_m26","d0_m25","d0_m24","d0_m23","d0_m22","d0_m21","d0_m20","d0_m19","d0_m18","d0_m17","d0_m16","d0_m15","d0_m14",
                                          "d0_m13","d0_m12","d0_m11","d0_m10","d0_m9","d0_m8","d0_m7","d0_m6","d0_m5","d0_m4","d0_m3","d0_m2","d0_m1","d0_0", 
                                          "d0_p1","d0_p2","d0_p3","d0_p4","d0_p5","d0_p6","d0_p7","d0_p8","d0_p9","d0_p10","d0_p11","d0_p12","d0_p13","d0_p14",
                                          "d0_p15","d0_p16","d0_p17","d0_p18","d0_p19","d0_p20", "d0_p21","d0_p22","d0_p23","d0_p24","d0_p25","d0_p26","d0_p27","d0_p28",
                                          "d0_p29","d0_p30","d0_p31","d0_p32","d0_p33","d0_p34","d0_p35","d0_p36","d0_p37","d0_p38","d0_p39"}; // Multiple selection cuts, applied independently in parallel
const std::vector<std::string> CAT_CUTS = {"pt_20_35", "pt_35_42", "pt_42_50", "pt_50_inf", "inclusive"};  // 5 muon pT bins: 20 - 35, 35 - 42.5, 42.5 - 50, > 50 GeV, and inclusive
const std::vector<std::string> ETA_CUTS = {"eta_0_0p9", "eta_0p9_1p7", "eta_1p7_inf", "eta_inc"};  // 5 muon pT bins: 20 - 35, 35 - 42.5, 42.5 - 50, > 50 GeV, and inclusive

// Command-line options for running in batch.  Running "root -b -l -q macros/GenRecoPtDiffVsD0VsPt.C" will use hard-coded options above.
void GenRecoPtDiffVsD0VsPt( TString sample = "", TString in_dir = "", TString out_dir = "",
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
    // float event_wgt = ( isData ? 1.0 : EventWeight(br, evt_wgt, verbose) );
    float event_wgt = 1.0;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Custom selection: require exactly one GEN Z --> mu-mu pair, matched to a pair of RECO muons //
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////////////////
    ///  Loop through alternate, optional selection cuts defined in src/SelectionCuts.cc  ///
    /////////////////////////////////////////////////////////////////////////////////////////

    for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
      std::string OPT_CUT = OPT_CUTS.at(iOpt);
      // Loop over possible central d0 values
      // double d0 = 0.0;
      // for (int iD0 = 0; iD0 < ((OPT_CUTS.size() + 1) / 2); iD0++) {
        // d0 = 0.00025 * iD0; // Central d0 value in the bin
        // if (iD0 < 10){
        //   if ( (iD0 == 0 && OPT_CUTS.at(iOpt).find("d0_0") != string::npos) || (iD0 != 0 && OPT_CUTS.at(iOpt).find( std::to_string(iD0) ) != string::npos) ) {
        //     // Switch central d0 value to negative if needed
        //     if (OPT_CUTS.at(iOpt).find("d0_m") != string::npos) d0 *= -1;
        //     // See if mu1 or mu2 passes the d0 cuts
        //     break;
        //     // if ( (mu1.d0_PV * mu1.charge > d0 - 0.0005) && (mu1.d0_PV * mu1.charge < d0 + 0.0005) )
        //     // if ( (mu2.d0_PV * mu2.charge > d0 - 0.0005) && (mu2.d0_PV * mu2.charge < d0 + 0.0005) )
        //   }
        // }
        // else{
        //   if ( (OPT_CUTS.at(iOpt).find( std::to_string(iD0) ) != string::npos) ) {
        //     // Switch central d0 value to negative if needed
        //     if (OPT_CUTS.at(iOpt).find("d0_m") != string::npos) d0 *= -1;
        //     // See if mu1 or mu2 passes the d0 cuts
        //     break;
        // }
      // } // End loop: for (int iD0 = 0; iD0 < ((OPT_CUTS.size() + 1) / 2); iD0++)

      double d0 = 0.0;
      if (OPT_CUT.find("d0_0") == string::npos){
        d0 =0.00025 * stof(OPT_CUT.substr(4));
      }
      if (OPT_CUT.substr(3).find("m") != string::npos) d0 *= -1;
      // std::cout << OPT_CUT << " " << d0 << std::endl;
    
      if ( SelectedMuPairs(obj_sel, br).size()!=1 or SelectedMuons(obj_sel, br).size()!=2 ) continue;
      MuPairInfo    dimu;
      dimu = SelectedCandPair(obj_sel, br);
      
      MuonInfo mu1, mu2, muP, muN;
      GenMuonInfo genMu1, genMu2, genMuDummy;
      mu1 = br.muons->at(dimu.iMu1);
      mu2 = br.muons->at(dimu.iMu2);


      TLorentzVector mu_vec1 = FourVec( mu1, "PF" );
      TLorentzVector mu_vec2 = FourVec( mu2, "PF" );

      // for MC samples, find gen match
      GenParentInfo    gen_dimu;
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
      if ( gen_dimu.mass == -999 and not isData) continue;
      // gen match found

      // float dimu_mass_KinRoch = dimu.mass_kinfit * dimu.mass_Roch / dimu.mass;
      // float dimu_pt_KinRoch = dimu.pt_kinfit * dimu.pt_Roch / dimu.pt;
      // float mu1_pt_KinRoch  = mu1.pt_kinfit * mu1.pt_Roch / mu1.pt;
      // float mu2_pt_KinRoch    = mu2.pt_kinfit * mu2.pt_Roch / mu2.pt;

      // float dimu_mass_KinKaMu = dimu.mass_kinfit * dimu.mass_KaMu / dimu.mass;
      // float dimu_pt_KinKaMu   = dimu.pt_kinfit * dimu.pt_KaMu / dimu.pt;
      // float mu1_pt_KinKaMu    = mu1.pt_kinfit * mu1.pt_KaMu / mu1.pt;
      // float mu2_pt_KinKaMu    = mu2.pt_kinfit * mu2.pt_KaMu / mu2.pt;

      // BookAndFill(h_map_1D, "h_"+OPT_CUTS.at(iOpt)+"_d0", 10000, -500,  500, 10000*mu.d0_PV*mu.charge );
      // BookAndFill(h_map_1D, "h_"+OPT_CUTS.at(iOpt)+"_pt", 10000,    0, 1000, mu.pt );

      // Loop through category cuts defined in src/CategoryCuts.cc

      for (int iCat = 0; iCat < CAT_CUTS.size(); iCat++) {
        TString optCatStr = OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat);
        // if ( not InCategory(obj_sel, br, CAT_CUTS.at(iCat), verbose) ) continue;
        // if (not InCategory(br, CAT_CUTS.at(iCat), verbose)) continue;  // Not using centrally-defined category cuts
        
        // Loop over pT ranges (each covers ~1/4 of the Z --> mu-mu phase space)
        double pt[2] = {20, 999};
        if        (CAT_CUTS.at(iCat).find("pt_20_35") != string::npos) {
          pt[0] = 20;
          pt[1] = 35;
        } 
        else if (CAT_CUTS.at(iCat).find("pt_35_42") != string::npos) {
          pt[0] = 35;
          pt[1] = 42.5;
        } 
        else if (CAT_CUTS.at(iCat).find("pt_42_50") != string::npos) {
          pt[0] = 42.5;
          pt[1] = 50;
        } 
        else if (CAT_CUTS.at(iCat).find("pt_50_inf") != string::npos) {
          pt[0] = 50;
          pt[1] = 999;
        }
        else if (CAT_CUTS.at(iCat).find("inclusive") != string::npos) {
          pt[0] = 20;
          pt[1] = 999;
        }
        // std::cout << "pt cut value: " << pt[0] << " - " << pt[1] << std::endl;
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
          // std::cout << "eta cut value: " << eta[0] << " - " << eta[1] << std::endl;
          std::string h_pre = "h_"+OPT_CUTS.at(iOpt)+"_"+CAT_CUTS.at(iCat)+"_"+ETA_CUTS.at(iEta)+"_";

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
            // See if muon passes the d0 cuts
            if ( (mu.d0_PV * mu.charge > d0 - 0.00025) && (mu.d0_PV * mu.charge < d0 + 0.00025) ){
              //See if muon passes the eta cuts
              if (eta[0] < fabs(mu.eta) && fabs(mu.eta) < eta[1]){
                // See if muon passes the pT cuts
                // PF first
                if (pt[0] < mu.pt && mu.pt < pt[1]){
                  if (verbose) std::cout << "Gen pT = " << genMu.pt << ", eta = " << genMu.eta << ", phi = " << genMu.phi << std::endl;
                  if (verbose) std::cout << "Mu pT = " << mu.pt << ", eta = " << mu.eta << ", phi = " << mu.phi << ", d0 = " << mu.d0_PV << ", charge = " << mu.charge << std::endl;
                  BookAndFill( h_map_1D, h_pre+"d0", 400, -0.01,  0.01, mu.d0_PV*mu.charge,event_wgt );
                  // BookAndFill( h_map_1D, h_pre+"pt", 10000,    0, 1000, mu.pt );
                  // BookAndFill( h_map_1D, h_pre+"eta", 10000,    -10, 10, mu.eta );
                  // BookAndFill( h_map_1D, h_pre+"dPt", 400,    -4, 4,  mu.pt      - genMu.pt );
                  BookAndFill( h_map_1D, h_pre+"dRelPt2p0", 4000,    -40, 40,  10000*(mu.pt      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
                  // BookAndFill( h_map_1D, h_pre+"dRelPt1p0", 400,    -40, 40,  100*(mu.pt      - genMu.pt) / genMu.pt );
                  int eta_sign = genMu.eta > 0 ? 1 : -1;
                  // BookAndFill( h_map_1D, h_pre+"dPhi", 400, -1, 1, 1000*(mu.phi - genMu.phi)*mu.charge );
                  // BookAndFill( h_map_1D, h_pre+"dEta", 400, -1, 1, 1000*(mu.eta - genMu.eta)*mu.charge*eta_sign );
                }
                // // Rochester
                if ((pt[0] < mu.pt_Roch && mu.pt_Roch < pt[1])){
                  // BookAndFill( h_map_1D, h_pre+"d0_Roch", 200, -0.01,  0.01, mu.d0_PV*mu.charge );
                  // BookAndFill( h_map_1D, h_pre+"dPt_Roch", 400,    -4, 4,  mu.pt_Roch      - genMu.pt );
                  BookAndFill( h_map_1D, h_pre+"dRelPt2p0_Roch", 4000,    -40, 40,  10000*(mu.pt_Roch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
                  // BookAndFill( h_map_1D, h_pre+"dRelPt1p0_Roch", 400,    -40, 40,  100*(mu.pt_Roch      - genMu.pt) / genMu.pt );
                }
                // // Kalman
                if ((pt[0] < mu.pt_KaMu && mu.pt_KaMu < pt[1])){
                  // BookAndFill( h_map_1D, h_pre+"d0_KaMu", 200, -0.01,  0.01, mu.d0_PV*mu.charge );
                  // BookAndFill( h_map_1D, h_pre+"dPt_KaMu", 400,    -4, 4,  mu.pt_KaMu      - genMu.pt );
                  BookAndFill( h_map_1D, h_pre+"dRelPt2p0_KaMu", 4000,    -40, 40,  10000*(mu.pt_KaMu      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
                  // BookAndFill( h_map_1D, h_pre+"dRelPt1p0_KaMu", 400,    -40, 40,  100*(mu.pt_KaMu      - genMu.pt) / genMu.pt );
                }
                // // Kinfit
                if ((pt[0] < mu.pt_kinfit && mu.pt_kinfit < pt[1])){
                  // BookAndFill( h_map_1D, h_pre+"d0_kinfit", 200, -0.01,  0.01, mu.d0_PV*mu.charge );
                  // BookAndFill( h_map_1D, h_pre+"dPt_kinfit", 400,    -4, 4,  mu.pt_kinfit      - genMu.pt );
                  BookAndFill( h_map_1D, h_pre+"dRelPt2p0_kinfit", 4000,    -40, 40,  10000*(mu.pt_kinfit      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
                  // BookAndFill( h_map_1D, h_pre+"dRelPt1p0_kinfit", 400,    -40, 40,  100*(mu.pt_kinfit      - genMu.pt) / genMu.pt  );
                }
                // // Roch+Kinfit
                float mu_pt_KinRoch  = mu.pt_kinfit * mu.pt_Roch / mu.pt;
                if ((pt[0] < mu_pt_KinRoch && mu_pt_KinRoch < pt[1])){
                  // BookAndFill( h_map_1D, h_pre+"d0_KinRoch", 200, -0.01,  0.01, mu.d0_PV*mu.charge );
                  // BookAndFill( h_map_1D, h_pre+"dPt_KinRoch", 400,    -4, 4,  mu_pt_KinRoch      - genMu.pt );
                  BookAndFill( h_map_1D, h_pre+"dRelPt2p0_KinRoch", 4000,    -40, 40,  10000*(mu_pt_KinRoch      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
                  // BookAndFill( h_map_1D, h_pre+"dRelPt1p0_KinRoch", 400,    -40, 40,  100*(mu_pt_KinRoch      - genMu.pt) / genMu.pt );
                }
                // // KaMu+Kinfit
                float mu_pt_KinKaMu  = mu.pt_kinfit * mu.pt_KaMu / mu.pt;
                if ((pt[0] < mu_pt_KinKaMu && mu_pt_KinKaMu < pt[1])){
                  // BookAndFill( h_map_1D, h_pre+"d0_KinKaMu", 200, -0.01,  0.01, mu.d0_PV*mu.charge );
                  // BookAndFill( h_map_1D, h_pre+"dPt_KinKaMu", 400,    -4, 4,  mu_pt_KinKaMu      - genMu.pt );
                  BookAndFill( h_map_1D, h_pre+"dRelPt2p0_KinKaMu", 4000,    -40, 40,  10000*(mu_pt_KinKaMu      - genMu.pt) / pow(genMu.pt, 2.0),event_wgt );
                  // BookAndFill( h_map_1D, h_pre+"dRelPt1p0_KinKaMu", 400,    -40, 40,  100*(mu_pt_KinKaMu      - genMu.pt) / genMu.pt );
                }
              } // End eta cuts
            } // End d0 cuts
          }// End iMuon loop
        } // End iEta loop
      } // End iCat loop
    } // End iOpt loop
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
  // for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++) {
    // std::string optStr = OPT_CUTS.at(iOpt);

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
      // std::string h_name = it->second->GetName();
      // if (h_name.find(optStr+"_") != std::string::npos) {
	// Remove optional selection and category cuts from histogram names
	// h_name.erase( h_name.find(optStr+"_"), optStr.length() + 1 );
	// it->second->SetName(h_name.c_str());
	// it->second->SetTitle(h_name.c_str());
	     std::cout << "  * Writing 1D histogram " << it->second->GetName() << std::endl;
	     it->second->Write();
      // }
    }
    // Write out 2D histograms
    for (std::map<TString, TH2*>::iterator it = h_map_2D.begin(); it != h_map_2D.end(); ++it) {
          // std::string h_name = it->second->GetName();
          // if (h_name.find(optStr+"_") != std::string::npos) {
    	// Remove optional selection and category cuts from histogram names
    	// h_name.erase( h_name.find(optStr+"_"), optStr.length() + 1 );
    	// it->second->SetName(h_name.c_str());
    	// it->second->SetTitle(h_name.c_str());
    	std::cout << "  * Writing 2D histogram " << it->second->GetName() << std::endl;
    	it->second->Write();
          // }
    }
    
    out_file->Write();
    std::cout << "Wrote output file " << out_file_name.Data() << std::endl;
    
  // } // End loop: for (int iOpt = 0; iOpt < OPT_CUTS.size(); iOpt++)
      
  std::cout << "\nExiting GenRecoPtDiffVsD0VsPt()\n";
  
} // End void GenRecoPtDiffVsD0VsPt()
