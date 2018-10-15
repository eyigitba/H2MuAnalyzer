
#ifndef LOAD_NTUPLE_BRANCHES
#define LOAD_NTUPLE_BRANCHES

#include <iostream>
#include <assert.h>

#include "TChain.h"

#include "Ntupliser/DiMuons/interface/NTupleBranches.h"

struct NTupleBranches {
  EventInfo * event = 0;
  VertexInfos * vertices = 0;
  
  MuonInfos * muons = 0;
  MuPairInfos * muPairs = 0;
  
  EleInfos * eles = 0;
  
  JetInfos * jets = 0;
  JetInfos * jets_JES_up = 0;
  JetInfos * jets_JES_down = 0;
  SlimJetInfos * slimJets = 0;
  SlimJetInfos * slimJets_JES_up = 0;
  SlimJetInfos * slimJets_JES_down = 0;
  
  JetPairInfos * jetPairs = 0;
  JetPairInfos * jetPairs_JES_up = 0;
  JetPairInfos * jetPairs_JES_down = 0;
  
  MetInfo * met = 0;
  MetInfo * met_JES_up = 0;
  MetInfo * met_JES_down = 0;
  
  MhtInfo * mht = 0;
  MhtInfo * mht_JES_up = 0;
  MhtInfo * mht_JES_down = 0;

  GenParentInfos * genParents = 0;
  GenMuonInfos * genMuons = 0;
  GenMuPairInfos * genMuPairs = 0;
  GenJetInfos * genJets = 0;
  
  int nVertices = -99;
  int nPU = -99;
  float LHE_HT = -99;
  int nMuons = -99;
  int nMuPairs = -99;
  int nEles = -99;
  int nJets = -99;
  int nJetPairs = -99;
  int nJetsCent = -99;
  int nJetsFwd = -99;
  int nBLoose = -99;
  int nBMed = -99;
  int nBTight = -99;
  int nJets_JES_up = -99;
  int nJetsCent_JES_up = -99;
  int nBLoose_JES_up = -99;
  int nBMed_JES_up = -99;
  int nBTight_JES_up = -99;
  int nJetsFwd_JES_up = -99;
  int nJets_JES_down = -99;
  int nJetsCent_JES_down = -99;
  int nJetsFws_JES_down = -99;
  int nBLoose_JES_down = -99;
  int nBMed_JES_down = -99;
  int nBTight_JES_down = -99;
  int nGenParents = -99;
  int nGenMuons = -99;
  int nGenMuPairs = -99;
  int nGenJets = -99;
  
  std::vector<std::string> hltPaths = {};
  std::string btagName = "INVALID";
  
  int Flag_all = -99;
  int Flag_badMu = -99;
  int Flag_dupMu = -99;
  int Flag_halo = -99;
  int Flag_PV = -99;
  int Flag_HBHE = -99;
  int Flag_HBHE_Iso = -99;
  int Flag_ECAL_TP = -99;
  int Flag_BadChCand = -99;
  int Flag_eeBadSc = -99;
  int Flag_ecalBadCalib = -99;
  
  float IsoMu_eff_3 = -99;
  float IsoMu_eff_3_up = -99;
  float IsoMu_eff_3_down = -99;
  float IsoMu_eff_4 = -99; // 2016 only
  float IsoMu_eff_4_up = -99; // 2016 only
  float IsoMu_eff_4_down = -99; // 2016 only
  float IsoMu_eff_bug = -99;
  float IsoMu_eff_bug_up = -99;
  float IsoMu_eff_bug_down = -99;

  float MuID_eff_3 = -99;
  float MuID_eff_3_up = -99;
  float MuID_eff_3_down = -99;
  float MuID_eff_4 = -99; // 2016 only
  float MuID_eff_4_up = -99; // 2016 only
  float MuID_eff_4_down = -99; // 2016 only

  float MuIso_eff_3 = -99;
  float MuIso_eff_3_up = -99;
  float MuIso_eff_3_down = -99;
  float MuIso_eff_4 = -99; // 2016 only
  float MuIso_eff_4_up = -99; // 2016 only
  float MuIso_eff_4_down = -99; // 2016 only

  float IsoMu_SF_3 = -99;
  float IsoMu_SF_3_up = -99;
  float IsoMu_SF_3_down = -99;
  float IsoMu_SF_4 = -99; // 2016 only
  float IsoMu_SF_4_up = -99; // 2016 only
  float IsoMu_SF_4_down = -99; // 2016 only
  float IsoMu_SF_bug = -99;
  float IsoMu_SF_bug_up = -99;
  float IsoMu_SF_bug_down = -99;

  float MuID_SF_3 = -99;
  float MuID_SF_3_up = -99;
  float MuID_SF_3_down = -99;
  float MuID_SF_4 = -99;
  float MuID_SF_4_up = -99;
  float MuID_SF_4_down = -99;

  float MuIso_SF_3 = -99;
  float MuIso_SF_3_up = -99;
  float MuIso_SF_3_down = -99;
  float MuIso_SF_4 = -99; // 2016 only
  float MuIso_SF_4_up = -99; // 2016 only
  float MuIso_SF_4_down = -99; // 2016 only

  float PU_wgt = -99;
  float PU_wgt_up = -99;
  float PU_wgt_down = -99;

  int GEN_wgt = -99;
};


void SetBranchAddresses(TChain & ch_, NTupleBranches & br_, std::vector<std::string> opts = {}, bool verbose = false);

JetInfos ConvertSlimJets(SlimJetInfos & _slimJets);


#endif  // #ifndef LOAD_NTUPLE_BRANCHES
