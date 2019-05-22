#ifndef OBJECT_HELPER
#define OBJECT_HELPER

#include <iostream>
#include <assert.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH2.h"

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"

bool  MuonID        ( const MuonInfo & muon, const std::string muon_ID );     // Return loose, medium, or tight muon ID
bool  MuonTrig      ( const MuonInfo & muon, const std::string year, const std::vector<std::string> trigNames ); // Return if muon fired HLT trigger
float MuonPt        ( const MuonInfo & muon, const std::string pt_corr );     // Return PF, Rochester, or Kalman corrected muon pT
float MuPairPt      ( const MuPairInfo & muPair, const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected dimuon pT
float MuPairMass    ( const MuPairInfo & muPair, const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected dimuon invariant mass
float MuPairMassErr ( const MuPairInfo & muPair, const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected dimuon invariant mass uncertainty
TH2F* LoadSFsLepMVA ( const std::string year, const std::string flavor, const std::string WP ); // Return a 2D histogram with LepMVA efficiency scale factors
float LepMVASF      ( const TH2F * h_SF, const float pt, const float eta );   // Return the LepMVA efficiency scale factor for a single lepton
bool  IsGenMatched  ( const MuPairInfo & muPair, const MuonInfos & muons, const GenMuonInfos & genMuons, const std::string gen_ID ); // Match di-muon pairs to GEN Higgs or Z

bool  EleID ( const EleInfo & ele, const std::string ele_ID ); // Return loose, medium, or tight electron ID

float GetLepMVASF(const std::string lep_type, float pt, float eta, float lepMVA_cut); // return lepMVA SF, lep_type is "muon" or "ele" 

bool JetPUID ( const JetInfo & jet, const std::string pu_ID, const std::string year ); // Return loose, medium, or tight jet PU ID from 2016 or 2017

JetPairInfo MakeJetPair( TLorentzVector jet1_vec, TLorentzVector jet2_vec );

TLorentzVector FourVec( const MuonInfo & muon, const std::string pt_corr, const std::string opt = "" );
TLorentzVector FourVec( const MuPairInfo & muPair, const std::string pt_corr, const std::string opt = "", const MuonInfos & muons = {} );
TLorentzVector FourVec( const EleInfo & ele, const std::string opt = "" );
TLorentzVector FourVec( const JetInfo & jet, const std::string opt = "" );
TLorentzVector FourVec( const MetInfo & met, const std::string opt = "" );
TLorentzVector FourVec( const MhtInfo & mht, const std::string opt = "" );
TLorentzVector FourVec( const GenParentInfo & genPar, const std::string opt = "" );
TLorentzVector FourVec( const GenMuonInfo & genMu, const std::string opt = "" );
TLorentzVector FourVec( const GenJetInfo & genJet, const std::string opt = "" );
TLorentzVector FourVec( const MuonInfos & muons, const std::string pt_corr, const EleInfos & eles, const JetInfos & jets, const std::string opt = "" );

float CosThetaStar( TLorentzVector vec1, TLorentzVector vec2 );
float Cos_CS_Theta( TLorentzVector vec1, TLorentzVector vec2 );
float Cos_CS_Phi  ( TLorentzVector vec1, TLorentzVector vec2 );
float Sin_CS_Phi  ( TLorentzVector vec1, TLorentzVector vec2 );
float SignedDEta  ( TLorentzVector vec1, TLorentzVector vec2 );


#endif  // #ifndef OBJECT_HELPER
