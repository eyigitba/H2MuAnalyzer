#ifndef OBJECT_HELPER
#define OBJECT_HELPER

#include <iostream>
#include <assert.h>

#include "TLorentzVector.h"

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"

bool  MuonID        ( const MuonInfo & muon,     const std::string muon_ID ); // Return loose, medium, or tight muon ID
bool  MuonTrig      ( const MuonInfo & muon,     const std::string year    ); // Return if muon fired HLT trigger
float MuonPt        ( const MuonInfo & muon,     const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected muon pT
float MuPairPt      ( const MuPairInfo & muPair, const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected dimuon pT
float MuPairMass    ( const MuPairInfo & muPair, const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected dimuon invariant mass
float MuPairMassErr ( const MuPairInfo & muPair, const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected dimuon invariant mass uncertainty
bool  IsGenMatched  ( const MuPairInfo & muPair, const MuonInfos & muons, const GenMuonInfos & genMuons, const std::string gen_ID ); // Match di-muon pairs to GEN Higgs or Z

bool  EleID ( const EleInfo & ele, const std::string ele_ID ); // Return loose, medium, or tight electron ID

bool JetPUID ( const JetInfo & jet, const std::string pu_ID, const std::string year ); // Return loose, medium, or tight jet PU ID from 2016 or 2017


TLorentzVector FourVec( const MuonInfo & muon, const std::string pt_corr, const std::string opt = "" );
TLorentzVector FourVec( const MuPairInfo & muPair, const std::string pt_corr, const std::string opt = "" );
TLorentzVector FourVec( const EleInfo & ele, const std::string opt = "" );
TLorentzVector FourVec( const JetInfo & jet, const std::string opt = "" );
TLorentzVector FourVec( const MetInfo & met, const std::string opt = "" );
TLorentzVector FourVec( const MhtInfo & mht, const std::string opt = "" );
TLorentzVector FourVec( const GenParentInfo & genPar, const std::string opt = "" );
TLorentzVector FourVec( const GenMuonInfo & genMu, const std::string opt = "" );
TLorentzVector FourVec( const GenJetInfo & genJet, const std::string opt = "" );


#endif  // #ifndef OBJECT_HELPER