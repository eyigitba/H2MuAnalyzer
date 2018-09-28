#ifndef OBJECT_HELPER
#define OBJECT_HELPER

#include <iostream>
#include <assert.h>

#include "TLorentzVector.h"

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"

bool  MuonID     ( const MuonInfo & muon,     const std::string muon_ID ); // Return loose, medium, or tight muon ID
bool  MuonTrig   ( const MuonInfo & muon,     const std::string year    ); // Return if muon fired HLT trigger
float MuonPt     ( const MuonInfo & muon,     const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected muon pT
float MuPairMass ( const MuPairInfo & muPair, const std::string pt_corr ); // Return PF, Rochester, or Kalman corrected dimuon invariant mass

bool JetPUID ( const JetInfo & jet, const std::string pu_ID, const std::string year ); // Return loose, medium, or tight jet PU ID from 2016 or 2017

TLorentzVector FourVec( const MuonInfo & muon, const float pT );
TLorentzVector FourVec( const EleInfo & ele );
TLorentzVector FourVec( const JetInfo & jet );
TLorentzVector FourVec( const MetInfo & met );
TLorentzVector FourVec( const MhtInfo & mht );
TLorentzVector FourVec( const GenParentInfo & genPar );
TLorentzVector FourVec( const GenMuonInfo & genMu );
TLorentzVector FourVec( const GenJetInfo & genJet );


#endif  // #ifndef OBJECT_HELPER
