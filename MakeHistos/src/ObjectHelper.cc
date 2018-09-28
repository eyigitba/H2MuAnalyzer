#include "H2MuAnalyzer/MakeHistos/interface/ObjectHelper.h"


////////////////////////
///  Muon functions  ///
////////////////////////

// Return loose, medium, or tight muon ID
bool MuonID ( const MuonInfo & muon, const std::string muon_ID ) {
  if ( muon_ID == "loose" ) return muon.isLooseID;
  if ( muon_ID == "medium") return muon.isMediumID;
  if ( muon_ID == "tight" ) return muon.isTightID;
  std::cout << "\n\nInside ObjectSelections.cc, invalid option muon_ID = " << muon_ID << std::endl;
  assert(false);
}

// Return if muon fired HLT trigger
bool MuonTrig ( const MuonInfo & muon, const std::string year) {
  if        ( year == "2016" ) {
    // [2] = HLT_IsoMu24, [3] = HLT_IsoTkMu24
    return ( muon.isHltMatched[2] || muon.isHltMatched[3] );
  } else if ( year == "2017" ) {
    // [2] = HLT_IsoMu27, [3] = HLT_IsoTkMu27
    return ( muon.isHltMatched[2] || muon.isHltMatched[3] );
  }
  std::cout << "\n\nInside ObjectSelections.cc, invalid option year = " << year << std::endl;
  assert(false);
}

// Return PF, Rochester, or Kalman corrected muon pT
float MuonPt ( const MuonInfo & muon, const std::string pt_corr ) {
  if ( pt_corr == "PF"   ) return muon.pt;
  if ( pt_corr == "Roch" ) return muon.pt_Roch;
  if ( pt_corr == "KaMu" ) return muon.pt_KaMu;
  std::cout << "\n\nInside ObjectSelections.cc, invalid option pt_corr = " << pt_corr << std::endl;
  assert(false);
}

// Return PF, Rochester, or Kalman corrected dimuon invariant mass
float MuPairMass ( const MuPairInfo & muPair, const std::string pt_corr ) {
  if ( pt_corr == "PF"   ) return muPair.mass;
  if ( pt_corr == "Roch" ) return muPair.mass_Roch;
  if ( pt_corr == "KaMu" ) return muPair.mass_KaMu;
  std::cout << "\n\nInside ObjectSelections.cc, invalid option pt_corr = " << pt_corr << std::endl;
  assert(false);
}


///////////////////////
///  Jet functions  ///
///////////////////////

 // Return loose, medium, or tight jet PU ID from 2016 or 2017
bool JetPUID ( const JetInfo & jet, const std::string PU_ID, const std::string year ) {

  float puID_cut = 999;  // Minimum puID value to pass cut

  if (year == "2016") return true;
  
  else if (year == "2017") {

    if (jet.pt >= 50) return true;

    // Where do these numbers come from?!?  No recommendation on the following twiki - AWB 24.09.2018
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
    if ( PU_ID == "loose" ) {
      if      ( abs(jet.eta) < 2.50 ) {
	puID_cut = (jet.pt < 30 ? -0.97 : -0.89); }
      else if ( abs(jet.eta) < 2.75 ) {
	puID_cut = (jet.pt < 30 ? -0.68 : -0.52); }
      else if ( abs(jet.eta) < 3.00 ) {
	puID_cut = (jet.pt < 30 ? -0.53 : -0.38); }
      else if ( abs(jet.eta) < 5.00 ) {
	puID_cut = (jet.pt < 30 ? -0.47 : -0.30); }
      else { std::cout << "Inside JetPUID, invalid jet eta = " << jet.eta << std::endl; assert(false); }
    }
    else if ( PU_ID == "medium" ) {
      if      ( abs(jet.eta) < 2.50 ) {
	puID_cut = (jet.pt < 30 ? +0.18 : +0.61); }
      else if ( abs(jet.eta) < 2.75 ) {
	puID_cut = (jet.pt < 30 ? -0.55 : -0.35); }
      else if ( abs(jet.eta) < 3.00 ) {
	puID_cut = (jet.pt < 30 ? -0.42 : -0.23); }
      else if ( abs(jet.eta) < 5.00 ) {
	puID_cut = (jet.pt < 30 ? -0.36 : -0.17); }
      else { std::cout << "Inside JetPUID, invalid jet eta = " << jet.eta << std::endl; assert(false); }
    }
    else if ( PU_ID == "tight" ) {
      if      ( abs(jet.eta) < 2.50 ) {
	puID_cut = (jet.pt < 30 ? +0.69 : +0.86); }
      else if ( abs(jet.eta) < 2.75 ) {
	puID_cut = (jet.pt < 30 ? -0.35 : -0.10); }
      else if ( abs(jet.eta) < 3.00 ) {
	puID_cut = (jet.pt < 30 ? -0.42 : -0.23); }
      else if ( abs(jet.eta) < 5.00 ) {
	puID_cut = (jet.pt < 30 ? -0.36 : -0.17); }
      else { std::cout << "Inside JetPUID, invalid jet eta = " << jet.eta << std::endl; assert(false); }
    }

  } // End conditional: if (year == "2017")

  else { std::cout << "Inside JetPUID, invalid year = " << year << std::endl; assert(false); }

  return (jet.puID >= puID_cut);

} // End function: bool JetPUID ()


/////////////////////////////
///  Kinematic functions  ///
/////////////////////////////

TLorentzVector FourVec( const MuonInfo & muon, const float pT ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(pT, muon.eta, muon.phi, 0.105658367 );
  return vec;
}
TLorentzVector FourVec( const EleInfo & ele ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(ele.pt, ele.eta, ele.phi, 0.000511 );
  return vec;
}
TLorentzVector FourVec( const JetInfo & jet ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass );
  return vec;
}
TLorentzVector FourVec( const MetInfo & met ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(met.pt, 0, met.phi, 0 );
  return vec;
}
TLorentzVector FourVec( const MhtInfo & mht ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(mht.pt, 0, mht.phi, 0 );
  return vec;
}
TLorentzVector FourVec( const GenParentInfo & genPar ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(genPar.pt, genPar.eta, genPar.phi, genPar.mass );
  return vec;
}
TLorentzVector FourVec( const GenMuonInfo & genMu ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(genMu.pt, genMu.eta, genMu.phi, 0.105658367 );
  return vec;
}
TLorentzVector FourVec( const GenJetInfo & genJet ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(genJet.pt, genJet.eta, genJet.phi, genJet.mass );
  return vec;
}




