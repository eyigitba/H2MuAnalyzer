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

// Return PF, Rochester, or Kalman corrected dimuon invariant pT
float MuPairPt ( const MuPairInfo & muPair, const std::string pt_corr ) {
  if ( pt_corr == "PF"   ) return muPair.pt;
  if ( pt_corr == "Roch" ) return muPair.pt_Roch;
  if ( pt_corr == "KaMu" ) return muPair.pt_KaMu;
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

// Return PF, Rochester, or Kalman corrected dimuon invariant mass uncertainty
float MuPairMassErr ( const MuPairInfo & muPair, const std::string pt_corr ) {
  if ( pt_corr == "PF"   ) return muPair.massErr;
  if ( pt_corr == "Roch" ) return muPair.massErr_Roch;
  if ( pt_corr == "KaMu" ) return muPair.massErr_KaMu;
  std::cout << "\n\nInside ObjectSelections.cc, invalid option pt_corr = " << pt_corr << std::endl;
  assert(false);
}

// Determine if dimuon pair is matched to GEN pair
bool IsGenMatched( const MuPairInfo & muPair, const MuonInfos & muons, const GenMuonInfos & genMuons, const std::string gen_ID ) {

  TLorentzVector mu_vec1 = FourVec( muons.at(muPair.iMu1), "PF" );
  TLorentzVector mu_vec2 = FourVec( muons.at(muPair.iMu2), "PF" );
  bool mu1_matched = false;
  bool mu2_matched = false;

  for (const auto & genMu : genMuons) {
    if      (gen_ID == "Z") { if (genMu.mother_ID != 23) continue; }
    else if (gen_ID == "H") { if (genMu.mother_ID != 25) continue; }
    else assert(gen_ID == "Z" || gen_ID == "H");

    TLorentzVector gen_vec = FourVec( genMu );
    if ( mu_vec1.DeltaR(gen_vec) < 0.05 ) mu1_matched = true;
    if ( mu_vec2.DeltaR(gen_vec) < 0.05 ) mu2_matched = true;
  }

  return (mu1_matched && mu2_matched);
} // End function: bool IsMatchedToGen()


////////////////////////////
///  Electron functions  ///
////////////////////////////

// Return loose, medium, or tight electron ID
bool EleID ( const EleInfo & ele, const std::string ele_ID ) {
  if ( ele_ID == "loose" ) return ele.isLooseID;
  if ( ele_ID == "medium") return ele.isMediumID;
  if ( ele_ID == "tight" ) return ele.isTightID;
  std::cout << "\n\nInside ObjectSelections.cc, invalid option ele_ID = " << ele_ID << std::endl;
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

TLorentzVector FourVec( const MuonInfo & muon, const std::string pt_corr, const std::string opt ) {
  TLorentzVector vec;
  if (opt == "T")
    vec.SetPtEtaPhiM(MuonPt(muon, pt_corr),        0, muon.phi, 0.105658367 );
  else
    vec.SetPtEtaPhiM(MuonPt(muon, pt_corr), muon.eta, muon.phi, 0.105658367 );
  return vec;
}
TLorentzVector FourVec( const MuPairInfo & muPair, const std::string pt_corr, const std::string opt, const MuonInfos & muons ) {
  TLorentzVector vec;
  if (opt == "T") {
    TLorentzVector mu1_vec = FourVec( muons.at(muPair.iMu1), pt_corr, "T");
    TLorentzVector mu2_vec = FourVec( muons.at(muPair.iMu2), pt_corr, "T");
    vec = mu1_vec + mu2_vec;
  }
  else
    vec.SetPtEtaPhiM(MuPairPt(muPair, pt_corr), muPair.eta, muPair.phi, MuPairMass(muPair, pt_corr) );
  return vec;
}
TLorentzVector FourVec( const EleInfo & ele, const std::string opt ) {
  TLorentzVector vec;
  if (opt == "T")
    vec.SetPtEtaPhiM(ele.pt,       0, ele.phi, 0.000511 );
  else
    vec.SetPtEtaPhiM(ele.pt, ele.eta, ele.phi, 0.000511 );
  return vec;
}
TLorentzVector FourVec( const JetInfo & jet, const std::string opt ) {
  TLorentzVector vec;
  if (opt == "T")
    vec.SetPtEtaPhiM(jet.pt,       0, jet.phi, jet.mass );
  else
    vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass );
  return vec;
}
TLorentzVector FourVec( const MetInfo & met, const std::string opt ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(met.pt, 0, met.phi, 0 );
  return vec;
}
TLorentzVector FourVec( const MhtInfo & mht, const std::string opt ) {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(mht.pt, 0, mht.phi, 0 );
  return vec;
}
TLorentzVector FourVec( const GenParentInfo & genPar, const std::string opt ) {
  TLorentzVector vec;
  if (opt == "T")
    vec.SetPtEtaPhiM(genPar.pt,          0, genPar.phi, genPar.mass );
  else
    vec.SetPtEtaPhiM(genPar.pt, genPar.eta, genPar.phi, genPar.mass );
  return vec;
}
TLorentzVector FourVec( const GenMuonInfo & genMu, const std::string opt ) {
  TLorentzVector vec;
  if (opt == "T")
    vec.SetPtEtaPhiM(genMu.pt,         0, genMu.phi, 0.105658367 );
  else
    vec.SetPtEtaPhiM(genMu.pt, genMu.eta, genMu.phi, 0.105658367 );
  return vec;
}
TLorentzVector FourVec( const GenJetInfo & genJet, const std::string opt ) {
  TLorentzVector vec;
  if (opt == "T")
    vec.SetPtEtaPhiM(genJet.pt,          0, genJet.phi, genJet.mass );
  else
    vec.SetPtEtaPhiM(genJet.pt, genJet.eta, genJet.phi, genJet.mass );
  return vec;
}


float CosThetaStar( TLorentzVector vec1, TLorentzVector vec2 ) {
  TLorentzVector parent_vec = vec1 + vec2;
  TVector3 parent_p = parent_vec.BoostVector(), p1, p2;
  float cos_theta_star; // default to return the cos_theta_star of vec1, cannot tell charge from TLorentzVector.  - XWZ 23.10.2018

  vec1.Boost( -parent_p );
  p1 = vec1.BoostVector();

  cos_theta_star = parent_p * p1 / (parent_p.Mag() * p1.Mag());
  return cos_theta_star;
}



