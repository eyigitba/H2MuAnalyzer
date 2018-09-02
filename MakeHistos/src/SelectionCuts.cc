#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelections.h"
#include "H2MuAnalyzer/MakeHistos/interface/SelectionCuts.h"

bool PassSelection(NTupleBranches & br, std::string sel, bool verbose) {

  if (verbose) std::cout << "\nInside PassSelection, checking if event passes " << sel << std::endl;

  bool pass = false;
  
  if (sel.compare("NONE") == 0) pass = true;
  
  else if (sel.compare("Presel2017") == 0) {
    if (verbose) std::cout << "  * Applying Presel2017 cuts" << std::endl;
    bool passMu = false;
    bool hltmatched  = true; // if not using hlts

    for (int i = 0; i < br.nMuPairs; i++) {
      MuonInfo & Mu1 = br.muons->at(br.muPairs->at(i).iMu1);
      MuonInfo & Mu2 = br.muons->at(br.muPairs->at(i).iMu2);      

      if ( MuonPass(Mu1,30) and MuonPass(Mu2,20))  passMu = true;
      if (verbose) std::cout << "    - Mu1 pT = " << Mu1.pt 
			     << ", Mu2 pT = " << Mu2.pt << ", passMu = " << passMu << std::endl;

    } // End loop for (int i = 0; i < nMuPairs; i++)

    for (int iMu =0; iMu < br.nMuons; iMu++) {
      MuonInfo & Mu = br.muons->at(iMu) ;
      if ( MuonPass(Mu) and (Mu.isHltMatched[2] == 1 or Mu.isHltMatched[3] == 1 or Mu.isHltMatched[4] == 1 or Mu.isHltMatched[5] == 1 or Mu.isHltMatched[6] == 1 or Mu.isHltMatched[7] == 1) )  hltmatched = true;
      if (verbose) std::cout << "hltmatched = " << hltmatched << std::endl;
    } // End loop for (int iMu =0; iMu < br.nMuons; iMu++)
    
    if (passMu and hltmatched) pass = true;
    
  } // End if (sel.compare("Presel2017") == 0)
 

  else if (sel.compare("Dimu_sel_2017") == 0) {
    if (verbose) std::cout << "  * Applying Dimu_sel_2017 cut" << std::endl;
    bool passDimu = false;

    for (int i = 0; i < br.nMuPairs; i++) {
      if( br.muPairs->at(i).iMu1 == 0 and br.muPairs->at(i).iMu2 == 1 ) {
	MuPairInfo & Dimu = br.muPairs->at(i);
	if( Dimu.mass_Roch > 60 and MuonPass(br.muons->at(0),30) and MuonPass(br.muons->at(1),20) ) passDimu = true;  
      }
    } // End loop for (int i = 0; i < br.nMuPairs; i++)
    if (verbose) std::cout << "passDimu = " << passDimu << std::endl;
    if (passDimu) pass = true;
  } // End if (sel.compare("Dimu_sel_2017") == 0)

 
  else {
    std::cout << "\nInside SelectionCuts.cc, don't recognize selection " << sel << " - returning 'false'" << std::endl;
  }
  
  return pass;
  
}
