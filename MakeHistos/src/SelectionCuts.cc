
#include "H2MuAnalyzer/MakeHistos/interface/SelectionCuts.h"

bool PassSelection(NTupleBranches & br, std::string sel, bool verbose) {

  if (verbose) std::cout << "\nInside PassSelection, checking if event passes " << sel << std::endl;

  bool pass = false;
  
  if (sel.compare("NONE") == 0) pass = true;
  
  else if (sel.compare("Presel2017") == 0) {
    if (verbose) std::cout << "  * Applying Presel2017 cuts" << std::endl;
    bool passMu = false;
    
    for (int i = 0; i < br.nMuPairs; i++) {
      int iMu1 = br.muPairs->at(i).iMu1;
      int iMu2 = br.muPairs->at(i).iMu2;
      
      if (br.muons->at(iMu1).pt > 30 && br.muons->at(iMu2).pt > 20) passMu = true;
      if (verbose) std::cout << "    - Mu1 pT = " << br.muons->at(iMu1).pt 
			     << ", Mu2 pT = " << br.muons->at(iMu2).pt << ", passMu = " << passMu << std::endl;
    } // End loop for (int i = 0; i < nMuPairs; i++)
    
    if (passMu) pass = true;
    
  } // End if (sel.compare("Presel2017") == 0)
  
  else {
    std::cout << "\nInside SelectionCuts.cc, don't recognize selection " << sel << " - returning 'false'" << std::endl;
  }
  
  return pass;
  
}
