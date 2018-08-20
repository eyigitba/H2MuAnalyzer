
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"

bool InCategory(NTupleBranches & br, std::string sel, bool verbose) {
  
  if (verbose) std::cout << "\nInside InCategory, checking if event passes " << sel << std::endl;
  
  bool pass = false;

  if (sel.compare("NONE") == 0) pass = true;
  
  else if (sel.compare("Mu1Pos") == 0) {
    if (verbose) std::cout << "  * Applying Mu1Pos cuts" << std::endl;

    for (int i = 0; i < br.nMuPairs; i++) {
      int iMu1 = br.muPairs->at(i).iMu1;

      if (br.muons->at(iMu1).charge == 1) pass = true;
      if (verbose) std::cout << "    - Mu1 charge = " << br.muons->at(iMu1).charge << ", pass = " << pass << std::endl;
      break; // Only look at the first muon pair
    } // End loop for (int i = 0; i < nMuPairs; i++)
  } // End if (sel.compare("Mu1Pos") == 0)
  
  else if (sel.compare("Mu1Neg") == 0) {
    if (verbose) std::cout << "  * Applying Mu1Neg cuts" << std::endl;
    
    for (int i = 0; i < br.nMuPairs; i++) {
      int iMu1 = br.muPairs->at(i).iMu1;
      
      if (br.muons->at(iMu1).charge == -1) pass = true;
      if (verbose) std::cout << "    - Mu1 charge = " << br.muons->at(iMu1).charge << ", pass = " << pass << std::endl;
      break; // Only look at the first muon pair
    } // End loop for (int i = 0; i < nMuPairs; i++)
  } // End if (sel.compare("Mu1Neg") == 0)

  else {
    std::cout << "\nInside CategoryCuts.cc, don't recognize category " << sel << " - returning 'false'" << std::endl;
  }

  return pass;

}
