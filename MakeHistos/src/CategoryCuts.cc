#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelections.h"
#include "H2MuAnalyzer/MakeHistos/interface/CategoryCuts.h"
#include "TLorentzVector.h"

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

  else if (sel.compare("WHlep") == 0) {  // real toy category, needs more study on cuts
    if (verbose) std::cout << "  * Applying WHlep cuts" << std::endl;
    
    if (br.nBTight == 0 and (br.met)->pt > 20) {
      if (br.nMuons == 3) {
	pass = true;
	for (int i = 0; i < br.nMuons; i++) {
	  if ( !MuonPass(br.muons->at(i),20) or abs(br.muons->at(i).d0_PV)>0.02 ) pass = false;
 	}
      }//end if (br.nMuons == 3)
      if (br.nMuons == 2 and br.nEles == 1) {
	if (br.eles->at(0).pt > 20 and abs(br.eles->at(0).d0_PV)<0.02 ) pass = true;
      }
    }//end if (br.nBTight == 0 and (br.met)->pt > 30)
  }//end if (sel.compare("WHlep") == 0)

  else if (sel.compare("WHhad2J") == 0) {
    if (verbose) std::cout << "  * Applying WHhad2J cuts" << std::endl;
    
    if (not InCategory(br, "WHlep", verbose) ) {
	for (int i = 0; i < br.nJetPairs; i++) {
	  JetPairInfo & dijet = br.jetPairs->at(i);
	  if ( DijetPass(br, dijet) and abs(dijet.mass-80)<20 and abs(dijet.eta)<4 and dijet.dR<4 ) pass = true;
	}//end for (int i = 0; i < br.nJetPairs; i++)
    }
  }//end if (sel.compare("WHhad2J") == 0) 

  else if (sel.compare("WHhad1J") == 0) {
    if (verbose) std::cout << "  * Applying WHhad1J cuts" << std::endl;

    if (not InCategory(br, "WHlep", verbose) ) {
      if (not InCategory(br, "WHhad2J", verbose) ) {
	for (int i = 0; i < br.nJets; i++) {
	   JetInfo & jet = br.jets->at(i);
	   if ( JetPass(jet) and abs(jet.eta)<2 and not BJetPass(jet) ) pass = true;
	}// end for (int i = 0; i < br.nJets; i++) {
      }
    }
  }//end if (sel.compare("WHhad1J") == 0)


  else if (sel.compare("ZHmu") == 0) {  // real toy category, needs more study on cuts
    if (verbose) std::cout << "  * Applying ZHmu cuts" << std::endl;

    if (br.nMuons == 4) {
	for (int i = 0; i < br.nMuPairs; i ++) {
	  if (br.muPairs->at(i).iMu1 != 0 and br.muPairs->at(i).iMu1 != 1 and br.muPairs->at(i).iMu2 != 0 and br.muPairs->at(i).iMu2 != 1) {
	    if (br.muPairs->at(i).mass_Roch > 80 and br.muPairs->at(i).mass_Roch < 100) pass = true;
	  }
	} 
    }
  }// end if (sel.compare("ZHmu") == 0)

  else if (sel.compare("ZHele") == 0) {   // real toy category, needs more study on cuts
    if (verbose) std::cout << "  * Applying ZHele cuts" << std::endl;

    if (br.nMuons == 2 and br.nEles ==2) {
      EleInfo & ele1 = br.eles->at(0);
      EleInfo & ele2 = br.eles->at(1);
      if(ele1.pt > 20 and ele2.pt > 20 and ele1.charge + ele2.charge == 0) {
    	TLorentzVector v1, v2, vz;
    	v1.SetPtEtaPhiM( ele1.pt, ele1.eta, ele1.phi, 0.0005);
    	v2.SetPtEtaPhiM( ele2.pt, ele2.eta, ele2.phi, 0.0005);
    	vz = v1 + v2;
        std::cout << vz.M() << std::endl;
    	if( vz.M()>80 and vz.M()<100 ) pass = true;
      }
    }
  } // end if (sel.compare("ZHele") == 0)

  else {
    std::cout << "\nInside CategoryCuts.cc, don't recognize category " << sel << " - returning 'false'" << std::endl;
  }

  return pass;

}
