#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelections.h"
#include "H2MuAnalyzer/MakeHistos/interface/SelectionCuts.h"

bool PassSelection(NTupleBranches & br, std::string sel, bool verbose) {

  if (verbose) std::cout << "\nInside PassSelection, checking if event passes " << sel << std::endl;

  bool pass = false;
  
  if (sel.compare("NONE") == 0) pass = true;
  
  // Added dimuon selection into preselection, instead of a stand alone selection - XWZ 20.09.2018
  else if (sel.compare("Presel2017") == 0 || sel.compare("Presel2018") == 0) {
    if (verbose) std::cout << "  * Applying Presel2017/2018 cuts" << std::endl;
    bool passDimu = false;
    bool hltmatched  = false; 
    bool passflag = false;

    for (int i = 0; i < br.nMuPairs; i++) {
      if ( DimuPass(br, br.muPairs->at(i)) ) passDimu = true;
      if (verbose) std::cout  << " passDimu = " << passDimu << std::endl;

    } // End loop for (int i = 0; i < nMuPairs; i++)

    std::string year = (sel.compare("Presel2017") == 0 ? "2017" : "2018");

    for (int iMu =0; iMu < br.nMuons; iMu++) {
      MuonInfo & Mu = br.muons->at(iMu) ;
      hltmatched = ( MuonPass(Mu) and MuonTrig(Mu, year) );
      if (verbose) std::cout << "hltmatched = " << hltmatched << std::endl;
    } // End loop for (int iMu =0; iMu < br.nMuons; iMu++)

    if (br.Flag_all == 1 ) passflag = true;
    if (passDimu and hltmatched and passflag) pass = true;

  } // End if (sel.compare("Presel2017") == 0 || sel.compare("Presel2018") == 0)
 

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
