#include "H2MuAnalyzer/MakeHistos/interface/EventHelper.h"


// Return the number of selected jets, including optional cuts
int NumJets ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel, const bool verbose ) {

  int nJets = 0;
  for (const auto & jet : (*br.jets)) {
    if (JetPass(cfg, jet, (*br.muons), sel)) nJets += 1;
  }
  return nJets;
}
