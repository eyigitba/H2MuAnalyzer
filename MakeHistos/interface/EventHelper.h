#ifndef EVENT_HELPER
#define EVENT_HELPER

#include <iostream>
#include <assert.h>

#include "TLorentzVector.h"

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"


int NumJets ( const ObjectSelectionConfig & cfg, const NTupleBranches & br, const std::string sel = "", const bool verbose = false );


#endif  // #ifndef EVENT_HELPER
