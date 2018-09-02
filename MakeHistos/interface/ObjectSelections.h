#ifndef OBJECT_SELECTIONS
#define OBJECT_SELECTIONS

#include <iostream>
#include <cmath>

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"

bool MuonPass	(MuonInfo & muon, int pt = 20, std::string id = "medium", bool verbose = false);
bool JetPass	(JetInfo & jet  , int pt = 30, std::string PU_ID = "medium", bool verbose = false);

#endif  // #ifndef OBJECT_SELECTIONS
