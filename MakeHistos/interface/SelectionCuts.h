
#ifndef SELECTION_CUTS
#define SELECTION_CUTS

#include <iostream>

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"

bool PassSelection(NTupleBranches & br, std::string sel, bool verbose = false);

#endif  // #ifndef SELECTION_CUTS
