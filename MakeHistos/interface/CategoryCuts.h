
#ifndef CATEGORY_CUTS
#define CATEGORY_CUTS

#include <iostream>

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"

bool InCategory(NTupleBranches & br, std::string sel, bool verbose = false);

#endif  // #ifndef CATEGORY_CUTS
