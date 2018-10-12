
#ifndef CATEGORY_CUTS
#define CATEGORY_CUTS

#include <iostream>

#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"
#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelection.h"

bool InCategory( const ObjectSelectionConfig & cfg, const NTupleBranches & br,
		 const std::string sel, const bool verbose = false);

#endif  // #ifndef CATEGORY_CUTS
