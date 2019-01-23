
#include "H2MuAnalyzer/MakeHistos/interface/MiniNTupleHelper.h"

TTree * PlantTree(TString tree_name, TString tree_title) {
  return new TTree(tree_name, tree_title);
}

