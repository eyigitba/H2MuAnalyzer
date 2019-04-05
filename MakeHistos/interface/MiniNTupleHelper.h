
#ifndef MININTUPLE_HELPER
#define MININTUPLE_HELPER

#include <iostream>
#include <assert.h>

#include "TTree.h"
#include "TBranch.h"

TTree * PlantTree(TString tree_name, TString tree_title);

void BookBranch(std::map<TString, float>       & b_map_flt, TTree * out_tree, const TString name );
void BookBranch(std::map<TString, int>         & b_map_int, TTree * out_tree, const TString name );
void BookBranch(std::map<TString, std::string> & b_map_str, TTree * out_tree, const TString name );

#endif  // #ifndef HISTO_HELPER
