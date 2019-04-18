
#include "H2MuAnalyzer/MakeHistos/interface/MiniNTupleHelper.h"

TTree * PlantTree(TString tree_name, TString tree_title) {
  return new TTree(tree_name, tree_title);
}

// Book a branch for floating point values ("F")
void BookBranch(std::map<TString, float> & b_map_flt, TTree * out_tree, const TString b_name ) {

  if (b_map_flt.find(b_name) == b_map_flt.end()) {
    b_map_flt[b_name] = -99.0;
  }
  else {
    std::cout << "Why is there already a floating point branch booked with name " << b_name << "?!?  QUITING!!!" << std::endl;
    assert(false);
  }

  TString branch_str;
  branch_str.Form("%s/F", b_name.Data());
  out_tree->Branch(b_name, &b_map_flt[b_name], branch_str);
}

// Book a branch for integer values ("I")
void BookBranch(std::map<TString, int> & b_map_int, TTree * out_tree, const TString b_name ) {

  if (b_map_int.find(b_name) == b_map_int.end()) {
    b_map_int[b_name] = -99;
  }
  else {
    std::cout << "Why is there already an integer branch booked with name " << b_name << "?!?  QUITING!!!" << std::endl;
    assert(false);
  }

  TString branch_str;
  branch_str.Form("%s/I", b_name.Data());
  out_tree->Branch(b_name, &b_map_int[b_name], branch_str);
}

// Book a branch for TString values ("C")
void BookBranch(std::map<TString, std::string> & b_map_str, TTree * out_tree, const TString b_name ) {

  if (b_map_str.find(b_name) == b_map_str.end()) {
    b_map_str[b_name] = "-99";
  }
  else {
    std::cout << "Why is there already a string branch booked with name " << b_name << "?!?  QUITING!!!" << std::endl;
    assert(false);
  }

  out_tree->Branch(b_name, &b_map_str[b_name]);
}
