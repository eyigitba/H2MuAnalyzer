
#ifndef READMVA
#define READMVA

#include "TXMLEngine.h"
#include "TMVA/Reader.h"

// Class to store all information needed by TMVA to evaluate MVAs via XMLs
class MVA {
  
 public:
  MVA();
  MVA(TString _dir, TString _xml, TString _method);
  ~MVA();

  // Input file
  TString file_name;
  // MVA method
  TString method;
  // Reader for input XML files
  TMVA::Reader * reader;
  // Training and spectator variables
  std::vector<TString> train_vars;
  std::vector<TString> spect_vars;
  // Values for all variables
  std::map<TString, Float_t> train_vals;
  std::map<TString, Float_t> spect_vals;

  // Evaluate MVA output value for a given event
  float Evaluate( const std::map<TString, float> val_map_flt,
		  const std::map<TString, int>   val_map_int );

}; // End class MVA

// Supporting function to get variable names from XML file
void GetNamesRecursive( TXMLEngine* xml_eng, XMLNodePointer_t node, TString node_name,
			TString node_att, std::vector<TString> & names );

#endif  // #ifndef READMVA
