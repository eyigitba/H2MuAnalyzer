
// Functions to read in XMLs for MVAs from TMVA
//   * Modeled on https://github.com/acarnes/UFDimuAnalysis/blob/master/bin/TMVAClassificationApplication_H2Mu.cxx
//            and https://github.com/acarnes/UFDimuAnalysis/blob/master/tools/TMVATools.cxx
#include "H2MuAnalyzer/MakeHistos/interface/ReadMVA.h"

bool VERBOSE = false;

// MVA class constructor
MVA::MVA() {
}

// Initialization of MVA class with input XML file names and MVA method 
MVA::MVA(TString _dir, TString _xml, TString _method) {

  file_name = _dir + _xml;
  method    = _method;

  std::cout << "\nLoading MVA XML from " << file_name << " with method " << method << std::endl;

  // ******************************************************************************** //
  std::cout << "  * Loading the input variable names from the XML file" << std::endl;
  // ******************************************************************************** //

  // First create the engine
  TXMLEngine* xml_eng = new TXMLEngine;
  // Now try to parse xml file
  XMLDocPointer_t xml_doc = xml_eng->ParseFile(file_name);
  if (xml_doc == 0) {
    delete xml_eng;
    std::cout << "\n\n*** FAILED TO PARSE " << _dir << _xml << " !!!  Exiting." << std::endl;
    return;
  }
  // Get access to main node of the xml file
  XMLNodePointer_t main_node = xml_eng->DocGetRootElement(xml_doc);
  // Recursively connect nodes together
  GetNamesRecursive(xml_eng, main_node, "Variable",  "Label", train_vars);
  GetNamesRecursive(xml_eng, main_node, "Spectator", "Label", spect_vars);
  // Release memory before exit
  xml_eng->FreeDoc(xml_doc);
  delete xml_eng;

  // ************************************************************************ //
  std::cout << "  * Booking the variables into the TMVA reader" << std::endl;
  // ************************************************************************ //

  // Book a new TMVA reader
  reader = new TMVA::Reader("!Color:!Silent");    
  // Load training and spectator variables into the reader
  // Must have the same names as in the XML, and be booked in the same order
  for (auto & var: train_vars) {
    if (VERBOSE) std::cout << "    - Booking input variable " << var << std::endl;
    train_vals[var] = -999;
    reader->AddVariable(var, &train_vals[var]);
  }
  for (auto & var: spect_vars) {
    // std::cout << "    - Booking spectator variable " << var << std::endl;
    spect_vals[var] = -999;
    reader->AddSpectator(var, &spect_vals[var]);
  }
  // Book BDT method into reader, now that the variables are set up
  std::cout << "  * Completing the MVA booking ... ";
  reader->BookMVA( method, file_name );  
  std::cout << "complete!" << std::endl;

} // End MVA class constructor MVA::MVA(TString _dir, TString _xml, TString _method)


// MVA class destructor
MVA::~MVA() {

  // Free pointed to memory!
  if (reader != 0) delete reader;
  
  train_vars.clear();
  spect_vars.clear();
  train_vals.clear();
  spect_vals.clear();

} // End MVA class destructor MVA::~MVA()


// Evaluate MVA output score
float MVA::Evaluate( const std::map<TString, float> val_map_flt,
		     const std::map<TString, int>   val_map_int ) {

  // Loop over variable [name, value] map
  for (auto & var : train_vals) {
    if ( val_map_flt.find(var.first) != val_map_flt.end() ) {
      train_vals[var.first] = val_map_flt.find(var.first)->second;
    } else if ( val_map_int.find(var.first) != val_map_int.end() ) {
      train_vals[var.first] = val_map_int.find(var.first)->second;
    } else {
      std::cout << "\n\n*** COULDN'T FIND VARIABLE " << var.first << " IN val_map_flt OR val_map_int!!!  Returning -999" << std::endl;
      return -999;
    }
  }

  return reader->EvaluateMVA(method);
} // End function: MVA::Evaluate()


// Supporting function to get variable names from XML file
void GetNamesRecursive( TXMLEngine* xml_eng, XMLNodePointer_t node, TString node_name,
			TString node_att, std::vector<TString> & names ) {

  TString nname = xml_eng->GetNodeName(node);
   
  // Display attributes
  XMLAttrPointer_t attr = xml_eng->GetFirstAttr(node);
  while (attr != 0) {
    TString att_string = xml_eng->GetAttrName (attr);
    TString val_string = xml_eng->GetAttrValue(attr);
    if (nname == node_name && att_string == node_att) {
      names.push_back(val_string);
    }
    attr = xml_eng->GetNextAttr(attr);
  }
   
  // Display all child nodes
  XMLNodePointer_t child = xml_eng->GetChild(node);
  while (child != 0) {
    GetNamesRecursive(xml_eng, child, node_name, node_att, names);
    child = xml_eng->GetNext(child);
  }

} // End function: GetNamesRecursive()
