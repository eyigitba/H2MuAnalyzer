#include "H2MuAnalyzer/MakeHistos/interface/Sample.h"

///////////////////////////////////
///  Constructors / Destructor  ///
///////////////////////////////////

// Null initialization
Sample::Sample() {
  // lumiWeights = 0;
  xsec = -999; 
  lumi = -999;
}

// Initialization with input file names, sample name, and sample type
Sample::Sample(TString _name, TString _sampleType, std::vector<TString> _filenames) {
  name = _name;
  if (_sampleType != "NONE")
    sampleType = _sampleType;
  filenames = _filenames;
  
  // treename = TString("dimuons/tree");
  // chain = new TChain(treename);
  // for (int i = 0; i < filenames.size(); i++) {
  //   chain->Add(filenames.at(i));
  //   //std::cout << i+1 << ", ";
  // }
  // N = chain->GetEntries();
  
  // lumiWeights = 0;
  // xsec = -999; 
  // lumi = -999;
  
  // setBranchAddresses();
  calculateNoriginal();
}

Sample::~Sample() {
  // Free pointed to memory!
  if (chain != 0) {
    delete chain;
  }
  if (files.size() !=0) {
    files.clear();
  }
  // if (lumiWeights !=0) {
  //   delete lumiWeights;
  // }
}


//////////////////////////
///  Useful functions  ///
//////////////////////////

// Calculate number of events processed in sample
void Sample::calculateNoriginal(const bool verbose) {

  if (verbose) std::cout << "\nInside calculateNoriginal for sample " << name << " of type " << sampleType << std::endl;
  gROOT->SetBatch(1); // Don't draw TCanvas
  // Calculate the number of original events using the meta data tree
  // Calculate the weighted number of original events as well
  TChain* metadata = new TChain("dimuons/metadata");
  for (auto f_name : filenames) {
    metadata->Add(f_name);
  }
  metadata->Draw("sumEventWeights>>eweights_"+name);
  TH1F* sumEventWeightsHist = (TH1F*) gDirectory->Get("eweights_"+name); 
  
  // There are many ttrees combined so the histogram has a numEvents entry for each
  // file. The total number is equal to the average number times the total number of entries.
  nOriginalWeighted = sumEventWeightsHist->GetEntries()*sumEventWeightsHist->GetMean();
  
  metadata->Draw("originalNumEvents>>nevents_"+name);
  TH1F* nEventsHist = (TH1F*) gDirectory->Get("nevents_"+name); 
  nOriginal = nEventsHist->GetEntries()*nEventsHist->GetMean();
  if (metadata !=0) delete metadata;
  if (verbose) std::cout << "Computed nOriginal = " << nOriginal << std::endl;
}

// Calculate the appropriate scale factor based on the luminosity
float Sample::getLumiScaleFactor(const float luminosity, const bool verbose) {
  // Scale the MC histograms based upon the data luminosity, the number of events
  // that the CMSSW analyzer looked at, and the xsec for the process
  if (verbose) std::cout << "\nInside getLumiScaleFactor, about to return " << (sampleType.EqualTo("data") ? 1.0 : (luminosity*xsec/nOriginalWeighted)) << std::endl;
  if(sampleType.EqualTo("data")) return 1.0;
  else return luminosity*xsec/nOriginalWeighted;
}
