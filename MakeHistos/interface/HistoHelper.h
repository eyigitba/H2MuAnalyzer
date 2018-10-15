
#ifndef HISTO_HELPER
#define HISTO_HELPER

#include <iostream>
#include <map>

#include "TH1.h"
#include "TH2.h"


// Book a 1D histogram (defaults to TH1D)
TH1 * BookHisto(TString h_name, int nBins, float min, float max, TString opt1 = "TH1D");
// Book a 2D histogram (defaults to TH2D)
TH2 * BookHisto(TString h_name, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, TString opt1 = "TH2D");

// Fill a 1D histogram 
void FillHisto(TH1 * hist, double val, float weight = 1.0, bool overflow = true) ;
// Fill a 2D histogram
void FillHisto(TH2 * hist, double valX, double valY, float weight = 1.0, bool overflow = true);

// Book a 1D histogram (if it does not already exist), then fill it, default weight 1.0
void BookAndFill( std::map<TString, TH1*> & h1_map, TString h_name,
		  int nBins, float min, float max,
		  double val, float weight = 1.0, bool overflow = true );
// Book a 1D histogram (if it does not already exist), then fill it, default weight 1.0
void BookAndFill( std::map<TString, TH2*> & h2_map, TString h_name,
                  int nBinsX, float minX, float maxX,
                  int nBinsY, float minY, float maxY,
                  double valX, double valY, float weight = 1.0, bool overflow = true );

//Initialize a default set of histos for data/MC comparison
void BookForMCvsData(std::map<TString, TH1*> & h1_map, std::string sample, std::vector<std::string> OPT_CUTS, std::vector<std::string> CAT_CUTS);




#endif  // #ifndef HISTO_HELPER
