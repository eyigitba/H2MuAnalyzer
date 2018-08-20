
#ifndef HISTO_HELPER
#define HISTO_HELPER

#include <iostream>
#include <map>

#include "TH1.h"
#include "TH2.h"

// Map of <name, hist> for already-booked histograms
std::map<TString, TH1*> h_map_1D;
std::map<TString, TH2*> h_map_2D;

// Book a 1D histogram (defaults to TH1D)
TH1 * BookHisto(TString h_name, int nBins, float min, float max, TString opt1 = "TH1D");
// Book a 2D histogram (defaults to TH2D)
TH2 * BookHisto(TString h_name, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, TString opt1 = "TH2D");

// Fill a 1D histogram
void FillHisto(TH1 * hist, double val);
// Fill a 2D histogram
void FillHisto(TH2 * hist, double valX, double valY);

// Book a 1D histogram (if it does not already exist), then fill it
void BookAndFill(TString h_name, int nBins, float min, float max, double val);
// Book a 1D histogram (if it does not already exist), then fill it
void BookAndFill(TString h_name, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, double valX, double valY);

// Return maps from HistoHelper
void RetreiveMap(std::map<TString, TH1*> & new_map);
void RetreiveMap(std::map<TString, TH2*> & new_map);

#endif  // #ifndef HISTO_HELPER
