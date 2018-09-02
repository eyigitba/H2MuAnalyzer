
#ifndef HISTO_HELPER
#define HISTO_HELPER

#include <iostream>
#include <map>

#include "TH1.h"
#include "TH2.h"


//struct Histo1 {
//  TH1 * histo = 0;
//  TString sample_name 	= "";
//  TString begin_file_name = "";
//  TString end_file_name = "";
//  TString variable_name = "";
//  TString sel_name 	= "";
//  TString opt_sel_name	= "";
//  TString cat_sel_name 	= "";
//}
//
//struct Histo2 {
//
//}
// wrapping histo may not make much sense, since information is passed on through the output root file, in which elements in this struct are lost // XWZ 31.08.2018


// Map of <name, hist> for already-booked histograms
//std::map<TString, TH1*> h_map_1D;
//std::map<TString, TH2*> h_map_2D;
// histo maps can be defined here or in macros. For clarity, they are defined in each macros // XWZ 31.08.2018

// Book a 1D histogram (defaults to TH1D)
TH1 * BookHisto(TString h_name, int nBins, float min, float max, TString opt1 = "TH1D");
// Book a 2D histogram (defaults to TH2D)
TH2 * BookHisto(TString h_name, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, TString opt1 = "TH2D");

// Fill a 1D histogram 
void FillHisto(TH1 * hist, double val, float weight = 1.0) ;
// Fill a 2D histogram
void FillHisto(TH2 * hist, double valX, double valY, float weight = 1.0);

// Book a 1D histogram (if it does not already exist), then fill it, default weight 1.0
void BookAndFill(std::map<TString, TH1*> & h1_map, TString h_name, int nBins, float min, float max, double val, float weight = 1.0);
// Book a 1D histogram (if it does not already exist), then fill it, default weight 1.0
void BookAndFill(std::map<TString, TH2*> & h2_map, TString h_name, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, double valX, double valY, float weight = 1.0);

// Return maps from HistoHelper
//void RetreiveMap(std::map<TString, TH1*> & new_map);
//void RetreiveMap(std::map<TString, TH2*> & new_map);
// retrieve can only work if h_map is defined here in the header // XWZ 31.08.2018

#endif  // #ifndef HISTO_HELPER
