
#ifndef HISTO_HELPER
#define HISTO_HELPER

#include <iostream>
#include <map>

#include "TH1.h"
#include "TH2.h"

#include "H2MuAnalyzer/MakeHistos/interface/MiniNTupleHelper.h"  // BookBranch function


// Book a 1D histogram (defaults to TH1D)
TH1 * BookHisto( const TString h_name, const int nBins, const float min, const float max, const TString opt1 = "TH1D");
// Book a 1D histogram with variable binning (defaults to TH1D)
TH1 * BookHisto( const TString h_name, const int nBins, const std::vector<float> binning, const TString opt1 = "TH1D" ); 
// Book a 2D histogram (defaults to TH2D)
TH2 * BookHisto( const TString h_name, const int nBinsX, const float minX, const float maxX,
		 const int nBinsY, const float minY, const float maxY, const TString opt1 = "TH2D" );

// Fill a 1D histogram 
void FillHisto( TH1 * hist, double val, const float weight = 1.0, const bool overflow = true );
// Fill a 2D histogram
void FillHisto( TH2 * hist, double valX, double valY, const float weight = 1.0, const bool overflow = true );

// Book a 1D histogram (if it does not already exist), then fill it, default weight 1.0
void BookAndFill( std::map<TString, TH1*> & h1_map, const TString h_name,
		  const int nBins, const float min, const float max,
		  const double val, const float weight = 1.0, const bool overflow = true );
// Book a 1D histogram with variable binning (if it does not already exist), then fill it, default weight 1.0
void BookAndFill( std::map<TString, TH1*> & h1_map, const TString h_name,
                  const int nBins, const std::vector<float> binning,
                  const double val, const float weight, const bool overflow = true);
// Book a 2D histogram (if it does not already exist), then fill it, default weight 1.0
void BookAndFill( std::map<TString, TH2*> & h2_map, const TString h_name,
                  const int nBinsX, const float minX, const float maxX,
                  const int nBinsY, const float minY, const float maxY,
                  const double valX, const double valY, const float weight = 1.0, const bool overflow = true );

// Book a floating point branch in a tree (if it does not already exist), then fill it
void BookAndFill( std::map<TString, float> & bf_map, TTree * tree,
		  const TString h_pre, const TString b_name, const double val );
// Book an integer branch in a tree (if it does not already exist), then fill it
void BookAndFill( std::map<TString, int> & bi_map, TTree * tree,
		  const TString h_pre, const TString b_name, const double val );
// Book a string branch in a tree (if it does not already exist), then fill it
void BookAndFill( std::map<TString, std::string> & bs_map, TTree * tree,
		  const TString h_pre, const TString b_name, const TString val );


// Book a 1D histogram (if it does not already exist), then fill it
// Book a floating point branch in a tree (if it does not already exist), then fill it
void BookAndFill( const TString opt, std::map<TString, float> & bf_map, TTree * tree,
		  std::map<TString, TH1*> & h1_map, const TString h_pre, const TString h_name,
		  const int nBins, const float min, const float max,
		  const double val, const float weight = 1.0, const bool overflow = true );
void BookAndFill( std::tuple< const TString, std::map<TString, float> &, TTree *, std::map<TString, TH1*> &, const TString > tup,
		  const TString h_name, const int nBins, const float min, const float max,
		  const double val, const float weight = 1.0, const bool overflow = true );
// Book a 1D histogram (if it does not already exist), then fill it
// Book an integer branch in a tree (if it does not already exist), then fill it
void BookAndFill( const TString opt, std::map<TString, int> & bi_map, TTree * tree,
		  std::map<TString, TH1*> & h1_map, const TString h_pre, const TString h_name,
		  const int nBins, const float min, const float max,
		  const double val, const float weight = 1.0, const bool overflow = true );
void BookAndFill( std::tuple< const TString, std::map<TString, int> &, TTree *, std::map<TString, TH1*> &, const TString > tup,
		  const TString h_name, const int nBins, const float min, const float max,
		  const double val, const float weight = 1.0, const bool overflow = true );


// Initialize a default set of histos for data/MC comparison
void BookForMCvsData(std::map<TString, TH1*> & h1_map, std::string sample, std::vector<std::string> OPT_CUTS, std::vector<std::string> CAT_CUTS);

// Test if a string starts or ends with a sub-string
bool StartsWith(const std::string str, const std::string beg);
bool EndsWith(const std::string str, const std::string end);



#endif  // #ifndef HISTO_HELPER
