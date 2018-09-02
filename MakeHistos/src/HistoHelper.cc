
#include "H2MuAnalyzer/MakeHistos/interface/HistoHelper.h"

// Book a 1D histogram (defaults to TH1D)
TH1 * BookHisto(TString h_name, int nBins, float min, float max, TString opt1) {
  // std::cout << "  * Inside BookHisto(" << h_name << ", " << nBins << ", " << min << ", " << max << ", " << opt1 << ")" << std::endl;
  if      (opt1.CompareTo("TH1D") == 0) return new TH1D(h_name, h_name, nBins, min, max);
  else if (opt1.CompareTo("TH1F") == 0) return new TH1F(h_name, h_name, nBins, min, max);
  else { std::cout << "BookHisto option " << opt1 << " is neither TH1D nor TH1F - exiting" << std::endl; return 0; }
  
}
// Book a 2D histogram (defaults to TH2D)
TH2 * BookHisto(TString h_name, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, TString opt1) {
  if      (opt1.CompareTo("TH2D") == 0) return new TH2D(h_name, h_name, nBinsX, minX, maxX, nBinsY, minY, maxY);
  else if (opt1.CompareTo("TH2F") == 0) return new TH2F(h_name, h_name, nBinsX, minX, maxX, nBinsY, minY, maxY);
  else { std::cout << "BookHisto option " << opt1 << " is neither TH2D nor TH2F - exiting" << std::endl; return 0; }

}

// Fill a 1D histogram
void FillHisto(TH1 * hist, double val, float weight) {
  // std::cout << "    * Filling " << hist->GetName() << " with value " << val << std::endl;
  hist->Fill(val, weight);
}
// Fill a 2D histogram
void FillHisto(TH2 * hist, double valX, double valY, float weight) { hist->Fill(valX, valY, weight); }

// Book a 1D histogram (if it does not already exist), then fill it
void BookAndFill(std::map<TString, TH1*> & h1_map, TString h_name, int nBins, float min, float max, double val, float weight) {
  // std::cout << "\nInside BookAndFill(" << h_name << ", " << nBins << ", " << min << ", " << max << ", " << val << ")" << std::endl;
  if (h1_map.find(h_name) == h1_map.end()) {
    // std::cout << "  * Did not find histogram in map, now booking ... " << std::endl;
    TH1 * hist = BookHisto(h_name, nBins, min, max);
    h1_map[h_name] = hist;
    // std::cout << "  * ... booked!" << std::endl;
  }
  FillHisto(h1_map.find(h_name)->second, val, weight);
  // for (const auto it : h_map_1D) std::cout << "Inside HistoHelper, map at " << it.first << " = " << it.second->Integral() << std::endl;
}
// Book a 2D histogram (if it does not already exist), then fill it
void BookAndFill(std::map<TString, TH2*> & h2_map, TString h_name, int nBinsX, float minX, float maxX, int nBinsY, float minY, float maxY, double valX, double valY, float weight) {
  if (h2_map.find(h_name) == h2_map.end()) {
    TH2 * hist = BookHisto(h_name, nBinsX, minX, maxX, nBinsY, minY, maxY);
    h2_map[h_name] = hist;
  }
  FillHisto(h2_map.find(h_name)->second, valX, valY, weight);
}

// Return maps from HistoHelper
//void RetrieveMap(std::map<TString, TH1*> & new_map) {
//  for (const auto it : h_map_1D) new_map[it.first] = it.second;
//}
//void RetrieveMap(std::map<TString, TH2*> & new_map) {
//  for (const auto it : h_map_2D) new_map[it.first] = it.second;
//}
// retrieve can only work if h_map is defined in the header // XWZ 31.08.2018 
