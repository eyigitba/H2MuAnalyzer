#ifndef Mass_Cal_Plots
#define Mass_Cal_Plots

#include <iostream>
#include <map>

#include "H2MuAnalyzer/MakeHistos/interface/HistoHelper.h" 

////////////////////////////////////////////
////     struct def of MassCalConfig    ////
////////////////////////////////////////////
struct MassCalConfig { 
  // default values
  TString 	peak 	= "None"; 	// options are: "Z", "H"
  TString 	nameX 	= "NA";    	// options are: "dimu_pt", "dimu_eta", "dimu_phi", "mu1_pt", "mu1_eta", "mu1_d0", "mu1_phi", "mu2...", "dimu_dEta", "dimu_dPhi", "diff_d0" ...
  TString 	nameY 	= "NA"; 	// "NA" means it is a 1D plot
  
  float 	minM 	= -99.0;
  float 	maxM 	= -99.0;
  int   	nbinsM 	= -99;            // mass has fixed binning

  float 	minX	= -99.0;
  float 	maxX	= -99.0;
  int   	nbinsX  = -99;
  float 	minY	= -99.0;
  float 	maxY	= -99.0;
  int   	nbinsY	= -99;
  std::vector<Float_t> binningX = {-99.0};
  std::vector<Float_t> binningY = {-99.0};
};

void ConfigureMassCal(MassCalConfig & cfg, const TString peak_name, const TString xname, const TString yname = "NA"); 



///////////////////////////////////////////
////    class def of MassCalPlots      ////
///////////////////////////////////////////
class MassCalPlots {
  public:
    MassCalConfig cfg;
    std::map<TString, TH1*>  mass_plots;
    TH1D* summary_plot_1D;
    TH2D* summary_plot_2D;

    MassCalPlots(const TString h_pre, MassCalConfig _cfg);

    void FillEvent(const TString h_pre, float m_value, float x_value, float weight, bool overflow = false);
    void FillEvent(const TString h_pre, float m_value, float x_value, float y_value, float weight, bool overflow = false);
}; 



#endif  // #ifndef Mass_Cal_Plots
