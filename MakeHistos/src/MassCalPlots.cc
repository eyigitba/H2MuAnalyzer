#include "H2MuAnalyzer/MakeHistos/interface/MassCalPlots.h"

/////////////////////////////////////////////////
////    functions for MassCalConfig struct   ////
/////////////////////////////////////////////////
void ConfigureMassCal(MassCalConfig & cfg, const TString peak_name, const TString xname, const TString yname) {

  // setting configurations for the mass peak
  cfg.peak = peak_name;
  if (peak_name == "Z") {
    cfg.minM    = 86.2;
    cfg.maxM    = 96.2;
    cfg.nbinsM  = 50;
  }
  else if (peak_name == "H") {
    cfg.minM    = 120.0;
    cfg.maxM    = 130.0;
    cfg.nbinsM  = 50;
  }
  else {
    std::cout << "invalid peak option: peak_name = " << peak_name << std::endl;
    assert(false);
  }
  
  // setting configurations for the x variable 
  cfg.nameX = xname;
  if (xname == "dimu_pt") {
    cfg.minX     = 0;
    cfg.maxX     = 500;
    cfg.binningX = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,110,120,130,140,150,200,500}; 
    cfg.nbinsX   = cfg.binningX.size()-1;
  }
  else if (xname == "dimu_eta") {
    cfg.minX	= -5;
    cfg.maxX	= 5;
    cfg.nbinsX	= 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "dimu_phi") {
    cfg.minX    = -3.2;
    cfg.maxX    = 3.2;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "dimu_dEta") {
    cfg.minX    = -5;
    cfg.maxX    = 5;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "dimu_dPhi") {
    cfg.minX    = -3.2;
    cfg.maxX    = 3.2;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }

  else if (xname == "d0_mean") {
    cfg.minX    = -0.01;
    cfg.maxX    = 0.01;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "d0_diff") {
    cfg.minX    = -0.01;
    cfg.maxX    = 0.01;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "muP_d0") {
    cfg.minX    = -0.01;
    cfg.maxX    = 0.01;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "muN_d0") {
    cfg.minX    = -0.01;
    cfg.maxX    = 0.01;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "muP_eta") {
    cfg.minX    = -2.5;
    cfg.maxX    = 2.5;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "muN_eta") {
    cfg.minX    = -2.5;
    cfg.maxX    = 2.5;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }

  else if (xname == "mu1_pt") {
    cfg.minX     = 0;
    cfg.maxX     = 500;
    cfg.binningX = {10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,110,120,130,150,200,500};
    cfg.nbinsX   = cfg.binningX.size()-1;
  } 
  else if (xname == "mu1_eta") {
    cfg.minX    = -2.5;
    cfg.maxX    = 2.5;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "mu1_phi") {
    cfg.minX    = -3.2;
    cfg.maxX    = 3.2;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "mu2_pt") {
    cfg.minX     = 0;
    cfg.maxX     = 500;
    cfg.binningX = {10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,110,120,130,150,200,500};
    cfg.nbinsX   = cfg.binningX.size()-1;
  }
  else if (xname == "mu2_eta") {
    cfg.minX    = -2.5;
    cfg.maxX    = 2.5;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }
  else if (xname == "mu2_phi") {
    cfg.minX    = -3.2;
    cfg.maxX    = 3.2;
    cfg.nbinsX  = 20;
    cfg.binningX.clear();
    float width = 1.0 * (cfg.maxX - cfg.minX) / cfg.nbinsX;
    for (int i=0; i<cfg.nbinsX+1; i++) cfg.binningX.push_back( cfg.minX + width * i );
  }


  else {
    std::cout << "invalid xname option: xname = " << xname << std::endl;
    assert(false);
  }

  // setting attributes for the y variable
  if (yname == "NA") {
  }
  else if (yname == "") {
  }
  else {
    std::cout << "invalid yname option: yname = " << yname << std::endl;
    assert(false);
  }  
} // end of void ConfigureMassCal(MassCalConfig & cfg, const TString peak_name, const TString xname, const TString yname)



//////////////////////////////////////////////////////
////       functions for MassCalPlots class      /////
//////////////////////////////////////////////////////

// book the plots needed for a mass cal study specified in cfg
MassCalPlots::MassCalPlots(const TString h_pre, MassCalConfig _cfg) {
  cfg = _cfg;

  if (cfg.nbinsX == -99) {
    std::cout << "weird case: nBinsX = -99" << std::endl;
    assert(false);
  }
  float edgesX[cfg.nbinsX+1];
  for (int i = 0; i < cfg.nbinsX+1; i++) edgesX[i] = cfg.binningX.at(i);

  //---------------Book histos for 1D MassCal study--------------
  if (cfg.nameY == "NA") {
    summary_plot_1D = new TH1D("summary_" + h_pre, "summary_" + h_pre, cfg.nbinsX, edgesX);
    for (int i = 0; i < cfg.nbinsX; i++) {  // only need to loop nbinsX times
      TString bin_name;
      bin_name.Form( "_%s_%8.4f_to_%8.4f", cfg.nameX.Data(), cfg.binningX.at(i),cfg.binningX.at(i+1) );
      bin_name.ReplaceAll(" ","");   //  "" for space
      bin_name.ReplaceAll("-","m");  // "m" for "minus"
      bin_name.ReplaceAll(".","p");  // "p" for "point"
      TString MP_name = h_pre + bin_name;
      BookAndFill( mass_plots, MP_name, cfg.nbinsM, cfg.minM, cfg.maxM, -99.0, 0, false);
      // need to have empty plots in order in the map, 
      // sending val = -99.0, weight = 0, overflow = false, so it does not mess with the actual events
    }
  } // end of if (cfg.nameY == "NA") 

  //---------------Book histos for 2D MassCal study--------------
  else {   // nameY is also assigned
    float edgesY[cfg.nbinsY+1];
    for (int i = 0; i < cfg.nbinsY+1; i++) edgesY[i] = cfg.binningY.at(i);
    summary_plot_2D = new TH2D("summary_" + h_pre, "summary_" + h_pre, cfg.nbinsX, edgesX, cfg.nbinsY, edgesY);
    for (int i = 0; i < cfg.nbinsX; i++) {    // only need to loop nbinsX times
      for (int j = 0; j < cfg.nbinsY; j++) {  // only need to loop nbinsY times
        TString bin_name;
        bin_name.Form( "_%s_%8.4f_to_%8.4f_%s_%8.4f_to_%8.4f", cfg.nameX.Data(), cfg.binningX.at(i), cfg.binningX.at(i+1), cfg.nameY.Data(), cfg.binningY.at(j), cfg.binningY.at(j+1) );
	bin_name.ReplaceAll(" ","");   //  "" for space
        bin_name.ReplaceAll("-","m");  // "m" for "minus"
        bin_name.ReplaceAll(".","p");  // "p" for "point"
      	TString MP_name = h_pre + bin_name;
	BookAndFill( mass_plots, MP_name, cfg.nbinsM, cfg.minM, cfg.maxM, -99.0, 0, false);
	// need to have empty plots in order in the map, 
	// sending val = -99.0, weight = 0, overflow = false, so it does not mess with the actual events
      }
    } // end of for (int i = 0; i < cfg.nbinsX+1; i++)
  } // end of else

} // end of MassCalPlots::MassCalPlots(const TString h_pre, MassCalConfig _cfg) 



// fill mass plots for 1D study
void MassCalPlots::FillEvent(const TString h_pre, float m_value, float x_value, float weight, bool overflow) {
  TString MP_name = ""; 
  // find the mass plot to fill
  for (int i = 0; i < cfg.nbinsX; i++) {  // only need to loop nbinsX times
    if ( x_value >= cfg.binningX.at(i) and x_value < cfg.binningX.at(i+1) ) {
	TString bin_name;
	bin_name.Form( "_%s_%8.4f_to_%8.4f", cfg.nameX.Data(), cfg.binningX.at(i), cfg.binningX.at(i+1) );
	bin_name.ReplaceAll(" ","");   //  "" for space
        bin_name.ReplaceAll("-","m");  // "m" for "minus"
        bin_name.ReplaceAll(".","p");  // "p" for "point"
        MP_name = h_pre + bin_name;
	break;
    }
  } // end of for (int i = 0; i < cfg.nbinsX; i++)

  // fill histograms
  if (MP_name.Length() != 0) {
    BookAndFill( mass_plots, MP_name, cfg.nbinsM, cfg.minM, cfg.maxM, m_value, weight, overflow);
    summary_plot_1D->Fill(x_value, weight);
  }
} // end of MassCalPlots::FillEvent(const TString h_pre, m_value, x_value, weight, overflow)



// fill mass plots for 2D study
void MassCalPlots::FillEvent(const TString h_pre, float m_value, float x_value, float y_value, float weight, bool overflow) {
  TString MP_name = "";
  // find the mass plot to fill
  for (int i = 0; i < cfg.nbinsX; i++) {   // only need to loop nbinsX times
    if ( x_value >= cfg.binningX.at(i) and x_value < cfg.binningX.at(i+1) ) {
      for (int j = 0; j < cfg.nbinsY; j++) { // only need to loop nbinsY times
	if ( y_value >= cfg.binningY.at(i) and y_value < cfg.binningY.at(j+1) ) {
	  TString bin_name;
          bin_name.Form( "_%s_%8.4f_to_%8.4f_%s_%8.4f_to_%8.4f", cfg.nameX.Data(), cfg.binningX.at(i), cfg.binningX.at(i+1), cfg.nameY.Data(), cfg.binningY.at(j), cfg.binningY.at(j+1) );
	  bin_name.ReplaceAll(" ","");   //  "" for space
          bin_name.ReplaceAll("-","m");  // "m" for "minus"
          bin_name.ReplaceAll(".","p");  // "p" for "point"
          TString MP_name = h_pre + bin_name; 
          break;
	}
      } // end of for (int j = 0; j < cfg.nbinsY; j++) 
    }
    if (MP_name.Length() != 0) break;
  } // end of for (int i = 0; i < cfg.nbinsX; i++)

  // fill histograms
  if (MP_name.Length() != 0) {
    BookAndFill( mass_plots, MP_name, cfg.nbinsM, cfg.minM, cfg.maxM, m_value, weight, overflow);
    summary_plot_2D->Fill(x_value, y_value, weight);
  }
} // end of MassCalPlots::FillEvent(const TString h_pre, m_value, x_value, y_value, weight, overflow)


