#include "H2MuAnalyzer/MakeHistos/interface/EventWeight.h"

// Configure constants related to event weight
void ConfigureEventWeight( EventWeightConfig & cfg, const std::string _year, const std::string _SYS ) {

  if (_year == "Legacy2016") {
    cfg.year = _year;
    cfg.SYS  = _SYS;

    // Weights to apply
    cfg.PU         = true; // PU_wgt
    cfg.muon_ID    = true; // 0.5 * ( MuID_SF_3 +  MuID_SF_4)
    cfg.muon_Iso   = true; // 0.5 * (MuIso_SF_3 + MuIso_SF_4)
    cfg.trig_IsoMu = true; // 0.5 * (IsoMu_SF_3 + IsoMu_SF_4)
    cfg.GEN        = true; // GEN_wgt
  } // End if (_year == "Legacy2016")

  else if (_year == "2016" || _year == "2017" || _year == "2018") {
    cfg.year = _year;
    cfg.SYS  = _SYS;

    // Weights to apply
    cfg.PU         = true; // PU_wgt
    cfg.muon_ID    = true; // MuID_SF_3
    cfg.muon_Iso   = true; // MuIso_SF_3
    cfg.trig_IsoMu = true; // IsoMu_SF_3
    cfg.GEN        = true; // GEN_wgt
  } // End if (_year == "2016" || _year == "2017" || _year == "2018")

  else {
    std::cout << "Inside ConfigureEventWeight, invalid year = " << _year << std::endl;
    assert(false);
  }

} // End function: ConfigureEventWeight()


// TODO: implement systematic up/down weights for PU, muon_ID, muon_Iso, and trig_IsoMu - AWB 27.09.2018
float MuonWeight( const NTupleBranches & br, const EventWeightConfig & cfg, const bool verbose ) {

  float mu_weight = 1.0;
  float MuID_SF   = 1.0;
  float MuIso_SF  = 1.0;
  float IsoMu_SF  = 1.0;

  if (cfg.year == "Legacy2016") {

    MuID_SF  = (0.5 * ( br.MuID_SF_3 +  br.MuID_SF_4));
    MuIso_SF = (0.5 * (br.MuIso_SF_3 + br.MuIso_SF_4));
    IsoMu_SF = (0.5 * (br.IsoMu_SF_3 + br.IsoMu_SF_4));

    if (cfg.SYS == "MuID_SF_up")      MuID_SF  = (0.5 * ( br.MuID_SF_3_up   +  br.MuID_SF_4_up  ));
    if (cfg.SYS == "MuID_SF_down")    MuID_SF  = (0.5 * ( br.MuID_SF_3_down +  br.MuID_SF_4_down));

    if (cfg.SYS == "MuIso_SF_up")     MuIso_SF = (0.5 * (br.MuIso_SF_3_up   + br.MuIso_SF_4_up  ));
    if (cfg.SYS == "MuIso_SF_down")   MuIso_SF = (0.5 * (br.MuIso_SF_3_down + br.MuIso_SF_4_down));    

    if (cfg.SYS == "IsoMu_SF_up")     IsoMu_SF = (0.5 * (br.IsoMu_SF_3_up   + br.IsoMu_SF_4_up  ));
    if (cfg.SYS == "IsoMu_SF_down")   IsoMu_SF = (0.5 * (br.IsoMu_SF_3_down + br.IsoMu_SF_4_down)); 
  } // End conditional: if (cfg.year == "Legacy2016")

  else if (cfg.year == "2016" || cfg.year == "2017" || cfg.year == "2018") {

    MuID_SF = br.MuID_SF_3;
    MuIso_SF = br.MuIso_SF_3;
    IsoMu_SF = br.IsoMu_SF_3;

    if (cfg.SYS == "MuID_SF_up")    MuID_SF = br.MuID_SF_3_up;
    if (cfg.SYS == "MuID_SF_down")  MuID_SF = br.MuID_SF_3_down;

    if (cfg.SYS == "MuIso_SF_up")   MuIso_SF = br.MuIso_SF_3_up;
    if (cfg.SYS == "MuIso_SF_down") MuIso_SF = br.MuIso_SF_3_down;

    if (cfg.SYS == "IsoMu_SF_up")   IsoMu_SF = br.IsoMu_SF_3_up;
    if (cfg.SYS == "IsoMu_SF_down") IsoMu_SF = br.IsoMu_SF_3_down;
  } // End conditional: if (cfg.year == "2016" || cfg.year == "2017" || cfg.year == "2018")


  if (cfg.muon_ID   ) mu_weight *= MuID_SF;
  if (cfg.muon_Iso  ) mu_weight *= MuIso_SF;
  if (cfg.trig_IsoMu) mu_weight *= IsoMu_SF;

  else {
    std::cout << "Inside MuonWeight.cc, invalid year = " << cfg.year << std::endl;
    assert(false);
  }

  return mu_weight;

} // End function: MuonWeight()


// TODO: implement systematic up/down weights for PU, muon_ID, muon_Iso, and trig_IsoMu - AWB 27.09.2018
float EventWeight( const NTupleBranches & br, const EventWeightConfig & cfg, const bool verbose ) {

  float evt_weight = MuonWeight(br, cfg, verbose);
  float PU_wgt = 1.0;

  PU_wgt = br.PU_wgt;
  if (cfg.SYS == "PU_wgt_up")   PU_wgt = br.PU_wgt_up;
  if (cfg.SYS == "PU_wgt_down") PU_wgt = br.PU_wgt_down;

  // if (cfg.PU && br.PU_wgt <= 0) {
  //   std::cout << "\n\nTruly bizzare case where PU_wgt = " << br.PU_wgt << "!!!" << std::endl;
  //   std::cout << "Check computation!\n\n" << std::endl;
  // }

  if (cfg.year == "2016" || cfg.year == "2017" || cfg.year == "2018") {

    if (cfg.PU) {
      if (PU_wgt >= 99) {
	std::cout << "\n\nTruly bizzare case where PU_wgt = " << PU_wgt << "!!!" << std::endl;
	std::cout << "Check computation! Setting to 1.\n\n" << std::endl;
      } else evt_weight *= PU_wgt;
    }
    if (cfg.GEN) evt_weight *= br.GEN_wgt;

  } // End conditional: if (cfg.year == "2016" || cfg.year == "2017" || cfg.year == "2018")

  else {
    std::cout << "Inside EventWeight.cc, invalid year = " << cfg.year << std::endl;
    assert(false);
  }

  return evt_weight;

} // End function: EventWeight()
