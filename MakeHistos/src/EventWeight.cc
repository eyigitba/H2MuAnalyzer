#include "H2MuAnalyzer/MakeHistos/interface/EventWeight.h"

// Configure constants related to event weight
void ConfigureEventWeight( EventWeightConfig & cfg, const std::string _year ) {

  if (_year == "2016") {
    cfg.year = _year;

    // Weights to apply
    cfg.PU         = true; // PU_wgt
    cfg.muon_ID    = true; // 0.5 * ( MuID_SF_3 +  MuID_SF_4)
    cfg.muon_Iso   = true; // 0.5 * (MuIso_SF_3 + MuIso_SF_4)
    cfg.trig_IsoMu = true; // 0.5 * (IsoMu_SF_3 + IsoMu_SF_4)
    cfg.GEN        = true; // GEN_wgt
  } // End if (_year == "2016")

  else if (_year == "2017") {
    cfg.year = _year;

    // Weights to apply
    cfg.PU         = true; // PU_wgt
    cfg.muon_ID    = true; // MuID_SF_3
    cfg.muon_Iso   = true; // MuIso_SF_3
    cfg.trig_IsoMu = true; // IsoMu_SF_3
    cfg.GEN        = true; // GEN_wgt
  } // End if (_year == "2017")

  else {
    std::cout << "Inside ConfigureEventWeight, invalid year = " << _year << std::endl;
    assert(false);
  }

} // End function: ConfigureEventWeight()


// TODO: implement systematic up/down weights for PU, muon_ID, muon_Iso, and trig_IsoMu - AWB 27.09.2018
float EventWeight( const NTupleBranches & br, const EventWeightConfig & cfg, const bool verbose ) {

  float weight = 1.0;

  if (cfg.year == "2016") {

    if (cfg.PU        ) weight *= br.PU_wgt;
    if (cfg.muon_ID   ) weight *= (0.5 * ( br.MuID_SF_3 +  br.MuID_SF_4));
    if (cfg.muon_Iso  ) weight *= (0.5 * (br.MuIso_SF_3 + br.MuIso_SF_4));
    if (cfg.trig_IsoMu) weight *= (0.5 * (br.IsoMu_SF_3 + br.IsoMu_SF_4));
    if (cfg.GEN       ) weight *= br.GEN_wgt;

  } // End conditional: if (cfg.year == "2016")

  else if (cfg.year == "2017") {

    if (cfg.PU        ) weight *= br.PU_wgt;
    if (cfg.muon_ID   ) weight *= br.MuID_SF_3;
    if (cfg.muon_Iso  ) weight *= br.MuIso_SF_3;
    if (cfg.trig_IsoMu) weight *= br.IsoMu_SF_3;
    if (cfg.GEN       ) weight *= br.GEN_wgt;

  } // End conditional: if (cfg.year == "2017")

  else {
    std::cout << "Inside EventWeight.cc, invalid year = " << cfg.year << std::endl;
    assert(false);
  }

  return weight;

} // End function: EventWeight()
