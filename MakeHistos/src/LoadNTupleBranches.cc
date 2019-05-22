
#include "H2MuAnalyzer/MakeHistos/interface/LoadNTupleBranches.h"

// Manually force seg-fault because ROOT doesn't handle exceptions
void ASSERT(bool condition, std::string message) {
  if (!condition) {
    std::cout << "\n\nASSERTION FAILED WITH MESSAGE: " << message << "\n" << std::endl;
    assert(false);
  }
}

void SetBranchAddresses(TChain & ch_, NTupleBranches & br, std::vector<std::string> opts, bool verbose) {

  if (verbose) std::cout << "\nInside SetBranchAddress" << std::endl;

  TChain * ch = (&ch_);

  // Configure special options
  bool is2016    = false;
  bool is2017    = false;
  bool is2018    = false;
  bool isSlim    = false;
  bool notSlim   = false;
  bool loadGEN   = false;
  bool loadJES   = false;
  bool loadFlags = false;
  bool loadEffs   = false;
  bool loadWgts  = false;

  for (uint i = 0; i < opts.size(); i++) {
    if (verbose) std::cout << "  * Using option " << opts.at(i) << std::endl;
    if (opts.at(i) == "2016")    is2016    = true;
    if (opts.at(i) == "2017")    is2017    = true;
    if (opts.at(i) == "2018")    is2018    = true;
    if (opts.at(i) == "Slim")    isSlim    = true;
    if (opts.at(i) == "notSlim") notSlim   = true;
    if (opts.at(i) == "GEN")     loadGEN   = true;
    if (opts.at(i) == "JES")     loadJES   = true;
    if (opts.at(i) == "Flags")   loadFlags = true;
    if (opts.at(i) == "Effs")    loadEffs  = true;
    if (opts.at(i) == "Wgts")    loadWgts  = true;
  }

  assert(is2016 || is2017 || is2018);
  assert(isSlim || notSlim);

  if (verbose) std::cout << "is2016 = " << is2016 << ", is2017 = " << is2017 << ", is2018 = " << is2018
		         << ", isSlim = " << isSlim << ", loadGen = " << loadGEN << ", loadJES = " << loadJES
			 << ", loadFlags = " << loadFlags << ", loadEffs = " << loadEffs 
			 << ", loadWgts = " << loadWgts << std::endl;

  ch->SetBranchAddress("event", &(br.event));
  ch->SetBranchAddress("vertices", &(br.vertices));

  ch->SetBranchAddress("muons", &(br.muons));
  ch->SetBranchAddress("muPairs", &(br.muPairs));
  ch->SetBranchAddress("eles", &(br.eles));
  if ( is2016            ||  isSlim) ch->SetBranchAddress("jets", &(br.slimJets));
  if ((is2017 || is2018) && !isSlim) ch->SetBranchAddress("jets", &(br.jets));
  ch->SetBranchAddress("jetPairs", &(br.jetPairs));
  ch->SetBranchAddress("met", &(br.met));
  ch->SetBranchAddress("mht", &(br.mht));

  ch->SetBranchAddress("nVertices", &(br.nVertices));
  ch->SetBranchAddress("nMuons", &(br.nMuons));
  ch->SetBranchAddress("nMuPairs", &(br.nMuPairs));
  ch->SetBranchAddress("nEles", &(br.nEles));
  ch->SetBranchAddress("nJets", &(br.nJets));
  ch->SetBranchAddress("nJetPairs", &(br.nJetPairs));
  ch->SetBranchAddress("nJetsCent", &(br.nJetsCent));
  ch->SetBranchAddress("nJetsFwd", &(br.nJetsFwd));
  ch->SetBranchAddress("nBLoose", &(br.nBLoose));
  ch->SetBranchAddress("nBMed", &(br.nBMed));
  ch->SetBranchAddress("nBTight", &(br.nBTight));

  if (loadGEN) {
    ch->SetBranchAddress("nPU", &(br.nPU));
    ch->SetBranchAddress("LHE_HT", &(br.LHE_HT));
    ch->SetBranchAddress("genParents", &(br.genParents));
    ch->SetBranchAddress("genMuons", &(br.genMuons));
    ch->SetBranchAddress("genMuPairs", &(br.genMuPairs));
    if ((is2017 || is2018) && !isSlim) ch->SetBranchAddress("genJets", &(br.genJets));

    ch->SetBranchAddress("nGenParents", &(br.nGenParents));
    ch->SetBranchAddress("nGenMuons", &(br.nGenMuons));
    ch->SetBranchAddress("nGenMuPairs", &(br.nGenMuPairs));
    if ((is2017 || is2018) && !isSlim) ch->SetBranchAddress("nGenJets", &(br.nGenJets));
  }

  // ch->SetBranchAddress("hltPaths", &(br.hltPaths));  // Causes error messages after exiting code, for some reason - AWB 15.08.2018
  // ch->SetBranchAddress("btagName", &(br.btagName));  // Causes segfault, for some reason - AWB 15.08.2018

  if (loadJES) {
    if (isSlim)  ch->SetBranchAddress("jets_JES_up", &(br.slimJets_JES_up));
    if (isSlim)  ch->SetBranchAddress("jets_JES_down", &(br.slimJets_JES_down));
    if (!isSlim) ch->SetBranchAddress("jets_JES_up", &(br.jets_JES_up));
    if (!isSlim) ch->SetBranchAddress("jets_JES_down", &(br.jets_JES_down));
    ch->SetBranchAddress("jetPairs_JES_up", &(br.jetPairs_JES_up));
    // ch->SetBranchAddress("jetPairs_JES_down", &(br.jetPairs_JES_down));  // Causes segfault, for some reason - AWB 15.08.2018
    ch->SetBranchAddress("met_JES_up", &(br.met_JES_up));
    ch->SetBranchAddress("met_JES_down", &(br.met_JES_down));
    ch->SetBranchAddress("mht_JES_up", &(br.mht_JES_up));
    ch->SetBranchAddress("mht_JES_down", &(br.mht_JES_down));

    ch->SetBranchAddress("nJets_JES_up", &(br.nJets_JES_up));
    ch->SetBranchAddress("nJetsCent_JES_up", &(br.nJetsCent_JES_up));
    ch->SetBranchAddress("nBLoose_JES_up", &(br.nBLoose_JES_up));
    ch->SetBranchAddress("nBMed_JES_up", &(br.nBMed_JES_up));
    ch->SetBranchAddress("nBTight_JES_up", &(br.nBTight_JES_up));
    ch->SetBranchAddress("nJetsFwd_JES_up", &(br.nJetsFwd_JES_up));
    ch->SetBranchAddress("nJets_JES_down", &(br.nJets_JES_down));
    ch->SetBranchAddress("nJetsCent_JES_down", &(br.nJetsCent_JES_down));
    ch->SetBranchAddress("nJetsFws_JES_down", &(br.nJetsFws_JES_down));
    ch->SetBranchAddress("nBLoose_JES_down", &(br.nBLoose_JES_down));
    ch->SetBranchAddress("nBMed_JES_down", &(br.nBMed_JES_down));
    ch->SetBranchAddress("nBTight_JES_down", &(br.nBTight_JES_down));
  }

  if (loadFlags) {
    ch->SetBranchAddress("Flag_all", &(br.Flag_all));
    ch->SetBranchAddress("Flag_badMu", &(br.Flag_badMu));
    ch->SetBranchAddress("Flag_dupMu", &(br.Flag_dupMu));
    ch->SetBranchAddress("Flag_halo", &(br.Flag_halo));
    ch->SetBranchAddress("Flag_PV", &(br.Flag_PV));
    ch->SetBranchAddress("Flag_HBHE", &(br.Flag_HBHE));
    ch->SetBranchAddress("Flag_HBHE_Iso", &(br.Flag_HBHE_Iso));
    ch->SetBranchAddress("Flag_ECAL_TP", &(br.Flag_ECAL_TP));
    ch->SetBranchAddress("Flag_BadChCand", &(br.Flag_BadChCand));
    ch->SetBranchAddress("Flag_eeBadSc", &(br.Flag_eeBadSc));
    ch->SetBranchAddress("Flag_ecalBadCalib", &(br.Flag_ecalBadCalib));
  }

  if (loadEffs) {
    ch->SetBranchAddress("IsoMu_eff_3", &(br.IsoMu_eff_3));
    ch->SetBranchAddress("IsoMu_eff_3_up", &(br.IsoMu_eff_3_up));
    ch->SetBranchAddress("IsoMu_eff_3_down", &(br.IsoMu_eff_3_down));
    ch->SetBranchAddress("MuID_eff_3", &(br.MuID_eff_3));
    ch->SetBranchAddress("MuID_eff_3_up", &(br.MuID_eff_3_up));
    ch->SetBranchAddress("MuID_eff_3_down", &(br.MuID_eff_3_down));
    ch->SetBranchAddress("MuIso_eff_3", &(br.MuIso_eff_3));
    ch->SetBranchAddress("MuIso_eff_3_up", &(br.MuIso_eff_3_up));
    ch->SetBranchAddress("MuIso_eff_3_down", &(br.MuIso_eff_3_down));

    if (is2016) {
      ch->SetBranchAddress("IsoMu_eff_4", &(br.IsoMu_eff_4));
      ch->SetBranchAddress("IsoMu_eff_4_up", &(br.IsoMu_eff_4_up));
      ch->SetBranchAddress("IsoMu_eff_4_down", &(br.IsoMu_eff_4_down));
      ch->SetBranchAddress("MuID_eff_4", &(br.MuID_eff_4));
      ch->SetBranchAddress("MuID_eff_4_up", &(br.MuID_eff_4_up));
      ch->SetBranchAddress("MuID_eff_4_down", &(br.MuID_eff_4_down));
      ch->SetBranchAddress("MuIso_eff_4", &(br.MuIso_eff_4));
      ch->SetBranchAddress("MuIso_eff_4_up", &(br.MuIso_eff_4_up));
      ch->SetBranchAddress("MuIso_eff_4_down", &(br.MuIso_eff_4_down));
      ch->SetBranchAddress("IsoMu_eff_bug", &(br.IsoMu_eff_bug));
      ch->SetBranchAddress("IsoMu_eff_bug_up", &(br.IsoMu_eff_bug_up));
      ch->SetBranchAddress("IsoMu_eff_bug_down", &(br.IsoMu_eff_bug_down));
    }

  }

  if (loadWgts) {

    ch->SetBranchAddress("PU_wgt", &(br.PU_wgt));
    ch->SetBranchAddress("PU_wgt_up", &(br.PU_wgt_up));
    ch->SetBranchAddress("PU_wgt_down", &(br.PU_wgt_down));
    ch->SetBranchAddress("GEN_wgt", &(br.GEN_wgt));
 
    ch->SetBranchAddress("IsoMu_SF_3", &(br.IsoMu_SF_3));
    ch->SetBranchAddress("IsoMu_SF_3_up", &(br.IsoMu_SF_3_up));
    ch->SetBranchAddress("IsoMu_SF_3_down", &(br.IsoMu_SF_3_down));
    ch->SetBranchAddress("MuID_SF_3", &(br.MuID_SF_3));
    ch->SetBranchAddress("MuID_SF_3_up", &(br.MuID_SF_3_up));
    ch->SetBranchAddress("MuID_SF_3_down", &(br.MuID_SF_3_down));
    ch->SetBranchAddress("MuIso_SF_3", &(br.MuIso_SF_3));
    ch->SetBranchAddress("MuIso_SF_3_up", &(br.MuIso_SF_3_up));
    ch->SetBranchAddress("MuIso_SF_3_down", &(br.MuIso_SF_3_down));

    if (is2016) {
      ch->SetBranchAddress("IsoMu_SF_4", &(br.IsoMu_SF_4));
      ch->SetBranchAddress("IsoMu_SF_4_up", &(br.IsoMu_SF_4_up));
      ch->SetBranchAddress("IsoMu_SF_4_down", &(br.IsoMu_SF_4_down));
      ch->SetBranchAddress("MuID_SF_4", &(br.MuID_SF_4));
      ch->SetBranchAddress("MuID_SF_4_up", &(br.MuID_SF_4_up));
      ch->SetBranchAddress("MuID_SF_4_down", &(br.MuID_SF_4_down));
      ch->SetBranchAddress("MuIso_SF_4", &(br.MuIso_SF_4));
      ch->SetBranchAddress("MuIso_SF_4_up", &(br.MuIso_SF_4_up));
      ch->SetBranchAddress("MuIso_SF_4_down", &(br.MuIso_SF_4_down));
      ch->SetBranchAddress("IsoMu_SF_bug", &(br.IsoMu_SF_bug));
      ch->SetBranchAddress("IsoMu_SF_bug_up", &(br.IsoMu_SF_bug_up));
      ch->SetBranchAddress("IsoMu_SF_bug_down", &(br.IsoMu_SF_bug_down));
    }

  }

  if (verbose) std::cout << "Exiting SetBranchAddress\n" << std::endl;

} // End function: void SetBranchAddresses(TChain * ch)


// Convert "slim jet" collection (used in 2016) into regular jet collection
JetInfos ConvertSlimJets(SlimJetInfos & _slimJets) {

  JetInfos _jets;
  for (const auto & _slimJet : _slimJets) {
    JetInfo _jet;
    _jet.pt        = _slimJet.pt;
    _jet.eta       = _slimJet.eta;
    _jet.phi       = _slimJet.phi;
    _jet.mass      = _slimJet.mass;
    _jet.partonID  = _slimJet.partonID;

    _jet.jecFactor  = _slimJet.jecFactor;
    _jet.jecUnc     = _slimJet.jecUnc;

    _jet.CSV     = _slimJet.CSV;
    _jet.deepCSV = _slimJet.deepCSV;
    _jet.puID    = _slimJet.puID;
    
    _jets.push_back(_jet);
  } // End loop: for (const auto & _slimJet : _slimJets)

  return _jets;
} // End function: JetInfos ConvertSlimJets()
