####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os
import sys

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from MNT_Helper import LinearStack, RatioPlot, FillHistTerm, GetSF
import Plot_Configs as PC

gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/Run2/Histograms'

#LABEL = 'miniNtuple_WH_2016_v5'  ## Sub-folder within PLOT_DIR containing histograms
#LABEL = 'VH_selection_2019april/pt10_iso04/ZH_mu_massBDT'
LABEL = 'ZH_lep_2019_08_14'
#LEP = "muon"

def InitHists(histos, terms, signals, bkgs):
    for sample in signals + bkgs + ["Data"]:
	histos["mu1_pt"][sample] 	= TH1F("mu1_pt" + "_" + sample, "mu1_pt" + "_" + sample,			50,0,400 )
	histos["mu2_pt"][sample] 	= TH1F("mu2_pt" + "_" + sample, "mu2_pt" + "_" + sample, 			50,0,250 )
	histos["mu1_abs_eta"][sample] 	= TH1F("mu1_abs_eta" + "_" + sample, "mu1_abs_eta" + "_" + sample, 		50,0,2.5)
        histos["mu2_abs_eta"][sample] 	= TH1F("mu2_abs_eta" + "_" + sample, "mu2_abs_eta" + "_" + sample, 		50,0,2.5)
	histos["mu1_lepMVA"][sample]	= TH1F("mu1_lepMVA" + "_" + sample, "mu1_lepMVA" + "_" + sample, 		50,-1,1)
	histos["mu2_lepMVA"][sample]    = TH1F("mu2_lepMVA" + "_" + sample, "mu2_lepMVA" + "_" + sample, 		50,-1,1)
	histos["mu1_charge"][sample]	= TH1F("mu1_charge" + "_" + sample, "mu1_charge" + "_" + sample, 		8,-2,2)
	histos["mu2_charge"][sample]	= TH1F("mu2_charge" + "_" + sample, "mu2_charge" + "_" + sample, 		8,-2,2)

	histos["dimu_gen_ID"][sample]   = TH1F("dimu_gen_ID" + "_" + sample, "dimu_gen_ID" + "_" + sample,              30,0,30)
	histos["dimu_mass"][sample] 	= TH1F("dimu_mass" + "_" + sample, "dimu_mass" + "_" + sample, 			30,100,160)
        histos["dimu_pt"][sample] 	= TH1F("dimu_pt" + "_" + sample, "dimu_pt" + "_" + sample, 			50,0,400)	
	histos["dimu_mass_err"][sample] = TH1F("dimu_mass_err" + "_" + sample, "dimu_mass_err" + "_" + sample, 		50,0,10)
	histos["dimu_abs_eta"][sample] 	= TH1F("dimu_abs_eta" + "_" + sample, "dimu_abs_eta" + "_" + sample, 		50,0,6)
        histos["dimu_abs_dEta"][sample] = TH1F("dimu_abs_dEta" + "_" + sample, "dimu_abs_dEta" + "_" + sample, 		50,0,4)
        histos["dimu_abs_dPhi"][sample] = TH1F("dimu_abs_dPhi" + "_" + sample, "dimu_abs_dPhi" + "_" + sample, 		50,0,4)
        histos["dimu_dR"][sample] 	= TH1F("dimu_dR" + "_" + sample, "dimu_dR" + "_" + sample, 			50,0,5)
        histos["cts_mu1"][sample] 	= TH1F("cts_mu1" + "_" + sample, "cts_mu1" + "_" + sample, 			50,-1,1)
        histos["cts_mu_pos"][sample] 	= TH1F("cts_mu_pos" + "_" + sample, "cts_mu_pos" + "_" + sample, 		50,-1,1)

	histos["lep1_pt"][sample] 	= TH1F("lep1_pt" + "_" + sample, "lep1_pt" + "_" + sample, 			50,0,300)
        histos["lep1_abs_eta"][sample] 	= TH1F("lep1_abs_eta" + "_" + sample, "lep1_abs_eta" + "_" + sample, 		50,0,2.5)
	histos["lep1_lepMVA"][sample]   = TH1F("lep1_lepMVA" + "_" + sample, "lep1_lepMVA" + "_" + sample,	 	50,-1,1)
	histos["lep1_charge"][sample]	= TH1F("lep1_charge" + "_" + sample, "lep1_charge" + "_" + sample, 		8,-2,2)
	histos["lep2_pt"][sample]       = TH1F("lep2_pt" + "_" + sample, "lep1_pt" + "_" + sample,                      50,0,200)
        histos["lep2_abs_eta"][sample]  = TH1F("lep2_abs_eta" + "_" + sample, "lep2_abs_eta" + "_" + sample,            50,0,2.5)
        histos["lep2_lepMVA"][sample]   = TH1F("lep2_lepMVA" + "_" + sample, "lep2_lepMVA" + "_" + sample,              50,-1,1)
	histos["lep2_charge"][sample]	= TH1F("lep2_charge" + "_" + sample, "lep2_charge" + "_" + sample, 		8,-2,2)

	histos["lep_ID"][sample]        = TH1F("lep_ID" + "_" + sample, "lep_ID" + "_" + sample,                        20,0,20)
	histos["dilep_gen_ID"][sample]  = TH1F("dilep_gen_ID" + "_" + sample, "dilep_gen_ID" + "_" + sample,            30,0,30)
	histos["dilep_mass"][sample]	= TH1F("dilep_mass" + "_" + sample, "dilep_mass" + "_" + sample, 		30,70,110)
	histos["dilep_pt"][sample]	= TH1F("dilep_pt" + "_" + sample, "dilep_pt" + "_" + sample, 			50,0,300)
	histos["dilep_abs_eta"][sample]	= TH1F("dilep_abs_eta" + "_" + sample, "dilep_abs_eta" + "_" + sample, 		50,0,6)
	histos["dilep_abs_dEta"][sample]= TH1F("dilep_abs_dEta" + "_" + sample, "dilep_abs_dEta" + "_" + sample, 	50,0,4)
	histos["dilep_abs_dPhi"][sample]= TH1F("dilep_abs_dPhi" + "_" + sample, "dilep_abs_dPhi" + "_" + sample, 	50,0,4)
	histos["dilep_dR"][sample]	= TH1F("dilep_dR" + "_" + sample, "dilep_dR" + "_" + sample, 			50,0,5)
	histos["cts_lep1"][sample]	= TH1F("cts_lep1" + "_" + sample, "cts_lep1" + "_" + sample, 			50,-1,1)
	histos["cts_lep_pos"][sample]	= TH1F("cts_lep_pos" + "_" + sample, "cts_lep_pos" + "_" + sample, 		50,-1,1)

	histos["quadlep_mass"][sample]	= TH1F("quadlep_mass" + "_" + sample, "quadlep_mass" + "_" + sample, 		60,100,700)
	histos["quadlep_pt"][sample]	= TH1F("quadlep_pt" + "_" + sample, "quadlep_pt" + "_" + sample, 		50,0,400)
	histos["quadlep_abs_eta"][sample]= TH1F("quadlep_abs_eta" + "_" + sample, "quadlep_abs_eta" + "_" + sample, 	50,-10,10)
	histos["dipair_dEta_H"][sample]	= TH1F("dipair_dEta_H" + "_" + sample, "dipair_dEta_H" + "_" + sample, 		50,-10,10)
	histos["dipair_dPhi_H"][sample]	= TH1F("dipair_dPhi_H" + "_" + sample, "dipair_dPhi_H" + "_" + sample, 		50,-4,4)
	histos["dipair_dR_H"][sample]	= TH1F("dipair_dR_H" + "_" + sample, "dipair_dR_H" + "_" + sample, 		50,0,10)
	histos["cts_dipair_H"][sample]	= TH1F("cts_dipair_H" + "_" + sample, "cts_dipair_H" + "_" + sample, 		50,-1,1)

	histos["cs_costheta"][sample] 	= TH1F("cs_costheta" + "_" + sample, "cs_costheta" + "_" + sample, 		50,-1,1)
 	histos["cs_cosphi"][sample]	= TH1F("cs_cosphi" + "_" + sample, "cs_cosphi" + "_" + sample,			50,-1,1)
	histos["cs_sinphi"][sample]     = TH1F("cs_sinphi" + "_" + sample, "cs_sinphi" + "_" + sample,    		50,-1,1)
	histos["cos_theta1"][sample]    = TH1F("cos_theta1" + "_" + sample, "cos_theta1" + "_" + sample,    		50,-1,1)
	histos["cos_phi1"][sample]      = TH1F("cos_phi1" + "_" + sample, "cos_phi1" + "_" + sample,    		50,-1,1)
	histos["cos_phiH"][sample]      = TH1F("cos_phiH" + "_" + sample, "cos_phiH" + "_" + sample,    		50,-1,1)

        histos["met_pt"][sample] 	= TH1F("met_pt" + "_" + sample, "met_pt" + "_" + sample, 			50,0,250)
	histos["abs_dPhi_4lPt_met"][sample] = TH1F("abs_dPhi_4lPt_met" + "_" + sample, "abs_dPhi_4lPt_met" + "_" + sample,      50,0,4)
	histos["met_Rtrans_4lPt"][sample]   = TH1F("met_Rtrans_4lPt" + "_" + sample, "met_Rtrans_4lPt" + "_" + sample,          50,-1.5,1.5)
	histos["met_Rlongi_4lPt"][sample]   = TH1F("met_Rlongi_4lPt" + "_" + sample, "met_Rlongi_4lPt" + "_" + sample,          50,-1.5,1.5)
        histos["mht_pt"][sample] 	= TH1F("mht_pt" + "_" + sample, "mht_pt" + "_" + sample, 			50,0,250)
        histos["mht_mass"][sample] 	= TH1F("mht_mass" + "_" + sample, "mht_mass" + "_" + sample, 			50,0,5000)
	histos["abs_dPhi_4lPt_mht"][sample] = TH1F("abs_dPhi_4lPt_mht" + "_" + sample, "abs_dPhi_4lPt_mht" + "_" + sample,      50,0,4)
        histos["mht_Rtrans_4lPt"][sample]   = TH1F("mht_Rtrans_4lPt" + "_" + sample, "mht_Rtrans_4lPt" + "_" + sample,          50,-1.5,2.5)
        histos["mht_Rlongi_4lPt"][sample]   = TH1F("mht_Rlongi_4lPt" + "_" + sample, "mht_Rlongi_4lPt" + "_" + sample,          50,-1.5,1.5)
	histos["abs_dPhi_met_mht"][sample]  = TH1F("abs_dPhi_met_mht" + "_" + sample, "abs_dPhi_met_mht" + "_" + sample,        50,0,4)

        histos["dijet_mass"][sample] 	= TH1F("dijet_mass" + "_" + sample, "dijet_mass" + "_" + sample, 		50,0,2000)
        histos["dijet_pt"][sample] 	= TH1F("dijet_pt" + "_" + sample, "dijet_pt" + "_" + sample, 			50,0,1000)
        histos["dijet_abs_eta"][sample] = TH1F("dijet_abs_eta" + "_" + sample, "dijet_abs_eta" + "_" + sample, 		50,0,20)
        histos["dijet_abs_dEta"][sample] = TH1F("dijet_abs_dEta" + "_" + sample, "dijet_abs_dEta" + "_" + sample, 	50,0,20)
        histos["dijet_abs_dPhi"][sample] = TH1F("dijet_abs_dPhi" + "_" + sample, "dijet_abs_dPhi" + "_" + sample, 	50,0,4)
	histos["dijet_dR"][sample] 	= TH1F("dijet_dR" + "_" + sample, "dijet_dR" + "_" + sample, 			50,0,20)
        histos["jet1_pt"][sample] 	= TH1F("jet1_pt" + "_" + sample, "jet1_pt" + "_" + sample, 			50,0,1000)
        histos["jet1_abs_eta"][sample] 	= TH1F("jet1_abs_eta" + "_" + sample, "jet1_abs_eta" + "_" + sample, 		50,0,5)
        histos["jet2_pt"][sample] 	= TH1F("jet2_pt" + "_" + sample, "jet2_pt" + "_" + sample, 			50,0,500)
        histos["jet2_abs_eta"][sample] 	= TH1F("jet2_abs_eta" + "_" + sample, "jet2_abs_eta" + "_" + sample, 		50,0,5)
        histos["jet0_pt"][sample] 	= TH1F("jet0_pt" + "_" + sample, "jet0_pt" + "_" + sample, 			50,0,1000)
        histos["jet0_abs_eta"][sample] 	= TH1F("jet0_abs_eta" + "_" + sample, "jet0_abs_eta" + "_" + sample, 		50,0,5)

        histos["nJets"][sample] 	= TH1F("nJets" + "_" + sample, "nJets" + "_" + sample, 				10,0,10)
        histos["nCentJets"][sample] 	= TH1F("nCentJets" + "_" + sample, "nCentJets" + "_" + sample, 			10,0,10)
        histos["nFwdJets"][sample] 	= TH1F("nFwdJets" + "_" + sample, "nFwdJets" + "_" + sample, 			10,0,10)
        histos["nBJets_Med"][sample] 	= TH1F("nBJets_Med" + "_" + sample, "nBJets_Med" + "_" + sample, 		5,0,5)
        histos["nBJets_Loose"][sample] 	= TH1F("nBJets_Loose" + "_" + sample, "nBJets_Loose" + "_" + sample, 		5,0,5)
        histos["nBJets_Tight"][sample] 	= TH1F("nBJets_Tight" + "_" + sample, "nBJets_Tight" + "_" + sample, 		5,0,5)
        histos["nMuons"][sample] 	= TH1F("nMuons" + "_" + sample, "nMuons" + "_" + sample, 			5,0,5)
        histos["nEles"][sample] 	= TH1F("nEles" + "_" + sample, "nEles" + "_" + sample, 				5,0,5)
	histos["BDT_noMass_v3"][sample] = TH1F("BDT_noMass_v3" + "_" + sample, "BDT_noMass_v3" + "_" + sample,          50,-1,1)
	histos["BDT_mass_more"][sample] = TH1F("BDT_mass_more" + "_" + sample, "BDT_mass_more" + "_" + sample,          50,-1,1)
	histos["BDT_mass_min"][sample]  = TH1F("BDT_mass_min"  + "_" + sample, "BDT_mass_min"  + "_" + sample,          50,-1,1)
	histos["BDT_and_mass"][sample]  = TH1F("BDT_and_mass"  + "_" + sample, "BDT_and_mass"  + "_" + sample,          50,-1,1)


def PrintEventYield(histos, signals, bkgs):
    sig_num = 0
    ZZ_num = 0
    bkgs_num = 0
    range_start = histos["dimu_mass"]["ZH"].FindBin(120.1)
    range_end   = histos["dimu_mass"]["ZH"].FindBin(129.9)
    
    for sample in signals:
	sig_num += histos["dimu_mass"][sample].Integral(range_start, range_end) 
    for sample in bkgs:
	if sample == "ZZ":
	    ZZ_num += histos["dimu_mass"][sample].Integral(range_start, range_end)
	else:
	    bkgs_num += histos["dimu_mass"][sample].Integral(range_start, range_end)
    print "sig_num = %f" %sig_num + "      ZZ_num = %f" %ZZ_num + "      bkgs_num = %f" %bkgs_num

def main():
    out_name = "stack_plots.root"

    file_dir = PLOT_DIR+"/"+LABEL+"/"
    out_file = TFile( file_dir + "plots/" + out_name , "RECREATE")
    file_chain = TChain("tree","chain");
    file_chain.Add(file_dir + "all_samples.root")
#    file_chain.Add( file_dir + "signal.root")
#    file_chain.Add( file_dir + "bkg.root")

    terms = ["mu1_pt", "mu2_pt", "mu1_lepMVA", "mu1_charge", "mu1_abs_eta", "mu2_abs_eta", "mu2_lepMVA", "mu2_charge", 	#mu_vars
	    "dimu_mass", "dimu_pt", "dimu_mass_err", "dimu_abs_eta", "dimu_abs_dEta", 	#mu_vars
	    "dimu_abs_dPhi", "dimu_dR", "cts_mu1", "cts_mu_pos", "dimu_gen_ID", 	#mu_vars
	    "lep1_pt", "lep1_abs_eta", "lep1_lepMVA", "lep1_charge", 
	    "lep2_pt", "lep2_abs_eta", "lep2_lepMVA", "lep2_charge", #lep_vars
	    "lep_ID", "dilep_gen_ID",
	    "dilep_mass", "dilep_pt", "dilep_abs_eta", "dilep_abs_dEta", #lep_vars
	    "dilep_abs_dPhi", "dilep_dR", "cts_lep1", "cts_lep_pos", #lep_vars
	    "quadlep_mass", "quadlep_pt", "quadlep_abs_eta", "dipair_dEta_H", #quadlep_vars
	    "dipair_dPhi_H", "dipair_dR_H", "cts_dipair_H", #quadlep_vars
	    "cs_costheta", "cs_cosphi", "cs_sinphi", "cos_theta1", "cos_phi1", "cos_phiH", #angles
	    "met_pt", "abs_dPhi_4lPt_met", "met_Rtrans_4lPt", "met_Rlongi_4lPt", # met vars
	    "mht_pt", "mht_mass", "abs_dPhi_4lPt_mht", "mht_Rtrans_4lPt", "mht_Rlongi_4lPt", "abs_dPhi_met_mht", 	#met vars
	    "dijet_mass", "dijet_pt", "dijet_abs_eta", "dijet_abs_dEta", "dijet_abs_dPhi", "dijet_dR", 	#jet vars
	    "jet1_pt", "jet1_abs_eta", "jet2_pt", "jet2_abs_eta", "jet0_pt", "jet0_abs_eta",		#jet vars
	    "nJets", "nCentJets", "nFwdJets", "nBJets_Med", "nBJets_Loose", "nBJets_Tight", "nMuons", "nEles", 	#evt vars
	    "BDT_noMass_v3", "BDT_mass_more", "BDT_mass_min", "BDT_and_mass",
	    ] #end

    cfg = PC.Plot_Config("ZH_4l")
    print cfg.signals
    print cfg.bkgs

    histos = {}
    stack_all = {}
    stack_sig = {}
    stack_data = {}
    ratios = {}
    legend = {}

    for term in terms:
        histos[term] = {}
        stack_all[term] = THStack("h_stack_"+term, term)
        stack_sig[term] = THStack("sig_stack_"+term, "sig_"+term)
        stack_data[term] = THStack("data_stack"+term, "data_"+term)
	
    
    MVA_str = "mu1_lepMVA > -0.4 && mu2_lepMVA > -0.4 && lep1_lepMVA > -0.4 && lep2_lepMVA > -0.4"
    wgt_str = "event_wgt * xsec_norm * all_lepMVA_SF"
    for bdt_cut in [ -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]:
      file_chain.Draw("dimu_mass >> sig_p(20,123,127)", "(Sample_ID > 0 && BDT_noMass_v3 > %f && %s ) * %s" %(bdt_cut, MVA_str, wgt_str) )
      file_chain.Draw("dimu_mass >> sig_n(20,123,127)", "(Sample_ID > 0 && BDT_noMass_v3 < %f && %s ) * %s" %(bdt_cut, MVA_str, wgt_str) )
      file_chain.Draw("dimu_mass >> bkg_p(20,123,127)", "(Sample_ID < 0 && BDT_noMass_v3 > %f && %s ) * %s" %(bdt_cut, MVA_str, wgt_str) )
      file_chain.Draw("dimu_mass >> bkg_n(20,123,127)", "(Sample_ID < 0 && BDT_noMass_v3 < %f && %s ) * %s" %(bdt_cut, MVA_str, wgt_str) )

      sig_p_temp = gDirectory.Get("sig_p")
      sig_n_temp = gDirectory.Get("sig_n")
      bkg_p_temp = gDirectory.Get("bkg_p")
      bkg_n_temp = gDirectory.Get("bkg_n")

      print sig_p_temp.Integral()
      sens_temp = ( (sig_p_temp.Integral())**2.0 / bkg_p_temp.Integral() + (sig_n_temp.Integral())**2.0 / bkg_n_temp.Integral() ) ** 0.5
      print bdt_cut
      print sens_temp


    InitHists(histos, terms, cfg.signals, cfg.bkgs)

    print file_chain.GetEntries()
    for iEvt in range( file_chain.GetEntries() ):
	file_chain.GetEvent(iEvt)
	if (iEvt % 10000 == 1):
	    print "looking at event %d" %iEvt

	MVA_cut = -0.4
	if ( file_chain.mu1_lepMVA < MVA_cut or file_chain.mu2_lepMVA < MVA_cut or file_chain.lep1_lepMVA < MVA_cut or file_chain.lep2_lepMVA < MVA_cut):
	    continue
#	mu1_SF = GetSF("muon", file_chain.mu1_pt, file_chain.mu1_abs_eta, MVA_cut)
#        mu2_SF = GetSF("muon", file_chain.mu2_pt, file_chain.mu2_abs_eta, MVA_cut)
#        lep1_SF = GetSF(LEP,    file_chain.lep1_pt, file_chain.lep1_abs_eta, MVA_cut)
#	lep2_SF = GetSF(LEP,    file_chain.lep2_pt, file_chain.lep2_abs_eta, MVA_cut)
#        MVA_SF = mu1_SF * mu2_SF * lep1_SF * lep2_SF
	MVA_SF = file_chain.all_lepMVA_SF

	FillHistTerm(histos, "mu1_pt"	      	, cfg, file_chain.mu1_pt		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu2_pt"		, cfg, file_chain.mu2_pt		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu1_abs_eta"   	, cfg, file_chain.mu1_abs_eta   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu2_abs_eta"   	, cfg, file_chain.mu2_abs_eta   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu1_lepMVA"	, cfg, file_chain.mu1_lepMVA		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu2_lepMVA"	, cfg, file_chain.mu2_lepMVA		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu1_charge"	, cfg, file_chain.mu1_charge 		, file_chain.Sample_ID 	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu2_charge"	, cfg, file_chain.mu2_charge	 	, file_chain.Sample_ID 	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_gen_ID"      , cfg, file_chain.dimu_gen_ID           , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_mass"   	, cfg, file_chain.dimu_mass   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_pt"   	, cfg, file_chain.dimu_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_mass_err"	, cfg, file_chain.dimu_mass_err 	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_abs_eta"   	, cfg, file_chain.dimu_abs_eta   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_abs_dEta"   	, cfg, file_chain.dimu_abs_dEta   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_abs_dPhi"   	, cfg, file_chain.dimu_abs_dPhi   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_dR"   	, cfg, file_chain.dimu_dR   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cts_mu1"   	, cfg, file_chain.cts_mu1   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "cts_mu_pos"   	, cfg, file_chain.cts_mu_pos   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lep1_pt"   	, cfg, file_chain.lep1_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lep1_abs_eta"   	, cfg, file_chain.lep1_abs_eta   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lep1_lepMVA"	, cfg, file_chain.lep1_lepMVA 		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lep1_charge"	, cfg, file_chain.lep1_charge 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lep2_pt"          , cfg, file_chain.lep2_pt             	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lep2_abs_eta"     , cfg, file_chain.lep2_abs_eta        	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lep2_lepMVA"      , cfg, file_chain.lep2_lepMVA         	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
 	FillHistTerm(histos, "lep2_charge"	, cfg, file_chain.lep2_charge 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lep_ID"           , cfg, file_chain.lep_ID              	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dilep_gen_ID"     , cfg, file_chain.dilep_gen_ID        	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dilep_mass"	, cfg, file_chain.dilep_mass 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dilep_pt"		, cfg, file_chain.dilep_pt 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dilep_abs_eta"	, cfg, file_chain.dilep_abs_eta 	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dilep_abs_dEta"	, cfg, file_chain.dilep_abs_dEta	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dilep_abs_dPhi"	, cfg, file_chain.dilep_abs_dPhi	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dilep_dR"		, cfg, file_chain.dilep_dR 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cts_lep1"		, cfg, file_chain.cts_lep1 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cts_lep_pos"	, cfg, file_chain.cts_lep_pos 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "quadlep_mass"	, cfg, file_chain.quadlep_mass 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "quadlep_pt"	, cfg, file_chain.quadlep_pt 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "quadlep_abs_eta"	, cfg, file_chain.quadlep_abs_eta 	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dipair_dEta_H"	, cfg, file_chain.dipair_dEta_H 	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dipair_dPhi_H"	, cfg, file_chain.dipair_dPhi_H 	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dipair_dR_H"	, cfg, file_chain.dipair_dR_H 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cts_dipair_H"	, cfg, file_chain.cts_dipair_H 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cs_costheta"      , cfg, file_chain.cs_costheta 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cs_cosphi" 	, cfg, file_chain.cs_cosphi 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cs_sinphi" 	, cfg, file_chain.cs_sinphi 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cos_theta1" 	, cfg, file_chain.cos_theta1 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cos_phi1" 	, cfg, file_chain.cos_phi1 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cos_phiH" 	, cfg, file_chain.cos_phiH 		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "met_pt"   	, cfg, file_chain.met_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "abs_dPhi_4lPt_met", cfg, file_chain.abs_dPhi_4lPt_met     , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "met_Rtrans_4lPt"  , cfg, file_chain.met_Rtrans_4lPt     	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "met_Rlongi_4lPt"  , cfg, file_chain.met_Rlongi_4lPt     	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mht_pt"   	, cfg, file_chain.mht_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mht_mass"   	, cfg, file_chain.mht_mass   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "abs_dPhi_4lPt_mht", cfg, file_chain.abs_dPhi_4lPt_mht   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mht_Rtrans_4lPt"  , cfg, file_chain.mht_Rtrans_4lPt     	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mht_Rlongi_4lPt"  , cfg, file_chain.mht_Rlongi_4lPt     	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "abs_dPhi_met_mht" , cfg, file_chain.abs_dPhi_met_mht    	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_mass"   	, cfg, file_chain.dijet_mass   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_pt"   	, cfg, file_chain.dijet_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_abs_eta"  	, cfg, file_chain.dijet_abs_eta  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_abs_dEta"   , cfg, file_chain.dijet_abs_dEta    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_abs_dPhi"   , cfg, file_chain.dijet_abs_dPhi    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_dR"   	, cfg, file_chain.dijet_dR   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet1_pt"  	, cfg, file_chain.jet1_pt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet1_abs_eta"  	, cfg, file_chain.jet1_abs_eta  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet2_pt"   	, cfg, file_chain.jet2_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet2_abs_eta"   	, cfg, file_chain.jet2_abs_eta   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet0_pt"   	, cfg, file_chain.jet0_pt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)	
        FillHistTerm(histos, "jet0_abs_eta"   	, cfg, file_chain.jet0_abs_eta   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nJets"   		, cfg, file_chain.nJets  		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nCentJets"   	, cfg, file_chain.nCentJets   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nFwdJets"   	, cfg, file_chain.nFwdJets  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nBJets_Med" 	, cfg, file_chain.nBJets_Med 	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nBJets_Loose"   	, cfg, file_chain.nBJets_Loose  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nBJets_Tight"   	, cfg, file_chain.nBJets_Tight   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nMuons"   	, cfg, file_chain.nMuons  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nEles"  		, cfg, file_chain.nEles  		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "BDT_noMass_v3"    , cfg, file_chain.BDT_noMass_v3         , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "BDT_mass_more"    , cfg, file_chain.BDT_mass_more         , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "BDT_mass_min"     , cfg, file_chain.BDT_mass_min          , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "BDT_and_mass"     , cfg, file_chain.BDT_and_mass          , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)

    out_file.cd()
    scaled_signal = {}

    #PrintEventYield(histos, signals, bkgs)

    for term in terms:
	histos[term]["Data"].SetMarkerStyle(20)
        for sample in cfg.signals:
            histos[term][sample].SetLineColor(cfg.colors[sample])
            histos[term][sample].SetLineWidth(2)
            stack_sig[term].Add(histos[term][sample])

	histos[term]["signal"] = stack_sig[term].GetStack().Last()
        stack_all[term].Add(histos[term]["signal"])
        histos[term]["signal"].SetFillColor( kGray )
        histos[term]["signal"].SetLineWidth(0)
        for sample in cfg.bkgs:
            histos[term][sample].SetFillColor(cfg.colors[sample])
            stack_all[term].Add(histos[term][sample])

	scaled_signal[term] =  histos[term]["signal"].Clone()
  	scaled_signal[term].Scale(20)
	scaled_signal[term].SetLineColor(kRed)
	scaled_signal[term].SetLineWidth(2)
	scaled_signal[term].SetFillStyle(0)

	legend[term] = TLegend(0.7,0.7,1,1)
	legend[term].AddEntry(histos[term]["Data"], "Data")
	legend[term].AddEntry(histos[term]["signal"], "signal sum")
	for sample in cfg.bkgs:
            legend[term].AddEntry(histos[term][sample], sample )
	legend[term].AddEntry(scaled_signal[term], "signal X20")
	LinearStack( term, stack_all[term], scaled_signal[term], histos[term]["Data"], legend[term], PLOT_DIR+"/"+LABEL+"/plots")


#	ratios[term] = TGraphAsymmErrors()
#        ratios[term].Divide(histos[term]["Data"], stack_all[term].GetStack().Last(), "pois")
#        ratios[term].SetName("ratiograph_"+term)
#	RatioPlot( term, stack_all[term], scaled_signal[term], histos[term]["Data"], ratios[term], legend[term], PLOT_DIR+"/"+LABEL+"/plots" )

    out_file.Close()
main()











	
