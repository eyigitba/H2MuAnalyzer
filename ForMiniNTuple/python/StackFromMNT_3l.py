####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os
import sys

sys.path.insert(0, '%s/lib' % os.getcwd() )
from ROOT import *
from MNT_Helper import LoadColors, LinearStack, RatioPlot, FillHistTerm, GetSF
gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

#LABEL = 'miniNtuple_WH_2016_v5'  ## Sub-folder within PLOT_DIR containing histograms
#LABEL = 'WH_ele_loose_ID_loose_iso_loose_mu_iso_v1'
#LABEL = 'WH_mu_med_ID_loose_iso_v1'

LABEL = 'VH_selection_2019april/pt10_iso04/WH_ele_high_dimu_pt_ID_fix'
LEP = "ele"

def InitHists(histos, terms, signals, bkgs):
    for sample in signals + bkgs + ["data"]:
	histos["mu1_pt"][sample] 	= TH1F("mu1_pt" + "_" + sample, "mu1_pt" + "_" + sample,			50,0,800 )
	histos["mu2_pt"][sample] 	= TH1F("mu2_pt" + "_" + sample, "mu2_pt" + "_" + sample, 			50,0,400 )
	histos["mu1_abs_eta"][sample] 	= TH1F("mu1_abs_eta" + "_" + sample, "mu1_abs_eta" + "_" + sample, 		50,0,2.5)
        histos["mu2_abs_eta"][sample] 	= TH1F("mu2_abs_eta" + "_" + sample, "mu2_abs_eta" + "_" + sample, 		50,0,2.5)
	histos["mu1_lepMVA"][sample]	= TH1F("mu1_lepMVA" + "_" + sample, "mu1_lepMVA" + "_" + sample, 		50,-1,1)
	histos["mu2_lepMVA"][sample]    = TH1F("mu2_lepMVA" + "_" + sample, "mu2_lepMVA" + "_" + sample, 		50,-1,1)

	histos["dimu_mass"][sample] 	= TH1F("dimu_mass" + "_" + sample, "dimu_mass" + "_" + sample, 			12,100,160)
#	histos["dimu_mass"][sample]     = TH1F("dimu_mass" + "_" + sample, "dimu_mass" + "_" + sample,                  40,70,110)
        histos["dimu_pt"][sample] 	= TH1F("dimu_pt" + "_" + sample, "dimu_pt" + "_" + sample, 			50,0,1000)	
	histos["dimu_mass_err"][sample] = TH1F("dimu_mass_err" + "_" + sample, "dimu_mass_err" + "_" + sample, 		50,-10,10)
	histos["dimu_abs_eta"][sample] 	= TH1F("dimu_abs_eta" + "_" + sample, "dimu_abs_eta" + "_" + sample, 		50,0,6)
        histos["dimu_abs_dEta"][sample] = TH1F("dimu_abs_dEta" + "_" + sample, "dimu_abs_dEta" + "_" + sample, 		50,0,4)
        histos["dimu_abs_dPhi"][sample] = TH1F("dimu_abs_dPhi" + "_" + sample, "dimu_abs_dPhi" + "_" + sample, 		50,0,4)
        histos["dimu_dR"][sample] 	= TH1F("dimu_dR" + "_" + sample, "dimu_dR" + "_" + sample, 			50,0,5)
	histos["dimu_gen_ID"][sample]   = TH1F("dimu_gen_ID" + "_" + sample, "dimu_gen_ID" + "_" + sample,              30,0,30)
        histos["cts_mu1"][sample] 	= TH1F("cts_mu1" + "_" + sample, "cts_mu1" + "_" + sample, 			50,-1,1)
        histos["cts_mu_pos"][sample] 	= TH1F("cts_mu_pos" + "_" + sample, "cts_mu_pos" + "_" + sample, 		50,-1,1)

	histos["lep_pt"][sample] 	= TH1F("lep_pt" + "_" + sample, "lep_pt" + "_" + sample, 			50,0,500)
        histos["lep_abs_eta"][sample] 	= TH1F("lep_abs_eta" + "_" + sample, "lep_abs_eta" + "_" + sample, 		50,0,2.5)
	histos["lep_lepMVA"][sample]    = TH1F("lep_lepMVA" + "_" + sample, "lep_lepMVA" + "_" + sample,	 	50,-1,1)
        histos["cts_ldimu"][sample] 	= TH1F("cts_ldimu" + "_" + sample, "cts_ldimu" + "_" + sample, 			50,-1,1)
        histos["ldimu_mass"][sample] 	= TH1F("ldimu_mass" + "_" + sample, "ldimu_mass" + "_" + sample, 		50,0,1000)
        histos["ldimu_pt"][sample] 	= TH1F("ldimu_pt" + "_" + sample, "ldimu_pt" + "_" + sample, 			50,0,500)
        histos["ldimu_abs_eta"][sample] = TH1F("ldimu_abs_eta" + "_" + sample, "ldimu_abs_eta" + "_" + sample,  	50,0,8)
        histos["ldimu_abs_dEta"][sample] = TH1F("ldimu_abs_dEta" + "_" + sample, "ldimu_abs_dEta" + "_" + sample, 	50,0,6)
        histos["ldimu_abs_dPhi"][sample] = TH1F("ldimu_abs_dPhi" + "_" + sample, "ldimu_abs_dPhi" + "_" + sample, 	50,0,4)
        histos["ldimu_dR"][sample] 	= TH1F("ldimu_dR" + "_" + sample, "ldimu_dR" + "_" + sample, 			50,0,8)
        histos["cts_lmuSS"][sample] 	= TH1F("cts_lmuSS" + "_" + sample, "cts_lmuSS" + "_" + sample, 			50,-1,1)
        histos["cts_lmuOS"][sample] 	= TH1F("cts_lmuOS" + "_" + sample, "cts_lmuOS" + "_" + sample, 			50,-1,1)

        histos["lmuSS_pt"][sample] 	= TH1F("lmuSS_pt" + "_" + sample, "lmuSS_pt" + "_" + sample, 			50,0,500)
	histos["lmuSS_abs_eta"][sample] = TH1F("lmuSS_abs_eta" + "_" + sample, "lmuSS_abs_eta" + "_" + sample, 		50,0,6)
        histos["lmuSS_abs_dEta"][sample] = TH1F("lmuSS_abs_dEta" + "_" + sample, "lmuSS_abs_dEta" + "_" + sample, 	50,0,5)
        histos["lmuSS_abs_dPhi"][sample] = TH1F("lmuSS_abs_dPhi" + "_" + sample, "lmuSS_abs_dPhi" + "_" + sample, 	50,0,4)
        histos["lmuSS_dR"][sample] 	= TH1F("lmuSS_dR" + "_" + sample, "lmuSS_dR" + "_" + sample, 			50,0,6)
	histos["lmuSS_mass"][sample]    = TH1F("lmuSS_mass" + "_" + sample, "lmuSS_mass" + "_" + sample,                50,0,200)
        histos["lmuOS_pt"][sample] 	= TH1F("lmuOS_pt" + "_" + sample, "lmuOS_pt" + "_" + sample, 			50,0,500)
        histos["lmuOS_abs_eta"][sample] = TH1F("lmuOS_abs_eta" + "_" + sample, "lmuOS_abs_eta" + "_" + sample, 		50,0,6)
        histos["lmuOS_abs_dEta"][sample] = TH1F("lmuOS_abs_dEta" + "_" + sample, "lmuOS_abs_dEta" + "_" + sample, 	50,0,5)
        histos["lmuOS_abs_dPhi"][sample] = TH1F("lmuOS_abs_dPhi" + "_" + sample, "lmuOS_abs_dPhi" + "_" + sample, 	50,0,4)
        histos["lmuOS_dR"][sample] 	= TH1F("lmuOS_dR" + "_" + sample, "lmuOS_dR" + "_" + sample, 			50,0,6)
	histos["lmuOS_mass"][sample]    = TH1F("lmuOS_mass" + "_" + sample, "lmuOS_mass" + "_" + sample,                50,0,200)

        histos["met_pt"][sample] 	= TH1F("met_pt" + "_" + sample, "met_pt" + "_" + sample, 			50,0,500)
        histos["mt_lmet"][sample] 	= TH1F("mt_lmet" + "_" + sample, "mt_lmet" + "_" + sample, 			50,0,200)
        histos["abs_dPhi_lmet"][sample] = TH1F("abs_dPhi_lmet" + "_" + sample, "abs_dPhi_lmet" + "_" + sample, 		50,0,4)
        histos["mht_pt"][sample] 	= TH1F("mht_pt" + "_" + sample, "mht_pt" + "_" + sample, 			50,0,500)
        histos["mht_mass"][sample] 	= TH1F("mht_mass" + "_" + sample, "mht_mass" + "_" + sample, 			50,0,5000)
        histos["mt_lmht"][sample] 	= TH1F("mt_lmht" + "_" + sample, "mt_lmht" + "_" + sample, 			50,0,500)
        histos["abs_dPhi_lmht"][sample] = TH1F("abs_dPhi_lmht" + "_" + sample, "abs_dPhi_lmht" + "_" + sample, 		50,0,4)
        histos["mlt_pt"][sample] 	= TH1F("mlt_pt" + "_" + sample, "mlt_pt" + "_" + sample, 			50,0,500)
        histos["mt_lmlt"][sample] 	= TH1F("mt_lmlt" + "_" + sample, "mt_lmlt" + "_" + sample, 			50,-0.0001,0.0001)
        histos["abs_dPhi_lmlt"][sample] = TH1F("abs_dPhi_lmlt" + "_" + sample, "abs_dPhi_lmlt" + "_" + sample, 		50,0,4)

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


def main():
    out_name = "stack_plots.root"

    file_dir = PLOT_DIR+"/"+LABEL+"/"
    out_file = TFile( file_dir + "plots/" + out_name , "RECREATE")
    file_chain = TChain("tree","chain");
    file_chain.Add(file_dir + "all_samples.root")
#    file_chain.Add( file_dir + "signal.root")
#    file_chain.Add( file_dir + "bkg.root")

    terms = ["mu1_pt", "mu2_pt", "mu1_lepMVA", "mu1_abs_eta", "mu2_abs_eta", "mu2_lepMVA",  	#mu_vars
	    "dimu_mass", "dimu_pt", "dimu_mass_err", "dimu_abs_eta", "dimu_abs_dEta", 	#mu_vars
	    "dimu_abs_dPhi", "dimu_dR", "dimu_gen_ID", "cts_mu1", "cts_mu_pos",  	#mu_vars
	    "lep_pt", "lep_abs_eta", "lep_lepMVA", "cts_ldimu", "ldimu_mass", "ldimu_pt",  #lep_vars
	    "ldimu_abs_eta", "ldimu_abs_dEta", "ldimu_abs_dPhi", "ldimu_dR", #lep_vars
	    "cts_lmuSS", "cts_lmuOS", "lmuSS_pt", "lmuSS_abs_eta", 		#lmu vars
	    "lmuSS_abs_dEta", "lmuSS_abs_dPhi", "lmuSS_dR", "lmuSS_mass", "lmuOS_pt", 	#lmu vars
	    "lmuOS_abs_eta", "lmuOS_abs_dEta", "lmuOS_abs_dPhi", "lmuOS_dR", "lmuOS_mass", #lmu vars
	    "met_pt", "mt_lmet", "abs_dPhi_lmet", 		#met vars
	    "mht_pt", "mht_mass", "mt_lmht", "abs_dPhi_lmht",  	#met vars
	    "mlt_pt", "mt_lmlt", "abs_dPhi_lmlt",		#met vars
	    "dijet_mass", "dijet_pt", "dijet_abs_eta", "dijet_abs_dEta", "dijet_abs_dPhi", "dijet_dR", 	#jet vars
	    "jet1_pt", "jet1_abs_eta", "jet2_pt", "jet2_abs_eta", "jet0_pt", "jet0_abs_eta",		#jet vars
	    "nJets", "nCentJets", "nFwdJets", "nBJets_Med", "nBJets_Loose", "nBJets_Tight", "nMuons", "nEles", 	#evt vars
	    ] #end

    signals = ["ttH", "ZH", "WH", "VBF", "ggH"]
    bkgs = ["others", "triboson", "tZq", "tW", "ttZ", "ttbar", "WW", "ZZ", "WZ", "DY"]
    data = ["data"]

    color = {}
    LoadColors(color, "WH")

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
	

    InitHists(histos, terms, signals, bkgs)

    print file_chain.GetEntries()
    for iEvt in range( file_chain.GetEntries() ):
	file_chain.GetEvent(iEvt)
	if (iEvt % 10000 == 1):
	    print "looking at event %d" %iEvt

#	if file_chain.mu1_lepMVA < 0.4 or file_chain.mu2_lepMVA < 0.4 or file_chain.lep_lepMVA < 0.4:
#            continue
#        mu1_SF = GetSF("muon", file_chain.mu1_pt, file_chain.mu1_eta, 0.4)
#        mu2_SF = GetSF("muon", file_chain.mu2_pt, file_chain.mu2_eta, 0.4)
#        lep_SF = GetSF(LEP,    file_chain.lep_pt, file_chain.lep_eta, 0.4)
	MVA_SF = 1.0
#        MVA_SF = mu1_SF * mu2_SF * lep_SF

	FillHistTerm(histos, "mu1_pt"	      	, signals, bkgs, file_chain.mu1_pt		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu2_pt"		, signals, bkgs, file_chain.mu2_pt		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu1_abs_eta"   	, signals, bkgs, abs(file_chain.mu1_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu2_abs_eta"   	, signals, bkgs, abs(file_chain.mu2_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu1_lepMVA"	, signals, bkgs, file_chain.mu1_lepMVA		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mu2_lepMVA"	, signals, bkgs, file_chain.mu2_lepMVA		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_mass"   	, signals, bkgs, file_chain.dimu_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_pt"   	, signals, bkgs, file_chain.dimu_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_mass_err"	, signals, bkgs, file_chain.dimu_mass_err 	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_abs_eta"   	, signals, bkgs, abs(file_chain.dimu_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_abs_dEta"   	, signals, bkgs, abs(file_chain.dimu_dEta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_abs_dPhi"   	, signals, bkgs, abs(file_chain.dimu_dPhi)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dimu_dR"   	, signals, bkgs, file_chain.dimu_dR   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "dimu_gen_ID"      , signals, bkgs, file_chain.dimu_gen_ID         , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "cts_mu1"   	, signals, bkgs, file_chain.cts_mu1   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "cts_mu_pos"   	, signals, bkgs, file_chain.cts_mu_pos   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lep_pt"   	, signals, bkgs, file_chain.lep_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lep_abs_eta"   	, signals, bkgs, abs(file_chain.lep_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lep_lepMVA"	, signals, bkgs, file_chain.lep_lepMVA 		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "cts_ldimu"   	, signals, bkgs, file_chain.cts_ldimu   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "ldimu_mass"   	, signals, bkgs, file_chain.ldimu_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "ldimu_pt"   	, signals, bkgs, file_chain.ldimu_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "ldimu_abs_eta"   	, signals, bkgs, abs(file_chain.ldimu_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "ldimu_abs_dEta"   , signals, bkgs, abs(file_chain.ldimu_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "ldimu_abs_dPhi"   , signals, bkgs, abs(file_chain.ldimu_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "ldimu_dR"   	, signals, bkgs, file_chain.ldimu_dR  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "cts_lmuSS"   	, signals, bkgs, file_chain.cts_lmuSS   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "cts_lmuOS"   	, signals, bkgs, file_chain.cts_lmuOS   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuSS_pt"   	, signals, bkgs, file_chain.lmuSS_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuSS_abs_eta"   	, signals, bkgs, abs(file_chain.lmuSS_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lmuSS_abs_dEta"   , signals, bkgs, abs(file_chain.lmuSS_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuSS_abs_dPhi"   , signals, bkgs, abs(file_chain.lmuSS_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuSS_dR"   	, signals, bkgs, file_chain.lmuSS_dR   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lmuSS_mass"       , signals, bkgs, file_chain.lmuSS_mass          , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuOS_pt"   	, signals, bkgs, file_chain.lmuOS_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuOS_abs_eta"   	, signals, bkgs, abs(file_chain.lmuOS_eta)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lmuOS_abs_dEta"   , signals, bkgs, abs(file_chain.lmuOS_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuOS_abs_dPhi"   , signals, bkgs, abs(file_chain.lmuOS_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "lmuOS_dR"   	, signals, bkgs, file_chain.lmuOS_dR   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "lmuOS_mass"       , signals, bkgs, file_chain.lmuOS_mass          , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "met_pt"   	, signals, bkgs, file_chain.met_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mt_lmet"   	, signals, bkgs, file_chain.mt_lmet   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "abs_dPhi_lmet"   	, signals, bkgs, abs(file_chain.dPhi_lmet)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mht_pt"   	, signals, bkgs, file_chain.mht_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mht_mass"   	, signals, bkgs, file_chain.mht_mass   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mt_lmht"   	, signals, bkgs, file_chain.mt_lmht   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "abs_dPhi_lmht"  	, signals, bkgs, abs(file_chain.dPhi_lmht)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
	FillHistTerm(histos, "mlt_pt"   	, signals, bkgs, file_chain.mlt_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "mt_lmlt"  	, signals, bkgs, file_chain.mt_lmlt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "abs_dPhi_lmlt"   	, signals, bkgs, abs(file_chain.dPhi_lmlt)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_mass"   	, signals, bkgs, file_chain.dijet_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_pt"   	, signals, bkgs, file_chain.dijet_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_abs_eta"  	, signals, bkgs, abs(file_chain.dijet_eta)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_abs_dEta"   , signals, bkgs, abs(file_chain.dijet_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_abs_dPhi"   , signals, bkgs, abs(file_chain.dijet_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "dijet_dR"   	, signals, bkgs, file_chain.dijet_dR   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet1_pt"  	, signals, bkgs, file_chain.jet1_pt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet1_abs_eta"  	, signals, bkgs, abs(file_chain.jet1_eta)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet2_pt"   	, signals, bkgs, file_chain.jet2_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet2_abs_eta"   	, signals, bkgs, abs(file_chain.jet2_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "jet0_pt"   	, signals, bkgs, file_chain.jet0_pt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)	
        FillHistTerm(histos, "jet0_abs_eta"   	, signals, bkgs, abs(file_chain.jet0_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nJets"   		, signals, bkgs, file_chain.nJets  		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nCentJets"   	, signals, bkgs, file_chain.nCentJets   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nFwdJets"   	, signals, bkgs, file_chain.nFwdJets  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nBJets_Med" 	, signals, bkgs, file_chain.nBJets_Med 	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nBJets_Loose"   	, signals, bkgs, file_chain.nBJets_Loose  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nBJets_Tight"   	, signals, bkgs, file_chain.nBJets_Tight   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nMuons"   	, signals, bkgs, file_chain.nMuons  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
        FillHistTerm(histos, "nEles"  		, signals, bkgs, file_chain.nEles  		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)


    out_file.cd()
    scaled_signal = {}
    ratios = {}
    for term in terms:
	histos[term]["data"].SetMarkerStyle(20)
        for sample in signals:
            histos[term][sample].SetLineColor(color[sample])
            histos[term][sample].SetLineWidth(2)
            stack_sig[term].Add(histos[term][sample])

	histos[term]["signal"] = stack_sig[term].GetStack().Last()
        stack_all[term].Add(histos[term]["signal"])
        histos[term]["signal"].SetFillColor( kGray )
        histos[term]["signal"].SetLineWidth(0)
        for sample in bkgs:
            histos[term][sample].SetFillColor(color[sample])
            stack_all[term].Add(histos[term][sample])

	scaled_signal[term] =  histos[term]["signal"].Clone()
  	scaled_signal[term].Scale(100)
	scaled_signal[term].SetLineColor(kRed)
	scaled_signal[term].SetLineWidth(4)
	scaled_signal[term].SetFillStyle(0)


	legend[term] = TLegend(0.7,0.7,1,1)
	legend[term].AddEntry(histos[term]["data"], "data")
	legend[term].AddEntry(histos[term]["signal"], "signal sum")
	for sample in bkgs:
            legend[term].AddEntry(histos[term][sample], sample )
	legend[term].AddEntry(scaled_signal[term], "signal X100")
	LinearStack( term, stack_all[term], scaled_signal[term], histos[term]["data"], legend[term], PLOT_DIR+"/"+LABEL+"/plots")

#	ratios[term] = TGraphAsymmErrors()
#        ratios[term].Divide(histos[term]["data"], stack_all[term].GetStack().Last(), "pois")
#        ratios[term].SetName("ratiograph_"+term)
#	RatioPlot( term, stack_all[term], scaled_signal[term], histos[term]["data"], ratios[term], legend[term], PLOT_DIR+"/"+LABEL+"/plots" )

    out_file.Close()
main()











	
