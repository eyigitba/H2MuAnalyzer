####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from MNT_Helper import LinearStack, RatioPlot, PassLepMVA, FillHistTerm
#R.gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

#LABEL = 'miniNtuple_WH_2016_v5'  ## Sub-folder within PLOT_DIR containing histograms
LABEL = 'ZH_4l_ele_angle_test_v1'


def InitHists(histos, terms, signals, bkgs):
    for sample in signals + bkgs + ["data"]:
	histos["mu1_pt"][sample] 	= TH1F("mu1_pt" + "_" + sample, "mu1_pt" + "_" + sample,			50,0,400 )
	histos["mu2_pt"][sample] 	= TH1F("mu2_pt" + "_" + sample, "mu2_pt" + "_" + sample, 			50,0,250 )
	histos["mu1_abs_eta"][sample] 	= TH1F("mu1_abs_eta" + "_" + sample, "mu1_abs_eta" + "_" + sample, 		50,0,2.5)
        histos["mu2_abs_eta"][sample] 	= TH1F("mu2_abs_eta" + "_" + sample, "mu2_abs_eta" + "_" + sample, 		50,0,2.5)
	histos["mu1_lepMVA"][sample]	= TH1F("mu1_lepMVA" + "_" + sample, "mu1_lepMVA" + "_" + sample, 		50,-1,1)
	histos["mu2_lepMVA"][sample]    = TH1F("mu2_lepMVA" + "_" + sample, "mu2_lepMVA" + "_" + sample, 		50,-1,1)
	histos["mu1_charge"][sample]	= TH1F("mu1_charge" + "_" + sample, "mu1_charge" + "_" + sample, 		8,-2,2)
	histos["mu2_charge"][sample]	= TH1F("mu2_charge" + "_" + sample, "mu2_charge" + "_" + sample, 		8,-2,2)

	histos["dimu_mass"][sample] 	= TH1F("dimu_mass" + "_" + sample, "dimu_mass" + "_" + sample, 			12,100,160)
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

	histos["dilep_mass"][sample]	= TH1F("dilep_mass" + "_" + sample, "dilep_mass" + "_" + sample, 		30,70,100)
	histos["dilep_pt"][sample]	= TH1F("dilep_pt" + "_" + sample, "dilep_pt" + "_" + sample, 			50,0,300)
	histos["dilep_abs_eta"][sample]	= TH1F("dilep_abs_eta" + "_" + sample, "dilep_abs_eta" + "_" + sample, 		50,0,6)
	histos["dilep_abs_dEta"][sample]= TH1F("dilep_abs_dEta" + "_" + sample, "dilep_abs_dEta" + "_" + sample, 	50,0,4)
	histos["dilep_abs_dPhi"][sample]= TH1F("dilep_abs_dPhi" + "_" + sample, "dilep_abs_dPhi" + "_" + sample, 	50,0,4)
	histos["dilep_dR"][sample]	= TH1F("dilep_dR" + "_" + sample, "dilep_dR" + "_" + sample, 			50,0,5)
	histos["cts_lep1"][sample]	= TH1F("cts_lep1" + "_" + sample, "cts_lep1" + "_" + sample, 			50,-1,1)
	histos["cts_lep_pos"][sample]	= TH1F("cts_lep_pos" + "_" + sample, "cts_lep_pos" + "_" + sample, 		50,-1,1)

	histos["quadlep_mass"][sample]	= TH1F("quadlep_mass" + "_" + sample, "quadlep_mass" + "_" + sample, 		50,150,650)
	histos["quadlep_pt"][sample]	= TH1F("quadlep_pt" + "_" + sample, "quadlep_pt" + "_" + sample, 		50,0,400)
	histos["quadlep_eta"][sample]	= TH1F("quadlep_eta" + "_" + sample, "quadlep_eta" + "_" + sample, 		50,-10,10)
	histos["dipair_dEta_H"][sample]	= TH1F("dipair_dEta_H" + "_" + sample, "dipair_dEta_H" + "_" + sample, 		50,-10,10)
	histos["dipair_dPhi_H"][sample]	= TH1F("dipair_dPhi_H" + "_" + sample, "dipair_dPhi_H" + "_" + sample, 		50,-4,4)
	histos["dipair_dR_H"][sample]	= TH1F("dipair_dR_H" + "_" + sample, "dipair_dR_H" + "_" + sample, 		50,0,10)
	histos["dipair_dEta_pt"][sample]= TH1F("dipair_dEta_pt" + "_" + sample, "dipair_dEta_pt" + "_" + sample, 	50,-10,10)
	histos["dipair_dPhi_pt"][sample]= TH1F("dipair_dPhi_pt" + "_" + sample, "dipair_dPhi_pt" + "_" + sample, 	50,-4,4)
	histos["dipair_dR_pt"][sample]	= TH1F("dipair_dR_pt" + "_" + sample, "dipair_dR_pt" + "_" + sample, 		50,0,10)
	histos["cts_dipair_H"][sample]	= TH1F("cts_dipair_H" + "_" + sample, "cts_dipair_H" + "_" + sample, 		50,-1,1)
	histos["cts_dipair_pt"][sample]	= TH1F("cts_dipair_pt" + "_" + sample, "cts_dipair_pt" + "_" + sample, 		50,-1,1)

	histos["cs_costheta"][sample] 	= TH1F("cs_costheta" + "_" + sample, "cs_costheta" + "_" + sample, 		50,-1,1)
 	histos["cs_cosphi"][sample]	= TH1F("cs_cosphi" + "_" + sample, "cs_cosphi" + "_" + sample,			50,-1,1)
	histos["cs_sinphi"][sample]     = TH1F("cs_sinphi" + "_" + sample, "cs_sinphi" + "_" + sample,    		50,-1,1)
	histos["cos_theta1"][sample]    = TH1F("cos_theta1" + "_" + sample, "cos_theta1" + "_" + sample,    		50,-1,1)
	histos["cos_phi1"][sample]      = TH1F("cos_phi1" + "_" + sample, "cos_phi1" + "_" + sample,    		50,-1,1)
	histos["cos_phiH"][sample]      = TH1F("cos_phiH" + "_" + sample, "cos_phiH" + "_" + sample,    		50,-1,1)

        histos["met_pt"][sample] 	= TH1F("met_pt" + "_" + sample, "met_pt" + "_" + sample, 			50,0,250)
        histos["mht_pt"][sample] 	= TH1F("mht_pt" + "_" + sample, "mht_pt" + "_" + sample, 			50,0,250)
        histos["mht_mass"][sample] 	= TH1F("mht_mass" + "_" + sample, "mht_mass" + "_" + sample, 			50,0,5000)
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
	    "dimu_abs_dPhi", "dimu_dR", "cts_mu1", "cts_mu_pos",  	#mu_vars
	    "lep1_pt", "lep1_abs_eta", "lep1_lepMVA", "lep1_charge", "lep2_pt", "lep2_abs_eta", "lep2_lepMVA", "lep2_charge", #lep_vars
	    "dilep_mass", "dilep_pt", "dilep_abs_eta", "dilep_abs_dEta", #lep_vars
	    "dilep_abs_dPhi", "dilep_dR", "cts_lep1", "cts_lep_pos", #lep_vars
	    "quadlep_mass", "quadlep_pt", "quadlep_eta", "dipair_dEta_H", "dipair_dPhi_H", "dipair_dR_H", #quadlep_vars
	    "dipair_dEta_pt", "dipair_dPhi_pt", "dipair_dR_pt", "cts_dipair_H", "cts_dipair_pt", #quadlep_vars
	    "cs_costheta", "cs_cosphi", "cs_sinphi", "cos_theta1", "cos_phi1", "cos_phiH", #angles
	    "met_pt", "mht_pt", "mht_mass",  	#met vars
	    "dijet_mass", "dijet_pt", "dijet_abs_eta", "dijet_abs_dEta", "dijet_abs_dPhi", "dijet_dR", 	#jet vars
	    "jet1_pt", "jet1_abs_eta", "jet2_pt", "jet2_abs_eta", "jet0_pt", "jet0_abs_eta",		#jet vars
	    "nJets", "nCentJets", "nFwdJets", "nBJets_Med", "nBJets_Loose", "nBJets_Tight", "nMuons", "nEles", 	#evt vars
	    ] #end

    signals = ["ttH", "ZH", "WH", "VBF", "ggH"]
    bkgs = ["triboson", "ttZ", "tZq", "tX+ttX", "ttbar", "WW", "ZZ", "WZ", "DY"]
    data = ["data"]

    color = {}
    color["data"] = kBlack

    color["ggH"] = kRed
    color["VBF"] =  kBlue + 1
    color["ZH"] =  kOrange + 7
    color["WH"] = kGreen + 2
#    color["WH_neg"] = kViolet + 1
    color["ttH"]  = kPink + 6

    color["DY"] =  kAzure + 7
    color["WZ"] =       kGreen - 9
    color["ZZ"] =  kCyan - 7
    color["WW"] =  kBlue - 8
    color["tX+ttX"] = kPink + 6 #
    color["ttbar"] = kYellow - 9
    color["tZq"] = kViolet -9
    color["ttZ"] = kRed -7
    color["triboson"] = kOrange + 6 #

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

	if not (PassLepMVA("0.4", file_chain.mu1_lepMVA) and PassLepMVA("0.4", file_chain.mu2_lepMVA) and PassLepMVA("0.4", file_chain.lep1_lepMVA) and PassLepMVA("0.4", file_chain.lep2_lepMVA)): 
		continue

	FillHistTerm(histos, "mu1_pt"	      	, signals, bkgs, file_chain.mu1_pt		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu2_pt"		, signals, bkgs, file_chain.mu2_pt		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu1_abs_eta"   	, signals, bkgs, abs(file_chain.mu1_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu2_abs_eta"   	, signals, bkgs, abs(file_chain.mu2_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu1_lepMVA"	, signals, bkgs, file_chain.mu1_lepMVA		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu2_lepMVA"	, signals, bkgs, file_chain.mu2_lepMVA		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu1_charge"	, signals, bkgs, file_chain.mu1_charge 		, file_chain.Sample_ID 	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu2_charge"	, signals, bkgs, file_chain.mu2_charge	 	, file_chain.Sample_ID 	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dimu_mass"   	, signals, bkgs, file_chain.dimu_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dimu_pt"   	, signals, bkgs, file_chain.dimu_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dimu_mass_err"	, signals, bkgs, file_chain.dimu_mass_err 	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_abs_eta"   	, signals, bkgs, abs(file_chain.dimu_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_abs_dEta"   	, signals, bkgs, abs(file_chain.dimu_dEta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_abs_dPhi"   	, signals, bkgs, abs(file_chain.dimu_dPhi)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_dR"   	, signals, bkgs, file_chain.dimu_dR   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cts_mu1"   	, signals, bkgs, file_chain.cts_mu1   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "cts_mu_pos"   	, signals, bkgs, file_chain.cts_mu_pos   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "lep1_pt"   	, signals, bkgs, file_chain.lep1_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "lep1_abs_eta"   	, signals, bkgs, abs(file_chain.lep1_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "lep1_lepMVA"	, signals, bkgs, file_chain.lep1_lepMVA 	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "lep1_charge"	, signals, bkgs, file_chain.lep1_charge 	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "lep2_pt"          , signals, bkgs, file_chain.lep2_pt             , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "lep2_abs_eta"     , signals, bkgs, abs(file_chain.lep2_eta)       , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "lep2_lepMVA"      , signals, bkgs, file_chain.lep2_lepMVA         , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
 	FillHistTerm(histos, "lep2_charge"	, signals, bkgs, file_chain.lep2_charge 	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dilep_mass"	, signals, bkgs, file_chain.dilep_mass 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dilep_pt"		, signals, bkgs, file_chain.dilep_pt 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dilep_abs_eta"	, signals, bkgs, abs(file_chain.dilep_eta) 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dilep_abs_dEta"	, signals, bkgs, abs(file_chain.dilep_dEta)	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dilep_abs_dPhi"	, signals, bkgs, abs(file_chain.dilep_dPhi)	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dilep_dR"		, signals, bkgs, file_chain.dilep_dR 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cts_lep1"		, signals, bkgs, file_chain.cts_lep1 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cts_lep_pos"	, signals, bkgs, file_chain.cts_lep_pos 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "quadlep_mass"	, signals, bkgs, file_chain.quadlep_mass 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "quadlep_pt"	, signals, bkgs, file_chain.quadlep_pt 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "quadlep_eta"	, signals, bkgs, file_chain.quadlep_eta 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dipair_dEta_H"	, signals, bkgs, file_chain.dipair_dEta_H 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dipair_dPhi_H"	, signals, bkgs, file_chain.dipair_dPhi_H 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dipair_dR_H"	, signals, bkgs, file_chain.dipair_dR_H 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dipair_dEta_pt"	, signals, bkgs, file_chain.dipair_dEta_pt 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dipair_dPhi_pt"	, signals, bkgs, file_chain.dipair_dPhi_pt 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dipair_dR_pt"	, signals, bkgs, file_chain.dipair_dR_pt 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cts_dipair_H"	, signals, bkgs, file_chain.cts_dipair_H 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cts_dipair_pt"	, signals, bkgs, file_chain.cts_dipair_pt 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cs_costheta"      , signals, bkgs, file_chain.cs_costheta 	, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cs_cosphi" 	, signals, bkgs, file_chain.cs_cosphi 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cs_sinphi" 	, signals, bkgs, file_chain.cs_sinphi 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cos_theta1" 	, signals, bkgs, file_chain.cos_theta1 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cos_phi1" 	, signals, bkgs, file_chain.cos_phi1 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cos_phiH" 	, signals, bkgs, file_chain.cos_phiH 		, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "met_pt"   	, signals, bkgs, file_chain.met_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "mht_pt"   	, signals, bkgs, file_chain.mht_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "mht_mass"   	, signals, bkgs, file_chain.mht_mass   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dijet_mass"   	, signals, bkgs, file_chain.dijet_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dijet_pt"   	, signals, bkgs, file_chain.dijet_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dijet_abs_eta"  	, signals, bkgs, abs(file_chain.dijet_eta)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dijet_abs_dEta"   , signals, bkgs, abs(file_chain.dijet_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dijet_abs_dPhi"   , signals, bkgs, abs(file_chain.dijet_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dijet_dR"   	, signals, bkgs, file_chain.dijet_dR   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "jet1_pt"  	, signals, bkgs, file_chain.jet1_pt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "jet1_abs_eta"  	, signals, bkgs, abs(file_chain.jet1_eta)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "jet2_pt"   	, signals, bkgs, file_chain.jet2_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "jet2_abs_eta"   	, signals, bkgs, abs(file_chain.jet2_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "jet0_pt"   	, signals, bkgs, file_chain.jet0_pt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)	
        FillHistTerm(histos, "jet0_abs_eta"   	, signals, bkgs, abs(file_chain.jet0_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nJets"   		, signals, bkgs, file_chain.nJets  		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nCentJets"   	, signals, bkgs, file_chain.nCentJets   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nFwdJets"   	, signals, bkgs, file_chain.nFwdJets  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nBJets_Med" 	, signals, bkgs, file_chain.nBJets_Med 	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nBJets_Loose"   	, signals, bkgs, file_chain.nBJets_Loose  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nBJets_Tight"   	, signals, bkgs, file_chain.nBJets_Tight   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nMuons"   	, signals, bkgs, file_chain.nMuons  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "nEles"  		, signals, bkgs, file_chain.nEles  		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)


    out_file.cd()
    scaled_signal = {}

    PrintEventYield(histos, signals, bkgs)

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
  	scaled_signal[term].Scale(200)
	scaled_signal[term].SetLineColor(kRed)
	scaled_signal[term].SetLineWidth(2)
	scaled_signal[term].SetFillStyle(0)


	legend[term] = TLegend(0.7,0.7,1,1)
	legend[term].AddEntry(histos[term]["data"], "data")
	legend[term].AddEntry(histos[term]["signal"], "signal sum")
	for sample in bkgs:
            legend[term].AddEntry(histos[term][sample], sample )
	legend[term].AddEntry(scaled_signal[term], "signal X200")
	LinearStack( term, stack_all[term], scaled_signal[term], histos[term]["data"], legend[term], PLOT_DIR+"/"+LABEL+"/plots")


    out_file.Close()
main()











	
