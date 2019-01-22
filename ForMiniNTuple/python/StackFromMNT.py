####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os

from ROOT import *
#R.gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

LABEL = 'miniNtuple_WH_2016_v5'  ## Sub-folder within PLOT_DIR containing histograms



def LinearStack( term, all_stack, scaled_signal, legend):
    canv = TCanvas("Stack_" + term, "Stack_" + term, 600,600)
    canv.Clear()

    canv.cd()
    all_stack.Draw("HIST")
    scaled_signal.Draw("HISTSAME")
    legend.Draw()
    canv.Update()
    canv.Write()
    canv.SaveAs(PLOT_DIR+"/"+LABEL+"/files" + "/plots" + "/stack_" + term + ".png")


def ratioplot( term, sig_stack, all_stack, h_data, ratio_graph, legend ):   # Do not use TRatioPlot! It is a devil!    -XWZ 19.09.2018
    canv = TCanvas("ratio_"+term, "ratio_"+term, 600,600)
    canv.Clear()

    upper_pad = TPad("upperpad_"+term, "upperpad_"+term, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.Draw()
    upper_pad.cd()
    upper_pad.SetLogy()
#    all_stack.SetMinimum(1e-5)
#    all_stack.SetMaximum(1e6)
    all_stack.SetMinimum(1e-3)
    all_stack.SetMaximum(1e8)
    all_stack.Draw("HIST")
    h_data.SetMarkerStyle(20)
    h_data.Draw("SAME")
    for histo in sig_stack:
        histo.Draw("SAMEHIST")
    legend.Draw()

    canv.cd()
    lower_pad = TPad("lowerpad_"+term, "lowerpad_"+term, 0,0.05, 1,0.2)
    lower_pad.SetTopMargin(0.05)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    ratio_graph.SetMinimum(0.5)
    ratio_graph.SetMaximum(1.5)
#    ratio_graph.SetMinimum(0.5)
#    ratio_graph.SetMaximum(2.5)

    ratio_graph.GetXaxis().SetRangeUser( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_graph.SetMarkerStyle(20)
#    ratio_graph.GetYaxis().SetNdivisions(510)
    ratio_graph.GetYaxis().SetNdivisions(505)
    ratio_graph.GetXaxis().SetLabelSize(0.15)
    ratio_graph.GetYaxis().SetLabelSize(0.13)
    ratio_graph.GetYaxis().SetTitle("data/MC")
    ratio_graph.GetYaxis().SetTitleSize(0.20)
    ratio_graph.GetYaxis().SetTitleOffset(0.2)
    ratio_graph.Draw()

    canv.Update()
    canv.Write()
    canv.SaveAs(PLOT_DIR+"/"+LABEL+"/files/sum" + "/plots_none" + "/ratio_" + term + ".png")



def InitHists(histos, terms, signals, bkgs):
    for sample in signals + bkgs:
	histos["mu1_pt"][sample] 	= TH1F("mu1_pt" + "_" + sample, "mu1_pt" + "_" + sample,			50,0,800 )
	histos["mu2_pt"][sample] 	= TH1F("mu2_pt" + "_" + sample, "mu2_pt" + "_" + sample, 			50,0,400 )
	histos["mu1_abs_eta"][sample] 	= TH1F("mu1_abs_eta" + "_" + sample, "mu1_abs_eta" + "_" + sample, 		50,0,2.5)
        histos["mu2_abs_eta"][sample] 	= TH1F("mu2_abs_eta" + "_" + sample, "mu2_abs_eta" + "_" + sample, 		50,0,2.5)
	histos["dimu_mass"][sample] 	= TH1F("dimu_mass" + "_" + sample, "dimu_mass" + "_" + sample, 			12,100,160)
        histos["dimu_pt"][sample] 	= TH1F("dimu_pt" + "_" + sample, "dimu_pt" + "_" + sample, 			50,0,1000)	
	histos["dimu_abs_eta"][sample] 	= TH1F("dimu_abs_eta" + "_" + sample, "dimu_abs_eta" + "_" + sample, 		50,0,6)
        histos["dimu_abs_dEta"][sample] = TH1F("dimu_abs_dEta" + "_" + sample, "dimu_abs_dEta" + "_" + sample, 		50,0,4)
        histos["dimu_abs_dPhi"][sample] = TH1F("dimu_abs_dPhi" + "_" + sample, "dimu_abs_dPhi" + "_" + sample, 		50,0,4)
        histos["dimu_dR"][sample] 	= TH1F("dimu_dR" + "_" + sample, "dimu_dR" + "_" + sample, 			50,0,5)
        histos["cts_mu1"][sample] 	= TH1F("cts_mu1" + "_" + sample, "cts_mu1" + "_" + sample, 			50,-1,1)
        histos["cts_mu_pos"][sample] 	= TH1F("cts_mu_pos" + "_" + sample, "cts_mu_pos" + "_" + sample, 		50,-1,1)
	histos["ele_pt"][sample] 	= TH1F("ele_pt" + "_" + sample, "ele_pt" + "_" + sample, 			50,0,500)
        histos["ele_abs_eta"][sample] 	= TH1F("ele_abs_eta" + "_" + sample, "ele_abs_eta" + "_" + sample, 		50,0,2.5)
        histos["cts_edimu"][sample] 	= TH1F("cts_edimu" + "_" + sample, "cts_edimu" + "_" + sample, 			50,-1,1)
        histos["edimu_mass"][sample] 	= TH1F("edimu_mass" + "_" + sample, "edimu_mass" + "_" + sample, 		50,0,1000)
        histos["edimu_pt"][sample] 	= TH1F("edimu_pt" + "_" + sample, "edimu_pt" + "_" + sample, 			50,0,500)
        histos["edimu_abs_eta"][sample] = TH1F("edimu_abs_eta" + "_" + sample, "edimu_abs_eta" + "_" + sample,  	50,0,8)
        histos["edimu_abs_dEta"][sample] = TH1F("edimu_abs_dEta" + "_" + sample, "edimu_abs_dEta" + "_" + sample, 	50,0,6)
        histos["edimu_abs_dPhi"][sample] = TH1F("edimu_abs_dPhi" + "_" + sample, "edimu_abs_dPhi" + "_" + sample, 	50,0,4)
        histos["edimu_dR"][sample] 	= TH1F("edimu_dR" + "_" + sample, "edimu_dR" + "_" + sample, 			50,0,8)
        histos["cts_emuSS"][sample] 	= TH1F("cts_emuSS" + "_" + sample, "cts_emuSS" + "_" + sample, 			50,-1,1)
        histos["cts_emuOS"][sample] 	= TH1F("cts_emuOS" + "_" + sample, "cts_emuOS" + "_" + sample, 			50,-1,1)
        histos["emuSS_pt"][sample] 	= TH1F("emuSS_pt" + "_" + sample, "emuSS_pt" + "_" + sample, 			50,0,500)
	histos["emuSS_abs_eta"][sample] = TH1F("emuSS_abs_eta" + "_" + sample, "emuSS_abs_eta" + "_" + sample, 		50,0,6)
        histos["emuSS_abs_dEta"][sample] = TH1F("emuSS_abs_dEta" + "_" + sample, "emuSS_abs_dEta" + "_" + sample, 	50,0,5)
        histos["emuSS_abs_dPhi"][sample] = TH1F("emuSS_abs_dPhi" + "_" + sample, "emuSS_abs_dPhi" + "_" + sample, 	50,0,4)
        histos["emuSS_dR"][sample] 	= TH1F("emuSS_dR" + "_" + sample, "emuSS_dR" + "_" + sample, 			50,0,6)
        histos["emuOS_pt"][sample] 	= TH1F("emuOS_pt" + "_" + sample, "emuOS_pt" + "_" + sample, 			50,0,500)
        histos["emuOS_abs_eta"][sample] = TH1F("emuOS_abs_eta" + "_" + sample, "emuOS_abs_eta" + "_" + sample, 		50,0,6)
        histos["emuOS_abs_dEta"][sample] = TH1F("emuOS_abs_dEta" + "_" + sample, "emuOS_abs_dEta" + "_" + sample, 	50,0,5)
        histos["emuOS_abs_dPhi"][sample] = TH1F("emuOS_abs_dPhi" + "_" + sample, "emuOS_abs_dPhi" + "_" + sample, 	50,0,4)
        histos["emuOS_dR"][sample] 	= TH1F("emuOS_dR" + "_" + sample, "emuOS_dR" + "_" + sample, 			50,0,6)
        histos["met_pt"][sample] 	= TH1F("met_pt" + "_" + sample, "met_pt" + "_" + sample, 			50,0,500)
        histos["mt_emet"][sample] 	= TH1F("mt_emet" + "_" + sample, "mt_emet" + "_" + sample, 			50,0,200)
        histos["abs_dPhi_emet"][sample] = TH1F("abs_dPhi_emet" + "_" + sample, "abs_dPhi_emet" + "_" + sample, 		50,0,4)
        histos["mht_pt"][sample] 	= TH1F("mht_pt" + "_" + sample, "mht_pt" + "_" + sample, 			50,0,500)
        histos["mht_mass"][sample] 	= TH1F("mht_mass" + "_" + sample, "mht_mass" + "_" + sample, 			50,0,5000)
        histos["mt_emht"][sample] 	= TH1F("mt_emht" + "_" + sample, "mt_emht" + "_" + sample, 			50,0,500)
        histos["abs_dPhi_emht"][sample] = TH1F("abs_dPhi_emht" + "_" + sample, "abs_dPhi_emht" + "_" + sample, 		50,0,4)
        histos["mlt_pt"][sample] 	= TH1F("mlt_pt" + "_" + sample, "mlt_pt" + "_" + sample, 			50,0,500)
        histos["mt_emlt"][sample] 	= TH1F("mt_emlt" + "_" + sample, "mt_emlt" + "_" + sample, 			50,-0.0001,0.0001)
        histos["abs_dPhi_emlt"][sample] = TH1F("abs_dPhi_emlt" + "_" + sample, "abs_dPhi_emlt" + "_" + sample, 		50,0,4)
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


def FillHistTerm(histos, term, signals, bkgs, value, Sample_ID, event_wgt):
    #signal
    if Sample_ID == 25:
	histos[term]["ggH"].Fill(value, event_wgt /2)
    elif Sample_ID == 010225:
	histos[term]["VBF"].Fill(value, event_wgt /2)
    elif Sample_ID == 2425:
	histos[term]["WH"].Fill(value, event_wgt /2)
    elif Sample_ID == 2325:
	histos[term]["ZH"].Fill(value, event_wgt /2)
    elif Sample_ID == 060625:
	histos[term]["ttH"].Fill(value, event_wgt /2)

    #background
    elif Sample_ID == -23:
	histos[term]["DY"].Fill(value, event_wgt / 4)  # used two DY samples
    elif Sample_ID == -0606:
	histos[term]["ttbar"].Fill(value, event_wgt / 3)  # used two ttbar samples
    elif Sample_ID == -2423:
	histos[term]["WZ"].Fill(value, event_wgt)
    elif Sample_ID == -2323:  #ZZ_4l and ZZ_2l
	histos[term]["ZZ"].Fill(value, event_wgt)
    elif Sample_ID == -2424 or Sample_ID/10000 < -23:   # WW and triboson
	histos[term]["triboson"].Fill(value, event_wgt)
    elif Sample_ID == -062300:
	histos[term]["tZq"].Fill(value, event_wgt)
    else:
	histos[term]["tX"].Fill(value, event_wgt)


def main():
    out_name = "test.root"

    file_dir = PLOT_DIR+"/"+LABEL+"/files/"
    out_file = TFile( file_dir + "plots/" + out_name , "RECREATE")
    file_chain = TChain("tree","chain");
    file_chain.Add( file_dir + "signal.root")
    file_chain.Add( file_dir + "bkg.root")

    terms = ["mu1_pt", "mu2_pt", "mu1_abs_eta", "mu2_abs_eta",   	#mu_vars
	    "dimu_mass", "dimu_pt", "dimu_abs_eta", "dimu_abs_dEta", 	#mu_vars
	    "dimu_abs_dPhi", "dimu_dR", "cts_mu1", "cts_mu_pos",  	#mu_vars
	    "ele_pt", "ele_abs_eta", "cts_edimu", "edimu_mass", "edimu_pt",  #ele_vars
	    "edimu_abs_eta", "edimu_abs_dEta", "edimu_abs_dPhi", "edimu_dR", #ele_vars
	    "cts_emuSS", "cts_emuOS", "emuSS_pt", "emuSS_abs_eta", 		#emu vars
	    "emuSS_abs_dEta", "emuSS_abs_dPhi", "emuSS_dR", "emuOS_pt", 	#emu vars
	    "emuOS_abs_eta", "emuOS_abs_dEta", "emuOS_abs_dPhi", "emuOS_dR", 	#emu vars
	    "met_pt", "mt_emet", "abs_dPhi_emet", 		#met vars
	    "mht_pt", "mht_mass", "mt_emht", "abs_dPhi_emht",  	#met vars
	    "mlt_pt", "mt_emlt", "abs_dPhi_emlt",		#met vars
	    "dijet_mass", "dijet_pt", "dijet_abs_eta", "dijet_abs_dEta", "dijet_abs_dPhi", "dijet_dR", 	#jet vars
	    "jet1_pt", "jet1_abs_eta", "jet2_pt", "jet2_abs_eta", "jet0_pt", "jet0_abs_eta",		#jet vars
	    "nJets", "nCentJets", "nFwdJets", "nBJets_Med", "nBJets_Loose", "nBJets_Tight", "nMuons", "nEles", 	#evt vars
	    ] #end

    signals = ["ttH", "ZH", "WH", "VBF", "ggH"]
    bkgs = ["triboson", "tZq", "tX", "ttbar", "ZZ", "WZ", "DY"]

    color = {}
    color["ggH"] = kRed
    color["VBF"] =  kBlue + 1
    color["ZH"] =  kOrange + 7
    color["WH"] = kGreen + 2
#    color["WH_neg"] = kViolet + 1
    color["ttH"]  = kPink + 6

    color["DY"] =  kAzure + 7
    color["WZ"] =       kGreen - 9
    color["ZZ"] =  kCyan - 7
    color["tX"] = kPink + 6 #
    color["ttbar"] = kYellow - 9
    color["tZq"] = kViolet -9
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

	FillHistTerm(histos, "mu1_pt"	      	, signals, bkgs, file_chain.mu1_pt		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu2_pt"		, signals, bkgs, file_chain.mu2_pt		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu1_abs_eta"   	, signals, bkgs, abs(file_chain.mu1_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mu2_abs_eta"   	, signals, bkgs, abs(file_chain.mu2_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dimu_mass"   	, signals, bkgs, file_chain.dimu_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "dimu_pt"   	, signals, bkgs, file_chain.dimu_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_abs_eta"   	, signals, bkgs, abs(file_chain.dimu_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_abs_dEta"   	, signals, bkgs, abs(file_chain.dimu_dEta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_abs_dPhi"   	, signals, bkgs, abs(file_chain.dimu_dPhi)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "dimu_dR"   	, signals, bkgs, file_chain.dimu_dR   		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "cts_mu1"   	, signals, bkgs, file_chain.cts_mu1   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "cts_mu_pos"   	, signals, bkgs, file_chain.cts_mu_pos   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "ele_pt"   	, signals, bkgs, file_chain.ele_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "ele_abs_eta"   	, signals, bkgs, abs(file_chain.ele_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "cts_edimu"   	, signals, bkgs, file_chain.cts_edimu   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "edimu_mass"   	, signals, bkgs, file_chain.edimu_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "edimu_pt"   	, signals, bkgs, file_chain.edimu_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "edimu_abs_eta"   	, signals, bkgs, abs(file_chain.edimu_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "edimu_abs_dEta"   , signals, bkgs, abs(file_chain.edimu_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "edimu_abs_dPhi"   , signals, bkgs, abs(file_chain.edimu_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "edimu_dR"   	, signals, bkgs, file_chain.edimu_dR  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "cts_emuSS"   	, signals, bkgs, file_chain.cts_emuSS   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "cts_emuOS"   	, signals, bkgs, file_chain.cts_emuOS   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuSS_pt"   	, signals, bkgs, file_chain.emuSS_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuSS_abs_eta"   	, signals, bkgs, abs(file_chain.emuSS_eta)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "emuSS_abs_dEta"   , signals, bkgs, abs(file_chain.emuSS_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuSS_abs_dPhi"   , signals, bkgs, abs(file_chain.emuSS_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuSS_dR"   	, signals, bkgs, file_chain.emuSS_dR   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuOS_pt"   	, signals, bkgs, file_chain.emuOS_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuOS_abs_eta"   	, signals, bkgs, abs(file_chain.emuOS_eta)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "emuOS_abs_dEta"   , signals, bkgs, abs(file_chain.emuOS_dEta)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuOS_abs_dPhi"   , signals, bkgs, abs(file_chain.emuOS_dPhi)    	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "emuOS_dR"   	, signals, bkgs, file_chain.emuOS_dR   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "met_pt"   	, signals, bkgs, file_chain.met_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "mt_emet"   	, signals, bkgs, file_chain.mt_emet   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "abs_dPhi_emet"   	, signals, bkgs, abs(file_chain.dPhi_emet)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "mht_pt"   	, signals, bkgs, file_chain.mht_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "mht_mass"   	, signals, bkgs, file_chain.mht_mass   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "mt_emht"   	, signals, bkgs, file_chain.mt_emht   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "abs_dPhi_emht"  	, signals, bkgs, abs(file_chain.dPhi_emht)  	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "mlt_pt"   	, signals, bkgs, file_chain.mlt_pt   	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "mt_emlt"  	, signals, bkgs, file_chain.mt_emlt  	  	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
        FillHistTerm(histos, "abs_dPhi_emlt"   	, signals, bkgs, abs(file_chain.dPhi_emlt)   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
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
    for term in terms:
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
	legend[term].AddEntry(histos[term]["signal"], "signal sum")
	for sample in bkgs:
            legend[term].AddEntry(histos[term][sample], sample )
	legend[term].AddEntry(scaled_signal[term], "signal X200")
	LinearStack( term, stack_all[term], scaled_signal[term], legend[term])


    out_file.Close()
main()











	
