####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os
import sys

sys.path.insert(0, '%s/lib' % os.getcwd() )
from ROOT import *
from MNT_Helper import LoadColors, LinearStack, RatioPlot, FillHistTerm
#R.gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

#LABEL = 'miniNtuple_WH_2016_v5'  ## Sub-folder within PLOT_DIR containing histograms
LABEL = 'WH_ele_loose_ID_loose_iso_loose_mu_iso_v1'
#LABEL = 'WH_mu_med_ID_loose_iso_v1'

def InitHists(histos, terms, signals, bkgs):
    for sample in signals + bkgs:
	histos["dimu_mass"][sample] 	= TH1F("dimu_mass" + "_Net_" + sample, "dimu_mass" + "_Net_" + sample, 			12,100,160)
	histos["mu2_pt"][sample]        = TH1F("mu2_pt" + "_" + sample, "mu2_pt" + "_" + sample,                        50,0,400 )
	histos["lep_pt"][sample]        = TH1F("lep_pt" + "_" + sample, "lep_pt" + "_" + sample,                        50,0,500)
    histos["dimu_mass"]["data"]     = TH1F("dimu_mass" + "_Net_" + "Data", "dimu_mass" + "_Net_" + "Data",                  12,100,160)
    histos["mu2_pt"]["data"]        = TH1F("mu2_pt" + "_" + "Data", "mu2_pt" + "_" + "Data",                        50,0,400 )
    histos["lep_pt"]["data"]        = TH1F("lep_pt" + "_" + "Data", "lep_pt" + "_" + "Data",                        50,0,500)

def main():
    out_name = "mass_hists_cut_4n44.root"

    file_dir = PLOT_DIR+"/"+LABEL+"/"
    out_file = TFile( file_dir + "plots/" + out_name , "RECREATE")
    file_chain = TChain("tree","chain");
    file_chain.Add(file_dir + "all_samples.root")

    terms = ["dimu_mass", "mu2_pt", "lep_pt"]

    signals = ["ttH", "ZH", "WH", "VBF", "ggH"]
    bkgs = ["others", "triboson", "tZq", "tW", "ttZ", "ttbar", "WW", "ZZ", "WZ", "DY"]
    data = ["data"]

    color = {}
    LoadColors(color, "WH")

    histos = {}
    stack_all = {}
    stack_bkg = {}
    stack_sig = {}
    stack_data = {}
    ratios = {}
    legend = {}

    for term in terms:
        histos[term] = {}
        stack_all[term] = THStack("h_stack_"+term, term)
	stack_bkg[term] = THStack("bkg_stack_"+term, "bkg_"+term)
        stack_sig[term] = THStack("sig_stack_"+term, "sig_"+term)
        stack_data[term] = THStack("data_stack"+term, "data_"+term)
	

    InitHists(histos, terms, signals, bkgs)

    print file_chain.GetEntries()
    for iEvt in range( file_chain.GetEntries() ):
	file_chain.GetEvent(iEvt)
	if file_chain.mu1_lepMVA < 0.4 or file_chain.mu2_lepMVA < -0.4 or file_chain.lep_lepMVA < 0.4:
	    continue
	FillHistTerm(histos, "dimu_mass"   	, signals, bkgs, file_chain.dimu_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
 	FillHistTerm(histos, "mu2_pt"           , signals, bkgs, file_chain.mu2_pt              , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	FillHistTerm(histos, "lep_pt"           , signals, bkgs, file_chain.lep_pt              , file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)

    out_file.cd()
    scaled_signal = {}
    for term in terms:
	histos[term]["data"].SetMarkerStyle(20)
        for sample in signals:
            histos[term][sample].SetLineColor(color[sample])
            histos[term][sample].SetLineWidth(2)
            stack_sig[term].Add(histos[term][sample])
	histos[term]["signal"] = stack_sig[term].GetStack().Last()
	histos[term]["signal"].SetNameTitle(term + "_Net_Sig", term + "_Net_Sig")

        stack_all[term].Add(histos[term]["signal"])
        histos[term]["signal"].SetFillColor( kGray )
        histos[term]["signal"].SetLineWidth(0)
        for sample in bkgs:
            histos[term][sample].SetFillColor(color[sample])
	    stack_bkg[term].Add(histos[term][sample])
            stack_all[term].Add(histos[term][sample])
	histos[term]["bkg"] = stack_bkg[term].GetStack().Last()
	histos[term]["bkg"].SetNameTitle(term + "_Net_Bkg", term + "_Net_Bkg")

	scaled_signal[term] =  histos[term]["signal"].Clone()
  	scaled_signal[term].Scale(100)
	scaled_signal[term].SetLineColor(kRed)
	scaled_signal[term].SetLineWidth(2)
	scaled_signal[term].SetFillStyle(0)


	legend[term] = TLegend(0.7,0.7,1,1)
	legend[term].AddEntry(histos[term]["data"], "data")
	legend[term].AddEntry(histos[term]["signal"], "signal sum")
	for sample in bkgs:
            legend[term].AddEntry(histos[term][sample], sample )
	legend[term].AddEntry(scaled_signal[term], "signal X100")
	LinearStack( term, stack_all[term], scaled_signal[term], histos[term]["data"], legend[term], PLOT_DIR+"/"+LABEL+"/plots")

	for sample in ["data", "signal", "bkg"] + signals + bkgs:
	    histos[term][sample].Write()


    out_file.Close()
main()











	
