####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os

from ROOT import *
from MNT_Helper import LinearStack, RatioPlot
#R.gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

#LABEL = 'lepMVA_3l_mu_v2_miniNtuple_baseline_v1'  ## Sub-folder within PLOT_DIR containing histograms
LABEL = 'lepMVA_ttH_3l_ele_loose_v1'


def InitHists(histos, terms, lepMVA_cuts, signals, bkgs):
    for lepMVA_cut in lepMVA_cuts:
    	for sample in signals + bkgs + ["data"]:
            histos["mu1_pt"][lepMVA_cut][sample]        = TH1F("mu1_pt" + "_" + lepMVA_cut + "_" + sample, "mu1_pt" + "_" + lepMVA_cut + "_" + sample,                        50,0,800 )
	    histos["mu2_pt"][lepMVA_cut][sample]        = TH1F("mu2_pt" + "_" + lepMVA_cut + "_" + sample, "mu2_pt" + "_" + lepMVA_cut + "_" + sample,                        50,0,400 )
	    histos["dimu_mass"][lepMVA_cut][sample]     = TH1F("dimu_mass" + "_" + lepMVA_cut + "_" + sample, "dimu_mass" + "_" + lepMVA_cut + "_" + sample,                  12,100,160)
	    histos["lep_pt"][lepMVA_cut][sample]        = TH1F("lep_pt" + "_" + lepMVA_cut + "_" + sample, "lep_pt" + "_" + lepMVA_cut + "_" + sample,                        50,0,500)



def PassLepMVA(lepMVA_cut, mu1_value, mu2_value, lep_value):
    Pass = False
    cut_value = float(lepMVA_cut)
    if mu1_value > cut_value and mu2_value > cut_value and lep_value > cut_value:
    #if lep_value > cut_value:
	Pass = True
    return Pass

def FillHistTerm(histos, term, lepMVA_cut, signals, bkgs, value, Sample_ID, event_wgt):
    #blind
    if Sample_ID == 0 and term == "dimu_mass" and value > 120 and value < 130:
	return
    #data
    if Sample_ID == 0:
	histos[term][lepMVA_cut]["data"].Fill(value, 1)

    #signal
    elif Sample_ID == 25:
	histos[term][lepMVA_cut]["ggH"].Fill(value, event_wgt )
    elif Sample_ID == 010225:
	histos[term][lepMVA_cut]["VBF"].Fill(value, event_wgt )
    elif Sample_ID == 2425:
	histos[term][lepMVA_cut]["WH"].Fill(value, event_wgt )
    elif Sample_ID == 2325:
	histos[term][lepMVA_cut]["ZH"].Fill(value, event_wgt )
    elif Sample_ID == 060625:
	histos[term][lepMVA_cut]["ttH"].Fill(value, event_wgt )

    #background
    elif Sample_ID == -23:
	histos[term][lepMVA_cut]["DY"].Fill(value, event_wgt / 3)  # used two DY samples
    elif Sample_ID == -0606:
	histos[term][lepMVA_cut]["ttbar"].Fill(value, event_wgt / 2)  # used two ttbar samples
    elif Sample_ID == -2423:
	histos[term][lepMVA_cut]["WZ"].Fill(value, event_wgt)
    elif Sample_ID == -2323:  #ZZ_4l and ZZ_2l
	histos[term][lepMVA_cut]["ZZ"].Fill(value, event_wgt)
    elif Sample_ID == -2424 or Sample_ID/10000 < -23:   # WW and triboson
	histos[term][lepMVA_cut]["triboson"].Fill(value, event_wgt)
    elif Sample_ID == -062300:
	histos[term][lepMVA_cut]["tZq"].Fill(value, event_wgt)
    else:
	histos[term][lepMVA_cut]["tX"].Fill(value, event_wgt)


def main():
    out_name = "cut_efficiency_with_data.root"

    file_dir = PLOT_DIR+"/"+LABEL+"/"
    out_file = TFile( file_dir + "plots/" + out_name , "RECREATE")
    file_chain = TChain("tree","chain");
    file_chain.Add(file_dir + "all_samples.root")
#    file_chain.Add( file_dir + "signal.root")
#    file_chain.Add( file_dir + "bkg.root")

    terms = ["mu1_pt", "mu2_pt", #"mu1_abs_eta", "mu2_abs_eta",
	     "dimu_mass", "lep_pt", #"lep_abs_eta", "nMuons", "nEles", 
	    ] #end

    signals = ["ttH", "ZH", "WH", "VBF", "ggH"]
    bkgs = ["triboson", "tZq", "tX", "ttbar", "ZZ", "WZ", "DY"]
    data = ["data"]

    lepMVA_cuts = ["-1.0",
		"-0.98","-0.96","-0.94","-0.92","-0.9","-0.88","-0.86","-0.84","-0.82","-0.8","-0.78","-0.76","-0.74","-0.72","-0.7",
		   "-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0.0","0.1","0.2","0.3","0.4","0.5","0.6",
		   "0.7","0.72","0.74","0.76","0.78", 
		"0.8","0.82","0.84","0.86","0.88","0.9","0.92","0.94","0.96","0.98"]

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
	stack_all[term] = {}
	stack_sig[term] = {}
	stack_data[term] = {}
	for lepMVA_cut in lepMVA_cuts:
	    histos[term][lepMVA_cut] = {}
            stack_all[term][lepMVA_cut] = THStack("h_stack_"+term+">"+lepMVA_cut, term+">"+lepMVA_cut)
            stack_sig[term][lepMVA_cut] = THStack("sig_stack_"+term+">"+lepMVA_cut, "sig_"+term+">"+lepMVA_cut)
            stack_data[term][lepMVA_cut] = THStack("data_stack"+term+">"+lepMVA_cut, "data_"+term+">"+lepMVA_cut)
	

    InitHists(histos, terms, lepMVA_cuts, signals, bkgs)

    print file_chain.GetEntries()
    for iEvt in range( file_chain.GetEntries() ):
#	if iEvt % 100 != 1:
#	    continue
	file_chain.GetEvent(iEvt)

	for lepMVA_cut in lepMVA_cuts:
	    if not PassLepMVA(lepMVA_cut, file_chain.mu1_lepMVA, file_chain.mu2_lepMVA, file_chain.lep_lepMVA):
		continue

	    FillHistTerm(histos, "mu1_pt"	, lepMVA_cut, signals, bkgs, file_chain.mu1_pt		, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)
	    FillHistTerm(histos, "mu2_pt"	, lepMVA_cut, signals, bkgs, file_chain.mu2_pt		, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
	    FillHistTerm(histos, "dimu_mass"   	, lepMVA_cut, signals, bkgs, file_chain.dimu_mass   	, file_chain.Sample_ID  , file_chain.xsec_norm * file_chain.event_wgt)
            FillHistTerm(histos, "lep_pt"  , lepMVA_cut, signals, bkgs, file_chain.lep_pt   	, file_chain.Sample_ID	, file_chain.xsec_norm * file_chain.event_wgt)


    out_file.cd()
    eff_graphs = {}
    for sample in signals + bkgs + data:
        eff_graphs[sample] = TGraph()
	eff_graphs[sample].SetNameTitle("efficiency_" + sample, "efficiency_" + sample)
    eff_graphs["sig"] = TGraph()
    eff_graphs["sig"].SetNameTitle("efficiency_sum_signal", "efficiency_sum_signal")
    eff_graphs["bkg"] = TGraph()
    eff_graphs["bkg"].SetNameTitle("efficiency_sum_bkg", "efficiency_sum_bkg") 
    sens_graph = TGraph()
    sens_graph.SetNameTitle("sensitivity_estimate","sensitivity_estimate")

    for sample in signals + bkgs + data:
	point_num = 0
        for lepMVA_cut in lepMVA_cuts:
	    if histos["dimu_mass"]["-1.0"][sample].Integral() == 0:
		eff_graphs[sample].SetPoint(point_num, float(lepMVA_cut), 0)
		point_num += 1
	    else:
	        eff_value = histos["dimu_mass"][lepMVA_cut][sample].Integral() / histos["dimu_mass"]["-1.0"][sample].Integral()
	        eff_graphs[sample].SetPoint(point_num, float(lepMVA_cut), eff_value)
		point_num += 1
	eff_graphs[sample].SetLineColor(color[sample])
	eff_graphs[sample].SetLineWidth(2)
	eff_graphs[sample].Write()

    point_num = 0
    range_start = histos["dimu_mass"][lepMVA_cut]["WH"].FindBin(120.1)
    range_end   = histos["dimu_mass"][lepMVA_cut]["WH"].FindBin(129.9)
    bin_width   = histos["dimu_mass"][lepMVA_cut]["WH"].GetBinWidth(1)
    print "bin_width = %d" %bin_width
    for lepMVA_cut in lepMVA_cuts:
	sig_pass = 0
	sig_tot = 0
	for sample in signals:
	    sig_pass += histos["dimu_mass"][lepMVA_cut][sample].Integral(range_start, range_end) 
	    sig_tot += histos["dimu_mass"]["-1.0"][sample].Integral(range_start, range_end)
	print "lepMVA_cut = " + lepMVA_cut
	print "sig_pass = %f" %sig_pass + "        sig_tot = %f" %sig_tot
	if sig_tot == 0:
	    eff_graphs["sig"].SetPoint(point_num, float(lepMVA_cut), 0)
	else:
	    eff_graphs["sig"].SetPoint(point_num, float(lepMVA_cut), sig_pass/sig_tot)

	bkg_pass = 0
	bkg_tot = 0
	for sample in bkgs:
	    #bkg_pass += histos["dimu_mass"][lepMVA_cut][sample].Integral(range_start, range_end)
            bkg_tot += histos["dimu_mass"]["-1.0"][sample].Integral(range_start, range_end)
	    bkg_pass += histos["dimu_mass"][lepMVA_cut][sample].Integral() * histos["dimu_mass"]["-1.0"][sample].Integral(range_start, range_end) / histos["dimu_mass"]["-1.0"][sample].Integral()
	print "lepMVA_cut = " + lepMVA_cut
        print "bkg_pass = %f" %bkg_pass + "        bkg_tot = %f" %bkg_tot
	if bkg_tot == 0:
            eff_graphs["bkg"].SetPoint(point_num, float(lepMVA_cut), 0)
        else:
            eff_graphs["bkg"].SetPoint(point_num, float(lepMVA_cut), bkg_pass/bkg_tot)
	if bkg_pass == 0:
	    sens_graph.SetPoint(point_num, float(lepMVA_cut), 1)
	else:
	    sens_graph.SetPoint(point_num, float(lepMVA_cut), sig_pass / bkg_pass ** 0.5)
	point_num += 1
    eff_graphs["sig"].SetLineColor(kRed)
    eff_graphs["sig"].SetLineWidth(2)
    eff_graphs["bkg"].SetLineColor(kBlack)
    eff_graphs["bkg"].SetLineWidth(2)
    eff_graphs["sig"].Write()
    eff_graphs["bkg"].Write()

    sens_graph.SetLineWidth(2)
    sens_graph.Write()


    dimu_dir = out_file.mkdir("dimu_mass")
    dimu_dir.cd()
    scaled_signal = {}
    ratios = {}
    for term in ["dimu_mass"]:
	scaled_signal[term] = {}
	ratios[term] = {}
	for lepMVA_cut in lepMVA_cuts:
	    histos[term][lepMVA_cut]["data"].SetMarkerStyle(20)
            for sample in signals:
            	histos[term][lepMVA_cut][sample].SetLineColor(color[sample])
            	histos[term][lepMVA_cut][sample].SetLineWidth(2)
            	stack_sig[term][lepMVA_cut].Add(histos[term][lepMVA_cut][sample])

	    histos[term][lepMVA_cut]["signal"] = stack_sig[term][lepMVA_cut].GetStack().Last()
            stack_all[term][lepMVA_cut].Add(histos[term][lepMVA_cut]["signal"])
            histos[term][lepMVA_cut]["signal"].SetFillColor( kGray )
            histos[term][lepMVA_cut]["signal"].SetLineWidth(0)
            for sample in bkgs:
            	histos[term][lepMVA_cut][sample].SetFillColor(color[sample])
            	stack_all[term][lepMVA_cut].Add(histos[term][lepMVA_cut][sample])

	    scaled_signal[term][lepMVA_cut] =  histos[term][lepMVA_cut]["signal"].Clone()
  	    scaled_signal[term][lepMVA_cut].Scale(20)
	    scaled_signal[term][lepMVA_cut].SetLineColor(kRed)
	    scaled_signal[term][lepMVA_cut].SetLineWidth(2)
	    scaled_signal[term][lepMVA_cut].SetFillStyle(0)
	

	    legend[term] = TLegend(0.7,0.7,1,1)
	    legend[term].AddEntry(histos[term][lepMVA_cut]["data"], "data")
	    legend[term].AddEntry(histos[term][lepMVA_cut]["signal"], "signal sum")
	    for sample in bkgs:
                legend[term].AddEntry(histos[term][lepMVA_cut][sample], sample )
	    legend[term].AddEntry(scaled_signal[term][lepMVA_cut], "signal X20")

	    ratios[term][lepMVA_cut] = TGraphAsymmErrors()
            ratios[term][lepMVA_cut].Divide(histos[term][lepMVA_cut]["data"], stack_all[term][lepMVA_cut].GetStack().Last(), "pois")
            ratios[term][lepMVA_cut].SetName("ratiograph_"+term + lepMVA_cut)
	    #LinearStack( term + "_lepMVA>" + lepMVA_cut, stack_all[term][lepMVA_cut], scaled_signal[term][lepMVA_cut], histos[term][lepMVA_cut]["data"], legend[term], PLOT_DIR+"/"+LABEL+"/plots")
	    RatioPlot( term + "_lepMVA>" + lepMVA_cut, stack_all[term][lepMVA_cut], scaled_signal[term][lepMVA_cut], histos[term][lepMVA_cut]["data"], ratios[term][lepMVA_cut], legend[term], PLOT_DIR+"/"+LABEL+"/plots")


    out_file.Close()
main()



