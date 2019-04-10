#! /usr/bin/env python

##################################################
#####  functions for making plots from MNT   #####
##################################################
import array
from ROOT import *

def LoadColors(color, color_set):
    if color_set == "WH":
    	color["data"] = kBlack

    	color["ggH"] = kRed
    	color["VBF"] =  kBlue + 1
    	color["ZH"] =  kOrange + 7
    	color["WH"] = kGreen + 2
#   	 color["WH_neg"] = kViolet + 1
    	color["ttH"]  = kPink + 6

    	color["DY"] =  kAzure + 7
    	color["WZ"] =       kGreen - 9
    	color["ZZ"] =  kCyan - 7
    	color["WW"] = kBlack
    	color["ttbar"] = kYellow - 9
    	color["ttZ"] = kBlack
    	color["tW"]  = kBlack
    	color["tZq"] = kViolet -9
    	color["triboson"] = kOrange + 6 #
    	color["others"] = kPink + 6 #
    else:
	return

def LoadSampleNames(signals, bkgs, data, sample_set):
    if sample_set == "lep_3l":
	signals.append("ttH")
	signals.append("ZH")
	signals.append("WH")
	signals.append("VBF")
	signals.append("ggH")

	bkgs.append("others")
	bkgs.append("triboson")
	bkgs.append("tZq")
	bkgs.append("tW")
	bkgs.append("ttZ")
	bkgs.append("ttbar")
	bkgs.append("WW")
	bkgs.append("ZZ")
	bkgs.append("WZ")
	bkgs.append("DY")

    	data.append("data")
    else:
	return


def LinearStack( name, all_stack, scaled_signal, h_data, legend, plot_dir):
    canv = TCanvas("Stack_" + name, "Stack_" + name, 600,600)
    canv.Clear()

    canv.cd()
    all_stack.Draw("HIST")
    scaled_signal.Draw("HISTSAME")
    h_data.Draw("SAMEPE")
    legend.Draw()
    canv.Update()
    canv.Write()
#    canv.SaveAs( plot_dir + "/linear_stack_" + name + ".png")




def RatioPlot( name, all_stack, scaled_signal, h_data, ratio_graph, legend, plot_dir, log_plot = False ):   # Do not use TRatioPlot! It is a devil!    -XWZ 19.09.2018
    canv = TCanvas("ratio_"+name, "ratio_"+name, 600,600)
    canv.Clear()

    upper_pad = TPad("upperpad_"+name, "upperpad_"+name, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.Draw()
    upper_pad.cd()
    if log_plot:
        upper_pad.SetLogy()
#       all_stack.SetMinimum(1e-5)
#       all_stack.SetMaximum(1e6)
    	all_stack.SetMinimum(1e-3)
    	all_stack.SetMaximum(1e8)
    all_stack.Draw("HIST")
    scaled_signal.Draw("HISTSAME")
    h_data.SetMarkerStyle(20)
    h_data.Draw("SAMEPE")
    legend.Draw()

    canv.cd()
    lower_pad = TPad("lowerpad_"+name, "lowerpad_"+name, 0,0.05, 1,0.2)
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
#    canv.SaveAs( plot_dir + "/linear_stack_" + name + ".png")


def PassLepMVA(lepMVA_cut, lep_value):
    Pass = False
    cut_value = float(lepMVA_cut)
    if lep_value > cut_value:
        Pass = True
    return Pass


def FillHistTerm(histos, term, signals, bkgs, value, Sample_ID, event_wgt):
    #blind
    if Sample_ID == 0 and term == "dimu_mass" and value > 120 and value < 130:
        return
    #data
    if Sample_ID == 0:
        histos[term]["data"].Fill(value, 1)
    #signal
    elif Sample_ID == 25:
        histos[term]["ggH"].Fill(value, event_wgt )
    elif Sample_ID == 010225:
        histos[term]["VBF"].Fill(value, event_wgt )
    elif Sample_ID == 2425:
        histos[term]["WH"].Fill(value, event_wgt )
    elif Sample_ID == 2325:
        histos[term]["ZH"].Fill(value, event_wgt )
    elif Sample_ID == 060625:
        histos[term]["ttH"].Fill(value, event_wgt )

    #background
    elif Sample_ID == -23:
        histos[term]["DY"].Fill(value, event_wgt / 3)  # used two DY samples
    elif Sample_ID == -0606:
        histos[term]["ttbar"].Fill(value, event_wgt / 2)  # used two ttbar samples
    elif Sample_ID == -2423:
        histos[term]["WZ"].Fill(value, event_wgt)
    elif Sample_ID == -2323:  #ZZ_4l and ZZ_2l
        histos[term]["ZZ"].Fill(value, event_wgt)
    elif Sample_ID == -2424: # WW
	histos[term]["WW"].Fill(value, event_wgt)
    elif Sample_ID == -242424 or Sample_ID == -242423 or Sample_ID == -242323 or Sample_ID == -232323:   # triboson
        histos[term]["triboson"].Fill(value, event_wgt)
    elif Sample_ID == -062300:
        histos[term]["tZq"].Fill(value, event_wgt)
    elif Sample_ID == -060623:
        histos[term]["ttZ"].Fill(value, event_wgt)
    elif Sample_ID == -0624:
	histos[term]["tW"].Fill(value, event_wgt) 
    elif Sample_ID != -062500 and Sample_ID != -062524 and Sample_ID != -06062424 and Sample_ID != -999: #tHq, tHW, and ttWW 
        histos[term]["others"].Fill(value, event_wgt)


def FindBinNum(binning, value):
    for i in range(1,len(binning)):
	if value >= binning[i-1] and value < binning[i]:
	    return i
    #print len(binning)
    #print "no bin covers this value %f" %value
    return 0

def GetSF(lep_type, pt, eta, lepMVA_cut):
    if lepMVA_cut < -1 or lepMVA_cut > 1:
	print "lepMVA cut out of range"
	return 0.0

    MVA_bins    = 200
    MVA_low     = -1.0
    MVA_high    = 1.0
    eta = abs(eta)
    pt_binning = array.array('d',[10,20,30,40,50,100,200])
    eta_binning = array.array('d', [0, 0.9, 1.2, 1.9, 2.4])
    MVA_binning = array.array('d', [ round( MVA*(MVA_high-MVA_low)/MVA_bins + MVA_low, 3) for MVA in range(MVA_bins+1)])

    SF_mu = TH3D()
    SF_ele_loose = TH2D()
    SF_ele_medium = TH2D()
    SF_ele_tight = TH2D()
    
    scale_factor = 0.0

    if lep_type == "muon":
	in_file = TFile.Open("/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/lepMVA_SF_v1/SF_muon_variable_eta_bin.root", "READ")
	SF_mu = in_file.Get("SFs/" + "SF_3D" + lep_type).Clone()
	pt_bin = FindBinNum(pt_binning, pt) 
	eta_bin = FindBinNum(eta_binning, eta)
	MVA_bin = FindBinNum(MVA_binning, lepMVA_cut)
	scale_factor = SF_mu.GetBinContent(pt_bin, eta_bin, MVA_bin)
    elif lep_type == "ele":
	in_file = TFile.Open("/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/LepMVA_efficiency/python/scaleFactors_ele_2017.root", "READ")
    	SF_ele_loose = in_file.Get("EleToTTVLeptonMvattZ4l").Clone()
	SF_ele_medium = in_file.Get("EleToTTVLeptonMvattZ3l").Clone()
	SF_ele_tight = in_file.Get("EleToTTVLeptonMvatZq").Clone()
	pt_bin = FindBinNum(pt_binning, pt)
        eta_bin = FindBinNum(eta_binning, eta)
	if lepMVA_cut == -0.4:
	    scale_factor = SF_ele_loose.GetBinContent(pt_bin, eta_bin)
	elif lepMVA_cut == -0.4:
	    scale_factor = SF_ele_medium.GetBinContent(pt_bin, eta_bin)
	elif lepMVA_cut == 0.8:
	    scale_factor = SF_ele_tight.GetBinContent(pt_bin, eta_bin)
	else:
	    scale_factor = 1.0

    return scale_factor


