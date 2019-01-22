###############################################
###   ID efficiency plots for each sample   ###
###############################################

## Basic python includes for manipulating files
import os
import math

from ROOT import *

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

LABEL = 'lepMVA_efficiency_test_v6'  ## Sub-folder within PLOT_DIR containing histograms

def overlay_pt_segments(efficiencies, sample, lep_type, pt_segs ):  # for one sample, one lepton
    canv = TCanvas("overlay_" + sample +"_"+ lep_type, "overlay_" + sample +"_"+ lep_type, 600,600)
    canv.Clear()

    canv.cd()
    legend = TLegend(0,0.8,1,1)
    legend.SetNColumns(3)
    efficiencies[sample][lep_type]["pt_10_20"].Draw()
    for pt_seg in pt_segs:
    	efficiencies[sample][lep_type][pt_seg].Draw("same")
	legend.AddEntry(efficiencies[sample][lep_type][pt_seg],pt_seg, "LP")
    legend.Draw()

    canv.Update()
    canv.Write()
    return canv   

def data_over_MC(efficiencies, SF_hists, sample, data, lep_type, pt_segs, colors): # return one ratio histo, for a given pt range
    canv_up = TCanvas()
    canv_up = overlay_pt_segments(efficiencies, sample, lep_type, pt_segs)

    canv_down = TCanvas("SF_data/" + sample + "_" + lep_type, "SF_data/" + sample + "_" + lep_type, 600,200)
    canv_down.Clear()
    for pt_seg in pt_segs:
	SF_hists[pt_seg] = TGraphAsymmErrors()
	SF_hists[pt_seg] = efficiencies[data][lep_type][pt_seg].Clone()
	SF_hists[pt_seg].SetName("SF_" + sample + "_" + pt_seg)
	x = efficiencies[data][lep_type][pt_seg].GetX()
	y_sam = efficiencies[sample][lep_type][pt_seg].GetY()
	y_data = efficiencies[data][lep_type][pt_seg].GetY()
	for ip in range( len(x) ):
	    SF_hists[pt_seg].SetPoint(ip, x[ip], y_data[ip] / y_sam[ip] )
	    err_y = math.sqrt( efficiencies[sample][lep_type][pt_seg].GetErrorY(ip) ** 2 + efficiencies[data][lep_type][pt_seg].GetErrorY(ip) ** 2)
	    SF_hists[pt_seg].SetPointEYhigh(ip, err_y)
	    SF_hists[pt_seg].SetPointEYlow (ip, err_y)
	SF_hists[pt_seg].SetMaximum(1.5)
	SF_hists[pt_seg].SetMinimum(0.5)
	SF_hists[pt_seg].SetLineColor(colors[pt_seg])
	SF_hists[pt_seg].SetLineWidth(2)

    canv_down.cd()
    SF_hists["pt_10_20"].Draw()
    for pt_seg in pt_segs:
	SF_hists[pt_seg].Draw("same") 
    canv_down.SetGridy()
    canv_down.Update()
    canv_down.Write()
   

    canv = TCanvas("ID_and_SF_" + sample + "_" + lep_type, "ID_and_SF_" + sample + "_" + lep_type, 600, 800)
    canv.Clear()

    upper_pad = TPad("upperpad_"+sample, "upperpad_"+sample, 0,0.25, 1,1)
    upper_pad.Draw()
    upper_pad.cd()
    canv_up.DrawClonePad()
 
    canv.cd()
    lower_pad = TPad("lowerpad_"+sample, "lowerpad_"+sample, 0,0, 1,0.25)
    lower_pad.Draw()
    lower_pad.cd()
    canv_down.DrawClonePad()

    canv.Update()
    canv.Write()



def main():
    file_dir=PLOT_DIR+"/"+LABEL
    out_file = TFile( file_dir + "/eff_all" + ".root", "RECREATE")
    in_file = TFile.Open( file_dir + "/all_raw.root", "READ")

    overlay_dir = out_file.mkdir("overlays")
    SF_dir = out_file.mkdir("SFs")

    lep_types = ["mu", "ele"]
    pt_segs = ["pt_10_20", "pt_20_30", "pt_30_40", "pt_40_50", "pt_50_100", "pt_100_200", "pt_200_more"]
    lep_IDs = ["all", "tight"]

    data = ["SingleMu_2017B", "SingleMu_2017C", "SingleMu_2017D", "SingleMu_2017E", "SingleMu_2017F"]
#    bkgs = ["ZJets_AMC", "tt_ll_MG", "tt_ll_POW", "tt_ljj_POW_1", "tt_ljj_POW_2", "tW_neg", "tW_pos", "WZ_3l", "ZZ_4l"]
    bkgs = ["ZJets_AMC", "tt_ll_MG", "tt_ll_POW", "tt_ljj_POW_1", "tt_ljj_POW_2", "tW_neg", "tW_pos"]
#    signals = ["H2Mu_gg", "H2Mu_ZH_120", "H2Mu_ZH_125", "H2Mu_ZH_130", "H2Mu_WH_pos_120", 
#	       "H2Mu_WH_pos_125", "H2Mu_WH_pos_130", "H2Mu_WH_neg_125", "H2Mu_WH_neg_130", "H2Mu_ttH_125"]
    samples = data + bkgs

    colors = {}
    colors["pt_10_20"] = kBlack
    colors["pt_20_30"] = kRed +1
    colors["pt_30_40"] = kBlue +1
    colors["pt_40_50"] = kGreen +2
    colors["pt_50_100"] = kYellow + 2
    colors["pt_100_200"] = kOrange - 2
    colors["pt_200_more"] = kViolet

    #copy all histos from in_file
    histos = {}
    for sample in samples:
	histos[sample] = {}
    	for lep_type in lep_types:
	    histos[sample][lep_type] = {}
	    for pt_seg in pt_segs:
		histos[sample][lep_type][pt_seg] = {}
	    	for lep_ID in lep_IDs:
		    hist_name = sample + "_" + lep_type + "_" + pt_seg + "_" + lep_ID
		    if not in_file.GetListOfKeys().Contains(hist_name):
			print hist_name + "  not in in_file\n"
			histos[sample][lep_type][pt_seg][lep_ID] = TH1D(hist_name, hist_name, 16,-2.4,2.4)
		    else:
			histos[sample][lep_type][pt_seg][lep_ID] = in_file.Get(hist_name).Clone()

    #get sum of data
    histos["sum_data"] = {}
    for lep_type in lep_types:
	histos["sum_data"][lep_type] = {}
	for pt_seg in pt_segs:
	    histos["sum_data"][lep_type][pt_seg] = {}
	    for lep_ID in lep_IDs:
		histos["sum_data"][lep_type][pt_seg][lep_ID] = None
		for sample in data:
		    if histos["sum_data"][lep_type][pt_seg][lep_ID] is None:
			histos["sum_data"][lep_type][pt_seg][lep_ID] = histos[sample][lep_type][pt_seg][lep_ID].Clone()
			histos["sum_data"][lep_type][pt_seg][lep_ID].SetName("sum_data" + "_" + lep_type + "_" + pt_seg + "_" + lep_ID)
		    else:
			histos["sum_data"][lep_type][pt_seg][lep_ID].Add(histos[sample][lep_type][pt_seg][lep_ID])

    #divide tight by all and make efficiency plots
    out_file.cd()
    efficiencies = {}
    for sample in samples + ["sum_data"]:
	efficiencies[sample] = {}
	for lep_type in lep_types:
	    efficiencies[sample][lep_type] = {}
	    for pt_seg in pt_segs:
		efficiencies[sample][lep_type][pt_seg] = TGraphAsymmErrors()
		print sample + "_" + lep_type + "_" + pt_seg
		efficiencies[sample][lep_type][pt_seg].Divide(histos[sample][lep_type][pt_seg]["tight"], histos[sample][lep_type][pt_seg]["all"], "pois")   
 		efficiencies[sample][lep_type][pt_seg].SetName("efficiency" + "_" + sample + "_" + lep_type + "_" + pt_seg)
		efficiencies[sample][lep_type][pt_seg].SetMaximum(1.4)
		efficiencies[sample][lep_type][pt_seg].SetMinimum(0)
		efficiencies[sample][lep_type][pt_seg].SetLineColor(colors[pt_seg])
		efficiencies[sample][lep_type][pt_seg].SetLineWidth(2)
		efficiencies[sample][lep_type][pt_seg].Write()

    #make overlay plots for each sample and lepton
    overlay_dir.cd()
    for sample in samples + ["sum_data"]:
	for lep_type in lep_types:
	    overlay_pt_segments(efficiencies, sample, lep_type, pt_segs )

    SF_dir.cd()
    SF_histos = {}
    SF_histos["tt_ll_MG"] = {}
    SF_histos["tt_ll_POW"] = {}
    SF_histos["ZJets_AMC"] = {}
    data_over_MC(efficiencies, SF_histos["ZJets_AMC"], "ZJets_AMC", "sum_data", "mu", pt_segs, colors)
    data_over_MC(efficiencies, SF_histos["tt_ll_MG"], "tt_ll_MG", "sum_data", "ele", pt_segs, colors)
    data_over_MC(efficiencies, SF_histos["tt_ll_POW"], "tt_ll_POW", "sum_data", "ele", pt_segs, colors)


    out_file.Close()
main()

