###############################################
###   ID efficiency plots for each sample   ###
###############################################

## Basic python includes for manipulating files
import os
import math
import array

from ROOT import *

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

LABEL = 'lepMVA_SF_v1'  ## Sub-folder within PLOT_DIR containing histograms

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
    out_file = TFile( file_dir + "/SF_ele_variable_eta_bin" + ".root", "RECREATE")
    in_file = TFile.Open( file_dir + "/all_raw.root", "READ")

    #overlay_dir = out_file.mkdir("overlays")
    eff_3D_dir = out_file.mkdir("eff_3Ds")
    eff_2D_dir = out_file.mkdir("eff_2Ds")
    SF_dir = out_file.mkdir("SFs")
    SF_scan_dir = out_file.mkdir("SF_scans")

    # binning setting from input file
    pt_bins     = 50
    pt_low      = 0
    pt_high     = 500
    eta_bins    = 16
    eta_low     = -2.4
    eta_high    = 2.4
    MVA_bins    = 200
    MVA_low     = -1.0
    MVA_high    = 1.0

    lep_type = "ele"
    pt_segs = ["pt_10_20", "pt_20_30", "pt_30_40", "pt_40_50", "pt_50_100", "pt_100_200", "pt_200_more"]
    # binnings that are used for output file
    pt_binning = array.array('d',[10,20,30,40,50,100,200,500])
#    eta_binning = array.array('d', [ round( eta*(eta_high-eta_low)/eta_bins + eta_low, 3) for eta in range(eta_bins+1)])
    eta_binning = array.array('d', [0, 0.9, 1.2, 1.9, 2.4])
    MVA_binning = array.array('d', [ round( MVA*(MVA_high-MVA_low)/MVA_bins + MVA_low, 3) for MVA in range(MVA_bins+1)])
    print "pt_binning is" 
    print  pt_binning
    print "eta_binning is" 
    print  eta_binning
    print "MVA_binning is" 
    print  MVA_binning

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

    # get original hists from input
    hist3D = {}
    for sample in samples:
	hist_name = sample + "_" + lep_type
	if not in_file.GetListOfKeys().Contains(hist_name):
	    hist3D[sample] = TH3D(hist_name, hist_name, pt_bins,pt_low,pt_high, eta_bins,eta_low,eta_high, MVA_bins,MVA_low,MVA_high)
	else:
	    hist3D[sample] = in_file.Get(hist_name).Clone()
    hist3D["sum_data"] = None
    for sample in data:
	if hist3D["sum_data"] is None:	
	    hist3D["sum_data"] = hist3D[sample].Clone()
	    hist3D["sum_data"].SetName("sum_data" + "_" + lep_type)
	else:
	    hist3D["sum_data"].Add( hist3D[sample] )

    # calculate efficiency vs 3D
    eff_3D = {}
    for sample in ["sum_data"] + bkgs:
	hist_name = "eff_" + sample + "_" + lep_type
	eff_3D[sample] = TH3D(hist_name, hist_name, len(pt_binning)-1,pt_binning, len(eta_binning)-1,eta_binning, len(MVA_binning)-1,MVA_binning)
	for pt_bin in range(1,len(pt_binning)):  # [1,2, ...., len(bin)-1]
	    for eta_bin in range(1,len(eta_binning)):
		for MVA_bin in range(1,len(MVA_binning)):
		    pt_width = round( (1.0*pt_high-pt_low)/pt_bins, 2 )
		    pt_integral_low  = int( (pt_binning[pt_bin-1] - pt_low)/pt_width ) + 1 # first element is array[0], first bin in histogram is bin=1
		    pt_integral_high = int( (pt_binning[pt_bin]   - pt_low)/pt_width ) 
		    # bin coverage includes lower bound, excludes upper bound
		    eta_width = round( (1.0*eta_high-eta_low)/eta_bins, 2 )
		    eta_neg_low  = int( (-eta_binning[eta_bin] - eta_low)/eta_width ) + 1
		    eta_neg_high = int( (-eta_binning[eta_bin-1] - eta_low)/eta_width )
		    eta_pos_low  = int( (eta_binning[eta_bin-1] - eta_low)/eta_width ) + 1
		    eta_pos_high = int( (eta_binning[eta_bin] - eta_low)/eta_width )
		    # integral includes both bin_low and bin_high
		    if ( hist3D[sample].Integral(pt_integral_low,pt_integral_high, eta_neg_low,eta_neg_high, 1,len(MVA_binning)-1) + hist3D[sample].Integral(pt_integral_low,pt_integral_high, eta_pos_low,eta_pos_high, 1,len(MVA_binning)-1) ) == 0:
			eff_value = 0
		    else:
		    	eff_value = ( hist3D[sample].Integral(pt_integral_low,pt_integral_high, eta_neg_low,eta_neg_high, MVA_bin,len(MVA_binning)-1) + hist3D[sample].Integral(pt_integral_low,pt_integral_high, eta_pos_low,eta_pos_high, MVA_bin,len(MVA_binning)-1) ) / ( hist3D[sample].Integral(pt_integral_low,pt_integral_high, eta_neg_low,eta_neg_high, 1,len(MVA_binning)-1) + hist3D[sample].Integral(pt_integral_low,pt_integral_high, eta_pos_low,eta_pos_high, 1,len(MVA_binning)-1) )
		    eff_3D[sample].SetBinContent(pt_bin, eta_bin, MVA_bin, eff_value)	

    # calculate SF vs 3D
    SF_3D = TH3D("SF_3D" + lep_type, "SF_3D" + lep_type, len(pt_binning)-1,pt_binning, len(eta_binning)-1,eta_binning, len(MVA_binning)-1,MVA_binning)
    for pt_bin in range(1,len(pt_binning)):
	for eta_bin in range(1,len(eta_binning)):
	    for MVA_bin in range(1,len(MVA_binning)):
		if eff_3D["tt_ll_POW"].GetBinContent(pt_bin, eta_bin, MVA_bin) == 0:
		    SF_value = 0
	 	else:
		    SF_value = eff_3D["sum_data"].GetBinContent(pt_bin, eta_bin, MVA_bin) / eff_3D["tt_ll_POW"].GetBinContent(pt_bin, eta_bin, MVA_bin)
		SF_3D.SetBinContent(pt_bin, eta_bin, MVA_bin, SF_value)

    # get value and fill efficiency vs 2D
    pt_eta_eff_2D = {}
    for sample in ["sum_data"] + bkgs:
	pt_eta_eff_2D[sample] = {}
	for MVA_bin in range(1,len(MVA_binning)):
	    wp_str = "%.2f"%((MVA_bin-1.0)/(len(MVA_binning)-1)*2.0 - 1.0)
# 	    print wp_str + "\n"
	    hist_name = "eff_2D_" + sample + "_" + wp_str + "_" + lep_type
	    pt_eta_eff_2D[sample][wp_str] = TH2D(hist_name, hist_name, len(pt_binning)-1,pt_binning, len(eta_binning)-1,eta_binning)
	    for pt_bin in range(1,len(pt_binning)):
		for eta_bin in range(1,len(eta_binning)):
		    eff_value = eff_3D[sample].GetBinContent(pt_bin, eta_bin, MVA_bin)
		    pt_eta_eff_2D[sample][wp_str].SetBinContent(pt_bin, eta_bin, eff_value)

    # get value and fill SF vs 2D
    pt_eta_SF_2D = {}
    for MVA_bin in range(1,len(MVA_binning)):
	wp_str = "%.2f"%((MVA_bin-1.0)/(len(MVA_binning)-1)*2.0 - 1.0)
	hist_name = "SF_2D_" + wp_str + "_" + lep_type
	pt_eta_SF_2D[wp_str] = TH2D(hist_name, hist_name, len(pt_binning)-1,pt_binning, len(eta_binning)-1,eta_binning)
      	for pt_bin in range(1,len(pt_binning)):
            for eta_bin in range(1,len(eta_binning)):
                SF_value = SF_3D.GetBinContent(pt_bin, eta_bin, MVA_bin)
		pt_eta_SF_2D[wp_str].SetBinContent(pt_bin, eta_bin, SF_value)

#    SF_scan_1D = {}
#    for pt_seg in pt_segs:
#	for eta_seg in eta_segs:


    eff_3D_dir.cd()
    for sample in ["sum_data"] + bkgs:
	eff_3D[sample].Write()
    eff_2D_dir.cd()
    for sample in ["sum_data"] + bkgs:
	for MVA_bin in range(1,len(MVA_binning)):
            wp_str = "%.2f"%((MVA_bin-1.0)/(len(MVA_binning)-1)*2.0 - 1.0)
	    pt_eta_eff_2D[sample][wp_str].Write()
	
    SF_dir.cd()
    SF_3D.Write()
    for MVA_bin in range(1,len(MVA_binning)):
	wp_str = "%.2f"%((MVA_bin-1.0)/(len(MVA_binning)-1)*2.0 - 1.0)
	pt_eta_SF_2D[wp_str].Write()

#    #divide tight by all and make efficiency plots
#    out_file.cd()
#    efficiencies = {}
#    for sample in samples + ["sum_data"]:
#	efficiencies[sample] = {}
#	for lep_type in lep_types:
#	    efficiencies[sample][lep_type] = {}
#	    for pt_seg in pt_segs:
#		efficiencies[sample][lep_type][pt_seg] = TGraphAsymmErrors()
#		print sample + "_" + lep_type + "_" + pt_seg
#		efficiencies[sample][lep_type][pt_seg].Divide(histos[sample][lep_type][pt_seg]["tight"], histos[sample][lep_type][pt_seg]["all"], "pois")   
# 		efficiencies[sample][lep_type][pt_seg].SetName("efficiency" + "_" + sample + "_" + lep_type + "_" + pt_seg)
#		efficiencies[sample][lep_type][pt_seg].SetMaximum(1.4)
#		efficiencies[sample][lep_type][pt_seg].SetMinimum(0)
#		efficiencies[sample][lep_type][pt_seg].SetLineColor(colors[pt_seg])
#		efficiencies[sample][lep_type][pt_seg].SetLineWidth(2)
#		efficiencies[sample][lep_type][pt_seg].Write()
#
#    #make overlay plots for each sample and lepton
#    overlay_dir.cd()
#    for sample in samples + ["sum_data"]:
#	for lep_type in lep_types:
#	    overlay_pt_segments(efficiencies, sample, lep_type, pt_segs )
#
#    SF_dir.cd()
#    SF_histos = {}
#    SF_histos["tt_ll_MG"] = {}
#    SF_histos["tt_ll_POW"] = {}
#    SF_histos["ZJets_AMC"] = {}
#    data_over_MC(efficiencies, SF_histos["ZJets_AMC"], "ZJets_AMC", "sum_data", "mu", pt_segs, colors)
#    data_over_MC(efficiencies, SF_histos["tt_ll_MG"], "tt_ll_MG", "sum_data", "ele", pt_segs, colors)
#    data_over_MC(efficiencies, SF_histos["tt_ll_POW"], "tt_ll_POW", "sum_data", "ele", pt_segs, colors)


    out_file.Close()
main()

