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



def main():
    file_dir=PLOT_DIR+"/"+LABEL
    out_file = TFile( file_dir + "/compare_with_ref_ele" + ".root", "RECREATE")
    in_file_ref = TFile.Open( "scaleFactors_ele_2017" + ".root", "READ")
    in_file_calc = TFile.Open( file_dir + "/SF_ele_variable_eta_bin" + ".root", "READ")

    pull_dir = out_file.mkdir("SF_pull")
    scan_dir = out_file.mkdir("SF_scan")

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
    # binnings that are used for output file
    pt_binning = array.array('d',[10,20,30,40,50,100,200])
    #eta_binning = array.array('d', [ round( eta*(eta_high-eta_low)/eta_bins + eta_low, 3) for eta in range(eta_bins+1)])
    eta_binning = array.array('d', [0, 0.9, 1.2, 1.9, 2.4])
    MVA_binning = array.array('d', [ round( MVA*(MVA_high-MVA_low)/MVA_bins + MVA_low, 3) for MVA in range(MVA_bins+1)])
    print "pt_binning is" 
    print  pt_binning
    print "eta_binning is" 
    print  eta_binning
    print "MVA_binning is" 
    print  MVA_binning

    eta_ref = array.array('d', [0, 0.9, 1.2, 1.9, 2.4])
    print "eta_binning in ref is"
    print eta_ref

    colors = {}
    colors["ref"] = kRed
    colors["calc"] = kBlack

    wps = ["tight", "medium", "loose"]

    ref_name = {}
#    ref_name["tight"] = "MuonToTTVLeptonMvatZq"
#    ref_name["medium"] = "MuonToTTVLeptonMvattZ3l"
#    ref_name["loose"] = "MuonToTTVLeptonMvattZ4l"

    ref_name["tight"] = "EleToTTVLeptonMvatZq"
    ref_name["medium"] = "EleToTTVLeptonMvattZ3l"
    ref_name["loose"] = "EleToTTVLeptonMvattZ4l"

    wp_value = {}
    wp_value["tight"] = 0.8
    wp_value["medium"] = 0.4
    wp_value["loose"] = -0.4

    # get original hists from input
    SF_ref = {}
    for wp in wps:
	if not in_file_ref.GetListOfKeys().Contains( ref_name[wp] ):
	    print "cannot find " + ref_name[wp]
	else:
	    SF_ref[wp] = in_file_ref.Get(ref_name[wp]).Clone()

    SF_name = "SFs/" + "SF_3D" + lep_type
    SF_3D = in_file_calc.Get(SF_name).Clone()

    # make pull plots
    pull_plots = {}
    for wp in wps:
	hist_name = "pull_" + wp + "_" + lep_type
	pull_plots[wp] = TH2D(hist_name, hist_name, len(pt_binning)-1,pt_binning, len(eta_binning)-1,eta_binning)
	for pt_bin in range(1,len(pt_binning)):
	    for eta_bin in range(1,len(eta_binning)):
		MVA_bin = int( (wp_value[wp] - MVA_low)/0.01 ) + 1
		calc_value = SF_3D.GetBinContent(pt_bin, eta_bin, MVA_bin) 
		ref_value = SF_ref[wp].GetBinContent(pt_bin, eta_bin)
		ref_error = SF_ref[wp].GetBinError(pt_bin, eta_bin) 
		pull_plots[wp].SetBinContent( pt_bin, eta_bin, 100.0*(calc_value - ref_value)/ref_value )

    # make 1D SF_scan and fill
    SF_scan = {}
    ref_scan = {}
    for eta_bin in range(1,len(eta_ref)):
	eta_name = "eta_" + "%.1f"%eta_ref[eta_bin-1] + "_" + "%.1f"%eta_ref[eta_bin]
	SF_scan[eta_name] = {}
 	ref_scan[eta_name] = {}
#	eta_bin_calc = int( (eta_ref[eta_bin-1] + 2.4)/0.3 ) + 1   # currently taking the bin at lower edge
	eta_bin_calc = eta_bin
	for pt_bin in range(1,len(pt_binning)):
	    pt_name = "pt_" + "%d"%pt_binning[pt_bin-1] + "_" + "%d"%pt_binning[pt_bin]
	    SF_scan[eta_name][pt_name] = TGraphErrors()
	    SF_scan[eta_name][pt_name].SetNameTitle( "scan_" + eta_name+ "_" + pt_name, "scan_" + eta_name+ "_" + pt_name)
	    ref_scan[eta_name][pt_name] = TGraphErrors()
	    ref_scan[eta_name][pt_name].SetNameTitle( "ref_" + eta_name+ "_" + pt_name, "ref_" + eta_name+ "_" + pt_name)
	    for MVA_bin in range(1,len(MVA_binning)):
		SF_value = SF_3D.GetBinContent(pt_bin, eta_bin_calc, MVA_bin)
		SF_scan[eta_name][pt_name].SetPoint( MVA_bin-1, MVA_binning[MVA_bin-1], SF_value)
	    count = 0
	    for wp in wps:
	    	SF_value = SF_ref[wp].GetBinContent(pt_bin, eta_bin)
		SF_error = SF_ref[wp].GetBinError(pt_bin, eta_bin)
		ref_scan[eta_name][pt_name].SetPoint( count, wp_value[wp], SF_value)
		ref_scan[eta_name][pt_name].SetPointError( count, 0, SF_error)
		count += 1

    pull_dir.cd()
#    gStyle.SetPaintTextFormat("2.2f%%")
    for wp in wps:
	pull_plots[wp].SetMarkerSize(1.5)  # somehow this takes care of text size, default 1
	pull_plots[wp].Write()

    scan_dir.cd()
    for eta_bin in range(1,len(eta_ref)):
	eta_name = "eta_" + "%.1f"%eta_ref[eta_bin-1] + "_" + "%.1f"%eta_ref[eta_bin]
	for pt_bin in range(1,len(pt_binning)):
	    pt_name = "pt_" + "%d"%pt_binning[pt_bin-1] + "_" + "%d"%pt_binning[pt_bin]
	    SF_scan[eta_name][pt_name].SetLineColor(colors["calc"])
	    SF_scan[eta_name][pt_name].SetLineWidth(2)
	    SF_scan[eta_name][pt_name].SetMarkerStyle(10)
	    SF_scan[eta_name][pt_name].Write()
	    ref_scan[eta_name][pt_name].SetMarkerStyle(20)
	    ref_scan[eta_name][pt_name].SetMarkerColor(colors["ref"])
	    ref_scan[eta_name][pt_name].Write()

    out_file.Close()
main()

