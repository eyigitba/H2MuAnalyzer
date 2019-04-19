#! /usr/bin/env python

####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os
import array

from ROOT import *
gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

#LABEL = 'lepMVA_3l_mu_v2_miniNtuple_baseline_v1'  ## Sub-folder within PLOT_DIR containing histograms
#LABEL = 'lepMVA_ttH_3l_ele_loose_v1'

LABEL = 'WH_ele_loose_ID_loose_iso_loose_mu_iso_v1'
#LABEL = 'WH_mu_med_ID_loose_iso_v1'

def GrabHist( in_file, sample, mu1_cut, mu2_cut, lep_cut):
    hist_name = sample + "mu1_%.2f_"%mu1_cut + "mu2_%.2f_"%mu2_cut + "lep_%.2f"%lep_cut
    histo = None
    if not in_file.GetListOfKeys().Contains( hist_name ):
	print "do not have " + hist_name
	histo = TH1D(hist_name, hist_name, 12,100,160)
    else:
	histo = in_file.Get(hist_name).Clone()
    return histo
    
def GetYield( in_file, sample, mu1_cut, mu2_cut, lep_cut, is_sig): # assuming mass binning 12,100,160
    if is_sig:
	histo = GrabHist(in_file, sample, mu1_cut, mu2_cut, lep_cut)
	return histo.Integral(5,6) # bin 5 and bin 6 is mass[120,130]
    else:
	histo_tot = GrabHist(in_file, sample, -1.0, -1.0, -1.0)
	histo_it  = GrabHist(in_file, sample, mu1_cut, mu2_cut, lep_cut)
	bkg_tot_all    = histo_tot.Integral(3,10) # bin 3 to 10 is mass[110,150]
	if bkg_tot_all == 0:
	    return 0.0
	bkg_tot_window = histo_tot.Integral(5,6)
	bkg_it_all     = histo_it.Integral(3,10)
	return bkg_it_all * bkg_tot_window / bkg_tot_all

def main():
    out_name = "scan_results.root"

    file_dir = PLOT_DIR+"/"+LABEL+"/"
    out_file = TFile( file_dir + "plots/lepMVA_scan_round_2/" + out_name , "RECREATE")
    in_file  = TFile.Open(file_dir + "plots/lepMVA_scan_round_2/" + "all_out.root", "READ")

    signals = ["ttH", "ZH", "WH", "VBF", "ggH"]
    bkgs = ["others", "triboson", "tZq", "tW", "ttZ", "ttbar", "WW", "ZZ", "WZ", "DY"]
    data = ["data"]

#    mu1_MVAcuts = [ -1.0, -0.4, 0.4, 0.8, 1.0 ] # for standard cuts scanning
#    mu2_MVAcuts = [ -1.0, -0.4, 0.4, 0.8, 1.0 ]
#    lep_MVAcuts = [ -1.0, -0.4, 0.4, 0.8, 1.0 ]
#
#    mu1_cuts = array.array('d', mu1_MVAcuts)
#    mu2_cuts = array.array('d', mu2_MVAcuts)
#    lep_cuts = array.array('d', lep_MVAcuts)   # for standard cuts scanning

#    mu1_MVAcuts =  [ round(1.0*MVA/10 + 0.4,2) for MVA in range(7)]
#    mu2_MVAcuts =  [ round(1.0*MVA/10 + 0.4,2) for MVA in range(7)]
#    lep_MVAcuts =  [ round(1.0*MVA/20 + 0.0,2) for MVA in range(18)]

    mu1_MVAcuts =  [ round(1.0*MVA/5 - 0.6,2) for MVA in range(9)]
    mu2_MVAcuts =  [ round(1.0*MVA/5 - 0.6,2) for MVA in range(9)]
    lep_MVAcuts =  [ 0.4, 1.0 ]

    mu1_cuts = array.array('d', [-1.0] + mu1_MVAcuts)
    mu2_cuts = array.array('d', [-1.0] + mu2_MVAcuts)
    lep_cuts = array.array('d', [-1.0] + lep_MVAcuts)       # for custom cuts scanning


    sens_hist = TH3D("sensitivity_3D", "sensitivity_3D", len(mu1_cuts)-1, mu1_cuts, len(mu2_cuts)-1, mu2_cuts, len(lep_cuts)-1, lep_cuts)

    for mu1_bin in range(1,len(mu1_cuts)):
	mu1_cut = mu1_cuts[mu1_bin-1]
        for mu2_bin in range(1,len(mu2_cuts)):
	    mu2_cut = mu2_cuts[mu2_bin-1]
            for lep_bin in range(1,len(lep_cuts)):
		lep_cut = lep_cuts[lep_bin-1]
		sig_yield = 0
		bkg_yield = 0
		print "doing mu1_%.2f_"%mu1_cut + "mu2_%.2f_"%mu2_cut + "lep_%.2f"%lep_cut
		for sample in signals:
		    sig_yield += GetYield( in_file, sample, mu1_cut, mu2_cut, lep_cut, True)
		for sample in bkgs:
		    bkg_yield += GetYield( in_file, sample, mu1_cut, mu2_cut, lep_cut, False)
		if bkg_yield == 0:
		    sens_hist.SetBinContent(mu1_bin, mu2_bin, lep_bin, 0)
		else:
		    sens_hist.SetBinContent(mu1_bin, mu2_bin, lep_bin, sig_yield / bkg_yield ** 0.5)

    projections = {}
    for lep_cut in lep_MVAcuts:
        find_bin = 0
        for lep_bin in range(1,len(lep_cuts)):
	    print lep_bin-1
	    print lep_cuts[lep_bin-1]
            if lep_cut >= lep_cuts[lep_bin-1] and lep_cut < lep_cuts[lep_bin]:
                find_bin = lep_bin
        lep_name = "%.2f"%lep_cut
	print "cut is " + lep_name + "    bin is %d" %find_bin
        projections[lep_name] = TH2D("projection_lep_"+lep_name, "projection_lep_"+lep_name, len(mu1_cuts)-1, mu1_cuts, len(mu2_cuts)-1, mu2_cuts)
        for mu1_bin in range(1,len(mu1_cuts)):
            for mu2_bin in range(1,len(mu2_cuts)):
                sens = sens_hist.GetBinContent(mu1_bin, mu2_bin, find_bin)
                projections[lep_name].SetBinContent(mu1_bin, mu2_bin, sens)



    out_file.cd()
    sens_hist.Write()
    for key in projections.keys():
	projections[key].SetMarkerSize(1.5)
#	projections[key].GetXaxis().SetRange(2,len(mu1_cuts)-1)
#	projections[key].GetYaxis().SetRange(2,len(mu2_cuts)-1)
	projections[key].SetStats(False)
        projections[key].Write()

    out_file.Close()
main()



