#! /usr/bin/env python

####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os
import array

from ROOT import *
from MNT_Helper import LinearStack, RatioPlot, GetSF
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
LEP = 'ele'

#mu1_cuts = array.array('d', [ round(1.0*MVA/5 - 1.0 , 2) for MVA in range(10) ])
#mu2_cuts = array.array('d', [ round(1.0*MVA/5 - 1.0 , 2) for MVA in range(10) ])
mu1_cuts = array.array('d', [ -1.0, -0.4, 0.4, 0.8, 1.0])
mu2_cuts = array.array('d', [ -1.0, -0.4, 0.4, 0.8, 1.0])
lep_cuts = array.array('d', [ -1.0, -0.4, 0.4, 0.8, 1.0])


def main():
    out_name = "lepMVA_scan.root"

    file_dir = PLOT_DIR+"/"+LABEL+"/"
    out_file = TFile( file_dir + "plots/" + out_name , "RECREATE")
    file_chain = TChain("tree","chain");
    file_chain.Add(file_dir + "all_samples.root")
#    file_chain.Add( file_dir + "signal.root")
#    file_chain.Add( file_dir + "bkg.root")

    project_dir = out_file.mkdir("projections")

    cut_hist_sig = TH3D(LEP + "_sig_yield", LEP + "_sig_yield", len(mu1_cuts)-1, mu1_cuts, len(mu2_cuts)-1, mu2_cuts, len(lep_cuts)-1, lep_cuts)
    cut_hist_bkg = TH3D(LEP + "_bkg_yield", LEP + "_bkg_yield", len(mu1_cuts)-1, mu1_cuts, len(mu2_cuts)-1, mu2_cuts, len(lep_cuts)-1, lep_cuts)
    sens_hist = TH3D(LEP + "_sensitivity", LEP + "_sensitivity", len(mu1_cuts)-1, mu1_cuts, len(mu2_cuts)-1, mu2_cuts, len(lep_cuts)-1, lep_cuts)



    print file_chain.GetEntries()
    for iEvt in range( file_chain.GetEntries() ):
	if iEvt % 1000 == 1:
	    print "event %d"%iEvt
#	    continue
	file_chain.GetEvent(iEvt)
	if file_chain.Sample_ID == 0 :
            continue
	if file_chain.dimu_mass < 120 or file_chain.dimu_mass > 130:
	    continue

	for mu1_bin in range(1,len(mu1_cuts)):
	    mu1_cut = mu1_cuts[mu1_bin-1]
	    if file_chain.mu1_lepMVA < mu1_cut:
		continue
	    mu1_SF = GetSF("muon", file_chain.mu1_pt, file_chain.mu1_eta, mu1_cut)

	    for mu2_bin in range(1,len(mu2_cuts)):
	    	mu2_cut = mu2_cuts[mu2_bin-1]
		if file_chain.mu2_lepMVA < mu2_cut:
		    continue
		mu2_SF = GetSF("muon", file_chain.mu2_pt, file_chain.mu2_eta, mu2_cut)

		for lep_bin in range(1,len(lep_cuts)):
		    lep_cut = lep_cuts[lep_bin-1]
		    if file_chain.lep_lepMVA < lep_cut:
		 	continue
		    lep_SF = GetSF(LEP, file_chain.lep_pt, file_chain.lep_eta, lep_cut)
		    MVA_SF = mu1_SF * mu2_SF * lep_SF
		    if file_chain.Sample_ID > 0:
			cut_hist_sig.Fill(mu1_cut, mu2_cut, lep_cut, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
		    if file_chain.Sample_ID < 0:
			cut_hist_bkg.Fill(mu1_cut, mu2_cut, lep_cut, file_chain.xsec_norm * file_chain.event_wgt * MVA_SF)
    # looping over events finished
    for mu1_bin in range(1,len(mu1_cuts)):
	for mu2_bin in range(1,len(mu2_cuts)):
	    for lep_bin in range(1,len(lep_cuts)):
		sig_yield = cut_hist_sig.GetBinContent(mu1_bin, mu2_bin, lep_bin)
		bkg_yield = cut_hist_bkg.GetBinContent(mu1_bin, mu2_bin, lep_bin)
		if bkg_yield == 0 :
		    sens_hist.SetBinContent(mu1_bin, mu2_bin, lep_bin, 0)
		else:
		    sens = sig_yield / bkg_yield ** 0.5 
		    sens_hist.SetBinContent(mu1_bin, mu2_bin, lep_bin, sens)


    projections = {}
    for lep_cut in [-1.0, -0.4, 0.4, 0.8]:
	find_bin = 0
	for lep_bin in range(1,len(lep_cuts)):
	    if lep_cut >= lep_cuts[lep_bin-1] and lep_cut < lep_cuts[lep_bin]:
		fine_bin = lep_bin
	lep_name = "%.2f"%lep_cut
	projections[lep_name] = TH2D("projection_ele_"+lep_name, "projection_ele_"+lep_name, len(mu1_cuts)-1, mu1_cuts, len(mu2_cuts)-1, mu2_cuts)
	for mu1_bin in range(1,len(mu1_cuts)):
	    for mu2_bin in range(1,len(mu2_cuts)):
		sens = sens_hist.GetBinContent(mu1_bin, mu2_bin, find_bin)
		projections[lep_name].SetBinContent(mu1_bin, mu2_bin, sens)


    out_file.cd()
    cut_hist_sig.Write()
    cut_hist_bkg.Write()
    sens_hist.Write()

    project_dir.cd()
    for key in projections.keys():
	projections[key].Write()

    out_file.Close()
main()



