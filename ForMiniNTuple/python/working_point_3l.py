#! /usr/bin/env python

####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os

from ROOT import *
from MNT_Helper import LoadColors, LoadSampleNames, FillHistTerm, LinearStack, RatioPlot, GetSF
gROOT.SetBatch(True)

def MakeMassStack(in_file_name, out_file_dir, out_file_name, color_set, sample_set, mu1_MVAcut, mu2_MVAcut, lep_MVAcut):
    file_chain = TChain("tree","chain")
    file_chain.Add(in_file_name)
    out_file = TFile( out_file_dir + out_file_name, "RECREATE")

    signals = []
    bkgs = []
    data = []
    LoadSampleNames(signals, bkgs, data, sample_set)
    colors = {}
    LoadColors(colors, color_set)

    histos = {}
    histos["dimu_mass"] = {}
    for sample in signals + bkgs + ["data"]:
	hist_name = sample + "mu1_%.2f_"%mu1_MVAcut + "mu2_%.2f_"%mu2_MVAcut + "lep_%.2f"%lep_MVAcut
	print hist_name
	histos["dimu_mass"][sample] = TH1D(hist_name, hist_name, 12,100,160)

    for iEvt in range( file_chain.GetEntries() ):
	if iEvt % 1000 == 1:
	    print "event %d" %iEvt
	file_chain.GetEvent(iEvt)
	if file_chain.mu1_lepMVA < mu1_MVAcut or file_chain.mu2_lepMVA < mu2_MVAcut or file_chain.lep_lepMVA < lep_MVAcut:
	    continue
	lepMVA_SF = GetSF("muon", file_chain.mu1_pt, file_chain.mu1_eta, mu1_MVAcut) * GetSF("muon", file_chain.mu2_pt, file_chain.mu2_eta, mu2_MVAcut) * GetSF("muon", file_chain.lep_pt, file_chain.lep_eta, lep_MVAcut)# "muon" or "ele"
        if file_chain.Sample_ID==0:
            lepMVA_SF = 1.0

	FillHistTerm(histos, "dimu_mass", signals, bkgs, file_chain.dimu_mass, file_chain.Sample_ID, file_chain.xsec_norm * file_chain.event_wgt * lepMVA_SF)

    out_file.cd()
    for sample in signals + bkgs + ["data"]:
	histos["dimu_mass"][sample].Write()
    out_file.Close()


