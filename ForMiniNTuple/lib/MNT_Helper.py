#! /usr/bin/env python

##################################################
#####  functions for making plots from MNT   #####
##################################################
import array
import Plot_Configs as PC
from ROOT import *

#####################################################################

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
    ratio_graph.SetMinimum(0.4)
    ratio_graph.SetMaximum(2.0)
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

###########################################################

def Sample_Name(Sample_ID):
    if Sample_ID == 0:	return "Data"

    elif Sample_ID > 0:
      if Sample_ID == 25:	return "ggH"
      elif Sample_ID == 10225: 	return "VBF"
      elif Sample_ID == 2425: 	return "WH"
      elif Sample_ID == 2325: 	return "ZH"
      elif Sample_ID == 60625: 	return "ttH"
      else: 	return "unknown_sig"

    elif Sample_ID < 0:
      # main
      if Sample_ID == -23: 		return "DY"
      elif Sample_ID == -606: 		return "ttbar"
      # tX
      elif Sample_ID == -624: 		return "tW"
      elif Sample_ID == -62500: 	return "tHq"
      elif Sample_ID == -62524: 	return "tHW"
      elif Sample_ID == -62300: 	return "tZq"
      elif Sample_ID == -62324: 	return "tZW"
      # ttX
      elif Sample_ID == -60623: 	return "ttZ"
      elif Sample_ID == -60625: 	return "ttH_bkg"
      elif Sample_ID == -60624:         return "ttW"
      elif Sample_ID == -6062424: 	return "ttWW"
      # multi boson
      elif Sample_ID == -2323: 		return "ZZ"
      elif Sample_ID == -23230000:      return "ggZZ"
      elif Sample_ID == -2423:       	return "WZ"
      elif Sample_ID == -2424:       	return "WW"
      elif Sample_ID == -232323:       	return "ZZZ"
      elif Sample_ID == -242323:       	return "WZZ"
      elif Sample_ID == -242423:       	return "WWZ"
      elif Sample_ID == -242424:        return "WWW"
      else:	return "unknown_bkg"

    else:
	print "what happened"
	return "unknown_samp"
# end def Sample_Name(Sample_ID):



def FillHistTerm(histos, term, cfg, value, Sample_ID, event_wgt):
    samp = Sample_Name(Sample_ID)

    # unknown sample
    if samp == "unknown_sig" or samp == "unknown_bkg": return

    #blind
    if samp == "Data" and term == "dimu_mass" and value > 120 and value < 130:
        return
    if samp == "Data" and ("BDT_mass" in term) and value > 0.4:
	return

    # fill data
    if samp == "Data":
        histos[term]["Data"].Fill(value, 1)
    # fill signal
    elif samp in cfg.signals:
        histos[term][samp].Fill(value, event_wgt)

    # fill background
    elif samp in cfg.bkgs:
      if samp == "DY":
	histos[term][samp].Fill(value, event_wgt / 2)  # used two high mass DY samples
#       histos[term][samp].Fill(value, event_wgt / 3)  # used three regular mass DY samples
      elif samp == "ttbar":
	histos[term][samp].Fill(value, event_wgt / 2)  # used two ttbar samples
      else:
	histos[term][samp].Fill(value, event_wgt)

    # if ZZ and ggZZ are combined
    elif (samp == "ZZ") and ("ZZ+ggZZ" in cfg.bkgs):
        histos[term]["ZZ+ggZZ"].Fill(value, event_wgt)
    elif (samp == "ggZZ") and ("ZZ+ggZZ" in cfg.bkgs):
	histos[term]["ZZ+ggZZ"].Fill(value, event_wgt)

    # if WZ and WW are combined
    elif (samp == "WZ") and ("WZ+WW" in cfg.bkgs):
	histos[term]["WZ+WW"].Fill(value, event_wgt)
    elif (samp == "WW") and ("WZ+WW" in cfg.bkgs):
        histos[term]["WZ+WW"].Fill(value, event_wgt)

    # if tribosons are combined
    elif (samp == "ZZZ") and ("triboson" in cfg.bkgs):
	histos[term]["triboson"].Fill(value, event_wgt)
    elif (samp == "ZZW") and ("triboson" in cfg.bkgs):
        histos[term]["triboson"].Fill(value, event_wgt)
    elif (samp == "ZWW") and ("triboson" in cfg.bkgs):
        histos[term]["triboson"].Fill(value, event_wgt)
    elif (samp == "WWW") and ("triboson" in cfg.bkgs):
        histos[term]["triboson"].Fill(value, event_wgt)

    else:
	histos[term]["others"].Fill(value, event_wgt)
    return
# end def FillHistTerm(histos, term, cfg, value, Sample_ID, event_wgt):
#############################################################

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
	in_file = TFile.Open("/afs/cern.ch/work/x/xzuo/public/H2Mu/2017/Histograms/lepMVA_SF_v1/SF_muon_variable_eta_bin.root", "READ")
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


