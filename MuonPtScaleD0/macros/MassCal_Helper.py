#! /usr/bin/env python

##################################################
#####  functions for making plots from MNT   #####
##################################################
import array
import ctypes 
import math
from ROOT import *
import ROOT.RooFit as RF

from ROOT import gInterpreter, gSystem
#gInterpreter.ProcessLine('#include "H2MuAnalyzer/FitBackground/python/Custom_DSCB/My_double_CB.h"')
# gSystem.Load("/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/FitBackground/python/Custom_DSCB/HZZ2L2QRooPdfs_cc.so")
# from ROOT import RooDoubleCB
#####################################################################

def GetHistInfo(hist):
    mean_val = hist.GetMean()
    mean_err = 0.0
    reso_val = hist.GetRMS()
    reso_err = 0.0
    return mean_val, mean_err, reso_val, reso_err    

def FitVoigtian(hist, exp = False):
    if hist.Integral() < 1000.0:
	     return 91.0, 91.0, 10.0, 10.0

    minM     = round(hist.GetMean() - min(hist.GetRMS()*2, 3.5), 2)
    maxM     = round(hist.GetMean() + min(hist.GetRMS()*2, 3.5), 2)

    F_voigt = TF1("voigtian_exp", "[0]*TMath::Voigt(x-[1],[2],[3]) + [4]*TMath::Exp([5]*x + [6])", minM, maxM)
    F_voigt.SetParNames(   "Norm_voigt",        "Mean",           "Sigma",                "Gamma",  "Norm_exp",	         "coef_exp",  "const_exp")
    F_voigt.SetParameters( hist.GetMaximum(),   hist.GetMean(),   hist.GetRMS() / 2.0,    2.5,      hist.GetMaximum(),   0.00001,         0.0     )
    F_voigt.SetParLimits(0, 0,    hist.GetEntries()) ## Norm_voigt range
    F_voigt.SetParLimits(2, 0,    hist.GetRMS() * 5.0 )  ## Sigma range
    F_voigt.SetParLimits(1, minM, maxM			    )  ## Mean  range
#    F_voigt.SetParLimits(3, 2.0, 3.0) ## test floating Z width
    F_voigt.FixParameter(3, 2.5)  ## fix Z natural width

    if exp:
    	F_voigt.SetParLimits(4, 0, hist.GetMaximum()) ## norm_exp range
    	F_voigt.SetParLimits(5, 0, 0.0001) ## coef_exp range
    	#F_voigt.SetParLimits(6, 0, 80) ## const_exp range
    	F_voigt.FixParameter(6, 0.0) ## test fixing const_exp
    else:
    	F_voigt.FixParameter(4, 0.0) ## fix exp norm to 0
    	F_voigt.FixParameter(5, 1.0) ## fix coef_exp. value does not matter
    	F_voigt.FixParameter(6, 0.0) ## fix const_exp. value does not matter


    for i in range(10):
        hist.Fit("voigtian_exp", "QR")
        F_voigt.SetParameters( F_voigt.GetParameter(0), F_voigt.GetParameter(1), F_voigt.GetParameter(2), F_voigt.GetParameter(3), F_voigt.GetParameter(4), F_voigt.GetParameter(5), F_voigt.GetParameter(6) )

    gStyle.SetOptFit(0111)
    mean_val = F_voigt.GetParameter(1)
    mean_err = F_voigt.GetParError(1)
    reso_val = F_voigt.GetParameter(2)
    reso_err = F_voigt.GetParError(2)

    if hist.GetTitle() == "color blue": hist.GetFunction("voigtian_exp").SetLineColor(kBlue)

    print "fit hist %s, "  %hist.GetName()
    print "chi^2 = %f, \n\n" %hist.GetFunction("voigtian_exp").GetChisquare()
    print "norm  = %f +- %f, " %( F_voigt.GetParameter(0), F_voigt.GetParError(0) )
    print "mean  = %f +- %f, " %( F_voigt.GetParameter(1), F_voigt.GetParError(1) )
    print "sigma = %f +- %f, " %( F_voigt.GetParameter(2), F_voigt.GetParError(2) )
    print "gamma = %f +- %f, " %( F_voigt.GetParameter(3), F_voigt.GetParError(3) )

    return mean_val, mean_err, reso_val, reso_err


def FitGaussian(hist):
    if hist.Integral() < 1000.0:
        return 125.0, 125.0, 10.0, 10.0

    minM     = round(hist.GetXaxis().GetXmin(), 2) + 1.0
    maxM     = round(hist.GetXaxis().GetXmax(), 2) - 1.0

    F_gaus = TF1("gauss", "[0]*TMath::Gaus(x,[1],[2])", minM, maxM)
    F_gaus.SetParNames(   "Norm_gaus",	     "Mean",         "Sigma")
    F_gaus.SetParameters( hist.GetMaximum(), hist.GetMean(), hist.GetRMS())

    F_gaus.SetParLimits(1, minM, maxM) ## Mean range
    F_gaus.SetParLimits(2, 0, hist.GetRMS()*5.0 ) ## Sigma range

    for i in range(10):
	     hist.Fit("gauss", "QR")
	     F_gaus.SetParameters( F_gaus.GetParameter(0), F_gaus.GetParameter(1), F_gaus.GetParameter(2) )

    gStyle.SetOptFit(0111)
    mean_val = F_gaus.GetParameter(1)
    mean_err = F_gaus.GetParError(1)
    reso_val = F_gaus.GetParameter(2)
    reso_err = F_gaus.GetParError(2)

    if hist.GetTitle() == "color blue": hist.GetFunction("gauss").SetLineColor(kBlue)

    print "fit hist %s, "  %hist.GetName()
    print "chi^2 = %f, \n\n" %hist.GetFunction("gauss").GetChisquare()
    print "norm  = %f +- %f, " %( F_gaus.GetParameter(0), F_gaus.GetParError(0) )
    print "mean  = %f +- %f, " %( F_gaus.GetParameter(1), F_gaus.GetParError(1) )
    print "sigma = %f +- %f, " %( F_gaus.GetParameter(2), F_gaus.GetParError(2) )

    return mean_val, mean_err, reso_val, reso_err


def FitCrystalBall(hist):
    if hist.Integral() < 1000.0:
        return 125.0, 125.0, 10.0, 10.0

#    minM     = round(hist.GetXaxis().GetXmin(), 2) + 1.0
#    maxM     = round(hist.GetXaxis().GetXmax(), 2) - 1.0

    minM  = round(hist.GetMean() - min(hist.GetRMS()*2, 3.5), 2)
    maxM  = round(hist.GetMean() + min(hist.GetRMS()*2, 3.5), 2)

    F_CB = TF1("CRB","crystalball", minM, maxM)
    F_CB.SetParNames(   "CB_Norm",         "Mean",         "Sigma", "Alpha", "N")
    F_CB.SetParameters( hist.GetMaximum(), hist.GetMean(), 1.3,     1.6,     1.3)

    F_CB.SetParLimits(1, minM, maxM) ## Mean range
    F_CB.SetParLimits(2, 0, hist.GetRMS()*5.0) ## Sigma range
    F_CB.SetParLimits(3, 0, 5) ## Alpha range
    F_CB.SetParLimits(4, 0, 100)

    for i in range(10):
        hist.Fit("CRB", "QR")
        F_CB.SetParameters( F_CB.GetParameter(0), F_CB.GetParameter(1), F_CB.GetParameter(2), F_CB.GetParameter(3), F_CB.GetParameter(4) )

    gStyle.SetOptFit(0111)
    mean_val = F_CB.GetParameter(1)
    mean_err = F_CB.GetParError(1)
    reso_val = F_CB.GetParameter(2)
    reso_err = F_CB.GetParError(2)

    if hist.GetTitle() == "color blue": hist.GetFunction("CRB").SetLineColor(kBlue)

    print "fit hist %s, "  %hist.GetName()
    print "chi^2 = %f, \n\n" %hist.GetFunction("CRB").GetChisquare()
    print "norm  = %f +- %f, " %( F_CB.GetParameter(0), F_CB.GetParError(0) )
    print "mean  = %f +- %f, " %( F_CB.GetParameter(1), F_CB.GetParError(1) )
    print "sigma = %f +- %f, " %( F_CB.GetParameter(2), F_CB.GetParError(2) )

    return mean_val, mean_err, reso_val, reso_err



# def FitDSCB(hist):
#     if hist.Integral() < 1000.0:
#         return (None, 125.0, 125.0, 10.0, 10.0)

# #    minM     = round(hist.GetXaxis().GetXmin(), 2) + 1.0
# #    maxM     = round(hist.GetXaxis().GetXmax(), 2) - 1.0

#     minM  = round(hist.GetMean() - min(hist.GetRMS()*2, 3.5), 2)
#     maxM  = round(hist.GetMean() + min(hist.GetRMS()*2, 3.5), 2)


#     h_mean = hist.GetMean()
#     h_rms  = hist.GetRMS()
#     var_name = 'mass'
#     x_var = RooRealVar(var_name, var_name, minM, maxM)

#     params = []
#     params.append( RooRealVar('mean',  'mean',  h_mean, h_mean - 1.5*h_rms, h_mean + 1.5*h_rms) )
#     params.append( RooRealVar('sigma', 'sigma', h_rms,  0.0,                5.0*h_rms) )
#     params.append( RooRealVar('a1',    'a1',    1.05,   1.05,               1.05) )
#     params.append( RooRealVar('n1',    'n2',    20.0,   0.0,                50.0) )
#     params.append( RooRealVar('a2',    'a2',    1.25,   1.25,               1.25) )
#     params.append( RooRealVar('n2',    'n2',    20.0,   0.0,                50.0) )

#     dat_hist = RooDataHist('dat_' + hist.GetName(), 'dat_'+hist.GetName(), RooArgList(RooArgSet(x_var)), hist)

#     func_DSCB = ( RooDoubleCB('DSCB', 'DSCB', x_var, params[0], params[1], params[2], params[3], params[4], params[5]) )
#     arg_sets  = RooArgSet(func_DSCB) 


#     func_DSCB.fitTo(dat_hist, RF.Save())
#     fit_hist = func_DSCB.createHistogram( var_name, hist.GetNbinsX() )

#     f_norm = fit_hist.Integral()
#     h_norm = hist.Integral()
#     fit_hist.Scale( h_norm / f_norm )
#     fit_hist.SetLineColor( kRed )


#     mean_val = params[0].getValV()
#     mean_err = params[0].getError()
#     reso_val = params[1].getValV()
#     reso_err = params[1].getError()


#     line_color = kBlack
#     if "GeoRoch"   in hist.GetName(): line_color = kBlue
#     if "GeoBSRoch" in hist.GetName(): line_color = kRed

#     frame_plot = x_var.frame()
#     frame_plot.SetNameTitle('RooPlot_' + hist.GetName(), 'RooPlot_' + hist.GetTitle())
#     dat_hist. plotOn(frame_plot, RF.MarkerStyle(1), RF.MarkerColor(line_color), RF.LineColor(line_color))
#     func_DSCB.plotOn(frame_plot, RF.LineColor(line_color), RF.Range('FULL'))
#     func_DSCB.paramOn(frame_plot, RF.Layout(0.45, 0.90, 0.90))

#     print "fit hist %s, "  %hist.GetName()
#     print "norm  = %f , "      %( h_norm / f_norm )
#     print "mean  = %f +- %f, " %( params[0].getValV(), params[0].getError() )
#     print "sigma = %f +- %f, " %( params[1].getValV(), params[1].getError() )
#     print "a1    = %f +- %f, " %( params[2].getValV(), params[2].getError() )
#     print "n1    = %f +- %f, " %( params[3].getValV(), params[3].getError() )
#     print "a2    = %f +- %f, " %( params[4].getValV(), params[4].getError() )
#     print "n2    = %f +- %f, " %( params[5].getValV(), params[5].getError() )

#     return (frame_plot, mean_val, mean_err, reso_val, reso_err)





def GetColor(sample, pt_cal):
    if sample == "ZJets_MG_1" or sample == "ZJets_AMC":
      if pt_cal == "PF" or pt_cal == "Kin_vs_d0kin_BB" or pt_cal == "Kin_vs_d0kin_d0PV_N50_N15":
        return kRed
      if pt_cal == "Roch" or pt_cal == "good_Kinfit" or pt_cal == "Kin_vs_d0kin_BE" or pt_cal == "Kin_vs_d0kin_d0PV_N15_N05":
        return kAzure
      if pt_cal == "Kinfit" or pt_cal == "Kin_vs_d0kin_EE" or pt_cal == "Kin_vs_d0kin_d0PV_P05_P15":
	      return kSpring
      if pt_cal == "KinRoch" or pt_cal == "Kin_vs_d0kin" or pt_cal == "Kin_vs_d0kin_d0PV_P15_P50":
	      return kOrange

    if sample == "data":
      if pt_cal == "PF" or pt_cal == "Kin_vs_d0kin_BB" or pt_cal == "Kin_vs_d0kin_d0PV_N50_N15":
        return kMagenta
      if pt_cal == "Roch" or pt_cal == "good_Kinfit" or pt_cal == "Kin_vs_d0kin_BE" or pt_cal == "Kin_vs_d0kin_d0PV_N15_N05":
        return kBlue + 2
      if pt_cal == "Kinfit" or pt_cal == "Kin_vs_d0kin_EE" or pt_cal == "Kin_vs_d0kin_d0PV_P05_P15":
	      return kGreen + 2
      if pt_cal == "KinRoch" or pt_cal == "Kin_vs_d0kin" or pt_cal == "Kin_vs_d0kin_d0PV_P15_P50":
	      return kOrange - 6
    else:
        return kBlack




def WriteOverlay(graphs, term, samples, pt_cals):
    for sample in samples:
      canv = TCanvas("overlay_" + term + "_" + sample, "overlay_" + term + "_" + sample, 600,600)
      canv.SetGridx(1)
      canv.SetGridy(1)
      legend = TLegend(0.7,0.7,1,1)
      canv.cd()

      for pt_cal in pt_cals:
        graphs[sample+pt_cal].Draw()  ## the first plot cannot have option "SAME", otherwise it is empty
      for pt_cal in pt_cals:
	      graphs[sample+pt_cal].Draw("SAME")
	      legend.AddEntry(graphs[sample+pt_cal], pt_cal.replace('KinRoch','Kinfit+Roch'), "LPE")
      legend.Draw()
      canv.Update()
      canv.Write()

    for pt_cal in pt_cals:
      canv = TCanvas("overlay_" + term + "_" + pt_cal, "overlay_" + term + "_" + pt_cal, 600,600)
      canv.SetGridx(1)
      canv.SetGridy(1)
      legend = TLegend(0.7,0.7,1,1)
      canv.cd()
      for sample in samples:
        graphs[sample+pt_cal].Draw()  ## the first plot cannot have option "SAME", otherwise it is empty
      for sample in samples:
        graphs[sample+pt_cal].Draw("SAME")
        legend.AddEntry(graphs[sample+pt_cal], sample, "LPE")

      legend.Draw()
      canv.Update()
      canv.Write()




def WriteSummary(graphs, term, nameX, samples, pt_cals, plot_dir):
    canv = TCanvas("summary_" + term, "summary_" + term, 600,600)
    
    ## upper pad
    upper_pad = TPad("UP_"+term, "UP_"+term, 0,0.3, 1,1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.SetGridx()
    upper_pad.SetGridy()
    upper_pad.Draw()
    upper_pad.cd()
    legend_U = TLegend(0.7,0.7,1,1)
    legend_U.SetHeader("mass_" + term + "_vs_" + nameX, "C")
    for pt_cal in pt_cals:
      for sample in samples:
        if term == "mean":
          graphs[sample+pt_cal].SetMaximum(92.5)
          graphs[sample+pt_cal].SetMinimum(90.0)
        else:
          graphs[sample+pt_cal].SetMaximum(2.8)
          graphs[sample+pt_cal].SetMinimum(0.8)
        graphs[sample+pt_cal].Draw()
    for pt_cal in pt_cals:
      for sample in samples:
        graphs[sample+pt_cal].Draw("SAME")
        legend_U.AddEntry(graphs[sample+pt_cal], pt_cal.replace('KinRoch','Kinfit+Roch') + "_" + sample, "LPE")
    legend_U.Draw()
   
    ## lower pad 
    canv.cd()
    lower_pad = TPad("LP_"+term, "LP_"+term, 0,0.05, 1,0.3)
    lower_pad.SetTopMargin(0.05)
    lower_pad.SetGridx()
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    ratios = {}
    legend_L = TLegend(0.8,0.6,1,1)
    legend_L.SetHeader("data/MC", "C")
    for pt_cal in pt_cals:
      ratios[pt_cal] = TGraphErrors()
      ratios[pt_cal].SetName("ratio_" + term + pt_cal)
      x_MC =   double(0.0)
      y_MC =   double(0.0)
      x_data = double(0.0)
      y_data = double(0.0)

      MC_sample = "ZJets_MG_1"
      if "ZJets_MG_1" not in samples and "ZJets_AMC" in samples:
	      MC_sample = "ZJets_AMC"
      for i in range( graphs[MC_sample+pt_cal].GetN() ):
        point_MC = graphs[MC_sample+pt_cal].GetPoint(i,  x_MC,  y_MC )
        point_data = graphs["data"+pt_cal].GetPoint(i,  x_data, y_data ) 
        x_err = graphs[MC_sample+pt_cal].GetErrorX(i)
        y_err_MC = graphs[MC_sample+pt_cal].GetErrorY(i)
        y_err_data = graphs["data"+pt_cal].GetErrorY(i)
        ratios[pt_cal].SetPoint(i, x_data, y_data/y_MC)
        ratios[pt_cal].SetPointError(i, x_err, y_data/y_MC * math.sqrt( (y_err_MC/y_MC) ** 2.0 + (y_err_data/y_data) ** 2.0 ) )

      ratios[pt_cal].SetLineColor( GetColor(MC_sample, pt_cal) )
      ratios[pt_cal].SetLineWidth(2)
      ratios[pt_cal].GetXaxis().SetLabelSize(0.12)
      ratios[pt_cal].GetYaxis().SetLabelSize(0.12)
      if term == "mean":
        ratios[pt_cal].SetMaximum(1.003)
        ratios[pt_cal].SetMinimum(0.997)
      else:
        ratios[pt_cal].SetMaximum(1.10)
        ratios[pt_cal].SetMinimum(0.90)
      ratios[pt_cal].Draw()

    for pt_cal in pt_cals:
      ratios[pt_cal].Draw("SAME")
      legend_L.AddEntry(ratios[pt_cal], pt_cal.replace('KinRoch','Kinfit+Roch'), "LPE")
    legend_L.Draw()

    canv.Update()
    canv.Write()
    canv.SaveAs( plot_dir + "/" + nameX + "_summary_" + term + ".png")





