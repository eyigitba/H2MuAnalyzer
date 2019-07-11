#! /usr/bin/env python

##################################################
#####  functions for making plots from MNT   #####
##################################################
import array
import ctypes 
import math
from ROOT import *

#####################################################################

def FitVoigtian(hist):
    if hist.Integral() < 1000.0:
	return 91.0, 91.0, 10.0, 10.0

    minM     = round(hist.GetXaxis().GetXmin(), 2) + 1.0
    maxM     = round(hist.GetXaxis().GetXmax(), 2) - 1.0

    F_voigt = TF1("voigtian", "[0]*TMath::Voigt(x-[1],[2],[3])", minM, maxM)
    F_voigt.SetParNames(   "Norm",              "Mean",           "Sigma", 		  "Gamma" )
    F_voigt.SetParameters( hist.GetMaximum(),   hist.GetMean(),   hist.GetRMS() / 2.0,    2.5     )
    F_voigt.SetParLimits(2, 0,    hist.GetRMS() * 5.0 )  ## Sigma range
    F_voigt.SetParLimits(1, minM, maxM			    )  ## Mean  range
    F_voigt.FixParameter(3, 2.5)  ## fix Z natural width

    gStyle.SetOptFit(0011)

    for i in range(10):
        hist.Fit("voigtian", "QR")
	F_voigt.SetParameters( F_voigt.GetParameter(0), F_voigt.GetParameter(1), F_voigt.GetParameter(2), F_voigt.GetParameter(3) )

    mean_val = F_voigt.GetParameter(1)
    mean_err = F_voigt.GetParError(1)
    reso_val = F_voigt.GetParameter(2)
    reso_err = F_voigt.GetParError(2)

    print "fit hist %s, "  %hist.GetName()
    print "chi^2 = %f, \n\n" %hist.GetFunction("voigtian").GetChisquare()
    print "norm  = %f +- %f, " %( F_voigt.GetParameter(0), F_voigt.GetParError(0) )
    print "mean  = %f +- %f, " %( F_voigt.GetParameter(1), F_voigt.GetParError(1) )
    print "sigma = %f +- %f, " %( F_voigt.GetParameter(2), F_voigt.GetParError(2) )
    print "gamma = %f +- %f, " %( F_voigt.GetParameter(3), F_voigt.GetParError(3) )

    return mean_val, mean_err, reso_val, reso_err




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
	legend.AddEntry(graphs[sample+pt_cal], pt_cal.replace('KinRoch','Kinfit+Roch'))
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
        legend.AddEntry(graphs[sample+pt_cal], sample)
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
    canv.SaveAs( plot_dir + "/" + nameX + "_kinfit_eta_" + term + ".png")





