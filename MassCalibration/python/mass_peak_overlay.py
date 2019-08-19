###################################################################
###    make 2D mass calibration plots from mass histograms      ###
###################################################################

## Basic python includes for manipulating files
import os
import sys
import array

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from MassCal_Helper import FitVoigtian, GetColor, WriteOverlay, WriteSummary

gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2016/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2016/Histograms'

LABEL = 'MassCal_KinRoch_approx/2D_muP_d0_muN_d0'
nameX = "muP_d0"
nameY = "muN_d0"
CAT   = "NONE_NONE"

def main():
    file_dir = PLOT_DIR+"/"+LABEL
    out_file = TFile( file_dir + "/mass_overlay" + ".root", "RECREATE")
    in_file  = TFile.Open( file_dir + "/all_samples" + ".root", "READ")

    gen_dir  = out_file.mkdir("pull_gen")
    cor_dir  = out_file.mkdir("data_over_MC")
    ovl_dir  = out_file.mkdir("overlay")
    cal_dir  = out_file.mkdir("individual_cal_plots")
    fit_dir  = out_file.mkdir("fits_specifics")

    samples = ["ZJets_AMC", "data"]
    pt_cals = ["gen", "PF", "Roch", "Kinfit", "KinRoch", "KaMu", "KinKaMu"]

    summary_info = in_file.Get( pt_cals[0] + "/summary_" + samples[0] + "_" + CAT + "_" + pt_cals[0]).Clone()
    binningX = summary_info.GetXaxis().GetXbins()  ## return type is TArrayD
    binningY = summary_info.GetYaxis().GetXbins()  ## return type is TArrayD

    print binningX

    mean_plots = {}
    reso_plots = {}
    pull_plots = {}
    mass_PP_NN = {}
    mass_PN_NP = {}
    canv = {}
    legend = {}
    for pt_cal in pt_cals:
      for sample in samples:
	mean_plots[sample+pt_cal] = summary_info.Clone("Mean_" + sample + "_" + pt_cal) 
	mean_plots[sample+pt_cal].SetTitle("Mean_" + sample + "_" + pt_cal) 
	reso_plots[sample+pt_cal] = summary_info.Clone("Reso_" + sample + "_" + pt_cal) 
	reso_plots[sample+pt_cal].SetTitle("Reso_" + sample + "_" + pt_cal)

	mass_PP_NN[sample+pt_cal] = None
	mass_PN_NP[sample+pt_cal] = None
	for i in range(binningX.GetSize()-1):
	  for j in range(binningY.GetSize()-1):
	    x_low  = binningX.GetAt(i)
	    x_high = binningX.GetAt(i+1)
	    y_low  = binningY.GetAt(j)
            y_high = binningY.GetAt(j+1)

	    print "looking at, %s,   %s,   %f to %f   " %(sample, pt_cal, x_low, x_high)
	    hist_name = pt_cal + "/" + sample + "_" + pt_cal + "_%s_%8.4f_to_%8.4f_%s_%8.4f_to_%8.4f" %( nameX, x_low, x_high, nameY, y_low, y_high)
	    hist_name = hist_name.replace(' ','').replace("-","m").replace(".","p")
	    print hist_name
	    mass_hist = in_file.Get(hist_name).Clone()

	    if x_low >= 0 and y_high <= 0:   #PP and NN
	      if mass_PP_NN[sample+pt_cal] == None: mass_PP_NN[sample+pt_cal] = mass_hist.Clone()
	      else: mass_PP_NN[sample+pt_cal].Add(mass_hist)
	    if x_high <= 0 and y_low >= 0:   #PN and NP
              if mass_PN_NP[sample+pt_cal] == None: mass_PN_NP[sample+pt_cal] = mass_hist.Clone()
              else: mass_PN_NP[sample+pt_cal].Add(mass_hist)

	mass_PP_NN[sample+pt_cal].SetTitle(sample + ", " + pt_cal + " mass")	
	mass_PP_NN[sample+pt_cal].SetLineColor(kRed)
	mass_PP_NN[sample+pt_cal].SetLineWidth(2)
	mean_PP, mean_err_PP, reso_PP, reso_err_PP = FitVoigtian(mass_PP_NN[sample+pt_cal])
	mass_PN_NP[sample+pt_cal].SetTitle("color blue")
	mass_PN_NP[sample+pt_cal].SetLineColor(kBlue)
        mass_PN_NP[sample+pt_cal].SetLineWidth(2)
	mean_PN, mean_err_PN, reso_PN, reso_err_PN = FitVoigtian(mass_PN_NP[sample+pt_cal])

	ovl_dir.cd() 

	canv[sample+pt_cal] = TCanvas("mass_" + sample + "_" + pt_cal, "mass_" + sample + "_" + pt_cal, 600,600)
 	canv[sample+pt_cal].cd()
	legend[sample+pt_cal] = TLegend(0.6,0.5,0.99,0.94)

	gStyle.SetOptFit(0000)
	mass_PP_NN[sample+pt_cal].Draw()
	mass_PN_NP[sample+pt_cal].Draw("SAME")	
	legend[sample+pt_cal].AddEntry( mass_PP_NN[sample+pt_cal], "d0_Pos > 0, d0_Neg < 0", "LPE")
	legend[sample+pt_cal].AddEntry( None, "mass_mean = %.2f +/- %.2f" %(mean_PP,mean_err_PP), "")
	legend[sample+pt_cal].AddEntry( None, "mass_reso = %.2f +/- %.2f" %(reso_PP,reso_err_PP), "")
	legend[sample+pt_cal].AddEntry( mass_PN_NP[sample+pt_cal], "d0_Pos < 0, d0_Neg > 0", "LPE")
        legend[sample+pt_cal].AddEntry( None, "mass_mean = %.2f +/- %.2f" %(mean_PN,mean_err_PN), "")
        legend[sample+pt_cal].AddEntry( None, "mass_reso = %.2f +/- %.2f" %(reso_PN,reso_err_PN), "")

	legend[sample+pt_cal].Draw()
	canv[sample+pt_cal].Write()



main()
