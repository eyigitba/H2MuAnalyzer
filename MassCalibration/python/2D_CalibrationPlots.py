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
    out_file = TFile( file_dir + "/mass_cal_plots" + ".root", "RECREATE")
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
    for pt_cal in pt_cals:
      for sample in samples:
	mean_plots[sample+pt_cal] = summary_info.Clone("Mean_" + sample + "_" + pt_cal) 
	mean_plots[sample+pt_cal].SetTitle("Mean_" + sample + "_" + pt_cal) 
	reso_plots[sample+pt_cal] = summary_info.Clone("Reso_" + sample + "_" + pt_cal) 
	reso_plots[sample+pt_cal].SetTitle("Reso_" + sample + "_" + pt_cal)
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
	    mean_val, mean_err, reso_val, reso_err = FitVoigtian(mass_hist)
	    if mean_err == 91:  # if not a good fit
		mean_plots[sample+pt_cal].SetBinContent(i+1,j+1, 0.0)
		mean_plots[sample+pt_cal].SetBinError  (i+1,j+1, 91.0)
	 	reso_plots[sample+pt_cal].SetBinContent(i+1,j+1, 0.0)
                reso_plots[sample+pt_cal].SetBinError  (i+1,j+1, 10.0)

	    else:
		mean_plots[sample+pt_cal].SetBinContent(i+1,j+1, mean_val)
                mean_plots[sample+pt_cal].SetBinError  (i+1,j+1, mean_err)
		reso_plots[sample+pt_cal].SetBinContent(i+1,j+1, reso_val)
                reso_plots[sample+pt_cal].SetBinError  (i+1,j+1, reso_err)

	    fit_dir.cd()
            mass_hist.Write()

	cal_dir.cd()
	gStyle.SetPaintTextFormat("4.2f%")
	mean_plots[sample+pt_cal].GetZaxis().SetRangeUser(89.5,92.5)
	mean_plots[sample+pt_cal].SetXTitle("d0_muPos")
	mean_plots[sample+pt_cal].SetYTitle("d0_muNeg")
	mean_plots[sample+pt_cal].Write()
	reso_plots[sample+pt_cal].GetZaxis().SetRangeUser(1,2.5)
	reso_plots[sample+pt_cal].SetXTitle("d0_muPos")
	reso_plots[sample+pt_cal].SetYTitle("d0_muNeg")
	reso_plots[sample+pt_cal].Write()

    for sample in samples:
      for pt_cal in pt_cals:
        if pt_cal == "gen": continue
        pull_plots[sample+pt_cal] = summary_info.Clone("pull_" + sample + "_" + pt_cal)
	pull_plots[sample+pt_cal].SetTitle("pull_" + sample + "_" + pt_cal)
   
        for i in range(binningX.GetSize()-1):
          for j in range(binningY.GetSize()-1):
  	    mass_gen = mean_plots[sample+"gen"]. GetBinContent(i+1, j+1)
  	    mass_cal = mean_plots[sample+pt_cal].GetBinContent(i+1, j+1)

	    if mass_gen == 0 or mass_cal == 0:
	      pull_plots[sample+pt_cal].SetBinContent(i+1,j+1, 0.0)
	    else:
	      pull_plots[sample+pt_cal].SetBinContent(i+1,j+1, 100.0 * (mass_cal - mass_gen)/mass_gen )

	gen_dir.cd()
	gStyle.SetPaintTextFormat("4.2f%%")
	pull_plots[sample+pt_cal].GetZaxis().SetRangeUser(-2.0,2.0)
	pull_plots[sample+pt_cal].SetXTitle("d0_muPos")
	pull_plots[sample+pt_cal].SetYTitle("d0_muNeg")
	pull_plots[sample+pt_cal].Write()


    mean_dataMC = {}
    reso_dataMC = {}
    for pt_cal in pt_cals:
      for sample in samples:
	if pt_cal == "gen":  continue
	if sample == "data": continue
	mean_dataMC[sample+pt_cal] = summary_info.Clone("Mean_data_over_MC_" + sample + "_" + pt_cal)
	mean_dataMC[sample+pt_cal].SetTitle("Mean_data_over_MC_" + sample + "_" + pt_cal)
	reso_dataMC[sample+pt_cal] = summary_info.Clone("Reso_data_over_MC_" + sample + "_" + pt_cal)
	reso_dataMC[sample+pt_cal].SetTitle("Reso_data_over_MC_" + sample + "_" + pt_cal)
	
	for i in range(binningX.GetSize()-1):
          for j in range(binningY.GetSize()-1):
	    mass_data = mean_plots["data"+pt_cal].GetBinContent(i+1, j+1)
	    mass_MC   = mean_plots[sample+pt_cal].GetBinContent(i+1, j+1)

	    reso_data = reso_plots["data"+pt_cal].GetBinContent(i+1, j+1)
            reso_MC   = reso_plots[sample+pt_cal].GetBinContent(i+1, j+1)

	    if mass_data == 0 or mass_MC == 0:
		mean_dataMC[sample+pt_cal].SetBinContent(i+1,j+1, 0.0)
	    else:
		mean_dataMC[sample+pt_cal].SetBinContent(i+1,j+1, 100.0 * (mass_data - mass_MC)/mass_MC)

	    if reso_data == 0 or reso_MC == 0:
		reso_dataMC[sample+pt_cal].SetBinContent(i+1,j+1, 0.0)
            else:
                reso_dataMC[sample+pt_cal].SetBinContent(i+1,j+1, 100.0 * (reso_data - reso_MC)/reso_MC)

	cor_dir.cd()
	gStyle.SetPaintTextFormat("4.2f%%")
	mean_dataMC[sample+pt_cal].GetZaxis().SetRangeUser(-0.5,0.5)
	mean_dataMC[sample+pt_cal].SetXTitle("d0_muPos")
	mean_dataMC[sample+pt_cal].SetYTitle("d0_muNeg")
	mean_dataMC[sample+pt_cal].Write()
	reso_dataMC[sample+pt_cal].GetZaxis().SetRangeUser(-10,10)
	reso_dataMC[sample+pt_cal].SetXTitle("d0_muPos")
	reso_dataMC[sample+pt_cal].SetYTitle("d0_muNeg")
	reso_dataMC[sample+pt_cal].Write()	






#    ovl_dir.cd()
#    WriteOverlay(mean_plots, "mean", samples, pt_cals)
#    WriteOverlay(reso_plots, "reso", samples, pt_cals)
#
#    sum_dir.cd()
#    WriteSummary(mean_plots, "mean", nameX, samples, pt_cals, file_dir)
#    WriteSummary(reso_plots, "reso", nameX, samples, pt_cals, file_dir)
#







main()
