################################################################
###    make mass calibration plots from mass histograms      ###
################################################################

## Basic python includes for manipulating files
import os
import sys

sys.path.insert(0, '%s/lib' % os.getcwd())
from ROOT import *
from MassCal_Helper import FitVoigtian, GetColor, WriteOverlay, WriteSummary

gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

LABEL = 'MassCal_KinRoch_approx/d0_diff_KinRoch_by_eta'
nameX = "d0_diff"
CAT   = "NONE_NONE"

def main():
    file_dir = PLOT_DIR+"/"+LABEL
    out_file = TFile( file_dir + "/mass_cal_plots" + ".root", "RECREATE")
    in_file  = TFile.Open( file_dir + "/all_samples" + ".root", "READ")

    sum_dir  = out_file.mkdir("summary")
    ovl_dir  = out_file.mkdir("overlay")
    cal_dir  = out_file.mkdir("individual_cal_plots")
    fit_dir  = out_file.mkdir("fits_specifics")

    samples = ["ZJets_MG_1", "data"]
#    pt_cals = ["PF", "Kinfit", "good_Kinfit", "Kin_vs_d0kin"]
    pt_cals = ["Kin_vs_d0kin", "Kin_vs_d0kin_BB", "Kin_vs_d0kin_BE", "Kin_vs_d0kin_EE"]
#    pt_cals = ["Kin_vs_d0kin_d0PV_N50_N15", "Kin_vs_d0kin_d0PV_N15_N05", "Kin_vs_d0kin_d0PV_P05_P15", "Kin_vs_d0kin_d0PV_P15_P50"]

    summary_info = in_file.Get( pt_cals[0] + "/summary_" + samples[0] + "_" + CAT + "_" + pt_cals[0]).Clone()
    binningX = summary_info.GetXaxis().GetXbins()  ## return type is TArrayD

    mean_plots = {}
    reso_plots = {}
    pull_plots = {}
    for pt_cal in pt_cals:
      for sample in samples:
	mean_plots[sample+pt_cal] = TGraphErrors()
	reso_plots[sample+pt_cal] = TGraphErrors()
	mean_plots[sample+pt_cal].SetName("Mean_" + sample + "_" + pt_cal)
	reso_plots[sample+pt_cal].SetName("Reso_" + sample + "_" + pt_cal)
	for i in range(binningX.GetSize()-1):
	  x_low  = binningX.GetAt(i)
	  x_high = binningX.GetAt(i+1)
	  x_val = ( x_high + x_low ) / 2.0
	  x_err = ( x_high - x_low ) / 2.0

	  print "looking at, %s,   %s,   %f to %f   " %(sample, pt_cal, x_low, x_high)
	  hist_name = pt_cal + "/" + sample + "_" + pt_cal + "_%s_%8.4f_to_%8.4f" %( nameX, x_low, x_high)
	  hist_name = hist_name.replace(' ','').replace("-","m").replace(".","p")
	  print hist_name
	  mass_hist = in_file.Get(hist_name).Clone()
	  mean_val, mean_err, reso_val, reso_err = FitVoigtian(mass_hist)
	  mean_plots[sample+pt_cal].SetPoint(i, x_val, mean_val)
	  mean_plots[sample+pt_cal].SetPointError(i, x_err, mean_err)

	  reso_plots[sample+pt_cal].SetPoint(i, x_val, reso_val)
          reso_plots[sample+pt_cal].SetPointError(i, x_err, reso_err)

	  mean_plots[sample+pt_cal].SetLineColor( GetColor(sample, pt_cal) )
	  mean_plots[sample+pt_cal].SetLineWidth(2)
 	  reso_plots[sample+pt_cal].SetLineColor( GetColor(sample, pt_cal) )
          reso_plots[sample+pt_cal].SetLineWidth(2)

	  fit_dir.cd()
          mass_hist.Write()

	cal_dir.cd()
	mean_plots[sample+pt_cal].Write()
	reso_plots[sample+pt_cal].Write()

    ovl_dir.cd()
    WriteOverlay(mean_plots, "mean", samples, pt_cals)
    WriteOverlay(reso_plots, "reso", samples, pt_cals)

    sum_dir.cd()
    WriteSummary(mean_plots, "mean", nameX, samples, pt_cals, file_dir)
    WriteSummary(reso_plots, "reso", nameX, samples, pt_cals, file_dir)








main()
