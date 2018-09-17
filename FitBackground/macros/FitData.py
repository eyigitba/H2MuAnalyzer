#! /usr/bin/env python

###############################################
###               FitData.py                ###
###                                         ###
###    Fit background data with different   ###
###    functions and compare their shapes   ###
###                                         ###
###           Andrew Brinkerhoff            ###
###               21.08.2018                ###
###############################################

## Basic python includes for manipulating files
import sys
import os

## ROOT includes
import ROOT as R
import ROOT.RooFit as RF

## FitFunction class from python/FitFunctions.py
sys.path.insert(0, '%s/python' % os.getcwd())
import FitFunctions as FF


INDIR  = '/afs/cern.ch/work/a/acarnes/public/h2mumu/rfiles'
INFILE = 'validate_UNBLINDED_dimu_mass_Roch_90_200_categories3_tree_bdt_res_mass_n16_mbg25_unc0_np1p0000_sf1_sd0_sb0_PoissonSignificance_36814_dyAMC-J_minpt20_b-4_sig-xlumi0.root'
SUBDIR = 'net_histos'

CATS = ['c8', 'c15']
# CATS = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c12', 'c13', 'c14', 'c15']

MASS_MIN = 110
MASS_MAX = 150
BLIND    = True

NUM_EXPO = 2

## Main function executed by ./batch/FitData.py
def main():

    print '\n\n*** Inside FitDatas.py ***'

    ## Import root file with data and MC invariant mass distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    if not os.path.exists('plots'):
        print '\nCreating output folders plots/png and plots/pdf'
        os.makedirs('plots/png')
        os.makedirs('plots/pdf')

    # print '\nCreating output root file plots/FitData.root'
    # out_file = R.TFile('plots/FitData.root', 'recreate')

    # ## Write histograms to output root file
    # out_file.cd()

    ## Input histograms and canvases to draw on
    h_data = {}
    c_data = {}

    ## Fits and models
    fits = {}

    ## Loop over categories
    for cat in CATS:
        print 'Processing category %s...' % cat
        c_data[cat] = R.TCanvas('c_%s_data' % cat, 'c_%s_data' % cat, 800, 600)
        c_data[cat].cd()

        ## Get histogram from input root file
        h_data = in_file.Get(SUBDIR+'/'+cat+'_Net_Data')

        ## Initialize an exponential function to fit the data
        ## The FitFunction class is defined in python/FitFunctions.py
        expo_str = 'exp%d_%s' % (NUM_EXPO, cat)
        fits[expo_str] = FF.FitFunction(expo_str, h_data, 'expo', NUM_EXPO, [MASS_MIN, MASS_MAX], [120, 130])

        ## Fit the function to the data
        FF.DoFit(fits[expo_str])

        ## Use a shorter name for drawing
        fit = fits[expo_str]

        ## Plot data into a frame
        fra = fit.var.frame()
        fit.dat.plotOn(fra)
        
        fit.model.plotOn(fra, RF.LineColor(R.kBlue), RF.Range('FULL'))
        for i in range(NUM_EXPO):
            fit.model.plotOn(fra, RF.Components(fit.arg_sets[i]), RF.LineStyle(R.kDashed), RF.LineColor(i+2), RF.Range('FULL'))

        ## Draw frame on canvas
        fra.SetMinimum(fit.hist.GetBinContent(300) / 3)
        fra.SetMaximum(fit.hist.GetBinContent(100) * 3)
        fra.Draw()

        c_data[cat].SetLogy()
        
        c_data[cat].Write()
        c_data[cat].SaveAs('plots/png/c_%s_data.png' % cat)
        c_data[cat].SaveAs('plots/pdf/c_%s_data.pdf' % cat)

    ## End loop: for cat in CATS


    ## Close output and input files
    print 'Finished, closing files ...'
    # out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
