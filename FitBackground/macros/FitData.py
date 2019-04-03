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
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

## FitFunction class from python/FitFunctions.py
sys.path.insert(0, '%s/python' % os.getcwd())
import FitFunctions as FF

# SOURCE = 'FEWZ'
# INDIR  = '/afs/cern.ch/user/a/abrinke1/HiggsToMuMu/CMSSW_9_4_9/src/H2MuAnalyzer/FitBackground/FEWZ/hists'
# INFILE = 'NNLO_Bourilkov_2017.root'
# CATS = ['full']
# # CATS = ['full', 'cc', 'ncnc', '1jet', '2jet']

SOURCE = 'acarnes_data'
# SOURCE = 'acarnes_MC'
INDIR  = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/acarnes/h2mumu/rfiles'
INFILE = 'validate_UNBLINDED_dimu_mass_Roch_90_200_categories3_tree_bdt_res_mass_n16_mbg25_unc0_np1p0000_sf1_sd0_sb0_PoissonSignificance_36814_dyAMC-J_minpt20_b-4_sig-xlumi0.root'
SUBDIR = 'net_histos'
# CATS = ['c15']
CATS = ['c3', 'c4', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c14', 'c15']
# CATS = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c12', 'c13', 'c14', 'c15']

# SOURCE = 'abrinke1'
# INDIR  = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/Limits/input_hists/Moriond17_Feb08/AWB_Mar21_test_v2'
# INFILE = 'MergedData.root'
# HISTNAME = 'DiMuonMass'
# CATS = ['NoCats', 'VBFTight']


MASS_MIN = 110
MASS_MAX = 150
BLIND    = [120, 130]  ## Set to "[]" for unblinded, "[120, 130]" to blind signal region

# MODELS = [['poly', 5]]
# MODELS = [['Bern', 5]]
MODELS = [['poly', 5], ['expo', 3], ['Bern', 5], ['RedBWZ', 1]]


## Main function executed by ./macros/FitData.py
def main():

    print '\n\n*** Inside FitData.py ***'

    ## Import root file with data and MC invariant mass distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    if not os.path.exists('plots/png/%s' % SOURCE):
        print '\nCreating output folders plots/png and plots/pdf'
        os.makedirs('plots/png/%s' % SOURCE)
        os.makedirs('plots/png/%s/diffs' % SOURCE)
        os.makedirs('plots/pdf/%s' % SOURCE)
        os.makedirs('plots/pdf/%s/diffs' % SOURCE)

    # print '\nCreating output root file plots/FitData.root'
    # out_file = R.TFile('plots/FitData.root', 'recreate')

    # ## Write histograms to output root file
    # out_file.cd()

    ## Input histograms and canvases to draw on
    h_data = {}
    c_data = {}
    c_fits = {}

    ## Fits and models
    fits = {}

    ## Loop over categories
    for cat in CATS:
        print 'Processing category %s...' % cat
        c_data[cat] = R.TCanvas('c_%s_data' % cat, 'c_%s_data' % cat, 800, 600)

        ## Get histogram from input root file
        if SOURCE == 'FEWZ':          h_data[cat] = in_file.Get(cat+'_36fb')
        if SOURCE == 'acarnes_data':  h_data[cat] = in_file.Get(SUBDIR+'/'+cat+'_Net_Data')
        if SOURCE == 'acarnes_MC':    h_data[cat] = in_file.Get(SUBDIR+'/'+cat+'_Net_Bkg')
        if SOURCE == 'abrinke1':      h_data[cat] = in_file.Get(cat+'/'+HISTNAME)


        ##################################
        ###  Fit each model in MODELS  ###
        ##################################
        for mod in MODELS:

            ## Initialize a function to fit the data
            ## The FitFunction class is defined in python/FitFunctions.py
            fit_str = '%s%d_%s' % (mod[0], mod[1], cat)
            fits[fit_str] = FF.FitFunction(fit_str, h_data[cat], mod[0], mod[1], [MASS_MIN, MASS_MAX], BLIND)

            ## Fit the function to the data
            FF.DoFit(fits[fit_str])

            ## Use a shorter name for drawing
            fit = fits[fit_str]

            ## Plot data into a frame
            c_data[cat].cd()
            fra = fit.var.frame()
            fit.dat.plotOn(fra)

            ## Plot fit onto the same frame
            fit.model.plotOn(fra, RF.LineColor(R.kBlue), RF.Range('FULL'))

            # for amp_var in fit.amp_vars:
            #     print 'Value of amp_var in %s = %.9f +/- %.9f' % (fit_str, amp_var.getValV(), amp_var.getError())
            # for param in fit.params[0]:
            #     print 'Value of param in %s = %.9f +/- %.9f' % (fit_str, param.getValV(), param.getError())


            ## Plot the components pdfs onto the same frame
            for i in range(mod[1]):
                fit.model.plotOn(fra, RF.Components(fit.arg_sets[i]), RF.LineStyle(R.kDashed), RF.LineColor(i+2), RF.Range('FULL'))

                ## Draw frame on canvas
                fra.SetMinimum( max(fit.hist.GetBinContent(fit.x_bins[2]) / 3, 0.1) )
                fra.SetMaximum(     fit.hist.GetBinContent(fit.x_bins[1]) * 3)
                fra.Draw()
                
            c_data[cat].SetLogy()
                
            # c_data[cat].Write()
            c_data[cat].SaveAs('plots/png/%s/c_%s_data.png' % (SOURCE, fit_str))
            # c_data[cat].SaveAs('plots/pdf/%s/c_%s_data.pdf' % (SOURCE, fit_str))


        ## End loop: for mod in MODELS:


        #############################################
        ###  Compare each fit to every other fit  ###
        #############################################

        ## Loop over all fits in category
        for key1 in fits.keys():
            if not cat in key1: continue

            max_diff = 0  ## Maximum difference between fit and other fits (usually a positive number)
            min_diff = 0  ## Minimum difference between fit and other fits (usually a negative number)

            c_fits[key1] = R.TCanvas('c_fits_vs_%s' % key1, 'c_fits_vs_%s' % key1, 800, 600)
            c_fits[key1].cd()

            h_diff = {}

            ## Subract "key1" fit histogram from data fit histogram
            keyD = 'data_%s' % cat
            h_diff[keyD] = R.TH1F('h_diff_data_%s' % cat, 'h_diff_data_%s' % cat, fits[key1].x_bins[0], MASS_MIN, MASS_MAX)

            off = fits[key1].x_bins[1] - 1  ## Offset for data hist, whose first bin may not be 1
            for i in range(1, fits[key1].x_bins[0]+1):
                if ( fits[key1].hist.GetBinContent(i+off) > 0 ):
                    h_diff[keyD].SetBinContent(i, fits[key1].hist.GetBinContent(i+off) - fits[key1].fit_hist.GetBinContent(i))
                else:
                    h_diff[keyD].SetBinContent(i, fits[key1].hist.GetBinContent(0))
                h_diff[keyD].SetBinError(i, fits[key1].hist.GetBinError(i+off))

            max_diff = max(max_diff, h_diff[keyD].GetMaximum())
            min_diff = min(min_diff, h_diff[keyD].GetMinimum())

            ## Loop over all other fits in category
            for key2 in fits.keys():
                if not cat in key2: continue
                if key2 == key1: continue

                ## Subract "key1" fit histogram from "key2" fit histogram
                h_diff[key2] = fits[key2].fit_hist.Clone('h_diff_%s' % key2)
                h_diff[key2].Add(fits[key1].fit_hist, -1)

                max_diff = max(max_diff, h_diff[key2].GetMaximum())
                min_diff = min(min_diff, h_diff[key2].GetMinimum())

            ## End loop: for key2 in fits.keys():

            ## Loop over fit differences and plot
            h_diff[keyD].SetMarkerStyle(10)
            # h_diff[keyD].SetMarkerSize(10)
            h_diff[keyD].SetMarkerColor(R.kBlack)
            h_diff[keyD].SetLineColor(R.kBlack)
            h_diff[keyD].SetLineWidth(1)
            h_diff[keyD].GetYaxis().SetRangeUser(min_diff*1.5, max_diff*1.5)
            h_diff[keyD].Draw('PE')

            iColor = 1
            for key2 in h_diff.keys():
                if key2 == keyD: continue
                iColor += 1
                h_diff[key2].SetLineColor(iColor)
                h_diff[key2].SetLineWidth(2)
                h_diff[key2].Draw('histsame')

            ## End loop: for key2 in h_diff.keys():

            if len(h_diff.keys()) > 0:
                c_fits[key1].SaveAs('plots/png/%s/diffs/c_fits_vs_%s.png' % (SOURCE, key1))
                # c_fits[key1].SaveAs('plots/pdf/%s/diffs/c_fits_vs_%s.pdf' % (SOURCE, key1))
 
        ## End loop: for key1 in fits.keys():

    ## End loop: for cat in CATS


    ## Close output and input files
    print 'Finished, closing files ...'
    # out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
