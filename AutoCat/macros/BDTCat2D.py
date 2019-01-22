#! /usr/bin/env python

###############################################
###              BDTCat2D.py                ###
###                                         ###
###      Take optimal ranges of BDT and     ###
###       test S/sqrt(B) with eta cuts      ###
###                                         ###
###           Andrew Brinkerhoff            ###
###               17.11.2018                ###
###############################################

## Basic python includes for manipulating files
import sys
import os

## Other python includes
import math

## ROOT includes
import ROOT as R
import ROOT.RooFit as RF


INDIR   = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2016/H2MuAnalyzer/TrainMVA/output'
INFILE  = '2018_11_16_TMVA_BDT_2016_massRes_ggH_vs_Z_AMC_MG.root'
# FACTORY = 'f_Opt_v1_etas_all_sig_all_bkg_ge0j'
FACTORY = 'f_Opt_v1_res_wgt_etas_all_sig_all_bkg_ge0j'
METHOD  = 'BDTG_UF_v1'

## MAX_EVT = -1  ## Must always run over all events - may be ordered sigal then background
PRT_EVT = 100000

LUMI     = 35.9
SIG_LOSS =    1  ## Percentage loss in significance tolerated from merging bins
MC_ERR   =  1.0  ## Scaling factor for MC error: 0 to ignore, 1 to use

MASS_MIN =  115  ## Minimum dimuon mass considered for signal peak
MASS_MAX =  135  ## Maximum dimuon mass considered for signal peak
NBINS    =   20  ## Number of mass bins

## Optimized BDT bins for f_Opt_v1_etas_all_sig_all_bkg_ge0j
# BDT_CUTS = [-1.0, -0.264, 0.032, 0.254, 0.438, 0.512, 0.660, 1.0]  ## 7 bins (assuming no statistical uncertainties)
# BDT_CUTS = [-1.0, -0.264, 0.069, 0.254,        0.512,        1.0]  ## 5 bins (assuming 50% of statistical uncertainties)
# BDT_CUTS = [-1.0, -0.264,        0.254,        0.512,        1.0]  ## 4 bins (assuming 100% of statistical uncertainties)

## Optimized BDT cuts for f_Opt_v1_res_wgt_etas_all_sig_all_bkg_ge0j
# BDT_CUTS = [-1.0, -0.137, 0.219, 0.378, 0.497, 0.576, 0.695, 1.0] ## 7 bins (assuming no statistical uncertainties)
BDT_CUTS = [-1.0, -0.137, 0.219, 0.338,        0.576,        1.0] ## 5 bins (assuming 50% of statistical uncertainties)
# BDT_CUTS = [-1.0, -0.137,        0.338,        0.576,        1.0] ## 4 bins (assuming 100% of statistical uncertainties)


ETA_CUTS = [0.0, 0.9, 1.9, 2.4]
# ETA_CUTS = [0.0, 2.4]

## Main function executed by ./macros/BDTCat2D.py
def main():

    print '\n\n*** Inside BDTCat2D.py ***'

    ## Import root file with data and MC BDT distributions and trees
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    if not os.path.exists('plots'):
        print '\nCreating output folders plots/png and plots/pdf'
        os.makedirs('plots/png')
        os.makedirs('plots/pdf')

    print '\nCreating output root file plots/BDTCat2D.root'
    out_file = R.TFile('plots/BDTCat2D.root', 'recreate')

    ## Mass distribution histograms
    h_sig = {}
    h_bkg = {}
    for i in range(len(BDT_CUTS) - 1):
        h_sig[i] = {}
        h_bkg[i] = {}
        for j in range(len(ETA_CUTS) - 1):
            h_sig[i][j] = R.TH1D('h_sig_%d_%d' % (i, j), '%.3f < BDT < %.3f, %.1f < |#eta| < %.1f' % (BDT_CUTS[i], BDT_CUTS[i+1], ETA_CUTS[j], ETA_CUTS[j+1]), NBINS, MASS_MIN, MASS_MAX)
            h_bkg[i][j] = R.TH1D('h_bkg_%d_%d' % (i, j), '%.3f < BDT < %.3f, %.1f < |#eta| < %.1f' % (BDT_CUTS[i], BDT_CUTS[i+1], ETA_CUTS[j], ETA_CUTS[j+1]), NBINS, MASS_MIN, MASS_MAX)


    ## Get TestTree
    tr = in_file.Get(FACTORY+'/TestTree')

    ## Number of signal and background events processed
    nSig    = 0
    nBkg    = 0
    nSigWgt = 0
    nBkgWgt = 0

    for iEvt in range(tr.GetEntries()):

        ## if iEvt > MAX_EVT and MAX_EVT > 0: break
        if iEvt % PRT_EVT is 0: print 'Event #', iEvt

        tr.GetEntry(iEvt)

        isSig = (tr.samp_ID == -1)  ## Mass 125 gluon-fusion Higgs
        isBkg = (tr.samp_ID > 0)    ## Any background sample

        if not (isSig or isBkg): continue

        ETA  = max(tr.mu1_abs_eta, tr.mu2_abs_eta)
        BDT  = tr.GetBranch(METHOD).GetLeaf(METHOD).GetValue()
        MASS = tr.dimu_mass

        if isSig:
            WGT      = tr.samp_wgt
            nSig    += 1
            nSigWgt += WGT
        if isBkg:
            WGT      = (tr.samp_wgt / 3)  ## Using 3 independent Drell-Yan samples
            nBkg    += 1
            nBkgWgt += WGT

        for i in range(len(BDT_CUTS) - 1):
            for j in range(len(ETA_CUTS) - 1):

                if (BDT < BDT_CUTS[i] or BDT > BDT_CUTS[i+1]): continue
                if (ETA < ETA_CUTS[j] or ETA > ETA_CUTS[j+1]): continue

                if isSig: h_sig[i][j].Fill(MASS, WGT)
                if isBkg: h_bkg[i][j].Fill(MASS, WGT)

            ## End loop: for j in range(len(ETA_CUTS) - 1):
        ## End loop: for i in range(len(BDT_CUTS) - 1):

    ## End loop: for iEvt in range(tr.GetEntries()):

    # print '\nNumber of signal events (weighted) = %d (%.1f)' % (nSig, nSigWgt)
    # print 'Number of background events (weighted) = %d (%.1f)' % (nBkg, nBkgWgt)

    ## Loop over histograms and extract significance
    net_sigma_sq_no_uncert = 0
    net_sigma_sq_uncert    = 0
    for i in range(len(BDT_CUTS) - 1):
        for j in range(len(ETA_CUTS) - 1):
            
            # print '\nLooking at BDT bin %d, |eta| bin %d' % (i, j)
            
            ## Estimate the significance of this BDT-eta category
            sigma_sq_no_uncert = 0
            sigma_sq_uncert    = 0
            for bin in range(1, NBINS+1):
                num = h_sig[i][j].GetBinContent(bin)
                den = h_bkg[i][j].GetBinContent(bin)
                err = pow(h_bkg[i][j].GetBinError(bin), 2)
                if den == 0:
                    # print '\n\nIn BDT cat %d, ETA cat %d, bin %d, num = %.4f, den = %.4f' % (i, j, bin, num, den)
                    den = 999999
                sigma_sq_no_uncert += ( pow(num, 2) / (den) )
                sigma_sq_uncert    += ( pow(num, 2) / (den + MC_ERR*math.sqrt(err)) )
                
            net_sigma_sq_no_uncert += sigma_sq_no_uncert
            net_sigma_sq_uncert    += sigma_sq_uncert

            # print 'Significance without statistical uncertainties = %.3f' % math.sqrt(sigma_sq_no_uncert)
            # print 'Significance with statistical uncertainties = %.3f (%.1f%% worse)' % ( math.sqrt(sigma_sq_uncert), 100 - 100*math.sqrt(sigma_sq_uncert/sigma_sq_no_uncert) )

        ## End loop: for j in range(len(ETA_CUTS) - 1):
    ## End loop: for i in range(len(BDT_CUTS) - 1):
    
    print '\n\nNet significance without statistical uncertainties = %.3f' % math.sqrt(net_sigma_sq_no_uncert)
    print 'Net significance with statistical uncertainties = %.3f (%.1f%% worse)' % ( math.sqrt(net_sigma_sq_uncert), 100 - 100*math.sqrt(net_sigma_sq_uncert/net_sigma_sq_no_uncert) )


    ## Write histograms to output root file
    out_file.cd()

    for i in range(len(BDT_CUTS) - 1):
        for j in range(len(ETA_CUTS) - 1):
            h_sig[i][j].Write()
            h_bkg[i][j].Write()

    ## Close output and input files
    print '\nFinished, closing files ...'
    out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
