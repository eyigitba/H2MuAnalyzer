#! /usr/bin/env python

###############################################
###              BDTCat1D.py                ###
###                                         ###
###    Slice BDT into optimal ranges for    ###
###       S/sqrt(B) with uncertainties      ###
###                                         ###
###           Andrew Brinkerhoff            ###
###               16.11.2018                ###
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

LUMI     = 35.9
SIG_LOSS =  1.0  ## Percentage loss in significance tolerated from merging bins
MC_ERR   =  0.0  ## Scaling factor for MC error: 0 to ignore, 1 to use


## Main function executed by ./macros/BDTCat1D.py
def main():

    print '\n\n*** Inside BDTCat1D.py ***'

    ## Import root file with data and MC BDT distributions and trees
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    if not os.path.exists('plots'):
        print '\nCreating output folders plots/png and plots/pdf'
        os.makedirs('plots/png')
        os.makedirs('plots/pdf')

    print '\nCreating output root file plots/BDTCat1D.root'
    out_file = R.TFile('plots/BDTCat1D.root', 'recreate')


    ## Get signal BDT histogram
    h_sig_BDT = in_file.Get(FACTORY+'/Method_'+METHOD+'/'+METHOD+'/MVA_BDTG_UF_v1_S')
    h_bkg_BDT = in_file.Get(FACTORY+'/Method_'+METHOD+'/'+METHOD+'/MVA_BDTG_UF_v1_B')

    ## Normalize to expected number of ggH signal events
    ## 9.62 fb cross-section, 59% selection efficiency
    print '\nScaling signal histogram by %.1f' % (LUMI * 9.62 * 0.59 / h_sig_BDT.Integral())
    h_sig_BDT.Scale( LUMI * 9.62 * 0.59 / h_sig_BDT.Integral() )
    ## Normalize to expected number of total background events in FWHM
    ## 13000 events/GeV at mass = 125 in 2016 (35.9 fb-1), FWHM = 3.9 GeV
    print 'Scaling background histogram by %.3f' % (13000 * 3.9 * (LUMI / 35.9) / h_bkg_BDT.Integral())
    h_bkg_BDT.Scale( 13000 * 3.9 * (LUMI / 35.9) / h_bkg_BDT.Integral() )

    ## Make a histogram of the signal^2
    h_sig_BDT_sq = h_sig_BDT.Clone('h_sig_BDT_sq')
    h_sig_BDT_sq.Multiply(h_sig_BDT)

    ## Print significance for a 1-bin experiment
    print '\nSignificance treating the whole BDT range as one bin = %.3f' % ( h_sig_BDT.Integral() / math.sqrt(h_bkg_BDT.Integral()) )

    ## Make a histogram of the signal^2 / background
    h_BDT_sigma_sq = h_sig_BDT_sq.Clone('h_BDT_sigma_sq')
    h_BDT_sigma_sq.Divide(h_bkg_BDT)
    sigma_integral = math.sqrt(h_BDT_sigma_sq.Integral())

    print '\nNet significance estimated directly from the S^2 / B integral = %.3f' % sigma_integral

    ## Estimate the significance excluding uncertainties
    sigma_no_uncert = 0
    nBins = h_sig_BDT.GetNbinsX()
    for bin in range(1, nBins+1):
        num = h_sig_BDT.GetBinContent(bin)
        den = h_bkg_BDT.GetBinContent(bin)
        sigma_no_uncert += ( pow(num, 2) / den )
    sigma_no_uncert = math.sqrt(sigma_no_uncert)
    
    print 'Net significance with original binning and no statistical uncertainties = %.3f\n' % sigma_no_uncert

    ## Estimate the significance including uncertainties
    sigma_with_uncert = 0
    nBins = h_sig_BDT.GetNbinsX()
    for bin in range(1, nBins+1):
        num = h_sig_BDT.GetBinContent(bin)
        den = h_bkg_BDT.GetBinContent(bin)
        err = pow(h_bkg_BDT.GetBinError(bin), 2)
        sigma_with_uncert += ( pow(num, 2) / (den + MC_ERR*math.sqrt(err)) )
        print '  * Bin %d has %.2f +/- %.2f (%.2f%%) background events' % (bin, den, math.sqrt(err), 100*math.sqrt(err) / den)
    sigma_with_uncert = math.sqrt(sigma_with_uncert)
    
    print 'Net significance with original binning and statistical uncertainties = %.3f' % sigma_with_uncert

    
    ##################################################################
    ###  Merge bins until there is a net XX% loss in significance  ###
    ##################################################################

    orig_binning = []
    for bin in range(1, nBins+1):
        orig_binning.append([bin])
    old_binning = orig_binning
    print '\n***  Original binning = '+str(orig_binning)+'  ***'

    sigma_merged_sq = pow(sigma_with_uncert, 2)  ## Significance after merging

    while (math.sqrt(sigma_merged_sq) > (1 - 0.01*SIG_LOSS)*sigma_with_uncert and len(old_binning) > 1):

        min_loss_sq = -999       ## Minimum loss in significance^2 for a given iteration
        merged_bins = [-99, -99] ## Indices of sets of bins merged in a given iteration 

        for iBins in range(len(old_binning) - 1):
            num_A  = 0
            num_B  = 0
            num_AB = 0
            den_A  = 0
            den_B  = 0
            den_AB = 0
            err_A  = 0
            err_B  = 0
            err_AB = 0
            for bin in old_binning[iBins]:
                num_A += h_sig_BDT.GetBinContent(bin)
                den_A += h_bkg_BDT.GetBinContent(bin)
                err_A += pow(h_bkg_BDT.GetBinError(bin), 2)
            for bin in old_binning[iBins+1]:
                num_B += h_sig_BDT.GetBinContent(bin)
                den_B += h_bkg_BDT.GetBinContent(bin)
                err_B += pow(h_bkg_BDT.GetBinError(bin), 2)

            sig_A_sq  = pow(num_A, 2) / (den_A + MC_ERR*math.sqrt(err_A))
            sig_B_sq  = pow(num_B, 2) / (den_B + MC_ERR*math.sqrt(err_B))
            sig_AB_sq = pow(num_A + num_B, 2) / (den_A + den_B + MC_ERR*math.sqrt(err_A + err_B))

            # print '\nBin range '+str(old_binning[iBins])+' has significance %.3f' % math.sqrt(sig_A_sq)
            # print 'Bin range '+str(old_binning[iBins+1])+' has significance %.3f' % math.sqrt(sig_B_sq)
            # print 'Combined range has significance %.3f, for a loss of %.4f' % ( math.sqrt(sig_AB_sq), math.sqrt(sig_AB_sq) - math.sqrt(sig_A_sq + sig_B_sq) )

            if ( (sig_AB_sq - sig_A_sq - sig_B_sq) > min_loss_sq ):
                min_loss_sq    = sig_AB_sq - sig_A_sq - sig_B_sq
                merged_bins[0] = iBins
                merged_bins[1] = iBins+1

        ## End loop: for iBins in range(len(old_binning) - 1):

        sigma_merged_sq += min_loss_sq
        print '  - Merged bins '+str(old_binning[merged_bins[0]])+' with '+str(old_binning[merged_bins[1]])
        print '    * '+str(len(old_binning) - 1)+' bins with a net significance = %.3f (%.3f%% loss)' % ( math.sqrt(sigma_merged_sq), 100*(1.0 - math.sqrt(sigma_merged_sq) / sigma_with_uncert) )

        new_binning = []
        for iBins in range(len(old_binning)):
            if   iBins != merged_bins[0] and iBins != merged_bins[1]:
                new_binning.append(old_binning[iBins])
            elif iBins == merged_bins[0]:
                new_binning.append(old_binning[iBins] + old_binning[iBins+1])
            elif iBins == merged_bins[1]:
                continue
            else: print '\nBIZZARE ERROR! iBins = %d, merged_bins[0] = %d, merged_bins[1] = %d' % (iBins, merged_bins[0], merged_bins[1])

        # print '***  New binning = '+str(new_binning)+'  ***'
        old_binning = new_binning

    ## End conditional: while (sqrt(sigma_merged_sq) > (1 - 0.01*SIG_LOSS)*sigma_with_uncert):
        
    print '\nReduced to '+str(len(new_binning))+' bins, with new binning = '+str(new_binning)
    print 'Net merged significance = %.3f (compare to original = %.3f)' % ( math.sqrt(sigma_merged_sq), sigma_with_uncert )

    sigma_merged_sq_check = 0
    for iBins in range(len(new_binning)):
        num = 0
        den = 0
        err = 0
        for bin in new_binning[iBins]:
            num += h_sig_BDT.GetBinContent(bin)
            den += h_bkg_BDT.GetBinContent(bin)
            err += pow(h_bkg_BDT.GetBinError(bin), 2)
        sigma_merged_sq_check += ( pow(num, 2) / (den + MC_ERR*math.sqrt(err)) )
        if (iBins > 0): print 'BDT cut #%d = %.3f' % (iBins, h_bkg_BDT.GetBinLowEdge(new_binning[iBins][0]))
    ## End loop: for iBins in range(len(new_binning)):
    print 'Explicit check on merged significance = %.3f' % math.sqrt(sigma_merged_sq_check)

    ## Write histograms to output root file
    out_file.cd()


    ## Close output and input files
    print '\nFinished, closing files ...'
    out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
