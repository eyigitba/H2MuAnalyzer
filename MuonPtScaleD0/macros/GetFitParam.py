#! /usr/bin/env python

###############################################
###              GetFitParam.py               ###
###                                         ###
###    Find peaks in dPt and dPhi plots     ###
###    (RECO - GEN) vs. d0 in ZJets MC.     ###
###                                         ###
###               Efe Yigitbasi             ###
###               16.09.2019                ###
###############################################

## Basic python includes for manipulating files
import sys
import os

## More python tools
import math
import numpy

## ROOT includes
import ROOT as R
import ROOT.RooFit as RF
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

YEAR = '2016'
# YEAR = '2017'
# YEAR = '2018'

INDIR = 'plots_fineBinned_fullD0Range__tanhFit_%s' % YEAR

INFILE = 'FindPeaks.root'

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Saving Fit parameters ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    out_file = open(INDIR + '/fit_param.txt','w') 

    # graph = {}


    for var in ['dRelPt2p0']: #, 'dPhi', 'dEta']:  ## Variables plotted vs. RECO d0 * charge
    # for var in ['dRelPt2']:  ## Variables plotted vs. RECO d0 * charge
        for corr in ['', 'Roch', 'KaMu', 'kinfit', 'KinRoch', 'KinKaMu']:
        # for corr in ['Roch']:

            ## Unique string identifying this canvas
            if (var == 'dPhi' or var == 'dEta'): p_str = var
            elif (corr == ''): p_str = '%s%s' % (var, corr)
            else: p_str = '%s_%s' % (var, corr)

            if ((var == 'dPhi' or var == 'dEta') and (corr == 'Roch' or corr == 'KaMu' or corr == 'kinfit' or corr == 'KinRoch' or corr == 'KinKaMu')): continue

            ## Unique string identifying this graph
            # c_str = c_str
            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf']: # |eta| binned
                c_str = eta + '_' + p_str 
                ## Canvas for the graphs of peak vs. d0
                # canv[c_str] = R.TCanvas(c_str, c_str, 1600, 1200)
                # iPt = -1
                for pt in ['inclusive']:
                # for pt in ['inclusive', 'pt_20_35', 'pt_35_42', 'pt_42_50', 'pt_50_inf']:
                    # iPt += 1
                    g_str = pt + '_' + c_str

                    print 'Getting data for %s' % g_str

                    ## Get the graph for g_str
                    graph = in_file.Get(g_str)

                    f = graph.GetFunction("f_%s" % g_str)

                    # print f.GetName()
                    # print graph.GetName()

                    out_file.write('f_%s : \n' % g_str)
                    out_file.write('Chi2 : ' + str(round(f.GetChisquare(),4)) + '\n')
                    out_file.write('NDF : ' + str(round(f.GetNDF(),4)) + '\n')
                    out_file.write('intercept : ' + str(round(f.GetParameter(0),4)) + '\n')
                    out_file.write('intercept_err : ' + str(round(f.GetParError(0),4)) + '\n')
                    out_file.write('amp : ' + str(round(f.GetParameter(1),4)) + '\n')
                    out_file.write('amp_err : ' + str(round(f.GetParError(1),4)) + '\n\n')
                    out_file.write('slope : ' + str(round(f.GetParameter(2),4)) + '\n')
                    out_file.write('slope_err : ' + str(round(f.GetParError(2),4)) + '\n\n')

                # End loop for pt
                out_file.write('\n\n')
            #End loop for eta
        ## End loop: for corr in ['PF', 'Roch']:
    ## End loop: for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:


    ## Write histograms to output root file
    # out_file.Write()

    ## Close output and input files
    print 'Finished, closing files ...'
    out_file.close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()