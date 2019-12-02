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

# YEAR = '2016'
YEAR = '2017'
# YEAR = '2018'

INDIR = 'plots_DY_d0BS_linear_%s' % YEAR

INFILE = 'FindPeaks.root'

OUTFILE = 'fit_param_d0BS_%s.txt' % YEAR

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Saving Fit parameters ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    out_file = open(INDIR + '/' + OUTFILE,'w') 

    for var in ['dRelPt2p0']:
        for corr in ['Roch']:

            ## Unique string identifying this canvas
            if (var == 'dPhi' or var == 'dEta'): p_str = var
            elif (corr == ''): p_str = '%s%s' % (var, corr)
            else: p_str = '%s_%s' % (var, corr)

            if ((var == 'dPhi' or var == 'dEta') and (corr == 'Roch' or corr == 'KaMu' or corr == 'kinfit' or corr == 'KinRoch' or corr == 'KinKaMu')): continue

            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf']: # |eta| binned
                c_str = eta + '_' + p_str 
                for pt in ['inclusive']:
                    g_str = pt + '_' + c_str

                    print 'Getting data for %s' % g_str

                    ## Get the graph for g_str
                    graph = in_file.Get(g_str)

                    f = graph.GetFunction("f_%s" % g_str)

                    out_file.write('f_%s : \n' % g_str)
                    out_file.write('chi2/ndf : ' + str(round(f.GetChisquare()/f.GetNDF(),2)) + ' \n')                    
                    out_file.write('intercept : ' + str(round(f.GetParameter(0),4)) + ' \n')
                    out_file.write('slope : ' + str(round(f.GetParameter(1),3)) + '\n')
                    

                # End loop for pt
                out_file.write('\n\n')
            #End loop for eta
        ## End loop: for corr in ['PF', 'Roch']:
    ## End loop: for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:
    
    ## Close output and input files
    print 'Finished, closing files ...'
    out_file.close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()