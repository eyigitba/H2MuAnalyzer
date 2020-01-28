#! /usr/bin/env python

###############################################
###            ANPlotMaking.py              ###
###                                         ###
###    Find peaks in dPt and dPhi plots     ###
###    (RECO - GEN) vs. d0 in ZJets MC.     ###
###                                         ###
###              Efe Yigitbasi              ###
###               24.10.2019                ###
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

# Bugfix DY samples
# INDIR = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots/plots_DY_d0BS_linear_%s' %YEAR
INDIR = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots/plots_DY_bugFix_%s' %YEAR

# ttbar samples
# INDIR_TT = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots_2D_ttbar_%s' %YEAR

# ttH samples
# INDIR_TT = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots_2D_ttH_%s' %YEAR

INFILE = 'FindPeaks.root'

OUTDIR = 'plots_DY_d0PV_AN_%s' % YEAR


D0_BINS = 19

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Inside FindPeaks.py ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    print '\nCreating output root file %s/FindPeaks.root' % OUTDIR
    if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

    ## Canvases to draw the graph on
    canv = {}

    # for var in ['dPt',  'dRelPt2p0', 'dRelPt1p0']: #, 'dPhi', 'dEta']:  ## Variables plotted vs. RECO d0 * charge
    for var in ['dRelPt2p0']: #, 'dPhi', 'dEta']:  ## Variables plotted vs. RECO d0 * charge
    # for var in ['dRelPt2']:  ## Variables plotted vs. RECO d0 * charge
        for corr in ['Roch']:
        # for corr in ['Roch']:

            ## Unique string identifying this canvas
            if (var == 'dPhi' or var == 'dEta'): p_str = var
            elif (corr == ''): p_str = '%s%s' % (var, corr)
            else: p_str = '%s_%s' % (var, corr)

            if ((var == 'dPhi' or var == 'dEta') and (corr == 'Roch' or corr == 'KaMu' or corr == 'kinfit' or corr == 'KinRoch' or corr == 'KinKaMu')): continue

            ## Unique string identifying this graph
            # c_str = c_str
            # for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf', 'eta_inc']: # |eta| binned
            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf']: # |eta| binned
                c_str = eta + '_' + p_str 
                ## Canvas for the graphs of peak vs. d0
                canv[c_str] = R.TCanvas(c_str, c_str, 1200, 1200)
                
                iPt = -1
                # for pt in ['pt_20_35']:
                # for pt in ['inclusive', 'pt_20_35', 'pt_35_42', 'pt_42_50', 'pt_50_inf']:
                for pt in ['inclusive']:
                # for pt in ['inclusive', 'nVtx_0_21', 'nVtx_22_26', 'nVtx_27_33', 'nVtx_34_inf']:
                    iPt += 1
                    g_str = pt + '_' + c_str

                    ## Create a graph with the peak values vs. d0
                    graph_DY = in_file.Get(g_str)


                    graph_DY.SetLineColor(1)
                    graph_DY.SetMarkerColor(1)
                    graph_DY.SetMarkerStyle(8)
                    graph_DY.SetMarkerSize(0.7)

                    graph_DY.SetTitle('')
                    graph_DY.GetYaxis().SetRangeUser(-15, 15)
                    graph_DY.GetYaxis().SetTitleOffset(1.1)
                    graph_DY.GetXaxis().SetTitleOffset(1.2)
                    # graph_DY.GetXaxis().SetLabelSize(0)
                    graph_DY.GetXaxis().SetTitle('d0_{PV} (cm) * charge')

                    graph_DY.Draw('AP')

                    if eta == 'eta_inc':
                        eta_str = 'inclusive\; \eta '
                    elif eta == 'eta_0_0p9':
                        eta_str = '0 < |\eta| < 0.9'
                    elif eta == 'eta_0p9_1p7':
                        eta_str = '0.9 < |\eta| < 1.7'
                    elif eta == 'eta_1p7_inf':
                        eta_str = '1.7 < |\eta|'


                    fit_DY = graph_DY.GetFunction('f_' + g_str)
                    # g_leg = R.TLegend(0.12,0.76,0.72,0.88)
                    # g_leg.AddEntry( graph_DY, 'DY : y = %.2f + %.2f \pm %.2f *d0' % (fit_DY.GetParameter(0), fit_DY.GetParameter(1), fit_DY.GetParError(1)), 'LAP' )
                    # g_leg.AddEntry( graph_DY, 'DY sample', 'LAP'  )
                    # g_leg.AddEntry( graph_tt, 'ttbar sample' , 'LAP' )
                    # g_leg.Draw('same')
                    g_latex = R.TLatex()
                    g_latex.SetTextAlign(10)
                    g_latex.SetTextSize(0.04)
                    g_latex.DrawLatex(-0.01, 12, '#font[42]{%s}' % eta_str)

                    cms_latex = R.TLatex()
                    cms_latex.SetTextAlign(11)
                    cms_latex.SetTextSize(0.03)
                    cms_latex.DrawLatexNDC(0.11, 0.93, '#scale[1.5]{CMS}')
                    cms_latex.DrawLatexNDC(0.21, 0.93,'#font[52]{#scale[1.2]{preliminary simulation}}')
                    cms_latex.DrawLatexNDC(0.82, 0.93,'#scale[1.3]{%s}' % YEAR)
                  
                # End loop for pt

                # canv[c_str].SaveAs('%s/%s.png' % (OUTDIR, c_str))
                canv[c_str].SaveAs('%s/%s.pdf' % (OUTDIR, c_str))
            #End loop for eta
        ## End loop: for corr in ['PF', 'Roch']:
    ## End loop: for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:


    ## Write histograms to output root file
    # out_file.Write()

    ## Close output and input files
    print 'Finished, closing files ...'
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
