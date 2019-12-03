#! /usr/bin/env python

###############################################
###              FindPeaks.py               ###
###                                         ###
###       Find peaks in dRelPt2p0 vs. d0    ###
###                                         ###
###              Efe Yigitbasi              ###
###               10.10.2019                ###
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
# YEAR = '2017'
YEAR = '2018'

# 2D DY samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_DY_beforeAfter/files/HADD/' % YEAR

# 2D ttbar samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ttbar/files/HADD/' % YEAR

# 2D ttH samples
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ttH_beforeAfter/files/HADD/' % YEAR

# 2D ggH samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ggH/files/HADD/' % YEAR


# INFILE = 'histos_ZJets_hadd.root'
# INFILE = 'histos_ttbar_hadd.root'
INFILE = 'histos_ttH_hadd.root'

OUTDIR = 'plots_ttH_pT_resolution_%s' % YEAR

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
    out_file = R.TFile('%s/FindPeaks.root' % OUTDIR, 'recreate')
    out_file.cd()

    if not os.path.exists('%s/png' % OUTDIR):
        print '\nCreating output folders %s/png'  % OUTDIR
        os.makedirs('%s/png' % OUTDIR)

    ## Input histograms and canvases to draw on
    hist = {}
    ## Graphs of delta(RECO, GEN) vs. RECO muon d0 * charge
    graph = {}
    ## Data to fill the graph
    data = {}
    ## Fits to the data
    fit = {}
    fit_left = {}
    fit_right = {}
    ## Canvases to draw the graph on
    canv = {}
    h_proj_tail_left = R.TH1D()
    h_proj_core = R.TH1D()
    h_proj_tail_right = R.TH1D()
    hh = R.TH1D()

    for var in ['dPt']:
        # for corr in ['', 'Roch', 'KaMu', 'kinfit', 'KinRoch', 'KinKaMu']:
        for corr in ['Roch','corr']:

            ## Unique string identifying this canvas
            if (corr == ''): p_str = '%s%s' % (var, corr)
            else: p_str = '%s_%s' % (var, corr)

            ## Unique string identifying this graph

            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf', 'eta_inc']: # |eta| binned
                c_str = eta + '_' + p_str 
                ## Canvas for the graphs of peak vs. d0
                for pt in ['inclusive']:
                    g_str = pt + '_' + c_str
                    for region in ['tail_left', 'core', 'tail_right']:
                        gg_str = g_str + '_' + region
                        canv[gg_str] = R.TCanvas(gg_str, gg_str, 1600, 1200)

                        h = in_file.Get('h_%s_vs_d0_BS' % (g_str)) # dPt2p0 vs d0_BS
                        h1 = h.Clone()

                        hh.Reset()

                        if region == 'tail_left':
                          hh = h1.ProjectionY('', 1, 140, 'e')
                        elif region == 'core':
                          hh = h1.ProjectionY('', 141, 260, 'e')  
                        elif region == 'tail_right':
                          hh = h1.ProjectionY('', 261, 400, 'e')

                        hh.SetName(gg_str)
                        # hh.Rebin(10)
                        hh.GetXaxis().SetRangeUser(-20,20)

                        # var_str = '10000 * (p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}'

                        if corr == '':
                            corr_str = '(PF)'
                        if corr == 'Roch':
                            corr_str = '(Rochester)'
                        if corr == 'KaMu':
                            corr_str = '(Kalman)'
                        if corr == 'kinfit':
                            corr_str = '(Kinfit)'                
                        if corr == 'KinRoch':
                            corr_str = '(Rochester + Kinfit)'
                        if corr == 'KinKaMu':
                            corr_str = '(Kalman + Kinfit)'

                        if eta == 'eta_inc':
                            eta_str = 'inclusive\; \eta '
                        elif eta == 'eta_0_0p9':
                            eta_str = '0 < |\eta| < 0.9'
                        elif eta == 'eta_0p9_1p7':
                            eta_str = '0.9 < |\eta| < 1.7'
                        elif eta == 'eta_1p7_inf':
                            eta_str = '1.7 < |\eta|'

                        # graph[g_str].SetTitle('%s  %s  vs.  d0' % (var_str, corr_str))
                        # graph[g_str].GetXaxis().SetTitle('RECO muon d0 (cm) * charge')
                        # graph[g_str].GetYaxis().SetTitle('%s %s' % (var_str, corr_str))

                        hh.Write()
                        hh.Draw()
                        canv[gg_str].SaveAs('%s/png/%s.png' % (OUTDIR, gg_str))

                    ## Draw the graph onto the canvas

                    # graph[g_str].SetLineColor(1)
                    # graph[g_str].SetMarkerColor(1)
                    # graph[g_str].SetMarkerStyle(8)
                    # graph[g_str].SetMarkerSize(0.7)
                    # fit[g_str].SetLineColor(1)
                    # graph[g_str].Fit('f_%s' % g_str,"R")
                    # print 'Writing: %s' % g_str
                    # out_file.cd()
                    # graph[g_str].Write()

                    # g_canv = R.TCanvas(g_str, g_str, 1600, 1200)
                    # g_canv.cd()
                    # graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                    # graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                    # graph[g_str].Draw('AP')
                    # g_leg = R.TLegend(0.12,0.76,0.72,0.88)

                    # # g_leg.AddEntry( fit[g_str], '\LARGE{y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)}' % (fit[g_str].GetParameter(0), fit[g_str].GetParameter(3), fit[g_str].GetParameter(1), fit[g_str].GetParameter(2)) )
                    # g_leg.AddEntry( fit[g_str], '\LARGE{y = %.2f + %.2f *d0}' % (fit[g_str].GetParameter(0), fit[g_str].GetParameter(1)) )
                    # g_leg.Draw('same')
                    # g_latex = R.TLatex()
                    # g_latex.SetTextAlign(12)
                    # g_latex.SetTextSize(0.04)
                    # g_latex.DrawLatex(-0.01, 8, eta_str)
                    # g_canv.SaveAs('%s/png/%s.png' % (OUTDIR, g_str))
                    
                    # canv[c_str].cd()
                    # graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                    # graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                    # graph[g_str].Draw('AP')

                    # latex = R.TLatex()
                    # latex.SetTextAlign(12)
                    # latex.SetTextSize(0.04)
                    # latex.DrawLatex(-0.01, 8, eta_str)
                    
                    


                # End loop for pt
            #End loop for eta
        ## End loop: for corr in ['PF', 'Roch']:
    ## End loop: for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:

    ## Close output and input files
    print 'Finished, closing files ...'
    out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
