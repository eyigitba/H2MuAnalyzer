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
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_DY_beforeAfter_withMass/files/HADD/' % YEAR

# 2D ttbar samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ttbar/files/HADD/' % YEAR

# 2D ttH samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ttH_beforeAfter_withMass/files/HADD/' % YEAR

# 2D ggH samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ggH_beforeAfter_withMass/files/HADD/' % YEAR

# 2D ttbar samples
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ttbar_beforeAfter_withMass/files/HADD/' % YEAR


# INFILE = 'histos_ZJets_hadd.root'
# INFILE = 'histos_ttbar_hadd.root'
# INFILE = 'histos_ttH_hadd.root'
# INFILE = 'histos_ggH_hadd.root'
# INFILE = 'histos_VBF_hadd.root'
INFILE = 'histos_ttbar_hadd.root'

OUTDIR = 'plots_ttbar_mass_resolution_%s' % YEAR

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Inside FindPeaks.py ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    print '\nCreating output root file %s' % OUTDIR
    if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)
    # out_file = R.TFile('%s/FindPeaks.root' % OUTDIR, 'recreate')
    # out_file.cd()

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

    for var in ['dimu_mass']:
        # for corr in ['', 'Roch', 'KaMu', 'kinfit', 'KinRoch', 'KinKaMu']:
        for corr in ['corr', 'corr_new', 'corr_ttH', 'corr_core', 'corr_tail', 'corr_lowPt']: #, 'corr_coreOnly', 'corr_tailOnly']:

            ## Unique string identifying this canvas
            p_str = '%s_%s' % (var, corr)
            p_str_base = '%s_Roch' % var

            ## Unique string identifying this graph

            for eta in ['eta_inc']: # |eta| binned
                c_str = eta + '_' + p_str 
                c_str_base = eta + '_' + p_str_base 
                ## Canvas for the graphs of peak vs. d0
                for pt in ['inclusive']:
                    g_str = pt + '_' + c_str
                    g_str_base = pt + '_' + c_str_base
                    canv[g_str] = R.TCanvas(g_str, g_str, 1600, 1200)

                    h = in_file.Get('h_%s' % (g_str)) # dPt2p0 vs d0_BS
                    h_base = in_file.Get('h_%s' % (g_str_base)) # dPt2p0 vs d0_BS
                    if h_base.GetNbinsX() == 4000:
                      h_base.Rebin(4)
                    h.Rebin(4)
                    if INFILE == 'histos_ZJets_hadd.root' :
                        h.GetXaxis().SetRangeUser(86.2,96.2)
                        h_base.GetXaxis().SetRangeUser(86.2,96.2)
                    elif INFILE == 'histos_ttbar_hadd.root':
                        h.GetXaxis().SetRangeUser(50,150)
                        h_base.GetXaxis().SetRangeUser(50,150)
                    else :
                        h.GetXaxis().SetRangeUser(121.2,128.4)
                        h_base.GetXaxis().SetRangeUser(121.2,128.4)

                    h_base.GetYaxis().SetRangeUser(0,h_base.GetMaximum()*1.2)

                    

                    # var_str = '10000 * (p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}'

                    if corr == 'corr':
                        corr_str = 'GeoFit'
                    if corr == 'corr_new':
                        corr_str = 'GeoFit_new'
                    if corr == 'corr_ttH':
                        corr_str = 'GeoFit_ttH'
                    if corr == 'corr_core':
                        corr_str = 'GeoFit_core'                
                    if corr == 'corr_tail':
                        corr_str = 'GeoFit_tail'
                    if corr == 'corr_lowPt':
                        corr_str = 'GeoFit_lowPt'
                    if corr == 'corr_coreOnly':
                        corr_str = 'GeoFit_coreOnly'
                    if corr == 'corr_tailOnly':
                        corr_str = 'GeoFit_tailOnly'

                    if INFILE == 'histos_ZJets_hadd.root':
                        sample_str = 'DY'
                    elif INFILE == 'histos_ggH_hadd.root':
                        sample_str = 'ggH'
                    elif INFILE == 'histos_ttH_hadd.root':
                        sample_str = 'ttH'
                    elif INFILE == 'histos_VBF_hadd.root':
                        sample_str = 'VBF'
                    elif INFILE == 'histos_ttbar_hadd.root':
                        sample_str = 'ttbar'

                    h.SetName('%s' % corr_str)
                    h_base.SetName('Baseline')


                    h_base.SetLineColor(1)
                    h_base.SetMarkerColor(1)
                    # h_base.SetMarkerStyle(8)
                    # h_base.SetMarkerSize(0.7)
                    h_base.SetTitle('%s %s sample %s vs. Baseline' % (YEAR, sample_str, corr_str))
                    h_base.GetXaxis().SetTitle('dimuon mass')
                    h_base.GetYaxis().SetTitle('Events')

                    h.SetLineColor(2)
                    h.SetMarkerColor(2)

                    # fit = R.TF1('fit_GeoFit', 'voigt', )
                    R.gStyle.SetOptStat(1101)

                    h_base.Draw()
                  
                    h.Draw('sames')
                    R.gPad.Update()
                    st = h.FindObject('stats')
                    st.SetY1NDC(0.68)
                    st.SetY2NDC(0.81)
                    # R.gPad.SetLogy()
                    canv[g_str].SaveAs('%s/png/%s.png' % (OUTDIR, g_str))

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
    # out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
