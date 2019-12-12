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
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ggH_beforeAfter_withMass/files/HADD/' % YEAR

# 2D VBF samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_VBF_beforeAfter_withMass/files/HADD/' % YEAR


# INFILE = 'histos_ZJets_hadd.root'
# INFILE = 'histos_ttbar_hadd.root'
# INFILE = 'histos_ttH_hadd.root'
INFILE = 'histos_ggH_hadd.root'
# INFILE = 'histos_VBF_hadd.root'

OUTDIR = 'plots_ggH_pT_resolution_%s' % YEAR

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
    hh_base = R.TH1D()

    for var in ['dPt', 'dRelPt2p0', 'dRelPt1p0']:
        # for corr in ['', 'Roch', 'KaMu', 'kinfit', 'KinRoch', 'KinKaMu']:
        for corr in ['corr', 'corr_new', 'corr_ttH', 'corr_lowPt']:
        # for corr in ['corr', 'corr_lowPt']:

            ## Unique string identifying this canvas
            p_str = '%s_%s' % (var, corr)
            p_str_base = '%s_Roch' % var

            ## Unique string identifying this graph

            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf', 'eta_inc']: # |eta| binned
                c_str = eta + '_' + p_str 
                c_str_base = eta + '_' + p_str_base 
                ## Canvas for the graphs of peak vs. d0
                for pt in ['inclusive']:
                    g_str = pt + '_' + c_str
                    g_str_base = pt + '_' + c_str_base
                    for region in ['tail_left', 'core', 'tail_right', 'tails', 'total']:
                        gg_str = g_str + '_' + region
                        gg_str_base = g_str_base + '_' + region
                        canv[gg_str] = R.TCanvas(gg_str, gg_str, 1600, 1200)

                        h = in_file.Get('h_%s_vs_d0_BS' % (g_str)) # dPt2p0 vs d0_BS
                        h_base = in_file.Get('h_%s_vs_d0_BS' % (g_str_base)) # dPt2p0 vs d0_BS

                        h1 = h.Clone()
                        h1_base = h_base.Clone()

                        hh.Reset()
                        hh_base.Reset()

                        if region == 'tail_left':
                          hh = h1.ProjectionY('', 1, 120, 'e').Clone()
                          hh_base = h1_base.ProjectionY('', 1, 120, 'e').Clone()
                        elif region == 'core':
                          hh = h1.ProjectionY('', 121, 280, 'e').Clone()  
                          hh_base = h1_base.ProjectionY('', 121, 280, 'e').Clone()  
                        elif region == 'tail_right':
                          hh = h1.ProjectionY('', 281, 400, 'e').Clone()
                          hh_base = h1_base.ProjectionY('', 281, 400, 'e').Clone()
                        elif region == 'tails':
                          hh = h1.ProjectionY('', 1, 120, 'e').Clone()
                          hh.Add( h1.ProjectionY('', 281, 400, 'e'))
                          hh_base = h1_base.ProjectionY('', 1, 120, 'e').Clone()
                          hh_base.Add( h1_base.ProjectionY('', 281, 400, 'e')) 
                        elif region == 'total':
                          hh = h1.ProjectionY('', 1, 400, 'e').Clone()
                          hh_base = h1_base.ProjectionY('', 1, 400, 'e').Clone()

                        # hh.SetNameTitle(gg_str, gg_str)
                        # hh_base.SetNameTitle(gg_str_base, gg_str_base)
                        # hh.Rebin(10)
                        if var == 'dPt':
                            hh_base.GetXaxis().SetRangeUser(-12,12)
                        if var == 'dRelPt1p0':
                            hh_base.GetXaxis().SetRangeUser(-40,40)
                        if var == 'dRelPt2p0':
                            hh_base.GetXaxis().SetRangeUser(-50,50)
                        
                        hh_base.GetYaxis().SetRangeUser(0,hh.GetMaximum()*1.05)

                        if var == 'dRelPt2p0':
                            var_str = '10000 * (p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}'
                        if var == 'dRelPt1p0':
                            var_str = '100 * (p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{GEN})'
                        if var == 'dPt':
                            var_str = 'p_{T}^{RECO} - p_{T}^{GEN}'

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

                        if region == 'tail_left':
                            reg_str = 'd0_BS < -0.04'
                        if region == 'tail_right':
                            reg_str = 'd0_BS > 0.04'
                        if region == 'tails':
                            reg_str = '|d0_BS| > 0.04'
                        if region == 'core':
                            reg_str = '|d0_BS| < 0.04'                
                        if region == 'total':
                            reg_str = 'all d0_BS'
                        
                        if eta == 'eta_inc':
                            eta_str = 'inclusive\; \eta '
                        elif eta == 'eta_0_0p9':
                            eta_str = '0 < |\eta| < 0.9'
                        elif eta == 'eta_0p9_1p7':
                            eta_str = '0.9 < |\eta| < 1.7'
                        elif eta == 'eta_1p7_inf':
                            eta_str = '1.7 < |\eta|'

                        if INFILE == 'histos_ZJets_hadd.root':
                            sample_str = 'DY'
                        elif INFILE == 'histos_ggH_hadd.root':
                            sample_str = 'ggH'
                        elif INFILE == 'histos_ttH_hadd.root':
                            sample_str = 'ttH'
                        elif INFILE == 'histos_VBF_hadd.root':
                            sample_str = 'VBF'

                        # graph[g_str].SetTitle('%s  %s  vs.  d0' % (var_str, corr_str))
                        # graph[g_str].GetXaxis().SetTitle('RECO muon d0 (cm) * charge')
                        # graph[g_str].GetYaxis().SetTitle('%s %s' % (var_str, corr_str))

                        # hh.Write()
                        # hh.Draw()

                        hh.SetNameTitle('%s' % corr_str,'%s' % corr_str)
                        hh_base.SetNameTitle('Roch','Roch' )


                        hh_base.SetLineColor(1)
                        hh_base.SetMarkerColor(1)
                        # h_base.SetMarkerStyle(8)
                        # h_base.SetMarkerSize(0.7)
                        hh_base.SetTitle('%s %s sample %s vs. Roch for %s' % (YEAR, sample_str, corr_str, reg_str))
                        hh_base.GetXaxis().SetTitle(var_str)
                        hh_base.GetXaxis().SetTitleOffset(1)
                        hh_base.GetYaxis().SetTitle('Events')

                        hh.SetLineColor(2)
                        hh.SetMarkerColor(2)



                        fit = R.TF1('fit_GeoFit', 'gaus', hh.GetBinCenter(hh.GetMaximumBin())-0.7, hh.GetBinCenter(hh.GetMaximumBin())+0.7 )
                        fit_base = R.TF1('fit_Roch', 'gaus', hh_base.GetBinCenter(hh_base.GetMaximumBin())-0.7, hh_base.GetBinCenter(hh_base.GetMaximumBin())+0.7 )
                        
                        fit_base.SetLineColor(1)
                        fit.SetLineColor(2)
                        hh.Fit('fit_GeoFit', 'R')
                        hh_base.Fit('fit_Roch', 'R')


                        R.gStyle.SetOptStat(0000)
                        # R.gStyle.SetOptFit(01)
                        hh_base.Draw()
                        R.gPad.Update()



                        # bin1 = hh.FindFirstBinAbove(hh.GetMaximum()/2);
                        # bin2 = hh.FindLastBinAbove(hh.GetMaximum()/2);
                        # fwhm = hh.GetBinCenter(bin2) - hh.GetBinCenter(bin1);

                        # bin1_base = hh_base.FindFirstBinAbove(hh_base.GetMaximum()/2);
                        # bin2_base = hh_base.FindLastBinAbove(hh_base.GetMaximum()/2);
                        # fwhm_base = hh_base.GetBinCenter(bin2_base) - hh_base.GetBinCenter(bin1_base);

                        # print fwhm_base
                        # print fwhm

                        hh.Draw('sames')
                        R.gPad.Update()
                        # st = hh.FindObject('stats')
                        # st.SetY1NDC(0.68)
                        # st.SetY2NDC(0.81)
                        # print hh.GetRMS()
                        # print hh_base.GetRMS()

                        latex = R.TLatex()
                        latex.SetTextAlign(12)
                        latex.SetTextSize(0.04)
                        latex.DrawLatex(-10, hh.GetMaximum()*0.7, eta_str)

                        corr_percent = (1 - fit.GetParameter(2)/fit_base.GetParameter(2))*100.0
                        corr_percent_str = 'imprv: %.1f'  % corr_percent
                        corr_percent_str += '%'


                        latex_corr = R.TLatex()
                        latex_corr.SetTextAlign(12)
                        latex_corr.SetTextSize(0.03)
                        latex_corr.DrawLatex(-10, hh.GetMaximum()*0.6, corr_percent_str)

                        leg = R.TLegend(0.12,0.70,0.44,0.88)
                        leg.AddEntry( fit_base, 'Roch : Mean = %.3f // Sigma = %.3f' % (fit_base.GetParameter(1), fit_base.GetParameter(2)) )
                        leg.AddEntry( fit, '%s : Mean = %.3f // Sigma = %.3f' % (corr_str, fit.GetParameter(1), fit.GetParameter(2)) )
                        leg.Draw('same')
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
    # out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
