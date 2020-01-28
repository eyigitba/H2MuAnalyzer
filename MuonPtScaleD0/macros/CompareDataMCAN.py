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
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn ZH_125 WH_neg_125 WH_pos_125 ttH_125

YEAR = '2016'
# YEAR = '2017'
# YEAR = '2018'

# Bugfix DY samples
# INDIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/%s/Histograms/MassCal_KinRoch_approx/GeoBSRoch_2D_muP_d0_rebin_muN_d0_rebin/' %YEAR
# INDIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/%s/Histograms/MassCal_KinRoch_approx/GeoBSRoch_2D_muP_d0_rebin_muP_phi/' %YEAR
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_2D_DY_eta_rebin/' % YEAR
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_DY_PVvsBS/files/HADD/' % YEAR


# ttbar samples
# INDIR_TT = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots_2D_ttbar_%s' %YEAR

# ttH samples
# INDIR_TT = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots_2D_ttH_%s' %YEAR

# INFILE = 'mass_cal_plots.root'
INFILE = 'histos_ZJets_hadd.root'

OUTDIR = 'plots_DY_PV_BS_AN_%s' % YEAR


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

    ccanv = R.TCanvas('WH', 'WH', 1800, 1200)
    ccanv.SetLeftMargin(0.15)
    ccanv.SetRightMargin(0.15)

    plot_PV = in_file.Get('h_inclusive_eta_inc_dRelPt2p0_Roch_vs_d0_PV').Clone()
    plot_BS = in_file.Get('h_inclusive_eta_inc_dRelPt2p0_Roch_vs_d0_BS').Clone()

    # mean_roch = float(str(dimu_mass_Roch.GetListOfFunctions().At(1).GetLine(5).GetTitle()).split()[2])
    # mean_geofit = float(str(dimu_mass_geofit.GetListOfFunctions().At(1).GetLine(5).GetTitle()).split()[2])

    # sigma_roch = float(str(dimu_mass_Roch.GetListOfFunctions().At(1).GetLine(6).GetTitle()).split()[2])
    # sigma_geofit = float(str(dimu_mass_geofit.GetListOfFunctions().At(1).GetLine(6).GetTitle()).split()[2])

    # plot_pull.SetStats(False)
    # plot_data.SetStats(False)

    plot_PV.SetStats(False)
    plot_BS.SetStats(False)
    # dimu_mass_Roch.GetFunction('voigtian_exp').SetLineColor(1)
    # dimu_mass_Roch.SetLineColor(1)
    # dimu_mass_geofit.SetLineColor(2)
    # dimu_mass_geofit.SetMarkerColor(2)
    # dimu_mass_Roch.remove('DSCB_paramBox')
    # dimu_mass_geofit.remove('DSCB_paramBox')
    R.gStyle.SetOptStat(0)

    R.gPad.Update()

    # graph_DY.SetLineColor(1)
    # graph_DY.SetMarkerColor(1)
    # graph_DY.SetMarkerStyle(8)
    # graph_DY.SetMarkerSize(0.7)

    # plot_pull.SetTitle('')
    # plot_data.SetTitle('')

    plot_PV.SetTitle('')
    plot_BS.SetTitle('')
    # dimu_mass_Roch.GetYaxis().SetRangeUser(0, dimu_mass_geofit.GetMaximum()*1.2)
    # dimu_mass_Roch.GetYaxis().SetTitleOffset(2.1)
    # dimu_mass_Roch.GetXaxis().SetTitleOffset(1.2)
    # dimu_mass_Roch.GetXaxis().SetLabelSize(0)
    # plot_pull.GetYaxis().SetTitle('d0_{BS,\mu^{-}} [cm]')
    # plot_pull.GetYaxis().SetTitle('\phi_{\mu^{+}}')
    # plot_pull.GetYaxis().SetTitle('\eta_{\mu^{+}}')
    # plot_pull.GetYaxis().SetTitleSize(0.05)
    # plot_pull.GetXaxis().SetTitle('d0_{BS,\mu^{+}} [cm]')
    # plot_pull.GetZaxis().SetTitle('(m_{\mu\mu}^{Roch} - m_{\mu\mu}^{Gen}) / m_{\mu\mu}^{Gen}')

    plot_PV.GetYaxis().SetTitle('10000* (p_{T}^{Roch}-p_{T}^{Gen})/(p_{T}^{Gen})^{2}')
    # plot_PV.GetYaxis().SetTitleSize(0.05)
    plot_PV.GetXaxis().SetTitle('d0_{PV} * muon charge [cm]')
    plot_PV.GetZaxis().SetTitle('Number of muons')
    # plot_data.GetYaxis().SetTitle('d0_{BS,\mu^{-}} [cm]')
    # plot_data.GetYaxis().SetTitle('\phi_{\mu^{+}}')
    # plot_data.GetYaxis().SetTitle('\eta_{\mu^{+}}')
    # plot_data.GetYaxis().SetTitleSize(0.05)
    # plot_data.GetXaxis().SetTitle('d0_{BS,\mu^{+}} [cm]')
    # plot_data.GetZaxis().SetTitle('(m_{\mu\mu}^{data} - m_{\mu\mu}^{MC}) / m_{\mu\mu}^{MC}')

    plot_BS.GetYaxis().SetTitle('10000* (p_{T}^{Roch}-p_{T}^{Gen})/(p_{T}^{Gen})^{2}')
    # plot_PV.GetYaxis().SetTitleSize(0.05)
    plot_BS.GetXaxis().SetTitle('d0_{BS} * muon charge [cm]')
    plot_BS.GetZaxis().SetTitle('Number of muons')
    # R.gStyle.SetPaintTextFormat(".1f %%");


    # plot_pull.Draw('colz text')
    # plot_data.Draw('colz text')
    R.gPad.SetLogz()
    # plot_PV.Draw('colz')
    plot_BS.Draw('colz')
    # dimu_mass_geofit.Draw('lpsame')


    # if eta == 'eta_inc':
    #     eta_str = 'inclusive\; \eta '
    # elif eta == 'eta_0_0p9':
    #     eta_str = '0 < |\eta| < 0.9'
    # elif eta == 'eta_0p9_1p7':
    #     eta_str = '0.9 < |\eta| < 1.7'
    # elif eta == 'eta_1p7_inf':
    #     eta_str = '1.7 < |\eta|'


    # fit_DY = dimu_mass_Roch.GetFunction('f_' + g_str)
    # g_leg = R.TLegend(0.12,0.76,0.72,0.88)
    # g_leg.AddEntry( dimu_mass_Roch, 'DY : y = %.2f + %.2f \pm %.2f *d0' % (fit_DY.GetParameter(0), fit_DY.GetParameter(1), fit_DY.GetParError(1)), 'LAP' )
    # g_leg.AddEntry( dimu_mass_Roch, 'DY sample', 'LAP'  )
    # g_leg.AddEntry( graph_tt, 'ttbar sample' , 'LAP' )
    # g_leg.Draw('same')
    # g_latex = R.TLatex()
    # g_latex.SetTextAlign(10)
    # g_latex.SetTextSize(0.04)
    # g_latex.DrawLatex(-0.01, 12, '#font[42]{%s}' % eta_str)

    cms_latex = R.TLatex()
    cms_latex.SetTextAlign(11)
    cms_latex.SetTextSize(0.03)
    cms_latex.DrawLatexNDC(0.15, 0.91, '#scale[1.5]{CMS}')
    cms_latex.DrawLatexNDC(0.22, 0.91,'#font[52]{#scale[1.2]{preliminary simulation}}')
    cms_latex.DrawLatexNDC(0.80, 0.91,'#scale[1.3]{%s}' % YEAR)
    # cms_latex.DrawLatexNDC(0.71, 0.91,'#font[42]{35.9 fb^{-1} (13 TeV)}') #35.9 41.5 59.7

    # cms_latex.DrawLatexNDC(0.18, 0.85, '#font[42]{#bf{Z+jets data Roch: }}')
    # cms_latex.DrawLatexNDC(0.40, 0.85, '#font[42]{mean: %.2f}, ' % mean_roch)
    # cms_latex.DrawLatexNDC(0.56, 0.85, '#font[42]{sigma: %.3f}' % sigma_roch)

    # cms_latex.DrawLatexNDC(0.18, 0.81, '#color[2]{#font[42]{#bf{Z+jets data GeoFit: }}}')
    # cms_latex.DrawLatexNDC(0.40, 0.81, '#color[2]{#font[42]{mean: %.2f}, }' % mean_geofit)
    # cms_latex.DrawLatexNDC(0.56, 0.81, '#color[2]{#font[42]{sigma: %.3f}}' % sigma_geofit)
    # cms_latex.DrawLatexNDC(0.18, 0.77, '#color[2]{#font[42]{#bf{%.1f %% improvement}}}' % (100.0*(1-sigma_geofit/sigma_roch)))

    ccanv.SaveAs('%s/%s_BS.pdf' % (OUTDIR, OUTDIR))
    
    # ccanv.SaveAs('%s/%s.pdf' % (OUTDIR, OUTDIR))


    ## Write histograms to output root file
    # out_file.Write()

    ## Close output and input files
    print 'Finished, closing files ...'
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
