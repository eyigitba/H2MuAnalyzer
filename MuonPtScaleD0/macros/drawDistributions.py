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
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_DY_beforeAfter_withMass/files/HADD/' % YEAR

# 2D ttbar samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ttbar/files/HADD/' % YEAR

# 2D ttH samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ttH_beforeAfter_withMass/files/HADD/' % YEAR

# 2D ggH samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_ggH_beforeAfter_withMass/files/HADD/' % YEAR

# 2D VBF samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_VBF_beforeAfter_withMass/files/HADD/' % YEAR


INFILE = 'histos_ZJets_hadd.root'
# INFILE = 'histos_ttbar_hadd.root'
# INFILE = 'histos_ttH_hadd.root'
# INFILE = 'histos_ggH_hadd.root'
# INFILE = 'histos_VBF_hadd.root'

OUTDIR = 'plots_2D_pT_d0_BS_%s' % YEAR

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
    # canv = {}
    canv = R.TCanvas('c', 'c', 1600, 1200)


    h = in_file.Get('h_inclusive_eta_inc_pt_Roch_vs_d0_BS')

    R.gPad.SetLogz()
    h.Scale(1/h.Integral())
    h.GetZaxis().SetRangeUser(0.00000001, 0.001)
    
    # if eta == 'eta_inc':
    #     eta_str = 'inclusive\; \eta '
    # elif eta == 'eta_0_0p9':
    #     eta_str = '0 < |\eta| < 0.9'
    # elif eta == 'eta_0p9_1p7':
    #     eta_str = '0.9 < |\eta| < 1.7'
    # elif eta == 'eta_1p7_inf':
    #     eta_str = '1.7 < |\eta|'

    if INFILE == 'histos_ZJets_hadd.root':
        sample_str = 'DY'
    elif INFILE == 'histos_ggH_hadd.root':
        sample_str = 'ggH'
    elif INFILE == 'histos_ttH_hadd.root':
        sample_str = 'ttH'
    elif INFILE == 'histos_VBF_hadd.root':
        sample_str = 'VBF'

    R.gStyle.SetOptStat(1101)

    h.SetNameTitle('%s pT vs d0_BS' % sample_str,'%s p_{T} vs d0_BS' % sample_str)

    h.Draw('colz')

    h.GetXaxis().SetTitle('d0_BS * muon charge')
    h.GetYaxis().SetTitle('muon p_{T}')


    canv.SaveAs('%s/png/%s.png' % (OUTDIR, sample_str))



    ## Close output and input files
    print 'Finished, closing files ...'
    # out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
