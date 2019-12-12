#! /usr/bin/env python

###############################################
###          CompareDYvsTTbar.py            ###
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

# YEAR = '2016'
# YEAR = '2017'
YEAR = '2018'

# Bugfix DY samples
INDIR = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots_2D_DY_%s' %YEAR

# ttbar samples
# INDIR_TT = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots_2D_ttbar_%s' %YEAR

# ttH samples
INDIR_TT = '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/plots_2D_ttH_%s' %YEAR

INFILE = 'FindPeaks.root'

OUTDIR = 'plots_DY_ttH_final_%s' % YEAR


D0_BINS = 19

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Inside FindPeaks.py ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)
    in_file_ttH = R.TFile(INDIR_TT+'/'+INFILE)

    print '\nCreating output root file %s/FindPeaks.root' % OUTDIR
    if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

    ## Canvases to draw the graph on
    canv = {}
    pad1 = {}
    pad2 = {}
    h_1 = R.TH1F('h_1', 'h_1', 22, -0.012, 0.012)
    h_2 = R.TH1F('h_2', 'h_2', 22, -0.012, 0.012)

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
            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf', 'eta_inc']: # |eta| binned
                c_str = eta + '_' + p_str 
                ## Canvas for the graphs of peak vs. d0
                canv[c_str] = R.TCanvas(c_str, c_str, 1600, 1200)
                pad1[c_str] = R.TPad(c_str+'_1', c_str+'_1', 0, 0.3, 1, 1.0)
                pad2[c_str] = R.TPad(c_str+'_2', c_str+'_2', 0, 0.05, 1, 0.3)

                pad1[c_str].SetBottomMargin(0)
                pad1[c_str].Draw()
                pad2[c_str].SetTopMargin(0);
                pad2[c_str].SetBottomMargin(0.2);
                pad2[c_str].Draw();
                iPt = -1
                # for pt in ['pt_20_35']:
                # for pt in ['inclusive', 'pt_20_35', 'pt_35_42', 'pt_42_50', 'pt_50_inf']:
                for pt in ['inclusive']:
                # for pt in ['inclusive', 'nVtx_0_21', 'nVtx_22_26', 'nVtx_27_33', 'nVtx_34_inf']:
                    iPt += 1
                    g_str = pt + '_' + c_str

                    ## Create a graph with the peak values vs. d0
                    graph_DY = in_file.Get(g_str)
                    graph_tt = in_file_ttH.Get(g_str)


                    graph_DY.SetLineColor(1)
                    graph_DY.SetMarkerColor(1)
                    graph_DY.SetMarkerStyle(8)
                    graph_DY.SetMarkerSize(0.7)

                    graph_tt.SetLineColor(2)
                    graph_tt.SetMarkerColor(2)
                    graph_tt.SetMarkerStyle(8)
                    graph_tt.SetMarkerSize(0.7)
                    graph_tt.GetFunction('f_'+g_str).SetLineColor(2)

                    pad1[c_str].cd()
                    graph_DY.GetYaxis().SetRangeUser(-15, 15)

                    graph_DY.Draw('AP')
                    graph_tt.Draw('Psame')

                    if eta == 'eta_inc':
                        eta_str = 'inclusive\; \eta '
                    elif eta == 'eta_0_0p9':
                        eta_str = '0 < |\eta| < 0.9'
                    elif eta == 'eta_0p9_1p7':
                        eta_str = '0.9 < |\eta| < 1.7'
                    elif eta == 'eta_1p7_inf':
                        eta_str = '1.7 < |\eta|'


                    fit_DY = graph_DY.GetFunction('f_' + g_str)
                    fit_tt = graph_tt.GetFunction('f_' + g_str)
                    g_leg = R.TLegend(0.12,0.76,0.72,0.88)
                    g_leg.AddEntry( graph_DY, 'DY : y = %.2f + %.2f \pm %.2f *d0' % (fit_DY.GetParameter(0), fit_DY.GetParameter(1), fit_DY.GetParError(1)), 'LAP' )
                    g_leg.AddEntry( graph_tt, 'ttbar : y = %.2f + %.2f \pm %.2f *d0' % (fit_tt.GetParameter(0), fit_tt.GetParameter(1), fit_tt.GetParError(1)), 'LAP' )
                    # g_leg.AddEntry( graph_DY, 'DY sample', 'LAP'  )
                    # g_leg.AddEntry( graph_tt, 'ttbar sample' , 'LAP' )
                    g_leg.Draw('same')
                    g_latex = R.TLatex()
                    g_latex.SetTextAlign(12)
                    g_latex.SetTextSize(0.04)
                    g_latex.DrawLatex(-0.01, 8, eta_str)
                  
                # End loop for pt

                # h_ratio = R.TH1F('h_ratio', 'h_ratio', 40, -0.01, 0.01)
                y_1 = graph_DY.GetY()
                yerr_1 = graph_DY.GetEY()
                y_2 = graph_tt.GetY()
                yerr_2 = graph_tt.GetEY()
                for i in range(D0_BINS):
                    # data['ratio']['x'][i] = data['nVtx_0_21_'+c_str]['x'][i]
                    # data['ratio']['xerr'][i] = data['nVtx_0_21_'+c_str]['xerr'][i]
                    # data['ratio']['y'][i] = data['nVtx_34_inf_'+c_str]['y'][i]/data['nVtx_0_21_'+c_str]['y'][i]
                    # data['ratio']['yerr'][i] = ((data['nVtx_34_inf_'+c_str]['yerr'][i]/data['nVtx_34_inf_'+c_str]['y'][i])+(data['nVtx_0_21_'+c_str]['yerr'][i]/data['nVtx_0_21_'+c_str]['y'][i]))*data['ratio']['y'][i]
                    h_1.SetBinContent(i+3, y_1[i])
                    h_1.SetBinError(i+3, yerr_1[i])
                    h_2.SetBinContent(i+3, y_2[i])
                    h_2.SetBinError(i+3, yerr_2[i])
                
                # pad1[c_str].cd()
                # h_1.Draw('ep')
                # h_2.Draw('ep same')

                pad2[c_str].cd()
                pad2[c_str].SetGridy()

                h_ratio = h_1.Clone()
                h_ratio.SetTitle('')
                h_ratio.SetLineColor(1)
                h_ratio.Sumw2()
                h_ratio.SetMinimum(0)
                h_ratio.SetMaximum(2)
                h_ratio.SetStats(0)
                h_ratio.SetMarkerStyle(20)
                h_ratio.Divide(h_2)

                h_ratio.GetYaxis().SetTitle('Ratio')
                h_ratio.GetYaxis().SetNdivisions(505)
                h_ratio.GetYaxis().SetTitleSize(20)
                h_ratio.GetYaxis().SetTitleFont(43)
                h_ratio.GetYaxis().SetTitleOffset(1.55)
                h_ratio.GetYaxis().SetLabelFont(43)
                h_ratio.GetYaxis().SetLabelSize(15)

                h_ratio.GetYaxis().SetTitle('d0_BS * muon charge')
                h_ratio.GetXaxis().SetTitleSize(20)
                h_ratio.GetXaxis().SetTitleFont(43)
                h_ratio.GetXaxis().SetTitleOffset(4.)
                h_ratio.GetXaxis().SetLabelFont(43)
                h_ratio.GetXaxis().SetLabelSize(15)
                # h_ratio.SetMarkerStyle(21)
                # r = R.TGraphErrors(len(data['ratio']['x']), data['ratio']['x'], data['ratio']['y'], data['ratio']['xerr'], data['ratio']['yerr'])
                # r.SetLineColor(1)
                # r.GetYaxis().SetRangeUser(0.5, 1.5)
                # r.Draw('AP')
                h_ratio.Draw('ep')
                # canv[c_str].SaveAs('%s/png/%s.png' % (OUTDIR, c_str))

                canv[c_str].SaveAs('%s/%s.png' % (OUTDIR, c_str))
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
