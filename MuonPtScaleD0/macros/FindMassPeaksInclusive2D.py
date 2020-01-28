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
sys.path.insert(0, '/afs/cern.ch/user/e/eyigitba/Hmm/CMSSW_9_4_13/src/H2MuAnalyzer/MuonPtScaleD0/macros')
import ROOT as R
import ROOT.RooFit as RF
from MassCal_Helper import FitVoigtian, GetColor, WriteOverlay, WriteSummary

R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

YEAR = '2016'
# YEAR = '2017'
# YEAR = '2018'

# 2D DY samples
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_dimu_mass_2D_DY_withData/files/HADD/' % YEAR

INFILE_DATA = 'histos_ZJets_data_hadd.root'
INFILE_MC = 'histos_ZJets_MC_hadd.root'

OUTDIR = 'plots_2D_DY_massPeaks_withData_%s' % YEAR

D0_BINS = 19

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Inside FindPeaks.py ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE_DATA)
    in_file_data = R.TFile(INDIR+'/'+INFILE_DATA)
    in_file_mc = R.TFile(INDIR+'/'+INFILE_MC)

    print '\nCreating output root file %s/%s.root' % (OUTDIR, OUTDIR)
    if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)
    out_file = R.TFile('%s/%s.root' % (OUTDIR,OUTDIR), 'recreate')
    out_file.cd()

    if not os.path.exists('%s/png' % OUTDIR):
        print '\nCreating output folders %s/png'  % OUTDIR
        os.makedirs('%s/png' % OUTDIR)
        os.makedirs('%s/pdf' % OUTDIR)

    ## Input histograms and canvases to draw on
    hist = {}
    ## Graphs of delta(RECO, GEN) vs. RECO muon d0 * charge
    graph = {}
    graph_data = {}
    graph_mc = {}
    ## Data to fill the graph
    data = {}
    mc = {}
    ## Fits to the data
    fit = {}
    ## Canvases to draw the graph on
    canv = {}
    pad1 = {}
    pad2 = {}
    h_1 = R.TH1F('h_1', 'h_1', 22, -0.012, 0.012)
    h_2 = R.TH1F('h_2', 'h_2', 22, -0.012, 0.012)

    for var in ['dimu_mass']:
        # for corr in ['', 'Roch', 'KaMu', 'kinfit', 'KinRoch', 'KinKaMu']:
        for corr in ['Roch','corr']:

            ## Unique string identifying this canvas
            if (corr == ''): p_str = '%s%s' % (var, corr)
            else: p_str = '%s_%s' % (var, corr)

            ## Unique string identifying this graph

            for eta in ['eta_inc']: # |eta| binned
                c_str = eta + '_' + p_str 
                ## Canvas for the graphs of peak vs. d0
                canv[c_str] = R.TCanvas(c_str, c_str, 1600, 1200)
                pad1[c_str] = R.TPad(c_str+'_1', c_str+'_1', 0, 0.31, 1, 1.0)
                pad2[c_str] = R.TPad(c_str+'_2', c_str+'_2', 0, 0.05, 1, 0.3)

                pad1[c_str].SetBottomMargin(0.02)
                pad1[c_str].Draw()
                pad2[c_str].SetTopMargin(0.02);
                pad2[c_str].SetBottomMargin(0.3);
                pad2[c_str].Draw();
                pad1[c_str].cd()
                for pt in ['inclusive']:
                    g_str = pt + '_' + c_str

                    # for 80 bins of d0
                    data[g_str] = {}
                    data[g_str]['x']    = numpy.empty([D0_BINS], dtype=float)
                    data[g_str]['y']    = numpy.empty([D0_BINS], dtype=float)
                    data[g_str]['xerr'] = numpy.empty([D0_BINS], dtype=float)
                    data[g_str]['yerr'] = numpy.empty([D0_BINS], dtype=float)

                    mc[g_str] = {}
                    mc[g_str]['x']    = numpy.empty([D0_BINS], dtype=float)
                    mc[g_str]['y']    = numpy.empty([D0_BINS], dtype=float)
                    mc[g_str]['xerr'] = numpy.empty([D0_BINS], dtype=float)
                    mc[g_str]['yerr'] = numpy.empty([D0_BINS], dtype=float)

                    # print 'Getting data for %s' % g_str
                    iDat = -1
                    for d0 in range(D0_BINS): # |d0| < 0.01 80 bins
                        iDat += 1
                     
                        h_data = in_file_data.Get('h_%s_vs_d0_BS' % (g_str)).Clone() # dPt2p0 vs d0_BS
                        h_mc = in_file_mc.Get('h_%s_vs_d0_BS' % (g_str)).Clone() # dPt2p0 vs d0_BS
                        # hRebin = h.Clone()
                        # hRebin.RebinY(10)

                        iX = [iDat*20+1, iDat*20+20]

                        h_data.GetXaxis().SetRange(iX[0], iX[1])

                        data[g_str]['x'][iDat] = h_data.GetMean(1)
                        data[g_str]['xerr'][iDat] = h_data.GetStdDev(1)
                        # h.GetXaxis().SetRange()
                        hh_data = h_data.ProjectionY('', iX[0], iX[1], 'e').Clone()

                        mean_val, mean_err, reso_val, reso_err = FitVoigtian(hh_data)


                        data[g_str]['y'][iDat] = mean_val
                        data[g_str]['yerr'][iDat] = mean_err

                        h_mc.GetXaxis().SetRange(iX[0], iX[1])

                        mc[g_str]['x'][iDat] = h_mc.GetMean(1)
                        mc[g_str]['xerr'][iDat] = h_mc.GetStdDev(1)
                        # h.GetXaxis().SetRange()
                        hh_mc = h_mc.ProjectionY('', iX[0], iX[1], 'e').Clone()

                        mean_val, mean_err, reso_val, reso_err = FitVoigtian(hh_mc)


                        mc[g_str]['y'][iDat] = mean_val
                        mc[g_str]['yerr'][iDat] = mean_err

                    ## End loop: for d0 

                    ## Create a graph with the peak values vs. d0
                    graph_data[g_str] = R.TGraphErrors(len(data[g_str]['x']), data[g_str]['x'], data[g_str]['y'], data[g_str]['xerr'], data[g_str]['yerr'])
                    graph_data[g_str].SetName(g_str)

                    graph_mc[g_str] = R.TGraphErrors(len(mc[g_str]['x']), mc[g_str]['x'], mc[g_str]['y'], mc[g_str]['xerr'], mc[g_str]['yerr'])
                    graph_mc[g_str].SetName(g_str + 'mc')

                    
                    var_str = 'm_{\mu\mu}'

                    if corr == 'corr':
                        corr_str = '(GeoFit)'
                    if corr == 'Roch':
                        corr_str = '(Rochester)'

                    if eta == 'eta_inc':
                        eta_str = 'inclusive \eta '

                    pad1[c_str].cd()

                    graph_data[g_str].SetTitle('')
                    graph_data[g_str].GetXaxis().SetTitle('d0_{BS,\mu^{+}} (cm)')
                    graph_data[g_str].GetYaxis().SetTitle('%s %s' % (var_str, corr_str))
                    graph_data[g_str].GetYaxis().SetTitleOffset(0.7)
                    print graph_data[g_str].GetYaxis().GetTitleSize()
                    graph_data[g_str].GetYaxis().SetTitleSize(0.05)
                    print graph_data[g_str].GetYaxis().GetTitleSize()

                    graph_mc[g_str].SetTitle('')
                    graph_mc[g_str].GetXaxis().SetTitle('d0_{BS,\mu^{+}} (cm)')
                    graph_mc[g_str].GetYaxis().SetTitle('#scale[2]{%s %s}' % (var_str, corr_str))
                    # graph_mc[g_str].GetYaxis().SetTitleOffset(1.1)
                    # graph_data[g_str].GetYaxis().SetTitleSize(20)
                    # graph[g_str].GetXaxis().SetLabelSize(0)

                    ## Draw the graph onto the canvas

                    graph_data[g_str].SetLineColor(1)
                    graph_data[g_str].SetMarkerColor(1)
                    graph_data[g_str].SetMarkerStyle(8)
                    graph_data[g_str].SetMarkerSize(0.7)
                    # graph_data[g_str].Fit('f_%s' % g_str,"R")
                    print 'Writing: %s' % g_str
                    out_file.cd()
                    graph_data[g_str].Write()

                    graph_mc[g_str].SetLineColor(2)
                    graph_mc[g_str].SetMarkerColor(2)
                    graph_mc[g_str].SetMarkerStyle(8)
                    graph_mc[g_str].SetMarkerSize(0.7)
                    # graph_data[g_str].Fit('f_%s' % g_str,"R")
                    print 'Writing: %s' % g_str
                    out_file.cd()
                    graph_mc[g_str].Write()

                    # g_canv = R.TCanvas(g_str, g_str, 1600, 1200)
                    pad1[c_str].cd()
                    graph_data[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                    graph_data[g_str].GetXaxis().SetLabelSize(0)
                    graph_data[g_str].GetYaxis().SetRangeUser(88, 94)
                    graph_mc[g_str].GetYaxis().SetRangeUser(88, 94)
                    graph_data[g_str].Draw('AP')
                    graph_mc[g_str].Draw('Psame')

                    # g_leg = R.TLegend(0.12,0.76,0.72,0.88)

                    # g_leg.Draw('same')
                    # g_latex = R.TLatex()
                    # g_latex.SetTextAlign(12)
                    # g_latex.SetTextSize(0.04)
                    # g_latex.DrawLatex(-0.01, 93.5, eta_str)

                    g_leg = R.TLegend(0.12,0.76,0.32,0.88)
                    g_leg.AddEntry( graph_data[g_str], 'Z+jets data' , 'LP' )
                    g_leg.AddEntry( graph_mc[g_str], 'Z+jets MC' , 'LP' )
                    # g_leg.AddEntry( graph_DY, 'DY sample', 'LAP'  )
                    # g_leg.AddEntry( graph_tt, 'ttbar sample' , 'LAP' )
                    g_leg.Draw('same')


                    cms_latex = R.TLatex()
                    cms_latex.SetTextAlign(11)
                    cms_latex.SetTextSize(0.03)
                    cms_latex.DrawLatexNDC(0.11, 0.91, '#scale[1.5]{CMS}')
                    cms_latex.DrawLatexNDC(0.17, 0.91,'#font[52]{#scale[1.2]{preliminary}}')
                    # cms_latex.DrawLatexNDC(0.84, 0.93,'#scale[1.3]{%s}' % YEAR)
                    cms_latex.DrawLatexNDC(0.79, 0.91,'#font[42]{35.9 fb^{-1} (13 TeV)}') #35.9 41.5 59.7

                    y_1 = graph_data[g_str].GetY()
                    yerr_1 = graph_data[g_str].GetEY()
                    y_2 = graph_mc[g_str].GetY()
                    yerr_2 = graph_mc[g_str].GetEY()
                    for i in range(D0_BINS):
                        # data['ratio']['x'][i] = data['nVtx_0_21_'+c_str]['x'][i]
                        # data['ratio']['xerr'][i] = data['nVtx_0_21_'+c_str]['xerr'][i]
                        # data['ratio']['y'][i] = data['nVtx_34_inf_'+c_str]['y'][i]/data['nVtx_0_21_'+c_str]['y'][i]
                        # data['ratio']['yerr'][i] = ((data['nVtx_34_inf_'+c_str]['yerr'][i]/data['nVtx_34_inf_'+c_str]['y'][i])+(data['nVtx_0_21_'+c_str]['yerr'][i]/data['nVtx_0_21_'+c_str]['y'][i]))*data['ratio']['y'][i]
                        h_1.SetBinContent(i+2, y_1[i])
                        h_1.SetBinError(i+2, yerr_1[i])
                        h_2.SetBinContent(i+2, y_2[i])
                        h_2.SetBinError(i+2, yerr_2[i])
                
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

                    h_ratio.GetYaxis().SetTitle('Data / MC')
                    h_ratio.GetYaxis().SetRangeUser(0.99,1.01)
                    h_ratio.GetYaxis().CenterTitle(1)
                    h_ratio.GetYaxis().SetNdivisions(505)
                    h_ratio.GetYaxis().SetTitleSize(40)
                    h_ratio.GetYaxis().SetTitleFont(43)
                    h_ratio.GetYaxis().SetTitleOffset(1.1)
                    h_ratio.GetYaxis().SetLabelFont(43)
                    h_ratio.GetYaxis().SetLabelSize(25)

                    h_ratio.GetXaxis().SetTitle('d0_{BS,\mu^{+}} (cm)')
                    h_ratio.GetXaxis().SetTitleSize(30)
                    h_ratio.GetXaxis().SetTitleFont(43)
                    h_ratio.GetXaxis().SetTitleOffset(3.6)
                    h_ratio.GetXaxis().SetLabelFont(43)
                    h_ratio.GetXaxis().SetLabelSize(25)
                    # h_ratio.SetMarkerStyle(21)
                    # r = R.TGraphErrors(len(data['ratio']['x']), data['ratio']['x'], data['ratio']['y'], data['ratio']['xerr'], data['ratio']['yerr'])
                    # r.SetLineColor(1)
                    # r.GetYaxis().SetRangeUser(0.5, 1.5)
                    # r.Draw('AP')
                    h_ratio.Draw('ep')

                    canv[c_str].SaveAs('%s/png/%s.png' % (OUTDIR, g_str))
                    canv[c_str].SaveAs('%s/pdf/%s.pdf' % (OUTDIR, g_str))
                    
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
    in_file_data.Close()
    in_file_mc.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
