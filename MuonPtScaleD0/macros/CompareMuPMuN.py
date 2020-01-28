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

YEAR = '2016'
# YEAR = '2017'
# YEAR = '2018'

# 2D DY samples 
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_2D_PUBinned/files/HADD' % YEAR

INFILE = 'histos_ZJets_hadd.root'
# INFILE = 'histos_ttbar_hadd.root'
# INFILE = 'histos_ttH_hadd.root'

OUTDIR = 'plots_muP_vs_muN_AN_%s' % YEAR

D0_BINS = 19

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Inside FindPeaks.py ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    print '\nCreating output root file %s/%s.root' % (OUTDIR,OUTDIR)
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
    ## Data to fill the graph
    data = {}
    ## Fits to the data
    fit = {}
    fit_left = {}
    fit_right = {}
    ## Canvases to draw the graph on
    canv = {}
    pad1 = {}
    pad2 = {}

    h_1 = R.TH1F('h_1', 'h_1', 22, -0.011, 0.011)
    h_2 = R.TH1F('h_2', 'h_2', 22, -0.011, 0.011)


    for var in ['dRelPt2p0']:
        for corr in ['Roch']:

            ## Unique string identifying this canvas
            if (corr == ''): p_str = '%s%s' % (var, corr)
            else: p_str = '%s_%s' % (var, corr)

            ## Unique string identifying this graph

            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf', 'eta_inc']: # |eta| binned
                c_str = eta + '_' + p_str 
                ## Canvas for the graphs of peak vs. d0
                canv[c_str] = R.TCanvas(c_str, c_str, 1200, 1200)
                pad1[c_str] = R.TPad(c_str+'_1', c_str+'_1', 0, 0.31, 1, 1.0)
                pad2[c_str] = R.TPad(c_str+'_2', c_str+'_2', 0, 0.05, 1, 0.3)

                pad1[c_str].SetBottomMargin(0.02)
                pad1[c_str].Draw()
                pad2[c_str].SetTopMargin(0.02);
                pad2[c_str].SetBottomMargin(0.2);
                pad2[c_str].Draw();

                iPt = -1
                for pt in ['inclusive']:
                    for charge in ['muP', 'muN']:
                        iPt += 1

                        g_str = pt + '_' + c_str + '_' + charge
                        gg_str = pt + '_' + c_str

                        # for bins of d0
                        data[g_str] = {}
                        data[g_str]['x']    = numpy.empty([D0_BINS], dtype=float)
                        data[g_str]['y']    = numpy.empty([D0_BINS], dtype=float)
                        data[g_str]['xerr'] = numpy.empty([D0_BINS], dtype=float)
                        data[g_str]['yerr'] = numpy.empty([D0_BINS], dtype=float)

                        # print 'Getting data for %s' % g_str
                        iDat = -1
                        for d0 in range(D0_BINS): # |d0| < 0.01 80 bins
                            iDat += 1
                            
                            h = in_file.Get('h_%s_vs_d0_BS_%s' % (gg_str, charge)) # dPt2p0 vs d0_BS
                            hRebin = h.Clone()
                            hRebin.RebinY(10)

                            iX = [iDat*20+1, iDat*20+20]

                            h.GetXaxis().SetRange(iX[0], iX[1])

                            data[g_str]['x'][iDat] = h.GetMean(1)
                            data[g_str]['xerr'][iDat] = h.GetStdDev(1)
                            h.GetXaxis().SetRange()

                            yMax = 0 ## Maximum integral over 101 bins  
                            iMax = 0 ## Central bin of maximum integral (i.e. the center of the peak) 
                            iLo = 0
                            iHi = 0
                            if (iDat > 7 and iDat < 13): 
                                for i in range(52, h.GetNbinsY() - 50):  ## Don't include over- or under-flow
                                    if (h.Integral(iX[0], iX[1], i-50, i+50) > yMax):
                                        yMax = h.Integral(iX[0], iX[1], i-50, i+50)
                                        iMax = i
                                iLo = h.GetNbinsY()    ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                                iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                                for i in range(52, h.GetNbinsY() - 50):  ## Don't include over- or under-flow
                                    if (h.Integral(iX[0], iX[1], i-50, i+50) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                                    if (h.Integral(iX[0], iX[1], i-50, i+50) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)
                                data[g_str]['y'][iDat] = h.GetYaxis().GetBinCenter(iMax)
                                data[g_str]['yerr'][iDat] = (h.GetYaxis().GetBinLowEdge(iHi+1) - h.GetYaxis().GetBinLowEdge(iLo)) / 2.
                            else:
                                for i in range(7, hRebin.GetNbinsY() - 5):  ## Don't include over- or under-flow
                                    if (hRebin.Integral(iX[0], iX[1], i-5, i+5) > yMax):
                                        yMax = hRebin.Integral(iX[0], iX[1], i-5, i+5)
                                        iMax = i
                                iLo = hRebin.GetNbinsY()  ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                                iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                                for i in range(7, hRebin.GetNbinsY() - 5):  ## Don't include over- or under-flow
                                    if (hRebin.Integral(iX[0], iX[1], i-5, i+5) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                                    if (hRebin.Integral(iX[0], iX[1], i-5, i+5) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)
                                data[g_str]['y'][iDat] = hRebin.GetYaxis().GetBinCenter(iMax)
                                data[g_str]['yerr'][iDat] = (hRebin.GetYaxis().GetBinLowEdge(iHi+1) - hRebin.GetYaxis().GetBinLowEdge(iLo)) / 2.



                        ## End loop: for d0 

                        ## Create a graph with the peak values vs. d0
                        graph[g_str] = R.TGraphErrors(len(data[g_str]['x']), data[g_str]['x'], data[g_str]['y'], data[g_str]['xerr'], data[g_str]['yerr'])
                        graph[g_str].SetName(g_str)

                        
                        var_str = '10000 * (p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}'

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


                        graph[g_str].SetTitle('')
                        graph[g_str].GetXaxis().SetTitle('d0_{BS} (cm) #times charge')
                        graph[g_str].GetYaxis().SetTitle('%s %s' % (var_str, corr_str))
                        graph[g_str].GetYaxis().SetTitleOffset(1.1)

                        ## Create a fit to the points on the graph
                        fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*x', -0.01, 0.01) # linear fit
                        fit[g_str].SetParameter(0, 0)

                        # Set fit parameters
                        #for linear
                        if eta == 'eta_0_0p9':
                          fit[g_str].SetParameter(1, 300) # linear
                        elif eta == 'eta_0p9_1p7':
                          fit[g_str].SetParameter(1, 600) # linear
                        elif eta == 'eta_1p7_inf':
                          fit[g_str].SetParameter(1, 900) # linear
                        elif eta == 'eta_inc':
                          fit[g_str].SetParameter(1, 300) # linear

                        ## Draw the graph onto the canvas

                        graph[g_str].SetLineColor(iPt+1)
                        graph[g_str].SetMarkerColor(iPt+1)
                        graph[g_str].SetMarkerStyle(8)
                        graph[g_str].SetMarkerSize(0.7)
                        fit[g_str].SetLineColor(iPt+1)
                        graph[g_str].Fit('f_%s' % g_str,"R")
                        print 'Writing: %s' % g_str
                        out_file.cd()
                        graph[g_str].Write()

                        g_canv = R.TCanvas(g_str, g_str, 1600, 1200)
                        g_canv.cd()
                        graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                        graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                        graph[g_str].GetXaxis().SetLabelSize(0)
                        graph[g_str].Draw('AP')
                        
                        g_leg = R.TLegend(0.12,0.76,0.72,0.88)
                        g_leg.AddEntry( fit[g_str], '\LARGE{y = %.2f + %.2f *d0}' % (fit[g_str].GetParameter(0), fit[g_str].GetParameter(1)) )
                        g_leg.Draw('same')
                        g_latex = R.TLatex()
                        g_latex.SetTextAlign(12)
                        g_latex.SetTextSize(0.04)
                        g_latex.DrawLatex(-0.01, 8, eta_str)
                        g_canv.SaveAs('%s/png/%s.png' % (OUTDIR, g_str))
                        
                        pad1[c_str].cd()
                        cms_latex = R.TLatex()
                        cms_latex.SetTextAlign(11)
                        cms_latex.SetTextSize(0.03)
                        cms_latex.DrawLatexNDC(0.11, 0.93, '#scale[1.5]{CMS}')
                        cms_latex.DrawLatexNDC(0.18, 0.93,'#font[52]{#scale[1.2]{preliminary simulation}}')
                        cms_latex.DrawLatexNDC(0.84, 0.93,'#scale[1.3]{%s}' % YEAR)
                        if (iPt == 0):
                            graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                            graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                            # graph[g_str].GetYaxis().SetLabelSize(0.);
                            graph[g_str].Draw('AP')
                            # axis = R.TGaxis( -5, 20, -5, 220, 20,220,510,"");
                            # axis.SetLabelFont(43);
                            # axis.SetLabelSize(15);
                            # axis.Draw();
                        else:
                            graph[g_str].Draw('Psame')


                        latex = R.TLatex()
                        latex.SetTextAlign(10)
                        latex.SetTextSize(0.04)
                        latex.DrawLatex(-0.01, 8, '#font[42]{%s}' % eta_str)

                    # End loop for charge

                # End loop for pt

                leg = R.TLegend(0.12,0.76,0.72,0.88)
                leg.AddEntry( fit['inclusive_'+c_str+'_muP'], 'muP : y = %.2f + %.2f \pm %.2f *d0' % (fit['inclusive_'+c_str+'_muP'].GetParameter(0), fit['inclusive_'+c_str+'_muP'].GetParameter(1), fit['inclusive_'+c_str+'_muP'].GetParError(1)) )
                leg.AddEntry( fit['inclusive_'+c_str+'_muN'], 'muN : y = %.2f + %.2f \pm %.2f *d0' % (fit['inclusive_'+c_str+'_muN'].GetParameter(0), fit['inclusive_'+c_str+'_muN'].GetParameter(1), fit['inclusive_'+c_str+'_muN'].GetParError(1)) )
                leg.Draw('same')
                
                # h_ratio = R.TH1F('h_ratio', 'h_ratio', 40, -0.01, 0.01)
                for i in range(D0_BINS):
                    # data['ratio']['x'][i] = data['nVtx_0_21_'+c_str]['x'][i]
                    # data['ratio']['xerr'][i] = data['nVtx_0_21_'+c_str]['xerr'][i]
                    # data['ratio']['y'][i] = data['nVtx_34_inf_'+c_str]['y'][i]/data['nVtx_0_21_'+c_str]['y'][i]
                    # data['ratio']['yerr'][i] = ((data['nVtx_34_inf_'+c_str]['yerr'][i]/data['nVtx_34_inf_'+c_str]['y'][i])+(data['nVtx_0_21_'+c_str]['yerr'][i]/data['nVtx_0_21_'+c_str]['y'][i]))*data['ratio']['y'][i]
                    h_1.SetBinContent(i+2, data['inclusive_'+c_str+'_muP']['y'][i])
                    h_1.SetBinError(i+2, data['inclusive_'+c_str+'_muP']['yerr'][i])
                    h_2.SetBinContent(i+2, data['inclusive_'+c_str+'_muN']['y'][i])
                    h_2.SetBinError(i+2, data['inclusive_'+c_str+'_muN']['yerr'][i])
                
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
                h_ratio.GetYaxis().CenterTitle(1)
                h_ratio.GetYaxis().SetNdivisions(505)
                h_ratio.GetYaxis().SetTitleSize(40)
                h_ratio.GetYaxis().SetTitleFont(43)
                h_ratio.GetYaxis().SetTitleOffset(1.1)
                h_ratio.GetYaxis().SetLabelFont(43)
                h_ratio.GetYaxis().SetLabelSize(25)

                h_ratio.GetXaxis().SetTitle('d0_{BS} (cm) * charge')
                h_ratio.GetXaxis().SetTitleSize(30)
                h_ratio.GetXaxis().SetTitleFont(43)
                h_ratio.GetXaxis().SetTitleOffset(3.2)
                h_ratio.GetXaxis().SetLabelFont(43)
                h_ratio.GetXaxis().SetLabelSize(25)
                # h_ratio.SetMarkerStyle(21)
                # r = R.TGraphErrors(len(data['ratio']['x']), data['ratio']['x'], data['ratio']['y'], data['ratio']['xerr'], data['ratio']['yerr'])
                # r.SetLineColor(1)
                # r.GetYaxis().SetRangeUser(0.5, 1.5)
                # r.Draw('AP')
                h_ratio.Draw('ep')
                canv[c_str].SaveAs('%s/png/%s.png' % (OUTDIR, c_str))
                canv[c_str].SaveAs('%s/pdf/%s.pdf' % (OUTDIR, c_str))

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
