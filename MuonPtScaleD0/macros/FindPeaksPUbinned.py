#! /usr/bin/env python

###############################################
###              FindPeaks.py               ###
###                                         ###
###    Find peaks in dPt and dPhi plots     ###
###    (RECO - GEN) vs. d0 in ZJets MC.     ###
###                                         ###
###           Andrew Brinkerhoff            ###
###               27.08.2018                ###
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
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_DY_PUbinned_bugFix/files/HADD' %YEAR

# DY samples with d0_BS
INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_d0BS/files/HADD' %YEAR

# ttbar samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_ttbar_PUbinned_bugFix/files/HADD' %YEAR

# ttH samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_ttH/files/HADD' %YEAR

INFILE = 'histos_ZJets_hadd.root'
# INFILE = 'histos_ttbar_hadd.root'
# INFILE = 'histos_ttH_hadd.root'

## Main function executed by ./macros/FindPeaks.py
def main():

    print '\n\n*** Inside FindPeaks.py ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    print '\nCreating output root file plots_DY_d0BS_linear_PUbinned_%s/FindPeaks.root' % YEAR
    if not os.path.exists('plots_DY_d0BS_linear_PUbinned_%s' % YEAR):
      os.makedirs('plots_DY_d0BS_linear_PUbinned_%s' % YEAR)
    out_file = R.TFile('plots_DY_d0BS_linear_PUbinned_%s/FindPeaks.root' % YEAR, 'recreate')
    out_file.cd()

    if not os.path.exists('plots_DY_d0BS_linear_PUbinned_%s/png' % YEAR):
        print '\nCreating output folders plots_DY_d0BS_linear_PUbinned_%s/png'  % YEAR
        os.makedirs('plots_DY_d0BS_linear_PUbinned_%s/png' % YEAR)

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

    for var in ['dRelPt2p0']: 
        for corr in ['', 'Roch', 'KaMu', 'kinfit', 'KinRoch', 'KinKaMu']:

            ## Unique string identifying this canvas
            if (corr == ''): p_str = '%s%s' % (var, corr)
            else: p_str = '%s_%s' % (var, corr)

            ## Unique string identifying this graph
            for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf']: # |eta| binned
                c_str = eta + '_' + p_str 
                ## Canvas for the graphs of peak vs. d0
                canv[c_str] = R.TCanvas(c_str, c_str, 1600, 1200)
                iVtx = -1
                for vtx in ['nVtx_0_21', 'nVtx_22_26', 'nVtx_27_33', 'nVtx_34_inf']:
                    iVtx += 1
                    g_str = vtx + '_' + c_str

                    # for 80 bins of d0
                    data[g_str] = {}
                    data[g_str]['x']    = numpy.empty([79], dtype=float)
                    data[g_str]['y']    = numpy.empty([79], dtype=float)
                    data[g_str]['xerr'] = numpy.empty([79], dtype=float)
                    data[g_str]['yerr'] = numpy.empty([79], dtype=float)

                    iDat = -1
                    for d0 in ['m39','m38','m37','m36','m35','m34','m33','m32','m31','m30','m29','m28','m27','m26','m25','m24','m23','m22','m21','m20','m19','m18','m17','m16','m15','m14','m13','m12','m11','m10','m9','m8','m7','m6','m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28','p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39']: # |d0| < 0.01 80 bins
                        iDat += 1

                        h_x = in_file.Get('h_d0_%s_%s_%s_d0' % (d0, vtx, eta)) #d0
                        h_y = in_file.Get('h_d0_%s_%s' % (d0, g_str)) #dPt

                        data[g_str]['x'][iDat] = h_x.GetMean()
                        data[g_str]['xerr'][iDat] = h_x.GetStdDev()

                        yMax = 0 ## Maximum integral over 101 bins  
                        iMax = 0 ## Central bin of maximum integral (i.e. the center of the peak) 
                        iLo = 0
                        iHi = 0
                        if (iDat > 29 and iDat < 49):
                            for i in range(52, h_y.GetNbinsX() - 50):  ## Don't include over- or under-flow
                                if (h_y.Integral(i-50, i+50) > yMax):
                                    yMax = h_y.Integral(i-50, i+50)
                                    iMax = i
                            iLo = h_y.GetNbinsX()  ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                            iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                            for i in range(52, h_y.GetNbinsX() - 50):  ## Don't include over- or under-flow
                                if (h_y.Integral(i-50, i+50) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                                if (h_y.Integral(i-50, i+50) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)
                        else:
                            if h_y.GetNbinsX() >= 4000:
                              h_y.Rebin(10) # Rebin the histograms at the tails of d0 distribution 
                            for i in range(7, h_y.GetNbinsX() - 5):  ## Don't include over- or under-flow
                                if (h_y.Integral(i-5, i+5) > yMax):
                                    yMax = h_y.Integral(i-5, i+5)
                                    iMax = i
                            iLo = h_y.GetNbinsX()  ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                            iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                            for i in range(7, h_y.GetNbinsX() - 5):  ## Don't include over- or under-flow
                                if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                                if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)
                        
                        

                        data[g_str]['y'][iDat] = h_y.GetBinCenter(iMax)
                        data[g_str]['yerr'][iDat] = (h_y.GetBinLowEdge(iHi+1) - h_y.GetBinLowEdge(iLo)) / 2.



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
                        eta_str = '0 < |\eta| '
                    elif eta == 'eta_0_0p9':
                        eta_str = '0 < |\eta| < 0.9'
                    elif eta == 'eta_0p9_1p7':
                        eta_str = '0.9 < |\eta| < 1.7'
                    elif eta == 'eta_1p7_inf':
                        eta_str = '1.7 < |\eta|'


                    graph[g_str].SetTitle('%s  %s  vs.  d0' % (var_str, corr_str))
                    graph[g_str].GetXaxis().SetTitle('RECO muon d0 (cm) * charge')
                    graph[g_str].GetYaxis().SetTitle('%s %s' % (var_str, corr_str))

                    ## Create a fit to the points on the graph
                    fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*x', -0.01, 0.01) # linear fit
                    # fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*tanh([2]*x) + [3]*x', -0.009, 0.009) # tanh + linear fit
                    
                    # Set fit parameters
                    #for linear
                    if eta == 'eta_0_0p9':
                      fit[g_str].SetParameter(1, 300) # linear
                    elif eta == 'eta_0p9_1p7':
                      fit[g_str].SetParameter(1, 600) # linear
                    elif eta == 'eta_1p7_inf':
                      fit[g_str].SetParameter(1, 900) # linear

                    # for tanh + linear
                    # if eta == 'eta_0_0p9':
                    #   fit[g_str].SetParameter(1, 3) # tanh
                    #   fit[g_str].SetParameter(2, 500) # tanh
                    #   fit[g_str].SetParameter(3, 100) # linear
                    # elif eta == 'eta_0p9_1p7':
                    #   fit[g_str].SetParameter(1, 3) # tanh
                    #   fit[g_str].SetParameter(2, 500) # tanh
                    #   fit[g_str].SetParameter(3, 400) # linear
                    # elif eta == 'eta_1p7_inf':
                    #   fit[g_str].SetParameter(1, 3) # tanh
                    #   fit[g_str].SetParameter(2, 500) # tanh
                    #   fit[g_str].SetParameter(3, 500) # linear


                    ## Draw the graph onto the canvas

                    graph[g_str].SetLineColor(iVtx+1)
                    graph[g_str].SetMarkerColor(iVtx+1)
                    graph[g_str].SetMarkerStyle(8)
                    graph[g_str].SetMarkerSize(0.7)
                    fit[g_str].SetLineColor(iVtx+1)
                    graph[g_str].Fit('f_%s' % g_str,"R")
                    print 'Writing: %s' % g_str
                    out_file.cd()
                    graph[g_str].Write()

                    # g_canv = R.TCanvas(g_str, g_str, 1600, 1200)
                    # g_canv.cd()
                    # graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                    # graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                    # graph[g_str].Draw('AP')
                    # g_leg = R.TLegend(0.12,0.76,0.72,0.88)

                    # g_leg.AddEntry( fit[g_str], '\LARGE{y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)}' % (fit[g_str].GetParameter(0), fit[g_str].GetParameter(3), fit[g_str].GetParameter(1), fit[g_str].GetParameter(2)) )
                    # g_leg.AddEntry( fit[g_str], '\LARGE{y = %.2f + %.2f *d0}' % (fit[g_str].GetParameter(0), fit[g_str].GetParameter(1)) )
                    # g_leg.Draw('same')
                    # g_latex = R.TLatex()
                    # g_latex.SetTextAlign(12)
                    # g_latex.SetTextSize(0.04)
                    # g_latex.DrawLatex(-0.01, 8, eta_str)
                    
                    canv[c_str].cd()
                    if (iVtx == 0):
                        graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                        graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                        graph[g_str].Draw('AP')
                    else:
                        graph[g_str].Draw('Psame')

                    latex = R.TLatex()
                    latex.SetTextAlign(12)
                    latex.SetTextSize(0.04)
                    latex.DrawLatex(-0.01, 8, eta_str)
                    
                    


                # End loop for pt

                # f5 = fit['inclusive_'+c_str]
                f1 = fit['nVtx_0_21_'+c_str]
                f2 = fit['nVtx_22_26_'+c_str]
                f3 = fit['nVtx_27_33_'+c_str]
                f4 = fit['nVtx_34_inf_'+c_str]



                if (var == 'dPhi' or var == 'dEta'): leg = R.TLegend(0.48,0.66,0.88,0.88)
                else:               leg = R.TLegend(0.12,0.66,0.52,0.88)

                # leg.AddEntry( f1, '0 <= nVtx <= 21: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f1.GetParameter(0), f1.GetParameter(3), f1.GetParameter(1), f1.GetParameter(2)) )
                # leg.AddEntry( f2, '22 <= nVtx <= 26: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f2.GetParameter(0), f2.GetParameter(3), f2.GetParameter(1), f2.GetParameter(2)) )
                # leg.AddEntry( f3, '27 <= nVtx <= 33: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f3.GetParameter(0), f3.GetParameter(3), f3.GetParameter(1), f3.GetParameter(2)) )
                # leg.AddEntry( f4, '34 <= nVtx: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f4.GetParameter(0), f4.GetParameter(3), f4.GetParameter(1), f4.GetParameter(2)) )

                leg.AddEntry( f1, '0 <= nVtx <= 21: y = %.2f + d0 * %.2f' % (f1.GetParameter(0), f1.GetParameter(1)) )
                leg.AddEntry( f2, '22 <= nVtx <= 26: y = %.2f + d0 * %.2f' % (f2.GetParameter(0), f2.GetParameter(1)) )
                leg.AddEntry( f3, '27 <= nVtx <= 33: y = %.2f + d0 * %.2f' % (f3.GetParameter(0), f3.GetParameter(1)) )
                leg.AddEntry( f4, '34 <= nVtx: y = %.2f + d0 * %.2f' % (f4.GetParameter(0), f4.GetParameter(1)) )
                leg.Draw('same')

                canv[c_str].SaveAs('plots_DY_d0BS_linear_PUbinned_%s/png/%s.png' % (YEAR, c_str))
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
