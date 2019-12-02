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
YEAR = '2017'
# YEAR = '2018'

# Bugfix DY samples
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_DY_PUbinned_bugFix/files/HADD' %YEAR

# DY samples with d0_BS
# INDIR = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/%s/Histograms/MassCal_ptVsd0_d0BS/files/HADD' %YEAR

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

    print '\nCreating output root file plots_DY_d0BS_linear_%s/FindPeaks.root' % YEAR
    if not os.path.exists('plots_DY_d0BS_linear_%s' % YEAR):
      os.makedirs('plots_DY_d0BS_linear_%s' % YEAR)
    out_file = R.TFile('plots_DY_d0BS_linear_%s/FindPeaks.root' % YEAR, 'recreate')
    out_file.cd()

    if not os.path.exists('plots_DY_d0BS_linear_%s/png' % YEAR):
        print '\nCreating output folders plots_DY_d0BS_linear_%s/png'  % YEAR
        os.makedirs('plots_DY_d0BS_linear_%s/png' % YEAR)
        os.makedirs('plots_DY_d0BS_linear_%s/png/histos' % YEAR)
        # os.makedirs('plots/pdf')

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
    ## Canvases for individual graphs
    # ccanv = {}

    # dPt     =        RECO - GEN pT
    # dRelPt  =   100*(RECO - GEN) / GEN pT
    # dRelPt2 = 10000*(RECO - GEN) / (GEN pT)^2 
    # dPhi    = (RECO - GEN phi) * GEN charge
    # for var in ['dPt',  'dRelPt2p0', 'dRelPt1p0']: #, 'dPhi', 'dEta']:  ## Variables plotted vs. RECO d0 * charge
    for var in ['dRelPt2p0']: #, 'dPhi', 'dEta']:  ## Variables plotted vs. RECO d0 * charge
    # for var in ['dRelPt2']:  ## Variables plotted vs. RECO d0 * charge
        for corr in ['', 'Roch', 'KaMu', 'kinfit', 'KinRoch', 'KinKaMu']:
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
                canv[c_str] = R.TCanvas(c_str, c_str, 1600, 1200)
                iPt = -1
                # for pt in ['pt_20_35']:
                # for pt in ['inclusive', 'pt_20_35', 'pt_35_42', 'pt_42_50', 'pt_50_inf']:
                for pt in ['inclusive']:
                # for pt in ['nVtx_0_21', 'nVtx_22_26', 'nVtx_34_inf']:
                    iPt += 1
                    g_str = pt + '_' + c_str

                    # for 20 bins
                    # data[g_str] = {}
                    # data[g_str]['x']    = numpy.empty([19], dtype=float)
                    # data[g_str]['y']    = numpy.empty([19], dtype=float)
                    # data[g_str]['xerr'] = numpy.empty([19], dtype=float)
                    # data[g_str]['yerr'] = numpy.empty([19], dtype=float)

                    # # for 40 bins
                    # data[g_str] = {}
                    # data[g_str]['x']    = numpy.empty([39], dtype=float)
                    # data[g_str]['y']    = numpy.empty([39], dtype=float)
                    # data[g_str]['xerr'] = numpy.empty([39], dtype=float)
                    # data[g_str]['yerr'] = numpy.empty([39], dtype=float)

                    # for 80 bins
                    data[g_str] = {}
                    data[g_str]['x']    = numpy.empty([79], dtype=float)
                    data[g_str]['y']    = numpy.empty([79], dtype=float)
                    data[g_str]['xerr'] = numpy.empty([79], dtype=float)
                    data[g_str]['yerr'] = numpy.empty([79], dtype=float)

                    # print 'Getting data for %s' % g_str
                    iDat = -1
                    for d0 in ['m39','m38','m37','m36','m35','m34','m33','m32','m31','m30','m29','m28','m27','m26','m25','m24','m23','m22','m21','m20','m19','m18','m17','m16','m15','m14','m13','m12','m11','m10','m9','m8','m7','m6','m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28','p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39']: # |d0| < 0.01 80 bins
                    # for d0 in ['m19','m18','m17','m16','m15','m14','m13','m12','m11','m10','m9','m8','m7','m6','m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19']: # |d0| < 0.005 40 bins
                    # for d0 in ['m9','m8','m7', 'm6','m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5','p6','p7','p8','p9']: # |d0| < 0.005 20 bins
                    # for d0 in ['m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5']: #only |d0| < 0.005
                    # for d0 in ['p8']:
                        iDat += 1
                        h_x = in_file.Get('h_d0_%s_%s_%s_d0' % (d0, pt, eta)) #d0
                        h_y = in_file.Get('h_d0_%s_%s' % (d0, g_str)) #dPt
                        # print 'h_d0_%s_%s' % (d0, g_str)

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
                              h_y.Rebin(10)
                            for i in range(7, h_y.GetNbinsX() - 5):  ## Don't include over- or under-flow
                                if (h_y.Integral(i-5, i+5) > yMax):
                                    yMax = h_y.Integral(i-5, i+5)
                                    iMax = i
                            iLo = h_y.GetNbinsX()  ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                            iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                            for i in range(7, h_y.GetNbinsX() - 5):  ## Don't include over- or under-flow
                                if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                                if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)
                        # if (iDat > 29 and iDat < 49):
                        #     for i in range(51, h_y.GetNbinsX() - 49):  ## Don't include over- or under-flow
                        #         if (h_y.Integral(i-50, i+50) > yMax):
                        #             yMax = h_y.Integral(i-50, i+50)
                        #             iMax = i
                        #     iLo = h_y.GetNbinsX()  ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                        #     iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                        #     for i in range(51, h_y.GetNbinsX() - 49):  ## Don't include over- or under-flow
                        #         if (h_y.Integral(i-50, i+50) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                        #         if (h_y.Integral(i-50, i+50) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)
                        # # elif (iDat < 9 and iDat > 69):
                        # #     if h_y.GetNbinsX() == 4000:
                        # #       h_y.Rebin(40)
                        # #     for i in range(6, h_y.GetNbinsX() - 4):  ## Don't include over- or under-flow
                        # #         if (h_y.Integral(i-5, i+5) > yMax):
                        # #             yMax = h_y.Integral(i-5, i+5)
                        # #             iMax = i
                        # #     iLo = h_y.GetNbinsX()  ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                        # #     iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                        # #     for i in range(6, h_y.GetNbinsX() - 4):  ## Don't include over- or under-flow
                        # #         if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                        # #         if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)
                        # else:
                        #     if h_y.GetNbinsX() == 4000:
                        #       h_y.Rebin(40)
                        #     for i in range(6, h_y.GetNbinsX() - 4):  ## Don't include over- or under-flow
                        #         if (h_y.Integral(i-5, i+5) > yMax):
                        #             yMax = h_y.Integral(i-5, i+5)
                        #             iMax = i
                        #     iLo = h_y.GetNbinsX()  ## Lowest bin with integral > yMax - 2*sqrt(yMax)
                        #     iHi = 0                ## Highest bin with integral > yMax - 2*sqrt(yMax)
                        #     for i in range(6, h_y.GetNbinsX() - 4):  ## Don't include over- or under-flow
                        #         if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iLo = min(iLo, i)
                        #         if (h_y.Integral(i-5, i+5) > yMax - 2*math.sqrt(yMax)): iHi = max(iHi, i)

                        

                        data[g_str]['y'][iDat] = h_y.GetBinCenter(iMax)
                        data[g_str]['yerr'][iDat] = (h_y.GetBinLowEdge(iHi+1) - h_y.GetBinLowEdge(iLo)) / 2.
                        # print h_y.GetNbinsX()

                        # # ## Canvas with the original histogram and its peak range
                        # h_can = R.TCanvas(g_str+'_'+d0, g_str+'_'+d0, 800, 600)
                        # h_can.cd()
                        # fit_y = R.TF1('fy_%s' % g_str, 'gaus', min(data[g_str]['y']), max(data[g_str]['y']) )
                        # fit_y.SetParameters(1,0,1)
                        # h_y.Fit('fy_%s' % g_str)
                        # h_y.SetLineColor(R.kBlack)
                        # h_y.Draw()
                        # fit_y.SetLineColor(R.kRed)
                        # fit_y.Draw('same')
                        # # h_y_peak = h_y.Clone()
                        # # h_y_peak.GetXaxis().SetRangeUser(h_y.GetBinLowEdge(iLo), h_y.GetBinLowEdge(iHi+1))
                        # # h_y_peak.SetLineWidth(2)
                        # # h_y_peak.SetLineColor(R.kBlue)
                        # # h_y_peak.Draw('same')
                        # h_can.SaveAs('plots_DY_d0BS_linear_%s/png/histos/%s_%s.png' % (YEAR, g_str, d0))
                        # # h_can.SaveAs('plots/pdf/%s_%s.pdf' % (g_str, d0))


                    ## End loop: for d0 in ['m9','m8','m7','m6','m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5','p6','p7','p8','p9']

                    ## Create a graph with the peak values vs. d0
                    graph[g_str] = R.TGraphErrors(len(data[g_str]['x']), data[g_str]['x'], data[g_str]['y'], data[g_str]['xerr'], data[g_str]['yerr'])
                    graph[g_str].SetName(g_str)

                    
                    # if var == 'dPt':
                    #     var_str = 'RECO - GEN p_{T}'
                    # if var == 'dRelPt':
                    #     var_str = '100 * (RECO - GEN) / GEN p_{T}'
                    # if var == 'dRelPt2p0':
                    #     var_str = '10000 * (RECO - GEN) / GEN p_{T}^{2}'
                    # if var == 'dRelPt1p0':
                    #     var_str = '100 * (RECO - GEN) / GEN p_{T}'
                    # if var == 'dPhi':
                    #     var_str = '(RECO - GEN #phi) * charge'
                    # if var == 'dEta':
                    #     var_str = '(RECO - GEN #eta) * charge'  
                    # if var == 'd0':
                    #     var_str = 'd0'

                    # if var == 'dPt':
                    #     var_str = 'RECO - GEN p_{T}'
                    # if var == 'dRelPt':
                    #     var_str = '100 * (RECO - GEN) / GEN p_{T}'
                    if var == 'dRelPt2p0':
                        var_str = '10000 * (p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}'
                    # if var == 'dRelPt1p0':
                    #     var_str = '100 * (RECO - GEN) / GEN p_{T}'
                    # if var == 'dPhi':
                    #     var_str = '(RECO - GEN #phi) * charge'
                    # if var == 'dEta':
                    #     var_str = '(RECO - GEN #eta) * charge'  
                    # if var == 'd0':
                    #     var_str = 'd0'

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

                    ## Create a linear fit to the points on the graph
                    # fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*x', min(data[g_str]['x']), max(data[g_str]['x']))
                    fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*x', -0.01, 0.01)
                    # fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*x + [2]*tanh([3]*x)', -0.008, 0.008)
                    # fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*tanh([2]*x) + [3]*x', -0.009, 0.009)
                    # fit_left[g_str] = R.TF1('fl_%s' % g_str, '[0] + [1]*x', -0.008, -0.0025)
                    # fit_right[g_str] = R.TF1('fr_%s' % g_str, '[0] + [1]*x', 0.0025, 0.008)
                    fit[g_str].SetParameter(0, 0)
                    # fit_left[g_str].SetParameter(0, 0)
                    # fit_right[g_str].SetParameter(0, 0)
                    # fit[g_str].SetParameter(1, (data[g_str]['y'][13] - data[g_str]['y'][4]) / (data[g_str]['x'][13] - data[g_str]['x'][4]) ) # for 20 bins
                    # fit[g_str].SetParameter(1, (data[g_str]['y'][25] - data[g_str]['y'][15]) / (data[g_str]['x'][25] - data[g_str]['x'][15]) ) # for 40 bins
                    # fit[g_str].SetParameter(1, (data[g_str]['y'][45] - data[g_str]['y'][35]) / (data[g_str]['x'][45] - data[g_str]['x'][35]) ) # for 80 bins
                    # fit_left[g_str].SetParameter(1, (data[g_str]['y'][30] - data[g_str]['y'][5]) / (data[g_str]['x'][30] - data[g_str]['x'][5]) ) # for 80 bins
                    # fit_right[g_str].SetParameter(1, (data[g_str]['y'][75] - data[g_str]['y'][50]) / (data[g_str]['x'][75] - data[g_str]['x'][50]) ) # for 80 bins
                    if eta == 'eta_0_0p9':
                      # fit[g_str].SetParameter(1, 3) # tanh
                      # fit[g_str].SetParameter(2, 500) # tanh
                      fit[g_str].SetParameter(1, 100) # linear
                    elif eta == 'eta_0p9_1p7':
                      # fit[g_str].SetParameter(1, 3) # tanh
                      # fit[g_str].SetParameter(2, 500) # tanh
                      fit[g_str].SetParameter(1, 400) # linear
                    elif eta == 'eta_1p7_inf':
                      # fit[g_str].SetParameter(1, 3) # tanh
                      # fit[g_str].SetParameter(2, 500) # tanh
                      fit[g_str].SetParameter(1, 500) # linear

                    ## Draw the graph onto the canvas
                    # canv[c_str].cd()
                    # if (muon == 'mu1'):
                    # graph[g_str].SetLineColor(R.kBlue)
                    # fit[g_str].SetLineColor(R.kBlue)
                    # graph[g_str].GetYaxis().SetRangeUser(min(data[g_str]['y'])*2.0, max(data[g_str]['y'])*2.0)
                    # graph[g_str].Fit('f_%s' % g_str)
                    # graph[g_str].Draw()
                    else:
                        graph[g_str].SetLineColor(R.kRed)
                        fit[g_str].SetLineColor(R.kRed)
                        graph[g_str].Fit('f_%s' % g_str)
                        # graph[g_str].Draw('same')
                    graph[g_str].SetLineColor(iPt+1)
                    graph[g_str].SetMarkerColor(iPt+1)
                    graph[g_str].SetMarkerStyle(8)
                    graph[g_str].SetMarkerSize(0.7)
                    fit[g_str].SetLineColor(iPt+1)
                    # fit_left[g_str].SetLineColor(iPt+1)
                    # fit_right[g_str].SetLineColor(iPt+1)
                    graph[g_str].Fit('f_%s' % g_str,"R")
                    # graph[g_str].Fit('fl_%s' % g_str,"R+")
                    # graph[g_str].Fit('fr_%s' % g_str,"R+")
                    print 'Writing: %s' % g_str
                    out_file.cd()
                    graph[g_str].Write()

                    g_canv = R.TCanvas(g_str, g_str, 1600, 1200)
                    g_canv.cd()
                    graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                    graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                    # graph[g_str].GetYaxis().SetRangeUser(-20, 20)
                    graph[g_str].Draw('AP')
                    g_leg = R.TLegend(0.12,0.76,0.72,0.88)
                    # f_str = g_str.replace(c_str, '')
                    # print f_str
                    # g_leg.AddEntry( fit[g_str], '%s = %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (g_str.replace(('_' + c_str), '').replace('_','-').replace('pt-','p_T ='), fit[g_str].GetParameter(0), fit[g_str].GetParError(0), fit[g_str].GetParameter(1), fit[g_str].GetParError(1)) )
                    # g_leg.AddEntry( fit_left[g_str], '%s = %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (g_str.replace(('_' + c_str), '').replace('_','-').replace('pt-','p_T ='), fit_left[g_str].GetParameter(0), fit_left[g_str].GetParError(0), fit_left[g_str].GetParameter(1), fit_left[g_str].GetParError(1)) )
                    # g_leg.AddEntry( fit_right[g_str], '%s = %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (g_str.replace(('_' + c_str), '').replace('_','-').replace('pt-','p_T ='), fit_right[g_str].GetParameter(0), fit_right[g_str].GetParError(0), fit_right[g_str].GetParameter(1), fit_right[g_str].GetParError(1)) )
                    # g_leg.AddEntry( fit_left[g_str], 'left = %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (fit_left[g_str].GetParameter(0), fit_left[g_str].GetParError(0), fit_left[g_str].GetParameter(1), fit_left[g_str].GetParError(1)) )
                    # g_leg.AddEntry( fit[g_str], 'middle = %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (fit[g_str].GetParameter(0), fit[g_str].GetParError(0), fit[g_str].GetParameter(1), fit[g_str].GetParError(1)) )
                    # g_leg.AddEntry( fit_right[g_str], 'right = %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (fit_right[g_str].GetParameter(0), fit_right[g_str].GetParError(0), fit_right[g_str].GetParameter(1), fit_right[g_str].GetParError(1)) )
                    # g_leg.AddEntry( fit[g_str], '\LARGE{y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)}' % (fit[g_str].GetParameter(0), fit[g_str].GetParameter(3), fit[g_str].GetParameter(1), fit[g_str].GetParameter(2)) )
                    g_leg.AddEntry( fit[g_str], '\LARGE{y = %.2f + %.2f *d0}' % (fit[g_str].GetParameter(0), fit[g_str].GetParameter(1)) )
                    g_leg.Draw('same')
                    g_latex = R.TLatex()
                    g_latex.SetTextAlign(12)
                    g_latex.SetTextSize(0.04)
                    g_latex.DrawLatex(-0.01, 8, eta_str)
                    g_canv.SaveAs('plots_DY_d0BS_linear_%s/png/histos/%s.png' % (YEAR, g_str))
                    
                    canv[c_str].cd()
                    if (iPt == 0):
                        # graph[g_str].GetXaxis().SetRangeUser(-0.06, 0.06)
                        graph[g_str].GetXaxis().SetLimits(-0.011, 0.011)
                        # graph[g_str].GetYaxis().SetRangeUser(min(data[g_str]['y'])*1.2, max(data[g_str]['y'])*1.2)
                        # graph[g_str].GetYaxis().SetRangeUser(-20, 20)
                        graph[g_str].GetYaxis().SetRangeUser(-15, 15)
                        # graph[g_str].GetYaxis().SetLimits(-5, 5)
                        graph[g_str].Draw('AP')
                    else:
                        graph[g_str].Draw('Psame')

                    latex = R.TLatex()
                    latex.SetTextAlign(12)
                    latex.SetTextSize(0.04)
                    latex.DrawLatex(-0.01, 8, eta_str)
                    
                    


                # End loop for pt

                # f1 = fit[g_str+'_mu1']
                # f2 = fit[g_str+'_mu2']
                # f1 = fit[c_str]
                # f5 = fit['inclusive_'+c_str]
                # f1 = fit['nVtx_0_21_'+c_str]
                # f2 = fit['nVtx_22_26_'+c_str]
                # f3 = fit['nVtx_27_33_'+c_str]
                # f4 = fit['nVtx_34_inf_'+c_str]

                # par_file = R.TFile('plots_DY_d0BS_linear_%s/fit_param.txt' % YEAR, 'recreate')
                # par_file



                if (var == 'dPhi' or var == 'dEta'): leg = R.TLegend(0.48,0.66,0.88,0.88)
                else:               leg = R.TLegend(0.12,0.66,0.52,0.88)
                # leg.AddEntry(graph[c_str+'_mu1'], 'Higher pT muon')
                # leg.AddEntry(graph[c_str+'_mu2'], 'Lower pT muon')
                # leg.AddEntry( f1, '%.2f \pm %.2f + d0 * %.0f \pm %.0f' % (f1.GetParameter(0), f1.GetParError(0), f1.GetParameter(1), f1.GetParError(1)) )
                # leg.AddEntry( f2, '%.2f \pm %.2f + d0 * %.0f \pm %.0f' % (f2.GetParameter(0), f2.GetParError(0), f2.GetParameter(1), f2.GetParError(1)) )
                # leg.AddEntry( f5, 'inclusive: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f5.GetParameter(0), f5.GetParameter(3), f5.GetParameter(1), f5.GetParameter(2)) )
                # leg.AddEntry( f1, '0 <= nVtx <= 21: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f1.GetParameter(0), f1.GetParameter(3), f1.GetParameter(1), f1.GetParameter(2)) )
                # leg.AddEntry( f2, '22 <= nVtx <= 26: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f2.GetParameter(0), f2.GetParameter(3), f2.GetParameter(1), f2.GetParameter(2)) )
                # leg.AddEntry( f3, '27 <= nVtx <= 33: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f3.GetParameter(0), f3.GetParameter(3), f3.GetParameter(1), f3.GetParameter(2)) )
                # leg.AddEntry( f4, '34 <= nVtx: y = %.2f + %.2f *d0 + %.2f * tanh(d0* %.2f)' % (f4.GetParameter(0), f4.GetParameter(3), f4.GetParameter(1), f4.GetParameter(2)) )

                # leg.AddEntry( f5, 'inclusive: %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (f5.GetParameter(0), f5.GetParError(0), f5.GetParameter(1), f5.GetParError(1)) )
                # leg.AddEntry( f1, '0 <= nVtx <= 21: %.2f \pm %.2f + d0 * %.2f \pm %.2f' % (f1.GetParameter(0), f1.GetParError(0), f1.GetParameter(1), f1.GetParError(1)) )
                # leg.AddEntry( f2, '22 <= nVtx <= 26:%.2f \pm %.2f + d0 * %.2f \pm %.2f' % (f2.GetParameter(0), f2.GetParError(0), f2.GetParameter(1), f2.GetParError(1)) )
                # leg.AddEntry( f3, '27 <= nVtx <= 33:%.2f \pm %.2f + d0 * %.2f \pm %.2f' % (f3.GetParameter(0), f3.GetParError(0), f3.GetParameter(1), f3.GetParError(1)) )
                # leg.AddEntry( f4, '34 <= nVtx:%.2f \pm %.2f + d0 * %.2f \pm %.2f' % (f4.GetParameter(0), f4.GetParError(0), f4.GetParameter(1), f4.GetParError(1)) )
                # leg.Draw('same')

                canv[c_str].SaveAs('plots_DY_d0BS_linear_%s/png/%s.png' % (YEAR, c_str))
                # canv[c_str].SaveAs('plots/pdf/%s.pdf' % c_str)
                # canv[c_str].Write()
            #End loop for eta
        ## End loop: for corr in ['PF', 'Roch']:
    ## End loop: for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:


    ## Write histograms to output root file
    # out_file.Write()

    ## Close output and input files
    print 'Finished, closing files ...'
    out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
