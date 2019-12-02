#! /usr/bin/env python

###############################################
###             FindPeaks2D.py              ###
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


# INDIR  = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/GenRecoPtDiffVsD0VsPt_2016_Sep01_v1/files'
INDIR  = '/afs/cern.ch/work/e/eyigitba/public/H2Mu/2016/Histograms/'
INFILE = 'histos_ZJets_AMC__NONE_1_20.root'
# INDIR  = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms/GenRecoPtDiffVsD0VsPt_Aug29_v2/files'
# INFILE = 'histos_ZJets_AMC__NONE.root'

D0_BINS = 20

## Main function executed by ./macros/FindPeaks2D.py
def main():

    print '\n\n*** Inside FindPeaks2Ds.py ***'

    ## Import root file with dPt and dPhi/dEta distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    if not os.path.exists('plots_2D'):
        print '\nCreating output folders plots_2D/png and plots_2D/png/histos'
        os.makedirs('plots_2D/png')
        os.makedirs('plots_2D/png/histos')
        # os.makedirs('plots_2D/pdf')

    print '\nCreating output root file plots_2D/FindPeaks2D.root'
    out_file = R.TFile('plots_2D/FindPeaks2D.root', 'recreate')
    out_file.cd()

    ## Input histograms and canvases to draw on
    hist = {}
    ## Graphs of delta(RECO, GEN) vs. RECO muon d0 * charge
    graph = {}
    ## Data to fill the graph
    data = {}
    ## Fits to the data
    fit = {}
    ## Canvases to draw the graph on
    canv = {}

    # dPt     =        RECO - GEN pT
    # dRelPt  =   100*(RECO - GEN) / GEN pT
    # dRelPt2 = 10000*(RECO - GEN) / (GEN pT)^2 
    # dPhi    = (RECO - GEN phi) * GEN charge
    # for var in ['dPt', 'dRelPt', 'dRelPt1p4', 'dRelPt1p6', 'dRelPt1p8', 'dRelPt2p0', 'dRelPt2p2', 'dPhi', 'dEta']:  ## Variables plotted vs. RECO d0 * charge
    # for var in ['dRelPt1p6', 'dRelPt1p8', 'dRelPt2p0', 'dRelPt2p2', 'dRelPt2p4', 'dPhi', 'dEta']:
    for var in ['dRelPt2p0']:
    # for var in ['dEta']:
        # for corr in ['PF', 'Roch']:
        for corr in ['Roch', 'RochKin']:

            ## Unique string identifying this canvas
            if (var == 'dPhi' or var == 'dEta'): c_str = var
            else: c_str = '%s_%s' % (var, corr)

            if ((var == 'dPhi' or var == 'dEta') and corr == 'Roch'): continue

            ## Canvas for the graphs of peak vs. d0
            canv[c_str] = R.TCanvas(c_str, c_str, 800, 600)

            iPt = -1
            # for pt in ['pt_20_35']:
            for pt in ['pt_20_35', 'pt_35_42', 'pt_42_50', 'pt_50_inf', 'inclusive']:
                iPt += 1

                ## Unique string identifying this graph
                g_str = pt+'_'+c_str

                data[g_str] = {}
                data[g_str]['x']    = numpy.empty([D0_BINS], dtype=float)
                data[g_str]['y']    = numpy.empty([D0_BINS], dtype=float)
                data[g_str]['xerr'] = numpy.empty([D0_BINS], dtype=float)
                data[g_str]['yerr'] = numpy.empty([D0_BINS], dtype=float)

                print 'Getting data for %s' % g_str
                h_2D = in_file.Get('h_%s_vs_d0' % g_str)

                binsX = h_2D.GetNbinsX() ## d0 axis
                binsY = h_2D.GetNbinsY() ## Variable axis (dPt, dRelPt, etc.) 
                integ = h_2D.Integral(1, binsX, 1, binsY)

                ## High and low bins covering 1/D0_BINS of the data
                iX = [0, 0]

                iDat = -1
                sumInteg = 0 
                for iD0 in range(D0_BINS):
                    iDat += 1

                    iInteg = 0  ## Integral for this range
                    iX = [iX[1]+1, iX[1]]  ## New "low" bin = old "high" bin + 1
                    while ( iInteg < (integ / (D0_BINS+1)) and iX[1] < binsX):
                        iX = [iX[0], iX[1]+1]  ## Extend the "high" bin by +1
                        iInteg = h_2D.Integral(iX[0], iX[1], 1, binsY)
                    sumInteg += iInteg

                    ## Get the mean and std. dev. for d0 in the bin range
                    h_2D.GetXaxis().SetRange(iX[0], iX[1])
                    data[g_str]['x'][iDat] = h_2D.GetMean(1)
                    data[g_str]['xerr'][iDat] = h_2D.GetStdDev(1)
                    h_2D.GetXaxis().SetRange()

                    # print '  * D0 range #%d: from bins %d through %d, integral = %d (sum = %d/%d)' % (iD0+1, iX[0], iX[1], iInteg, sumInteg, integ)
                    # print '    - Mean = %.6f, std. dev. = %.6f' % (data[g_str]['x'][iDat], data[g_str]['xerr'][iDat])

                    pWid = int(binsY / 20) ## Number of bins summed on each side of the peak 
                    yMax = 0 ## Maximum integral over 2*pWid + 1 bins
                    iMax = 0 ## Central bin of maximum integral (i.e. the center of the peak)
                    for i in range(pWid+1, binsY - (pWid-1)):  ## Don't include over- or under-flow
                        if (h_2D.Integral(iX[0], iX[1], i-pWid, i+pWid) > yMax):
                            yMax = h_2D.Integral(iX[0], iX[1], i-pWid, i+pWid)
                            iMax = i
                    iLo = binsY  ## Lowest bin with integral  > yMax - sqrt(yMax)
                    iHi = 0      ## Highest bin with integral > yMax - sqrt(yMax)
                    for i in range(pWid+1, binsY - (pWid-1)):  ## Don't include over- or under-flow
                        if (h_2D.Integral(iX[0], iX[1], i-pWid, i+pWid) > yMax - 1.5*math.sqrt(yMax)): iLo = min(iLo, i)
                        if (h_2D.Integral(iX[0], iX[1], i-pWid, i+pWid) > yMax - 1.5*math.sqrt(yMax)): iHi = max(iHi, i)

                    data[g_str]['y'][iDat] = h_2D.GetYaxis().GetBinCenter(iMax)
                    data[g_str]['yerr'][iDat] = (h_2D.GetYaxis().GetBinLowEdge(iHi+1) - h_2D.GetYaxis().GetBinLowEdge(iLo)) / 2.

                    ## Canvas with the original histogram and its peak range
                    h_can = R.TCanvas('%s_d0_%s' % (g_str, iD0), '%s_d0_%s' % (g_str, iD0), 800, 600)
                    h_can.cd()
                    h_y = h_2D.ProjectionY('h_y', iX[0], iX[1])
                    h_y.SetLineColor(R.kBlack)
                    h_y.Draw()
                    h_y_peak = h_y.Clone()
                    h_y_peak.GetXaxis().SetRangeUser(h_y.GetBinLowEdge(iLo), h_y.GetBinLowEdge(iHi+1))
                    h_y_peak.SetLineWidth(2)
                    h_y_peak.SetLineColor(R.kBlue)
                    h_y_peak.Draw('same')
                    h_can.SaveAs('plots_2D/png/histos/%s_d0_%s.png' % (g_str, iD0))
                    # h_can.SaveAs('plots_2D/pdf/histos/%s_d0_%s.pdf' % (g_str, iD0))
                    h_can.Write()


                ## End loop: for iD0 in range(D0_BINS):

                ## Create a graph with the peak values vs. d0
                graph[g_str] = R.TGraphErrors(len(data[g_str]['x']), data[g_str]['x'], data[g_str]['y'], data[g_str]['xerr'], data[g_str]['yerr'])

                if var == 'dPt':
                    var_str = 'RECO - GEN p_{T}'
                if var == 'dRelPt':
                    var_str = '100 * (RECO - GEN) / GEN p_{T}'
                if var == 'dRelPt1p4':
                    var_str = '10000 * (RECO - GEN) / GEN p_{T}^{1.4}'
                if var == 'dRelPt1p6':
                    var_str = '10000 * (RECO - GEN) / GEN p_{T}^{1.6}'
                if var == 'dRelPt1p8':
                    var_str = '10000 * (RECO - GEN) / GEN p_{T}^{1.8}'
                if var == 'dRelPt2p0':
                    var_str = '10000 * (RECO - GEN) / GEN p_{T}^{2.0}'
                if var == 'dRelPt2p2':
                    var_str = '10000 * (RECO - GEN) / GEN p_{T}^{2.2}'
                if var == 'dPhi':
                    var_str = '(RECO - GEN #phi) * charge'
                if var == 'dEta':
                    var_str = '(RECO - GEN #eta) * charge'

                if corr == 'PF':
                    corr_str = '(PF)'
                if corr == 'Roch':
                    corr_str = '(Rochester)'

                graph[g_str].SetTitle('%s %s vs. d0' % (var_str, corr_str))
                graph[g_str].GetXaxis().SetTitle('RECO muon d0 (#mum) * charge')
                graph[g_str].GetYaxis().SetTitle('%s %s' % (var_str, corr_str))

                ## Create a linear fit to the points on the graph
                fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*x', min(data[g_str]['x']), max(data[g_str]['x']))
                fit[g_str].SetParameter(0, 0)
                fit[g_str].SetParameter(1, (data[g_str]['y'][3*D0_BINS/4] - data[g_str]['y'][D0_BINS/4]) / (data[g_str]['x'][3*D0_BINS/4] - data[g_str]['x'][D0_BINS/4]) )

                ## Draw the graph onto the canvas
                canv[c_str].cd()
                
                graph[g_str].SetLineColor(iPt+1)
                graph[g_str].SetMarkerColor(iPt+1)
                graph[g_str].SetMarkerStyle(8)
                graph[g_str].SetMarkerSize(0.7)
                fit[g_str].SetLineColor(iPt+1)
                graph[g_str].Fit('f_%s' % g_str)
                if (iPt == 0):
                    graph[g_str].GetXaxis().SetRangeUser(-20, 20)
                    # graph[g_str].GetYaxis().SetRangeUser(min(data[g_str]['y'])*1.2, max(data[g_str]['y'])*1.2)
                    graph[g_str].GetYaxis().SetRangeUser(-5, 5)
                    graph[g_str].Draw('AP')
                else:
                    graph[g_str].Draw('Psame')

            ## End loop: for pt in ['pt_20_35', 'pt_35_42', 'pt_42_50', 'pt_50_inf']:

            # f1 = fit[c_str+'_mu1']
            # f2 = fit[c_str+'_mu2']
 
            f1 = fit['pt_20_35_'+c_str]
            f2 = fit['pt_35_42_'+c_str]
            f3 = fit['pt_42_50_'+c_str]
            f4 = fit['pt_50_inf_'+c_str]

            if (var == 'dPhi'): leg = R.TLegend(0.48,0.66,0.88,0.88)
            else:               leg = R.TLegend(0.12,0.66,0.52,0.88)
            # leg.AddEntry(graph['pt_20_35_'+c_str],  '20 < p_{T} < 35 GeV')
            # leg.AddEntry(graph['pt_35_42_'+c_str],  '35 < p_{T} < 42.5 GeV')
            # leg.AddEntry(graph['pt_42_50_'+c_str],  '42.5 < p_{T} < 50 GeV')
            # leg.AddEntry(graph['pt_50_inf_'+c_str], '50 < p_{T} GeV')
            leg.AddEntry( f1, '%.4f \pm %.4f + d0 * %.4f \pm %.4f' % (f1.GetParameter(0), f1.GetParError(0), f1.GetParameter(1), f1.GetParError(1)) )
            leg.AddEntry( f2, '%.4f \pm %.4f + d0 * %.4f \pm %.4f' % (f2.GetParameter(0), f2.GetParError(0), f2.GetParameter(1), f2.GetParError(1)) )
            leg.AddEntry( f3, '%.4f \pm %.4f + d0 * %.4f \pm %.4f' % (f3.GetParameter(0), f3.GetParError(0), f3.GetParameter(1), f3.GetParError(1)) )
            leg.AddEntry( f4, '%.4f \pm %.4f + d0 * %.4f \pm %.4f' % (f4.GetParameter(0), f4.GetParError(0), f4.GetParameter(1), f4.GetParError(1)) )
            leg.Draw('same')

            canv[c_str].SaveAs('plots_2D/png/%s.png' % c_str)
            # canv[c_str].SaveAs('plots_2D/pdf/%s.pdf' % c_str)
            canv[c_str].Write()

        ## End loop: for corr in ['PF', 'Roch']:
    ## End loop: for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:


    ## Write histograms to output root file
    out_file.Write()

    ## Close output and input files
    print 'Finished, closing files ...'
    out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
