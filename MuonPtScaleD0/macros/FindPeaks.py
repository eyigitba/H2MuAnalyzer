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




INDIR  = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms/GenRecoPtDiffVsD0_Aug24_v1/files'
INFILE = 'histos_ZJets_AMC__NONE.root'

## Main function executed by ./batch/FindPeaks.py
def main():

    print '\n\n*** Inside FindPeakss.py ***'

    ## Import root file with dPt and dPhi distributions
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    if not os.path.exists('plots'):
        print '\nCreating output folders plots/png and plots/pdf'
        os.makedirs('plots/png')
        os.makedirs('plots/pdf')

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
    for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:  ## Variables plotted vs. RECO d0 * charge
        for corr in ['PF', 'Roch']:

            ## Unique string identifying this canvas
            if (var == 'dPhi'): c_str = var
            else: c_str = '%s_%s' % (var, corr)

            if (var == 'dPhi' and corr == 'Roch'): continue

            canv[c_str] = R.TCanvas(c_str, c_str, 800, 600)
            canv[c_str].cd()

            for muon in ['mu1', 'mu2']:

                ## Unique string identifying this graph
                g_str = c_str+'_'+muon

                data[g_str] = {}
                data[g_str]['x']    = numpy.empty([19], dtype=float)
                data[g_str]['y']    = numpy.empty([19], dtype=float)
                data[g_str]['xerr'] = numpy.empty([19], dtype=float)
                data[g_str]['yerr'] = numpy.empty([19], dtype=float)

                print 'Getting data for %s' % g_str
                iDat = -1
                for d0 in ['m9','m8','m7','m6','m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5','p6','p7','p8','p9']:
                    iDat += 1

                    h_x = in_file.Get('h_d0_%s_d0_%s' % (d0, muon))
                    h_y = in_file.Get('h_d0_%s_%s' % (d0, g_str))
                    
                    data[g_str]['x'][iDat] = h_x.GetMean()
                    data[g_str]['xerr'][iDat] = h_x.GetStdDev()

                    yMax = 0 ## Maximum integral over 11 bins
                    iMax = 0 ## Central bin of maximum integral (i.e. the center of the peak)
                    for i in range(6, h_y.GetNbinsX() - 4):
                        if (h_y.Integral(i-5, i+5) > yMax):
                            yMax = h_y.Integral(i-5, i+5)
                            iMax = i
                    iLo = h_y.GetNbinsX()  ## Lowest bin with integral > yMax - sqrt(yMax)
                    iHi = 0                ## Highest bin with integral > yMax - sqrt(yMax)
                    for i in range(1, h_y.GetNbinsX() + 1):
                        if (h_y.Integral(i-5, i+5) > yMax - math.sqrt(yMax)): iLo = min(iLo, i)
                        if (h_y.Integral(i-5, i+5) > yMax - math.sqrt(yMax)): iHi = max(iHi, i)

                    data[g_str]['y'][iDat] = h_y.GetBinCenter(iMax)
                    data[g_str]['yerr'][iDat] = (h_y.GetBinCenter(iHi) - h_y.GetBinCenter(iLo)) / 2.

                ## End loop: for d0 in ['m9','m8','m7','m6','m5','m4','m3','m2','m1','0','p1','p2','p3','p4','p5','p6','p7','p8','p9']

                ## Create a graph with the peak values vs. d0
                graph[g_str] = R.TGraphErrors(len(data[g_str]['x']), data[g_str]['x'], data[g_str]['y'], data[g_str]['xerr'], data[g_str]['yerr'])

                if var == 'dPt':
                    var_str = 'RECO - GEN p_{T}'
                if var == 'dRelPt':
                    var_str = '100 * (RECO - GEN) / GEN p_{T}'
                if var == 'dRelPt2':
                    var_str = '10000 * (RECO - GEN) / GEN p_{T}^{2}'
                if var == 'dPhi':
                    var_str = '(RECO - GEN #phi) * charge'

                if corr == 'PF':
                    corr_str = '(PF)'
                if corr == 'Roch':
                    corr_str = '(Rochester)'

                graph[g_str].SetTitle('%s %s vs. d0' % (var_str, corr_str))
                graph[g_str].GetXaxis().SetTitle('RECO muon d0 (cm) * charge')
                graph[g_str].GetYaxis().SetTitle('%s %s' % (var_str, corr_str))

                ## Create a linear fit to the points on the graph
                fit[g_str] = R.TF1('f_%s' % g_str, '[0] + [1]*x', min(data[g_str]['x']), max(data[g_str]['x']))
                fit[g_str].SetParameter(0, 0)
                fit[g_str].SetParameter(1, (data[g_str]['y'][13] - data[g_str]['y'][4]) / (data[g_str]['x'][13] - data[g_str]['x'][4]) )

                if (muon == 'mu1'):
                    graph[g_str].SetLineColor(R.kBlue)
                    fit[g_str].SetLineColor(R.kBlue)
                    graph[g_str].GetYaxis().SetRangeUser(min(data[g_str]['y'])*2.0, max(data[g_str]['y'])*2.0)
                    graph[g_str].Fit('f_%s' % g_str)
                    graph[g_str].Draw()
                else:
                    graph[g_str].SetLineColor(R.kRed)
                    fit[g_str].SetLineColor(R.kRed)
                    graph[g_str].Fit('f_%s' % g_str)
                    graph[g_str].Draw('same')

            ## End loop: for muon in ['mu1', 'mu2']:

            f1 = fit[c_str+'_mu1']
            f2 = fit[c_str+'_mu2']

            if (var == 'dPhi'): leg = R.TLegend(0.48,0.66,0.88,0.88)
            else:               leg = R.TLegend(0.12,0.66,0.52,0.88)
            leg.AddEntry(graph[c_str+'_mu1'], 'Higher pT muon')
            leg.AddEntry(graph[c_str+'_mu2'], 'Lower pT muon')
            leg.AddEntry( f1, '%.3f \pm %.3f + d0 * %.0f \pm %.0f' % (f1.GetParameter(0), f1.GetParError(0), f1.GetParameter(1), f1.GetParError(1)) )
            leg.AddEntry( f2, '%.3f \pm %.3f + d0 * %.0f \pm %.0f' % (f2.GetParameter(0), f2.GetParError(0), f2.GetParameter(1), f2.GetParError(1)) )
            leg.Draw('same')

            canv[c_str].SaveAs('plots/png/%s.png' % c_str)
            canv[c_str].SaveAs('plots/pdf/%s.pdf' % c_str)

        ## End loop: for corr in ['PF', 'Roch']:
    ## End loop: for var in ['dPt', 'dRelPt', 'dRelPt2', 'dPhi']:


    # ## Write histograms to output root file
    # print '\nCreating output root file plots/FindPeaks.root'
    # out_file.cd()
    # out_file.Write()

    ## Close output and input files
    print 'Finished, closing files ...'
    # out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
