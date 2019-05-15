#! /usr/bin/env python 

######################################################
##                  AnalyzeToys.py                  ##
##                                                  ##
##  Analyzes workspace and FitDiagnostics output,   ##
##  including results from different toy datasets,  ##
##  producing best-fit, limit, and significance     ##
##  numbers as well as various plots.               ##
######################################################

#============================================
# Imports
#============================================

import array
import math
import ROOT as R

# LABELS = ['syst_0p3_toys_fit',
#           'syst_0p3_fit',
#           'noSyst']
# LABELS = ['T400_syst_0p3_fit_autoMCStats500',
#           'T400_syst_0p3_fit_autoMCStats0',
#           'T400_syst_0p3_fit']
# LABELS = ['test_BDT_mass_rebin_v1']
LABELS = ['test_XWZ_new_script_v1']

## Print additional debugging info
VERBOSE = False
## Small offset for edges of histograms
BIT = 0.001


def AnalyzeToys(label):

    print '\nInside AnalyzeToys.py for input %s' % label

    ## Don't print plots to screen while running (faster)
    R.gROOT.SetBatch(R.kTRUE)

    ## Input files
    in_FD = 'fitDiagnostics_%s.root' % label
    print '\nOpening file %s' % in_FD
    fitDiag   = R.TFile(in_FD)
    in_WS = 'higgsCombine_%s.FitDiagnostics.mH120.123456.root' % label
    print '\nOpening file %s' % in_WS
    workspace = R.TFile(in_WS)

    ## Fit diagnostic trees
    fit_tree = {}
    if VERBOSE: print '\nGetting tree_prefit from %s' % in_FD
    fit_tree['preFit']  = fitDiag.Get('tree_prefit')
    if VERBOSE: print '\nGetting tree_fit_b from %s' % in_FD
    fit_tree['bkgOnly'] = fitDiag.Get('tree_fit_b')
    if VERBOSE: print '\nGetting tree_fit_sb from %s' % in_FD
    fit_tree['sigBkg']  = fitDiag.Get('tree_fit_sb')
    ## Workspace tree
    ws_tree = workspace.Get('limit')

    ## Output file
    if VERBOSE: print '\nCreating output file AnalyzerToys.root\n'
    out_file = R.TFile('AnalyzerToys_%s.root' % label, 'recreate')


    ## ********************************* ##
    ##  Set binning and book histograms  ##
    ## ********************************* ##

    ## Have to get discriminant binning directly from a toy histogram
    ws_tree.GetEntry(0)
    h_tmp    = R.RooAbsData.createHistogram(workspace.Get('toys/toy_1'), 'CMS_th1x')
    bin_disc = [h_tmp.GetNbinsX(), h_tmp.GetBinLowEdge(1), h_tmp.GetBinLowEdge(h_tmp.GetNbinsX()+1)]
    
    bin_mu   = [200, -50,  50]
    bin_pull = [100,  -5,   5]
    bin_bkg  = [100,   0, 200]
    if 'BDT' in label:
        bin_ROI = [20,   0,  20]
    else:
        bin_ROI = [100, 0, 100]

    h_mu   = {}  ## Histogram of best-fit mu values
    h_pull = {}  ## Histogram of mu / muErr values

    for fit in fit_tree.keys():
        if VERBOSE: print '  * Booking histogram h_mu_%s' % fit
        h_mu  [fit] = R.TH1F('h_mu_%s'   % fit, 'h_mu_%s' % fit,   bin_mu[0],   bin_mu[1],   bin_mu[2])
        if VERBOSE: print '  * Booking histogram h_pull_%s' % fit
        h_pull[fit] = R.TH1F('h_pull_%s' % fit, 'h_pull_%s' % fit, bin_pull[0], bin_pull[1], bin_pull[2])

    ## Histogram of background yields
    if VERBOSE: print '  * Booking histogram h_bkg_toys'
    h_bkg_toys = R.TH1F('h_bkg_toys', 'h_bkg_toys', bin_bkg[0], bin_bkg[1], bin_bkg[2])
    ## Histogram of background yields in ROI (123 - 127 GeV or BDT in the last two bins)
    if VERBOSE: print '  * Booking histogram h_ROI_toys'
    h_ROI_toys = R.TH1F('h_ROI_toys', 'h_ROI_toys', bin_ROI[0], bin_ROI[1], bin_ROI[2])
    ## Histograms of background yields in each bin (110 - 160 GeV or BDT in the whole range)
    if VERBOSE: print '  * Booking histograms h_bins[50]'
    h_bins = []
    null = array.array('d')
    for iBin in range(bin_disc[0]):
        h_bins.append( R.TH1F('h_bin_%d' % iBin, 'h_bin_%d' % iBin, bin_ROI[0], bin_ROI[1], bin_ROI[2]) )
        null.append(0)

    ## Histograms of bin quantiles (median, 68% and 95%)
    ## Following: https://wiki.physik.uzh.ch/cms/limits:brazilianplotexample
    g_med_toys = R.TGraph(  bin_disc[0], null, null)
    g_q68_toys = R.TGraph(2*bin_disc[0], null, null)
    g_q95_toys = R.TGraph(2*bin_disc[0], null, null)


    ## ***************************************** ##
    ##  Extract information from FitDiagnostics  ##
    ## ***************************************** ##

    print '\nLooping over pre-fit, background-only, and signal+background fits'
    for fit in fit_tree.keys():

        ## Pre-fit does not contain signal strength errors
        if (fit == 'preFit'):  continue
        ## Background-only does not produce useful signal strength numbers
        if (fit == 'bkgOnly'): continue

        print 'Now looping over events in %s tree, each corresponding to a toy' % fit
        for iEvt in range(fit_tree[fit].GetEntries()):

            if VERBOSE: print '\nGetting entry %d in %s' % (iEvt, fit)
            fit_tree[fit].GetEntry(iEvt)
            
            ## Check the best-fit uncertainties
            fit_err    = fit_tree[fit].rErr
            fit_err_lo = fit_tree[fit].rLoErr
            fit_err_hi = fit_tree[fit].rHiErr

            ## Skip events with unphysical best-fit uncertainties
            if (fit_err_lo < 0.5 or fit_err_hi < 0.5 or fit_err_lo > 50 or fit_err_hi > 50):
                # print '\nIn event #%d, rErr = %.3f, rLoErr = %.3f, rHiErr = %.3f\n SKIPPING!!!\n' % (iEvt, fit_err, fit_err_lo, fit_err_hi)
                continue
            else:
                best_fit    = fit_tree[fit].r
                fit_err_avg = math.sqrt( 0.5*pow(fit_err_lo, 2) + 0.5*pow(fit_err_hi, 2) )

            ## Fill plots
            if VERBOSE: print '  * Filling histograms'
            h_mu  [fit].Fill( min( max(best_fit,               bin_mu[1]+BIT),   bin_mu[2]-BIT) )
            h_pull[fit].Fill( min( max(best_fit / fit_err_avg, bin_pull[1]+BIT), bin_pull[2]-BIT) )

        ## End loop: for iEvt in range(fit_tree[fit].GetEntries()):
    ## End loop:  for fit in fit_tree.keys():


    ## ************************************ ##
    ##  Extract information from Workspace  ##
    ## ************************************ ##

    ## Get the last entry in the "limit" tree
    ws_tree.GetEntry( ws_tree.GetEntries() - 1 )
    ## Number of toys equals last "iToy" in the tree
    nToys = ws_tree.iToy
    ## Create list of histograms extracted from toys
    h_toys = []

    print 'Looping over %d toys in workspace' % nToys
    for iToy in range(nToys):
        if VERBOSE: print '  * Extracting data for toy %d' % (iToy+1)
        ## PyROOT bug doesn't allow normal call to createHistogram returning a TH1 object:
        ## h_toys.append( workspace.Get('toys/toy_%d' % (iToy+1)).createHistogram('CMS_th1x') )
        ## Solution: https://root-forum.cern.ch/t/cant-access-th1f-constructors-of-a-roodataset-using-pyroot-only-the-th2f-constructors-work/26657
        h_toys.append( R.RooAbsData.createHistogram( workspace.Get('toys/toy_%d' % (iToy+1)), 'CMS_th1x') )
        h_toys[iToy].SetName('h_toy_%d' % (iToy+1))

        if VERBOSE: print '  * Filling histograms for toy %d' % (iToy+1)
        h_bkg_toys.Fill( min( max(h_toys[iToy].Integral(),       bin_bkg[1]+BIT), bin_bkg[2]-BIT) )

        ## Fill the yields in each discriminant bin (mass from 110 to 160, BDT from -1 to +1)
        for iBin in range(bin_disc[0]):
            h_bins[iBin].Fill( h_toys[iToy].GetBinContent(iBin+1) )
        if 'BDT' in label:
            ## Fill the yields for the last two BDT bins
            h_ROI_toys.Fill( min( max(h_toys[iToy].Integral(bin_disc[0]-1, bin_disc[0]), bin_ROI[1]+BIT), bin_ROI[2]-BIT) )
        else:
            ## Fill the yields for 123 - 127 GeV
            h_ROI_toys.Fill( min( max(h_toys[iToy].Integral(13, 16), bin_ROI[1]+BIT), bin_ROI[2]-BIT) )

    ## End loop: for iToy in range(nToys):
    

    ## ************************************************** ##
    ##  Fit and write out histograms from FitDiagnostics  ##
    ## ************************************************** ##

    print '\n\nWriting output file'
    out_file.cd()

    ## Loop over pre-fit, background-only, and signal+background fits
    for fit in fit_tree.keys():

        ## Only plot, fit, and write histograms with > 0 entries
        if (h_mu[fit].Integral() == 0): continue

        ## Set titles and axis labels
        h_mu[fit].SetTitle('Best-fit #mu with %d toys for %s' % (h_mu[fit].Integral(1, bin_mu[0]), label))
        h_mu[fit].GetXaxis().SetTitle('Best-fit signal strength #mu')
        h_mu[fit].GetYaxis().SetTitle('Number of toys')
        h_mu[fit].SetMarkerStyle(20)
        h_mu[fit].SetMarkerColor(R.kBlack)

        h_pull[fit].SetTitle('Pull value with %d toys for %s' % (h_pull[fit].Integral(1, bin_pull[0]), label))
        h_pull[fit].GetXaxis().SetTitle('Best-fit pull value #mu/#sigma(#mu)')
        h_pull[fit].GetYaxis().SetTitle('Number of toys')
        h_pull[fit].SetMarkerStyle(20)
        h_pull[fit].SetMarkerColor(R.kBlack)

        ## Save the original, un-fitted histograms
        h_mu  [fit].Write()
        h_pull[fit].Write()

        ## Get quantiles from histograms: 2.5%, 16%, 50%, 84%, 97.5%
        quants = array.array('d', [0.025, 0.16, 0.50, 0.84, 0.975])
        q_mu   = array.array('d', [0, 0, 0, 0, 0])
        q_pull = array.array('d', [0, 0, 0, 0, 0])
        h_mu  [fit].GetQuantiles(5, q_mu  , quants)
        h_pull[fit].GetQuantiles(5, q_pull, quants)

        ## Perform an asymmetric Gaussian fit to the distributions
        ##   *** Signal strength ***
        fit_mu = R.TF1('fit_mu', '[0]*exp(-0.5*pow((x-[1])/([2]*(x<[1])+[3]*(x>=[1])),2))', bin_mu[1], bin_mu[2])
        fit_mu.SetParameter( 0, h_mu[fit].GetMaximum() )
        fit_mu.SetParameter( 1, h_mu[fit].GetMean()    )
        fit_mu.SetParameter( 2, h_mu[fit].GetStdDev()  )
        fit_mu.SetParameter( 3, h_mu[fit].GetStdDev()  )
        h_mu[fit].Fit('fit_mu')
        ##   *** Pull value ***
        fit_pull = R.TF1('fit_pull', '[0]*exp(-0.5*pow((x-[1])/([2]*(x<[1])+[3]*(x>=[1])),2))', bin_pull[1], bin_pull[2])
        fit_pull.SetParameter( 0, h_pull[fit].GetMaximum() )
        fit_pull.SetParameter( 1, h_pull[fit].GetMean()    )
        fit_pull.SetParameter( 2, h_pull[fit].GetStdDev()  )
        fit_pull.SetParameter( 3, h_pull[fit].GetStdDev()  )
        h_pull[fit].Fit('fit_pull')

        ## Create canvases on which to draw the histograms and fits
        c_mu   = R.TCanvas('c_mu_%s'   % fit, 'c_mu_%s'   % fit, 800, 600)
        c_pull = R.TCanvas('c_pull_%s' % fit, 'c_pull_%s' % fit, 800, 600)
        
        ## Style options for histograms and legends
        R.gStyle.SetOptStat(0)
        R.gStyle.SetLegendTextSize(0.03)
        # ## Doesn't work, despite discussion here:
        # ## https://root-forum.cern.ch/t/setting-histogram-title-size-in-root/4468
        # R.gStyle.SetTitleFontSize(1.0)

        ## Draw histograms and fits
        c_mu.cd()
        h_mu[fit].SetLineWidth(2)
        h_mu[fit].Draw('e')

        c_pull.cd()
        h_pull[fit].SetLineWidth(2)
        h_pull[fit].Draw('e')
        
        ## Add legends
        c_mu.cd()
        l_mu = R.TLegend(0.68,0.55,0.895,0.895)
        l_mu.AddEntry(0, '2-side gaus #it{#Chi^{2}} = %.2f' % (fit_mu.GetChisquare()/fit_mu.GetNDF()), '')
        l_mu.AddEntry(0, 'Fit peak = %.3f' % fit_mu.GetParameter(1)         , '')
        l_mu.AddEntry(0, 'Median  = %.3f'  % q_mu[2]                        , '')
        l_mu.AddEntry(0, 'Fit #sigma = -%.2f/+%.2f' % (fit_mu.GetParameter(2),
                                                       fit_mu.GetParameter(3)), '')
        l_mu.AddEntry(0,       '68%% = %.2f/+%.2f'  % (q_mu[1], q_mu[3])      , '')
        l_mu.AddEntry(0,       '95%% = %.2f/+%.2f'  % (q_mu[0], q_mu[4])      , '')
        l_mu.SetMargin(0.05)
        l_mu.Draw('same')

        c_pull.cd()
        l_pull = R.TLegend(0.68,0.55,0.895,0.895)
        l_pull.AddEntry(0, '2-side gaus #it{#Chi^{2}} = %.2f' % (fit_pull.GetChisquare()/fit_pull.GetNDF()), '')
        l_pull.AddEntry(0, 'Fit peak = %.3f' % fit_pull.GetParameter(1)         , '')
        l_pull.AddEntry(0, 'Median  = %.3f'  % q_pull[2]                        , '')
        l_pull.AddEntry(0, 'Fit #sigma = -%.2f/+%.2f' % (fit_pull.GetParameter(2),
                                                         fit_pull.GetParameter(3)), '')
        l_pull.AddEntry(0,     '68%% = %.2f/+%.2f'    % (q_pull[1], q_pull[3])    , '')
        l_pull.AddEntry(0,     '95%% = %.2f/+%.2f'    % (q_pull[0], q_pull[4])    , '')
        l_pull.SetMargin(0.05)
        l_pull.Draw('same')

        ## Save canvases
        c_mu  .SaveAs('%s/c_mu_%s.png'   % (label, label))
        c_pull.SaveAs('%s/c_pull_%s.png' % (label, label))
        ## c_mu  .Write()  ## For some reason, can't open in output ROOT file: causes segfault
        ## c_pull.Write()  ## For some reason, can't open in output ROOT file: causes segfault

        ## Print the fit parameters
        print '\nHistogram mu mean = %.3f, std. dev. = %.3f' % (h_mu[fit].GetMean(), h_mu[fit].GetStdDev())
        print   'Gauss fit mu mean = %.3f, std. dev. = -%.3f/+%.3f, chi2 = %.3f' % (fit_mu.GetParameter(1), fit_mu.GetParameter(2),
                                                                                    fit_mu.GetParameter(3), fit_mu.GetChisquare()/fit_mu.GetNDF())

        print '\nHistogram pull mean = %.3f, std. dev. = %.3f' % (h_pull[fit].GetMean(), h_pull[fit].GetStdDev())
        print   'Gauss fit pull mean = %.3f, std. dev. = -%.3f/+%.3f, chi2 = %.3f' % (fit_pull.GetParameter(1), fit_pull.GetParameter(2),
                                                                                      fit_pull.GetParameter(3), fit_pull.GetChisquare()/fit_pull.GetNDF())
    

    ## ********************************************* ##
    ##  Fit and write out histograms from Workspace  ##
    ## ********************************************* ##

    h_bkg_toys.SetTitle('Background yield distribution in %d toys' % h_bkg_toys.Integral())
    h_bkg_toys.GetXaxis().SetTitle('Number of background events')
    h_bkg_toys.GetYaxis().SetTitle('Number of toys')
    h_bkg_toys.SetMarkerStyle(20)
    h_bkg_toys.SetMarkerColor(R.kBlack)
    h_bkg_toys.SetLineColor(R.kBlack)

    if 'BDT' in label:
        h_ROI_toys.SetTitle('Background yield in highest BDT bin in %d toys' % h_ROI_toys.Integral())
        h_ROI_toys.GetXaxis().SetTitle('Background yield in highest two BDT bins')
    else:
        h_ROI_toys.SetTitle('Background yield in 123 - 127 GeV in %d toys' % h_ROI_toys.Integral())
        h_ROI_toys.GetXaxis().SetTitle('Background yield in 123 - 127 GeV')
    h_ROI_toys.GetYaxis().SetTitle('Number of toys')
    h_ROI_toys.SetMarkerStyle(20)
    h_ROI_toys.SetMarkerColor(R.kBlue)
    h_ROI_toys.SetLineColor(R.kBlue)

    ## Save the original, un-fitted histograms
    h_bkg_toys.Write()
    h_ROI_toys.Write()
    
    ## Get quantiles from histograms: 2.5%, 16%, 50%, 84%, 97.5%
    quants = array.array('d', [0.025, 0.16, 0.50, 0.84, 0.975])
    q_bkg  = array.array('d', [0, 0, 0, 0, 0])
    q_ROI  = array.array('d', [0, 0, 0, 0, 0])
    h_bkg_toys.GetQuantiles(5, q_bkg, quants)
    h_ROI_toys.GetQuantiles(5, q_ROI, quants)
    ## Quadrature sum of shape variation, relative to median
    sum_bin_err_sq = 0
    for iBin in range(bin_disc[0]):
        q_bin = array.array('d', [0, 0, 0, 0, 0])
        h_bins[iBin].GetQuantiles(5, q_bin, quants)
        g_q95_toys.SetPoint(iBin,                    (iBin+1)*0.9-0.25, q_bin[0])
        g_q68_toys.SetPoint(iBin,                    (iBin+1)*0.9-0.25, q_bin[1])
        g_med_toys.SetPoint(iBin,                    (iBin+1)*0.9-0.25, q_bin[2])
        g_q68_toys.SetPoint(2*bin_disc[0] - iBin - 1, (iBin+1)*0.9-0.25, q_bin[3])
        g_q95_toys.SetPoint(2*bin_disc[0] - iBin - 1, (iBin+1)*0.9-0.25, q_bin[4])
        sum_bin_err_sq += pow((q_bin[1] - q_bin[2])/q_bin[2], 2)  ## Lower 68%
        sum_bin_err_sq += pow((q_bin[3] - q_bin[2])/q_bin[2], 2)  ## Upper 68%


    ## Perform an asymmetric Gaussian fit to the distributions
    ##   *** Total background yield ***
    fit_bkg = R.TF1('fit_bkg', '[0]*exp(-0.5*pow((x-[1])/([2]*(x<[1])+[3]*(x>=[1])),2))', bin_mu[1], bin_mu[2])
    fit_bkg.SetParameter( 0, h_bkg_toys.GetMaximum() )
    fit_bkg.SetParameter( 1, h_bkg_toys.GetMean()    )
    fit_bkg.SetParameter( 2, h_bkg_toys.GetStdDev()  )
    fit_bkg.SetParameter( 3, h_bkg_toys.GetStdDev()  )
    h_bkg_toys.Fit('fit_bkg')
    ##   *** Background yield in 123 - 127 GeV or the last BDT bin ***
    fit_ROI = R.TF1('fit_ROI', '[0]*exp(-0.5*pow((x-[1])/([2]*(x<[1])+[3]*(x>=[1])),2))', bin_pull[1], bin_pull[2])
    fit_ROI.SetParameter( 0, h_ROI_toys.GetMaximum() )
    fit_ROI.SetParameter( 1, h_ROI_toys.GetMean()    )
    fit_ROI.SetParameter( 2, h_ROI_toys.GetStdDev()  )
    fit_ROI.SetParameter( 3, h_ROI_toys.GetStdDev()  )
    h_ROI_toys.Fit('fit_ROI')


    ## Create canvases on which to draw the histograms and fits
    c_bkg_toys = R.TCanvas('c_bkg_toys', 'c_bkg_toys', 800, 600)
    c_ROI_toys = R.TCanvas('c_ROI_toys', 'c_ROI_toys', 800, 600)
    c_bin_toys = R.TCanvas('c_bin_toys', 'c_bin_toys', 800, 600)
        
    ## Style options for histograms and legends
    R.gStyle.SetOptStat(0)
    R.gStyle.SetLegendTextSize(0.03)
    # ## Doesn't work, despite discussion here:
    # ## https://root-forum.cern.ch/t/setting-histogram-title-size-in-root/4468
    # R.gStyle.SetTitleFontSize(1.0)

    ## Draw histograms and fits
    c_bkg_toys.cd()
    h_bkg_toys.SetLineWidth(2)
    h_bkg_toys.Draw('e')

    c_ROI_toys.cd()
    h_ROI_toys.SetLineWidth(2)
    h_ROI_toys.Draw('e')

    c_bin_toys.cd()
    f_bin_toys = c_bin_toys.DrawFrame(1.4,0.001, 4.1, 10)
    f_bin_toys.SetMinimum(0)
    f_bin_toys.SetMaximum( R.TMath.MaxElement(g_q95_toys.GetN(), g_q95_toys.GetY())*1.2 )
    if 'BDT' in label:
        f_bin_toys.GetXaxis().SetLimits(bin_disc[1], bin_disc[2])
        f_bin_toys.GetXaxis().SetTitle('Event BDT score')
    else:
        f_bin_toys.GetXaxis().SetLimits(bin_disc[1], bin_disc[2])
        f_bin_toys.GetXaxis().SetTitle('Candidate di-muon mass (1 GeV bins)')
    f_bin_toys.GetYaxis().SetTitle('Number of background events')
    f_bin_toys.SetTitle('Background yield from %d toys with %s' % (h_bkg_toys.Integral(), label))

    g_q95_toys.SetLineColor(R.kOrange)
    g_q95_toys.SetFillColor(R.kOrange)
    g_q95_toys.SetFillStyle(1001)
    g_q95_toys.Draw('F')

    g_q68_toys.SetLineColor(R.kGreen+1)
    g_q68_toys.SetFillColor(R.kGreen+1)
    g_q68_toys.SetFillStyle(1001)
    g_q68_toys.Draw('Fsame')

    g_med_toys.SetLineColor(R.kBlack)
    g_med_toys.SetLineWidth(3)
    g_med_toys.SetLineStyle(2)
    g_med_toys.Draw('L')
    f_bin_toys.Draw('sameaxis')

    ## Add legends
    c_bkg_toys.cd()
    l_bkg = R.TLegend(0.68,0.55,0.895,0.895)
    l_bkg.AddEntry(0, '2-side gaus #it{#Chi^{2}} = %.2f' % (fit_bkg.GetChisquare()/fit_bkg.GetNDF()), '')
    l_bkg.AddEntry(0, 'Fit peak = %.2f' % fit_bkg.GetParameter(1)         , '')
    l_bkg.AddEntry(0, 'Median  = %.2f'  % q_bkg[2]                        , '')
    l_bkg.AddEntry(0, 'Fit #sigma = -%.1f/+%.1f' % (fit_bkg.GetParameter(2),
                                                    fit_bkg.GetParameter(3)), '')
    l_bkg.AddEntry(0,       '68%% = %.1f/+%.1f'  % (q_bkg[1], q_bkg[3])     , '')
    l_bkg.AddEntry(0,       '95%% = %.1f/+%.1f'  % (q_bkg[0], q_bkg[4])     , '')
    l_bkg.SetMargin(0.05)
    l_bkg.Draw('same')

    c_ROI_toys.cd()
    l_ROI = R.TLegend(0.68,0.55,0.895,0.895)
    l_ROI.AddEntry(0, '2-side gaus #it{#Chi^{2}} = %.2f' % (fit_ROI.GetChisquare()/fit_ROI.GetNDF()), '')
    l_ROI.AddEntry(0, 'Fit peak = %.2f' % fit_ROI.GetParameter(1)         , '')
    l_ROI.AddEntry(0, 'Median  = %.2f'  % q_ROI[2]                        , '')
    l_ROI.AddEntry(0, 'Fit #sigma = -%.1f/+%.1f' % (fit_ROI.GetParameter(2),
                                                    fit_ROI.GetParameter(3)), '')
    l_ROI.AddEntry(0,       '68%% = %.1f/+%.1f'  % (q_ROI[1], q_ROI[3])     , '')
    l_ROI.AddEntry(0,       '95%% = %.1f/+%.1f'  % (q_ROI[0], q_ROI[4])     , '')
    l_ROI.SetMargin(0.05)
    l_ROI.Draw('same')

    c_bin_toys.cd()
    l_bin = R.TLegend(0.68,0.55,0.895,0.895)
    l_bin.AddEntry(0, 'Avg. bin err = %.1f%%' % (100*math.sqrt(sum_bin_err_sq/100)), '')
    l_bin.SetMargin(0.05)
    l_bin.Draw('same')


    ## Save canvases
    c_bkg_toys.SaveAs('%s/c_bkg_toys_%s.png' % (label, label))
    c_ROI_toys.SaveAs('%s/c_ROI_toys_%s.png' % (label, label))
    ## c_bkg_toys.Write()  ## For some reason, can't open in output ROOT file: causes segfault
    ## c_ROI_toys.Write()  ## For some reason, can't open in output ROOT file: causes segfault
    c_bin_toys.SaveAs('%s/c_bin_toys_%s.png' % (label, label))


    out_file.Close()
    print '\nDone analyzing input %s\n\n' % (label, label)

## End function: AnalyzeToys(label):


def main():

    for LABEL in LABELS:
        AnalyzeToys(LABEL)


main()  ## End of main function


