#! /usr/bin/env python

###############################################
###              RankVars.py                ###
###                                         ###
###  Rank signal-background discriminating  ###
###  variables by S/sqrt(B) + correlations  ###
###                                         ###
###           Andrew Brinkerhoff            ###
###               10.06.2018                ###
###############################################

## Basic python includes for manipulating files
import sys
import os

## Other python includes
import math
import array
import numpy
import operator

## ROOT includes
import ROOT as R
import ROOT.RooFit as RF

## Tools to compute significance
import H2MuAnalyzer.AutoCat.AutoCatTools as AC


INDIR   = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/Run2/TrainMVA/output'
# FACTORY = 'f_Opt_AWB_noMass_allVars_resWgt_all_sig_all_bkg_resWgt'
# FACTORY = 'f_Opt_AWB_noMass_allVars_ttH_vs_ttW_ttH_sig_ttW_bkg_ge0j'
# FACTORY = 'f_Opt_AWB_noMass_allVars_ttH_vs_ttbar_ttH_sig_ttbar_bkg_ge0j'

# INFILE  = 'Run2_ttH_3l_all_vs_all_2019_07_02_v1.root'
# # FACTORY = 'f_Opt_AWB_noMass_allVars_ttH_vs_ttbar_ttH_sig_ttbar_bkg_ge4j_btag'
# # FACTORY = 'f_Opt_AWB_noMass_allVars_ttH_vs_ttW_ttH_sig_ttW_bkg_ge4j_btag'
# # FACTORY = 'f_Opt_AWB_noMass_allVars_ttH_vs_ttZ_ttH_sig_ttZ_bkg_ge0j'
# FACTORY = 'f_Opt_AWB_noMass_allVars_all_sig_all_bkg_ge0j'

INFILE  = "Run2_ttH_3l_varStudy_ttH_vs_all_jetCuts_2019_07_13_test_v1.root"
# FACTORY = "f_Opt_AWB_noMass_allVars_ttH_sig_all_bkg_le3j_notag"
FACTORY = "f_Opt_AWB_noMass_allVars_ttH_sig_all_bkg_ge4j_btag"

## Script global settings
MAX_EVT =     -1  ## Total number of events in the tree to process
PRT_EVT =   1000  ## Print out every Nth event when processing
NBINS   =   1000  ## Number of bins for fine-grained version of each variable histogram
NQUANT  =     10  ## Number of quantile bins for re-binned signal and background distributions
NVARS   =     -1  ## Maximum number of variables to consider
RUN2D   =  False  ## Compute 2D significance values (takes ~2x as long to run)
 

## Main function executed by ./macros/RankVars.py
def main():

    print '\n\n*** Inside RankVars.py ***'

    ## Don't print plots to screen while running (faster)
    R.gROOT.SetBatch(R.kTRUE)
    ## Don't use the stats box in any plots
    R.gStyle.SetOptStat(0)

    ## Import root file with data and MC variable distributions and trees
    print '\nOpening file %s' % (INDIR+'/'+INFILE)
    in_file = R.TFile(INDIR+'/'+INFILE)

    out_dir = 'plots/%s/%s' % (INFILE.replace('.root',''), FACTORY)
    if not os.path.exists(out_dir):
        print '\nCreating output folders %s/png and pdf' % out_dir
        os.makedirs( '%s/png' % out_dir )
        os.makedirs( '%s/pdf' % out_dir )

    print '\nCreating output root file %s/RankVars.root' % out_dir
    out_file = R.TFile('%s/RankVars.root' % out_dir, 'recreate')


    ## Get list of input variable histograms and make new histograms with finer binning
    in_file.cd(FACTORY+'/InputVariables_Id')
    h_vars = {}
    for key in R.gDirectory.GetListOfKeys():
        hist   = in_file.Get(FACTORY+'/InputVariables_Id/'+key.GetName())
        nBins  = hist.GetNbinsX()
        x_low  = hist.GetBinLowEdge(0)
        x_high = hist.GetBinLowEdge(nBins+1)
        h_name = key.GetName().replace('__Signal_Id', '').replace('__Background_Id', '')
        if not h_name in h_vars.keys():
            h_vars[h_name] = {}

        if hist.GetName().endswith('__Signal_Id'):
            new_hist = R.TH1F(h_name+'_sig', h_name+'_sig', NBINS, x_low, x_high)
            h_vars[h_name]['sig'] = [hist, new_hist]
        if hist.GetName().endswith('__Background_Id'):
            new_hist = R.TH1F(h_name+'_bkg', h_name+'_bkg', NBINS, x_low, x_high)
            h_vars[h_name]['bkg'] = [hist, new_hist]
        ## Only process up to NVARS variables
        if NVARS > 0 and len(h_vars.keys()) >= NVARS and 'sig' in h_vars[h_name].keys() and 'bkg' in h_vars[h_name].keys(): break


    ## Loop over all events and fill re-binned histograms
    tr = in_file.Get(FACTORY+'/TrainTree')

    ## Number of signal and background events processed
    nSig    = 0
    nBkg    = 0
    nSigWgt = 0
    nBkgWgt = 0
    nEvt    = tr.GetEntries()
    nSkip   = int(nEvt/MAX_EVT) if MAX_EVT > 0 else 1

    print '\nAbout to process %d events (or up to %d) - will use 1/%d of the events' % (tr.GetEntries(), MAX_EVT, nSkip)
    for iEvt in range(tr.GetEntries()):

        if not iEvt % nSkip is 0: continue

        if iEvt % PRT_EVT is 0: print 'Event #', iEvt

        tr.GetEntry(iEvt)

        isSig = (tr.samp_ID < 0)  ## Any signal sample
        isBkg = (tr.samp_ID > 0)  ## Any background sample

        if not (isSig or isBkg): continue  ## Data samples

        WGT = tr.weight  ## Weight from BDT training, normalized to equal signal-background area
        ## WGT = tr.total_wgt  ## True weight, including proper Xsec * Lumi scaling

        for var in h_vars.keys():
            val = tr.GetBranch(var).GetLeaf(var).GetValue()
            if isSig: h_vars[var]['sig'][1].Fill(val, WGT)
            if isBkg: h_vars[var]['bkg'][1].Fill(val, WGT)

        nSig    +=   1*(isSig)
        nSigWgt += WGT*(isSig)
        nBkg    +=   1*(isBkg)
        nBkgWgt += WGT*(isBkg)

    print '\nFinished processing %d events' % iEvt
    print '  * %d signal     (%.2f weighted)' % (nSig, nSigWgt)
    print '  * %d background (%.2f weighted)' % (nBkg, nBkgWgt)


    ## Renormalize signal and background histograms to area = 1
    for var in h_vars.keys():
        for i in range(2):
            h_vars[var]['sig'][i].Scale(1.0 / h_vars[var]['sig'][i].Integral())
            h_vars[var]['bkg'][i].Scale(1.0 / h_vars[var]['bkg'][i].Integral())

    ## Rebin 1D histograms into 1/NQUANT quantiles in background
    for var in h_vars.keys():
        bins = [(1.0*bin/NQUANT) for bin in range(NQUANT+1)]
        quants = array.array('d', [(1.0*bin/NQUANT) for bin in range(NQUANT+1)])
        q_var  = array.array('d', [               0 for bin in range(NQUANT+1)])
        h_vars[var]['bkg'][1].GetQuantiles(NQUANT+1, q_var, quants)
        ## Reset lower and upper boundaries by hand: "quantile" sometimes truncates bins with negative weights
        q_var[0]  = h_vars[var]['sig'][1].GetBinLowEdge(1)
        q_var[-1] = h_vars[var]['sig'][1].GetBinLowEdge(NBINS+1)

        h_vars[var]['sig'].append( h_vars[var]['sig'][1].Rebin(NQUANT, var+'_sig_quant', q_var) )
        h_vars[var]['bkg'].append( h_vars[var]['bkg'][1].Rebin(NQUANT, var+'_bkg_quant', q_var) )
        h_vars[var]['bin'] = q_var

    if RUN2D:
        ## Book 2D histograms of varY vs. varX in quantile bins, and TGraphs of S/sqrt(B)
        for varX in h_vars.keys():
            h_vars[varX]['2D'] = {}
            for varY in h_vars.keys():
                if varX is varY: continue
                h_vars[varX]['2D'][varY] = {}
                h_vars[varX]['2D'][varY]['sig'] = R.TH2F(varX+'_vs_'+varY+'_sig', varX+'_vs_'+varY+'_sig',  NQUANT, h_vars[varX]['bin'], NQUANT, h_vars[varY]['bin'])
                h_vars[varX]['2D'][varY]['bkg'] = R.TH2F(varX+'_vs_'+varY+'_bkg', varX+'_vs_'+varY+'_bkg',  NQUANT, h_vars[varX]['bin'], NQUANT, h_vars[varY]['bin'])

        ## Loop over events again to fill 2D plots
        print '\nLoop #2: about to process %d events (or up to %d) - will use 1/%d of the events' % (tr.GetEntries(), MAX_EVT, nSkip)
        for jEvt in range(tr.GetEntries()):

            if not jEvt % nSkip is 0: continue

            if jEvt % PRT_EVT is 0: print 'Event #', jEvt

            tr.GetEntry(jEvt)

            isSig = (tr.samp_ID < 0)  ## Any signal sample
            isBkg = (tr.samp_ID > 0)  ## Any background sample

            if not (isSig or isBkg): continue  ## Data samples

            WGT = tr.weight  ## Weight from BDT training, normalized to equal signal-background area
            ## WGT = tr.total_wgt  ## True weight, including proper Xsec * Lumi scaling

            vals = {}
            for var in h_vars.keys():
                vals[var] = tr.GetBranch(var).GetLeaf(var).GetValue()

            for varX in h_vars.keys():
                for varY in h_vars.keys():
                    if varX is varY: continue
                    if isSig: h_vars[varX]['2D'][varY]['sig'].Fill(vals[varX], vals[varY], WGT)
                    if isBkg: h_vars[varX]['2D'][varY]['bkg'].Fill(vals[varX], vals[varY], WGT)


        ## Renormalize signal and background histograms to area = 1
        for varX in h_vars.keys():
            for varY in h_vars.keys():
                if varX is varY: continue
                h_vars[varX]['2D'][varY]['sig'].Scale(1.0 / h_vars[varX]['2D'][varY]['sig'].Integral())
                h_vars[varX]['2D'][varY]['bkg'].Scale(1.0 / h_vars[varX]['2D'][varY]['bkg'].Integral())

    ## End conditional: if RUN2D
            

    ## Store significance with conservative, nominal, and no systematics for each variable and each 2D combination of variables
    sig = {}
    sig['1D'] = {}
    sig['2D'] = {}
    for varX in h_vars.keys():

        ## Compute 1D significance ('nominal', 'no systematics', 'conservative')
        sig['1D'][varX] = [0, 0, 0]

        ## Re-bin for maximum significance under each scenario
        ## MergeBins(h_sig, h_bkg, syst = 'conserv', min_bkg = 0.0, max_loss = 0.02, max_bin_width = -1, verbose = False)
        bins_nominal = array.array('d', AC.MergeBins(h_vars[varX]['sig'][2], h_vars[varX]['bkg'][2], 'nominal', 0.0, 0.001))
        bins_noSyst  = array.array('d', AC.MergeBins(h_vars[varX]['sig'][2], h_vars[varX]['bkg'][2], 'none',    0.0, 0.001))
        bins_conserv = array.array('d', AC.MergeBins(h_vars[varX]['sig'][2], h_vars[varX]['bkg'][2], 'conserv', 0.0, 0.001))

        sig_nominal = h_vars[varX]['sig'][2].Rebin(len(bins_nominal)-1, 'sig_nominal', bins_nominal)
        sig_noSyst  = h_vars[varX]['sig'][2].Rebin(len(bins_noSyst )-1, 'sig_noSyst' , bins_noSyst )
        sig_conserv = h_vars[varX]['sig'][2].Rebin(len(bins_conserv)-1, 'sig_conserv', bins_conserv)
        bkg_nominal = h_vars[varX]['bkg'][2].Rebin(len(bins_nominal)-1, 'bkg_nominal', bins_nominal)
        bkg_noSyst  = h_vars[varX]['bkg'][2].Rebin(len(bins_noSyst )-1, 'bkg_noSyst' , bins_noSyst )
        bkg_conserv = h_vars[varX]['bkg'][2].Rebin(len(bins_conserv)-1, 'bkg_conserv', bins_conserv)

        ## Compute the significance under each scenario
        ## Scale the "conservative" and "nominal" approaches by the sensitivity of the un-binned integral
        err_sig = R.Double(-99)
        err_bkg = R.Double(-99)
        int_sig = sig_conserv.IntegralAndError(1, sig_conserv.GetNbinsX(), err_sig)
        int_bkg = bkg_conserv.IntegralAndError(1, bkg_conserv.GetNbinsX(), err_bkg)
        nominal_scale = pow(int_sig, 2)           / (int_bkg + pow(err_bkg, 2))
        conserv_scale = pow(int_sig - err_sig, 2) / (int_bkg + pow(err_bkg, 2) + err_bkg)
        ## ComputeSignificance(h_sig, h_bkg, syst = 'conserv', min_bkg = 0.0, binning = [], verbose = False):
        sig['1D'][varX][0] = AC.ComputeSignificance(sig_conserv, bkg_conserv, 'conserv', 0.0) / conserv_scale
        sig['1D'][varX][1] = AC.ComputeSignificance(sig_nominal, bkg_nominal, 'nominal', 0.0) / nominal_scale
        sig['1D'][varX][2] = AC.ComputeSignificance(sig_noSyst,  bkg_noSyst,  'none',    0.0)

        if RUN2D:
            ## Compute 2D significances
            for varY in h_vars.keys():
                if varX is varY: continue
                sig['2D'][varX+'_vs_'+varY]    = [0, 0, 0]
                sig['2D'][varX+'_vs_'+varY][0] = AC.ComputeSignificance2D(h_vars[varX]['2D'][varY]['sig'], h_vars[varX]['2D'][varY]['bkg'], 'conserv', 0.2/pow(NQUANT, 2)) / conserv_scale
                sig['2D'][varX+'_vs_'+varY][1] = AC.ComputeSignificance2D(h_vars[varX]['2D'][varY]['sig'], h_vars[varX]['2D'][varY]['bkg'], 'nominal', 0.2/pow(NQUANT, 2)) / nominal_scale
                sig['2D'][varX+'_vs_'+varY][2] = AC.ComputeSignificance2D(h_vars[varX]['2D'][varY]['sig'], h_vars[varX]['2D'][varY]['bkg'], 'none',    0.2/pow(NQUANT, 2))
                

    ## Sort by significance
    sig1D = sorted(sig['1D'].items(), key=operator.itemgetter(1), reverse=True)

    if RUN2D:
        ## Clean out duplicate 2D entries
        sig2D = sorted(sig['2D'].items(), key=operator.itemgetter(1), reverse=True)
        for i in range(len(sig1D)):
            for j in range(i+1, len(sig1D)):
                sig2D = [x for x in sig2D if x[0] != (sig1D[i][0]+'_vs_'+sig1D[j][0])]
    
    ## Print significances
    print '\nSingle variable significance values:\n'
    for var in sig1D:
        print '%20s = %.3f (conserv syst) / %.3f (nominal syst) / %.3f (no syst)' % ( var[0], var[1][0], var[1][1], var[1][2] )
    if RUN2D:
        print '\nCorrelated two-variable significance values:'
        for var in sig2D:
            var1 = var[0].split('_vs_')[0]
            var2 = var[0].split('_vs_')[1]
            gain = 100 * ((var[1][0] / max(sig['1D'][var1][0], sig['1D'][var2][0])) - 1)
            corr = var[1][0] / (sig['1D'][var1][0] * sig['1D'][var2][0])
            print '%20s vs. %20s = %.3f (conserv syst) / %.3f (nominal syst) / %.3f (no syst)  : %3d%% gain, %.3f significance from correlation' % (var1, var2, var[1][0], var[1][1], var[1][2], gain, corr)


    ## Write out histograms and draw plots
    out_file.cd()
    for var in h_vars.keys():
        for i in range(3):
            h_vars[var]['sig'][i].Write()
            h_vars[var]['bkg'][i].Write()

        c_1D = R.TCanvas('c_1D_%s' % var, 'c_1D_%s' % var, 800, 600)
        c_1D.cd()

        h_vars[var]['sig'][2].SetTitle('%s in %d quantiles' % (var, NQUANT))
        h_vars[var]['sig'][2].GetXaxis().SetTitle('%s' % var)
        h_vars[var]['sig'][2].GetYaxis().SetTitle('Fraction of events')
        h_vars[var]['sig'][2].GetYaxis().SetRangeUser(0, max(h_vars[var]['sig'][2].GetMaximum(), h_vars[var]['bkg'][2].GetMaximum())*1.5)

        h_vars[var]['sig'][2].SetLineColor(R.kBlue)
        h_vars[var]['sig'][2].SetLineWidth(3)
        h_vars[var]['sig'][2].SetFillColor(R.kBlue)
        h_vars[var]['sig'][2].SetFillStyle(3001)
        h_vars[var]['bkg'][2].SetLineColor(R.kRed)
        h_vars[var]['bkg'][2].SetLineWidth(3)
        h_vars[var]['bkg'][2].SetFillColor(R.kRed)
        h_vars[var]['bkg'][2].SetFillStyle(3004)

        h_vars[var]['sig'][2].Draw('histe')
        h_vars[var]['bkg'][2].Draw('histesame')

        l_1D = R.TLegend(0.50,0.75,0.895,0.895)
        l_1D.SetNColumns(2)
        l_1D.AddEntry(h_vars[var]['sig'][2], 'Sig.', 'lpfe')
        l_1D.AddEntry(0, 'Significance = %.3f' % sig['1D'][var][1], '')
        l_1D.AddEntry(h_vars[var]['bkg'][2], 'Bkg.', 'lpfe')
        l_1D.AddEntry(0, 'Conservative = %.3f' % sig['1D'][var][0], '')
        l_1D.SetMargin(0.05)
        l_1D.Draw('same')


        c_1D.SaveAs('%s/png/%s.png' % (out_dir, c_1D.GetName()))
        c_1D.SaveAs('%s/pdf/%s.pdf' % (out_dir, c_1D.GetName()))

    if RUN2D:
        for varX in h_vars.keys():
            for varY in h_vars.keys():
                if varX is varY: continue
                h_vars[varX]['2D'][varY]['sig'].Write()
                h_vars[varX]['2D'][varY]['bkg'].Write()

    ## Close output and input files
    print '\nFinished, closing files ...'
    out_file.Close()
    in_file.Close()
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
