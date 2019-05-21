
## Tools to help with auto-categorization

import sys
import math

import ROOT as R


## Compute significance of signal histogram w.r.t. background
def ComputeSignificance(h_sig, h_bkg, syst = 'conserv', min_bkg = 0.0, binning = [], verbose = False):

    if verbose: print '\nComputing significance for signal %s and background %s' % (h_sig.GetName(), h_bkg.GetName())
    if verbose: print '  * Options syst = '+syst+', binning = '+', '.join(':'.join(str(idx) for idx in bin) for bin in binning)

    if (h_sig.GetNbinsX() != h_bkg.GetNbinsX()):
        print '\n\nnBinsSig = %d, nBinsBkg = %d! What are you, nuts!!!' % (h_sig.GetNbinsX(), h_bkg.GetNbinsX())
        sys.exit()
    if (syst != 'none' and syst != 'nominal' and syst != 'conserv'):
        print '\n\nInvalid systematics option %s chosen: pick "none", "nominal", or "conserv"'
        sys.exit()

    nBins = h_sig.GetNbinsX()

    ## Use modified binning, if specified in input
    bins = []
    if len(binning) == 0:
        for i in range(1, nBins+1):
            ## Each modified bin specifies a range
            bins.append([i,i])
    else:
        bins = list(binning)

    if (verbose == 2 or syst == 'none'):
        ## Make a histogram of the signal^2
        h_sig_sq = h_sig.Clone('h_sig_sq')
        h_sig_sq.Multiply(h_sig)

        ## Print significance for a 1-bin experiment
        if verbose: print '\nSignificance treating the whole range as one bin = %.3f' % ( h_sig.Integral() / math.sqrt(h_bkg.Integral()) )

        ## Make a histogram of the signal^2 / background
        h_sigma_sq = h_sig_sq.Clone('h_sigma_sq')
        h_sigma_sq.Divide(h_bkg)
        sigma_integral = math.sqrt(h_sigma_sq.Integral())

        if verbose: print '\nNet significance estimated directly from the S^2 / B integral = %.3f' % sigma_integral

        ## Estimate the significance excluding uncertainties
        sigma_no_uncert = 0
        for bin in bins:
            num = h_sig.Integral(bin[0], bin[1])
            den = h_bkg.Integral(bin[0], bin[1])
            if den <= min_bkg:
                sigma_no_uncert += ( pow(num, 2) / den )
        sigma_no_uncert = math.sqrt(sigma_no_uncert)

        if verbose: print 'Net significance with original binning and no statistical uncertainties = %.3f\n' % sigma_no_uncert


    if (verbose == 2 or syst == 'nominal'):
        ## Estimate the significance including uncertainties
        sigma_with_uncert = 0
        for bin in bins:
            err = R.Double(-99)
            num = h_sig.Integral        (bin[0], bin[1])
            den = h_bkg.IntegralAndError(bin[0], bin[1], err)
            if den <= min_bkg:
                if verbose: print '  * Bins [%d, %d] have %.2f +/- %.2f background events (skipping)' % (bin[0], bin[1], den, math.sqrt(err))
                continue
            sigma_with_uncert += ( pow(num, 2) / (den + pow(err, 2)) )
            if verbose: print '  * Bins [%d, %d] have %.2f +/- %.2f (%.2f%%) background events' % (bin[0], bin[1], den, math.sqrt(err), 100*math.sqrt(err) / den)

        ## Don't return an estimate that is less than the integral of the whole range
        err = R.Double(-99)
        den = h_bkg.IntegralAndError(1, h_bkg.GetNbinsX(), err)
        sigma_with_uncert = max( sigma_with_uncert, pow(h_sig.Integral(), 2) / (den + pow(err, 2)) )
        ## Significance = sqrt(S^2 / B)
        sigma_with_uncert = math.sqrt(sigma_with_uncert)

        if verbose: print 'Net significance with original binning and nominal statistical uncertainties = %.3f\n' % sigma_with_uncert


    if (verbose == 2 or syst == 'conserv'):
        ## Estimate the significance with a conservative treatment of uncertainties
        ## Subtract signal uncertainty from signal yield, add background uncertainty to background yield
        sigma_conserv = 0
        for bin in bins:
            err_sig = R.Double(-99)
            err_bkg = R.Double(-99)
            num = h_sig.IntegralAndError(bin[0], bin[1], err_sig)
            den = h_bkg.IntegralAndError(bin[0], bin[1], err_bkg)
            if den <= min_bkg: continue
            sigma_conserv += ( pow(num - err_sig, 2) / (den + err_bkg + pow(err_bkg, 2)) )
            if verbose: print '  * Bins [%d, %d] have %.3f +/- %.3f signal events'              % (bin[0], bin[1], num, math.sqrt(err_sig))
            if verbose: print '  * Bins [%d, %d] have %.2f +/- %.2f (%.2f%%) background events' % (bin[0], bin[1], den, math.sqrt(err_bkg), 100*math.sqrt(err_bkg) / den)

        ## Don't return an estimate that is less than the integral of the whole range
        err = R.Double(-99)
        den = h_bkg.IntegralAndError(1, h_bkg.GetNbinsX(), err)
        sigma_conserv = max( sigma_conserv, pow(h_sig.Integral(), 2) / (den + pow(err, 2)) )
        ## Significance = sqrt(S^2 / B)
        sigma_conserv = math.sqrt(sigma_conserv)

        if verbose: print 'Net significance with original binning and conservative statistical uncertainties = %.3f\n' % sigma_conserv


    ## Return the final estimated significance
    if syst == 'none':    return sigma_no_uncert
    if syst == 'nominal': return sigma_with_uncert
    if syst == 'conserv': return sigma_conserv

## End function: def ComputeSignificance(h_sig, h_bkg, syst = 'conserv', min_bkg = 0.0, binning = [], verbose = False):


## Rebin a histogram to maximize the significance
## Continue rebinning as long as the relative significance loss compared to the maximum significance is less than "loss"
def MergeBins(h_sig, h_bkg, syst = 'conserv', min_bkg = 0.0, max_loss = 0.02, verbose = False):

    print '\nMerging %d bins for signal %s and background %s' % (h_sig.GetNbinsX(), h_sig.GetName(), h_bkg.GetName())
    print '  * Options syst = %s, min_bkg = %.2f, max_loss = %f' % (syst, min_bkg, max_loss)

    if (h_sig.GetNbinsX() != h_bkg.GetNbinsX()):
        print '\n\nnBinsSig = %d, nBinsBkg = %d! What are you, nuts!!!' % (h_sig.GetNbinsX(), h_bkg.GetNbinsX())
        sys.exit()
    if (syst != 'none' and syst != 'nominal' and syst != 'conserv'):
        print '\n\nInvalid systematics option %s chosen: pick "none", "nominal", or "conserv"'
        sys.exit()


    ## Compute the initial significance
    sigma_nom = ComputeSignificance(h_sig, h_bkg, 'nominal', 0.0, [], verbose)
    sigma_max = float(sigma_nom)
    sigma_new = float(sigma_nom)

    ## First attempt to simply merge by a factor of 2, to speed things up
    while( (sigma_new / sigma_max) > (1 - max_loss/10.) ):
        if h_sig.GetNbinsX() < 64: break
        h_sig_2 = (h_sig.Clone(h_sig.GetName())).Rebin(2)
        h_bkg_2 = (h_bkg.Clone(h_bkg.GetName())).Rebin(2)
        sigma_2 = ComputeSignificance(h_sig_2, h_bkg_2, 'nominal', 0.0, [], False)
        if (sigma_new / sigma_max) > (1 - max_loss/10.):
            print '*** Rebinning by a factor of 2 from %d bins to %d (significance %f to %f)' % (h_sig.GetNbinsX(), h_sig_2.GetNbinsX(), sigma_max, sigma_2)
            h_sig = h_sig_2.Clone(h_sig_2.GetName())
            h_bkg = h_bkg_2.Clone(h_bkg_2.GetName())
            sigma_max = float(sigma_2)
        else: break

    ## Reset "original" significance with rebin-by-2 histogram
    sigma_orig = ComputeSignificance(h_sig, h_bkg, syst, min_bkg, [], verbose)
    sigma_max  = float(sigma_orig)
    sigma_old  = float(sigma_orig)
    sigma_new  = float(sigma_orig)

    ## Start the regular, bin-by-bin merging sequence
    nBins = h_sig.GetNbinsX()

    ## Specify the original binning, i.e. one bin per bin "range"
    bins_orig = []
    for i in range(1, nBins+1):
        bins_orig.append([i,i])
    bins_max = list(bins_orig)  ## Make a hard copy, i.e. not a pointer
    bins_old = list(bins_orig)
    bins_new = list(bins_orig)

    while( (sigma_new / sigma_max) > (1 - max_loss) ):

        ## Reset "old" binning to "new" binning from previous iteration
        sigma_old = float(sigma_new)
        bins_old  = list(bins_new)

        ## Say whether we found adjacent low-stats background bins
        empty_bins = False
        ## Track the minimum loss [0], and the corresponding merged bins [1] and [2]
        min_loss = [99, -99, -99]
        for iBin in range(len(bins_old)-1):
            binA = list(bins_old[iBin])
            binB = list(bins_old[iBin+1])

            ## Automatically merge bins with low background stats
            if (h_bkg.Integral(binA[0], binA[1]) <= min_bkg and h_bkg.Integral(binB[0], binB[1]) <= min_bkg):
                min_loss   = [0, iBin, iBin+1]
                sigma_new  = float(sigma_old)
                empty_bins = True
                break

            ## Compute significance^2 for un-merged and merged bins
            sig_A_sq  = ComputeSignificance(h_sig, h_bkg, syst, min_bkg, [[binA[0], binA[1]]], False)
            sig_B_sq  = ComputeSignificance(h_sig, h_bkg, syst, min_bkg, [[binB[0], binB[1]]], False)
            sig_AB_sq = ComputeSignificance(h_sig, h_bkg, syst, min_bkg, [[binA[0], binB[1]]], False)

            loss = sig_A_sq + sig_B_sq - sig_AB_sq
            ## Found a new minimum loss
            if (loss < min_loss[0]):
                min_loss = [loss, iBin, iBin+1]
            ## If we actually gained by combining bins, simply merge these bins and continue
            if loss <= 0: break
        
        if verbose:
            if (min_loss[0] >= 0): print '\nMerging bins %d:%d with %d:%d produces a minimum loss of %f' % (bins_old[min_loss[1]][0], bins_old[min_loss[1]][1],
                                                                                                            bins_old[min_loss[2]][0], bins_old[min_loss[2]][1], min_loss[0])
            else:                  print '\nMerging bins %d:%d with %d:%d produces a maximum gain of %f' % (bins_old[min_loss[1]][0], bins_old[min_loss[1]][1],
                                                                                                            bins_old[min_loss[2]][0], bins_old[min_loss[2]][1], min_loss[0])
        ## Merge the bins that produce the minimum loss (or maximum gain)
        bins_new = list(bins_old)
        bins_new[min_loss[1]:min_loss[2]+1] = [[bins_old[min_loss[1]][0], bins_old[min_loss[2]][1]]]

        ## No need to re-compute significance if we were just merging empty bins
        if empty_bins: continue

        ## Compute the significance with the new binning
        sigma_new = ComputeSignificance(h_sig, h_bkg, syst, min_bkg, bins_new, False)

        if verbose:
            print 'Old binning was: '+','.join(':'.join(str(idx) for idx in bin) for bin in bins_old)
            print '  * Yielded significance = %f' % sigma_old
            print 'New binning is: '+','.join(':'.join(str(idx) for idx in bin) for bin in bins_new)
            print '  * Yields significance = %f (%f%%)' % (sigma_new, 100*(1 - (sigma_old/sigma_new)))

        ## If the new significance is a new maximum, reset the maximum significance value and binning
        if (sigma_new > sigma_max):
            if verbose: print '\nFound a new maximum significance!  %f (%d bins) vs. %f (%d bins)\n' % (sigma_new, len(bins_new), sigma_max, len(bins_max))
            sigma_max = float(sigma_new)
            bins_max  = list(bins_new)

    ## End conditional: while( (sigma_new / sigma_max) > (1 - max_loss) ):

    
    ## Use the last binning that passed the "while" loop
    print '\nFinal binning is '+','.join(':'.join(str(idx) for idx in bin) for bin in bins_old)
    print '  * Significance = %f, vs. original %f and maximum %f' % (sigma_old, sigma_orig, sigma_max)

    ## Convert bin indices to upper and lower edges per bin
    edges = []
    for bin in bins_old:
        edges.append( h_sig.GetBinLowEdge(bin[0]) )
    edges.append( h_sig.GetBinLowEdge(nBins+1) )
    print 'Bin edges are ['+','.join(('%.3f' % edge) for edge in edges)+']\n'

    return edges

## End function: def MergeBins(h_sig, h_bkg, min_bkg = 0.0, syst = 'conserv', verbose = False):

