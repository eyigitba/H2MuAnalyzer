
## Tools to help with auto-categorization

import sys
import math

import ROOT as R


## Compute significance of signal histogram w.r.t. background
def ComputeSignificance(h_sig, h_bkg, syst = 'conserv', binning = [], verbose = False):

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
        bins = binning

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
            if den == 0: continue
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
            if den == 0:
                if verbose: print '  * Bins [%d, %d] have %.2f +/- %.2f background events (skipping)' % (bin[0], bin[1], den, math.sqrt(err))
                continue
            sigma_with_uncert += ( pow(num, 2) / (den + pow(err, 2)) )
            if verbose: print '  * Bins [%d, %d] have %.2f +/- %.2f (%.2f%%) background events' % (bin[0], bin[1], den, math.sqrt(err), 100*math.sqrt(err) / den)
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
            # if den == 0: continue
            if den < 0.2: continue
            sigma_conserv += ( pow(num - err_sig, 2) / (den + err_bkg + pow(err_bkg, 2)) )
            if verbose: print '  * Bins [%d, %d] have %.3f +/- %.3f signal events'              % (bin[0], bin[1], num, math.sqrt(err_sig))
            if verbose: print '  * Bins [%d, %d] have %.2f +/- %.2f (%.2f%%) background events' % (bin[0], bin[1], den, math.sqrt(err_bkg), 100*math.sqrt(err_bkg) / den)
        sigma_conserv = math.sqrt(sigma_conserv)

        if verbose: print 'Net significance with original binning and conservative statistical uncertainties = %.3f\n' % sigma_conserv


    ## Return the final estimated significance
    if syst == 'none':    return sigma_no_uncert
    if syst == 'nominal': return sigma_with_uncert
    if syst == 'conserv': return sigma_conserv

## End function: def ComputeSignificance(h_sig, h_bkg, syst = 'conserv', binning = [], verbose = False):


## Rebin a histogram to maximize the significance
## Continue rebinning as long as the relative significance loss compared to the maximum significance is less than "loss"
def MergeBins(h_sig, h_bkg, syst = 'conserv', max_loss = 0.02, verbose = False):

    if verbose: print '\nMerging bins for signal %s and background %s' % (h_sig.GetName(), h_bkg.GetName())
    if verbose: print '  * Options syst = %s, max_loss = %f' % (syst, max_loss)

    if (h_sig.GetNbinsX() != h_bkg.GetNbinsX()):
        print '\n\nnBinsSig = %d, nBinsBkg = %d! What are you, nuts!!!' % (h_sig.GetNbinsX(), h_bkg.GetNbinsX())
        sys.exit()
    if (syst != 'none' and syst != 'nominal' and syst != 'conserv'):
        print '\n\nInvalid systematics option %s chosen: pick "none", "nominal", or "conserv"'
        sys.exit()

    nBins = h_sig.GetNbinsX()

    ## Specify the original binning, i.e. one bin per bin "range"
    bins_orig = []
    for i in range(1, nBins+1):
        bins_orig.append([i,i])
    bins_max = bins_orig
    bins_old = bins_orig
    bins_new = bins_orig

    ## Compute the initial significance
    sigma_orig = ComputeSignificance(h_sig, h_bkg, syst, [], verbose)
    sigma_max  = sigma_orig
    sigma_old  = sigma_orig
    sigma_new  = sigma_orig

    while( (sigma_new / sigma_max) > (1 - max_loss) ):

        ## Reset "old" binning to "new" binning from previous iteration
        sigma_old = sigma_new
        bins_old  = bins_new

        ## Track the minimum loss [0], and the corresponding merged bins [1] and [2]
        min_loss = [99, -99, -99]
        for iBin in range(len(bins_old)-1):
            binA = bins_old[iBin]
            binB = bins_old[iBin+1]

            sig_A_sq  = ComputeSignificance(h_sig, h_bkg, syst, [[binA[0], binA[1]]], False)
            sig_B_sq  = ComputeSignificance(h_sig, h_bkg, syst, [[binB[0], binB[1]]], False)
            sig_AB_sq = ComputeSignificance(h_sig, h_bkg, syst, [[binA[0], binB[1]]], False)

            loss = sig_A_sq + sig_B_sq - sig_AB_sq
            if (loss < min_loss[0]):
                min_loss = [loss, iBin, iBin+1]
        
        if verbose:
            if (min_loss[0] >= 0): print '\nMerging bins %d:%d with %d:%d produces a minimum loss of %f' % (bins_old[min_loss[1]][0], bins_old[min_loss[1]][1],
                                                                                                            bins_old[min_loss[2]][0], bins_old[min_loss[2]][1], min_loss[0])
            else:                  print '\nMerging bins %d:%d with %d:%d produces a maximum gain of %f' % (bins_old[min_loss[1]][0], bins_old[min_loss[1]][1],
                                                                                                            bins_old[min_loss[2]][0], bins_old[min_loss[2]][1], min_loss[0])
        ## Merge the bins that produce the minimum loss (or maximum gain)
        bins_new = []
        for iBin in range(len(bins_old)):
            if iBin == min_loss[1]:
                bins_new.append([bins_old[iBin][0], bins_old[iBin+1][1]])
            elif iBin == min_loss[2]:
                continue
            else:
                bins_new.append(bins_old[iBin])

        ## Compute the significance with the new binning
        sigma_new = ComputeSignificance(h_sig, h_bkg, syst, bins_new, False)

        if verbose:
            print 'Old binning was: '+','.join(':'.join(str(idx) for idx in bin) for bin in bins_old)
            print '  * Yielded significance = %f' % sigma_old
            print 'New binning is: '+','.join(':'.join(str(idx) for idx in bin) for bin in bins_new)
            print '  * Yields significance = %f (%f%%)' % (sigma_new, 100*(1 - (sigma_old/sigma_new)))

        ## If the new significance is a new maximum, reset the maximum significance value and binning
        if (sigma_new > sigma_max):
            if verbose: print '\nFound a new maximum significance!  %f vs. %f\n' % (sigma_new, sigma_max)
            sigma_max = sigma_new
            bins_max  = bins_new

    ## End conditional: while( (sigma_new / sigma_max) > (1 - max_loss) ):

    
    ## Use the last binning that passed the "while" loop
    print '\nFinal binning is '+','.join(':'.join(str(idx) for idx in bin) for bin in bins_old)
    print '  * Significance = %f, vs. original %f and maximum %f' % (sigma_old, sigma_orig, sigma_max)

    ## Convert bin indices to upper and lower edges per bin
    edges = []
    for bin in bins_old:
        edges.append( h_sig.GetBinLowEdge(bin[0]) )
    edges.append( h_sig.GetBinLowEdge(nBins+1) )
    print 'Bin edges are ['+','.join(('%.2f' % edge) for edge in edges)+']\n'

    return edges

## End function: def MergeBins(h_sig, h_bkg, syst = 'conserv', verbose = False):

