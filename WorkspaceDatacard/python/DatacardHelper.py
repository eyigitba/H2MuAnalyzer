
## Functions to write datacards

import os
import sys

import ROOT as R
import ROOT.RooFit as RF

sys.path.insert(1, '%s/../FitBackground/python' % os.getcwd())
import FitFunctions as FF  ## From FitBackground/python/FitFunctions.py


def WriteHeader(card, cat, out_dir, in_file):

    print '\nWriting datacard with input file '+out_dir+'/workspace/'+in_file+'.root'
    
    card.write('imax *\n')
    card.write('jmax *\n')
    card.write('kmax *\n')
    card.write('----------------------------------------------------------------------------------------------------------------------------------\n')

    if   ('template' in in_file or 'rebin' in in_file):
        card.write('shapes * * '+out_dir+'/workspace/'+in_file+'.root '+'$PROCESS\n')
    elif ('shape' in in_file):
        card.write('shapes * * '+out_dir+'/workspace/'+in_file+'.root '+in_file+':$PROCESS\n')
    else:
        print '\n\nWriting datacard in category %s, input file name %s does not match template, rebin, or shape' % (cat, in_file)
        sys.exit()

    card.write('----------------------------------------------------------------------------------------------------------------------------------\n')
    card.write('bin            '+cat+'\n')
    card.write('observation    -1.0\n')
    card.write('----------------------------------------------------------------------------------------------------------------------------------\n')

## End function: WriteHeader(card, cat, dist, fit):


def WriteSigBkgBody(card, cat, dist, fit, width, nSig, nBkg):

    ## Category name (same for each process)
    card.write('bin'.ljust(width))
    card.write((cat).ljust(width))
    card.write((cat).ljust(width))
    card.write('\n')

    ## Process names
    card.write('process'.ljust(width))

    if   (fit == 'template_stack' or fit == 'rebin_stack'):
        card.write(('Net_Sig').ljust(width))
        card.write(('Net_Bkg').ljust(width))
    elif (fit == 'shape_MC'):
        card.write(('sig_fit').ljust(width))
        card.write(('bkg_fit').ljust(width))
    elif (fit == 'shape_data'):
        card.write(('sig_fit').ljust(width))
        card.write(('data_fit').ljust(width))
    else:
        print '\n\nWriting datacard in category %s for %s, no valid fit type %s' % (cat, dist, fit)

    card.write('\n')

    ## Process numbers (<= 0 for signal, > 0 for background)
    card.write('process'.ljust(width))
    card.write('0'.ljust(width))
    card.write('1'.ljust(width))
    card.write('\n')

    ## Rate parameters (i.e. total integral) for each process
    card.write('rate'.ljust(width))
    card.write(('%f' % nSig).ljust(width))
    card.write(('%f' % nBkg).ljust(width))
    card.write('\n')

    card.write('----------------------------------------------------------------------------------------------------------------------------------\n')

    ## Systematic uncertainties: only on rate per process
    card.write(('bkg_norm').ljust(width-5))
    card.write(('lnN  ').ljust(2))
    card.write(('-').ljust(width))
    card.write(('9.99').ljust(width))
    card.write('\n')

    ## Final line to add bin-by-bin MC stats uncertainties
    if ('template' in fit or 'rebin' in fit):
        card.write('* autoMCStats 0 0 1\n')

    print 'Wrote out datacard %s\n' % card.name

## End function: def WriteSigBkgBody(card, cat, dist, fit, width):


def WriteGroupBody(card, cat, dist, fit, width, sig_hists, bkg_hists):

    ## Category name (same for each process)
    card.write('bin'.ljust(width))
    for i in range(2, len(sig_hists)+len(bkg_hists)):
        card.write((cat).ljust(width))
    card.write('\n')

    ## Process names
    card.write('process'.ljust(width))
    for i in range(1, len(sig_hists)):
        card.write((sig_hists[i].GetName()).ljust(width))
    for i in range(1, len(bkg_hists)):
        card.write((bkg_hists[i].GetName()).ljust(width))
    card.write('\n')

    ## Process numbers (<= 0 for signal, > 0 for background)
    card.write('process'.ljust(width))
    for i in range(1, len(sig_hists)):
        card.write('%d'.ljust(width+1) % (1 + i - len(sig_hists)))
    for i in range(1, len(bkg_hists)):
        card.write('%d'.ljust(width+1) % i)
    card.write('\n')

    ## Rate parameters (i.e. total integral) for each process
    card.write('rate'.ljust(width))
    for i in range(1, len(sig_hists)):
        card.write(('%f' % sig_hists[i].Integral()).ljust(width))
    for i in range(1, len(bkg_hists)):
        card.write(('%f' % bkg_hists[i].Integral()).ljust(width))
    card.write('\n')

    card.write('----------------------------------------------------------------------------------------------------------------------------------\n')

    ## Systematic uncertainties: only on rate per process
    card.write(('bkg_norm').ljust(width-5))
    card.write(('lnN  ').ljust(2))
    for i in range(1, len(sig_hists)):
        card.write(('-').ljust(width))
    for i in range(1, len(bkg_hists)):
        card.write(('9.99').ljust(width))
    card.write('\n')

    ## Final line to add bin-by-bin MC stats uncertainties
    card.write('* autoMCStats 0 0 1\n')

    print 'Wrote out datacard %s\n' % card.name

## End function: def WriteGroupBody(card, cat, dist, fit, width):
