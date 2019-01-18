#! /usr/bin/env python

## Basic python includes
import math

## ROOT includes
import ROOT as R
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

LUMI = 36  ## Luminosity (in fb-1) for scaling cross section

# author   = 'Regnery'
# in_dir   = 'NLO/data/Regnery'
# in_file  = 'NLO.mu13tev-higgs-110-160-XYZ-results.dat'
# in_cats  = ['full', 'cc', 'ncnc', '1jet']
# out_file = R.TFile('hists/NLO_Regnery_2017.root', 'recreate')

# author   = 'Bourilkov'
# in_dir   = 'NLO/data/Bourilkov'
# in_file  = 'h2mu-dsdm-13tev-xs-lux-XYZ-hs.dat'
# in_cats  = ['full', 'cc', 'ncnc', '1jet']
# out_file = R.TFile('hists/NLO_Bourilkov_2017.root', 'recreate')

author   = 'Bourilkov'
in_dir   = 'NNLO/data'
in_file  = 'h2mu-dsdm-13tev-xs-lux-XYZ-nnlo-hp.dat'
in_cats  = ['full', 'cc', 'ncnc', '1jet', '2jet']
out_file = R.TFile('hists/NNLO_Bourilkov_2017.root', 'recreate')

print '\nInside convertDatToTH1D.py'

for cat in in_cats:
    file_name = in_dir+'/'+in_file.replace('XYZ', cat)
    print 'Opening file %s' % file_name
    arrs = []

    ## In Brendan Regnery's files, only a portion corresponds to the Z mass section
    start_reading = (author != 'Regnery')
    end_reading   = False

    with open(file_name) as iFile:
        for line in iFile:
            ## In Brendan Regnery's files, the 'l-/lep. pT' section comes after the Z mass section
            if 'lep. pT' in line:
                end_reading = True

            if start_reading and not end_reading:
                arr = line.split()
                if len(arr) == 0: ## Skip empty lines
                    continue
                for i in range(len(arr)):
                    arr[i] = float(arr[i])
                arrs.append(arr)

            if not start_reading:
                ## In Brendan Regnery's files, the 'Q_ll Invaria' section is the Z mass section
                if 'Q_ll' in line: start_reading = True
            
            
    nBins = len(arrs)
    binW  = (arrs[-1][0] - arrs[0][0]) / (nBins - 1)
    xMin  = arrs[0][0]  - (binW/2)
    xMax  = arrs[-1][0] + (binW/2)

    out_file.cd()
    hist = R.TH1D('%s_xsec' % cat, '%s_xsec' % cat, nBins, xMin, xMax)
    for i in range(nBins):
        hist.SetBinContent(i+1, arrs[i][1])
        if author == 'Bourilkov':
            hist.SetBinError(i+1, arrs[i][2])
        if author == 'Regnery':
            hist.SetBinError(i+1, math.sqrt( pow(arrs[i][2], 2) + pow( (arrs[i][4]+arrs[i][3])/2, 2) ) )
    hist.Write()

    hist.Scale(LUMI*1000)  ## Scale xsec (in pb) by integrated luminosity (in fb-1)
    hist.SetName ('%s_%dfb' % (cat, LUMI))
    hist.SetTitle('%s_%dfb' % (cat, LUMI))
    hist.Write()


print '\nWriting output file %s' % out_file.GetName()
out_file.Write()
out_file.Close()
print '\nDone!\n'
