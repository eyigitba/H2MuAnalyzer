#! /usr/bin/env python

import ROOT as R

import os
import subprocess

## User-defined constants
MAX_FILE = 3  ## Maximum number of input files to process
MAX_EVT  = -1  ## Maximum number of events to process
PRT_EVT  = 1000  ## Print to screen every Nth event

def main():

###################
## Initialize files
###################

    ## Location of input files in eos on lxplus : SingleMuon data from 2017
    file_names = []
    store  = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/'
    in_dir = 'SingleMuon/SingleMu_2017F/180802_164117/0000/' ## 2017 v4 - https://github.com/UFLX2MuMu/Ntupliser/releases/tag/prod-v4.0

    ## List all files in directory

    for in_file_name in subprocess.check_output(['ls', store+in_dir]).splitlines():
        if not '.root' in in_file_name: continue
        if (len(file_names) >= MAX_FILE): break
        file_names.append(store+in_dir+in_file_name)
        print 'Opening file: '+store+in_dir+file_names[-1]

    ## Chain together trees from input files
    in_chains = []
    for i in range(len(file_names)):
        in_chains.append( R.TChain('dimuons/tree') )
        in_chains[i].Add( file_names[i] )

    ## Set output directories (create if they do not exist)
    if not os.path.exists('plots/png/ReadNTupleChain_py/'):
        os.makedirs('plots/png/ReadNTupleChain_py/')
        
    out_file = R.TFile('plots/ReadNTupleChain_py.root', 'recreate')
    png_dir  = 'plots/png/ReadNTupleChain_py/'

    
#############
## Histograms
#############

    ## Histogram bins: [# of bins, minimum, maximum]
    mass_bins = [80, 110, 150]

    ## Book 1D histogram
    ## Important to use '1D' instead of '1F' when dealing with large numbers of entries, and weighted events (higher precision)
    h_mass = R.TH1D('h_mass', 'h_mass', mass_bins[0], mass_bins[1], mass_bins[2])


#############
## Event loop
#############
                
    iEvt = -1
    for ch in in_chains:
        
        if iEvt > MAX_EVT and MAX_EVT > 0: break
                
        for jEvt in range(ch.GetEntries()):
            iEvt += 1
            
            if iEvt > MAX_EVT and MAX_EVT > 0: break
            if iEvt % PRT_EVT is 0: print 'Event #', iEvt

            ch.GetEntry(jEvt)

            nMuPairs = int(ch.nMuPairs)
            
            # print '\nIn event %d we find %d muon pairs'  % (iEvt, nMuPairs)

            ## Loop over RECO muon pairs
            for iPair in range(nMuPairs):

                ## Get variable quantities out of the tree
                mass = ch.muPairs[iPair].mass

                # print '  * mass(%d) = %.2f GeV' % (iPair, mass)

                ## Fill the histogram
                h_mass.Fill(mass)

            ## End loop over RECO muon pairs (iPair)
                
        ## End loop over events in chain (jEvt)
    ## End loop over chains (ch)


######################
## Save the histograms
######################

    out_file.cd()

    canv = R.TCanvas('canv')
    canv.cd()

    h_mass.SetLineWidth(2)
    h_mass.SetLineColor(R.kBlue)
    h_mass.Write()

    h_mass.Draw()
    canv.SaveAs(png_dir+h_mass.GetName()+'.png')

    del out_file


## Define 'main' function as primary executable
if __name__ == '__main__':
    main()
