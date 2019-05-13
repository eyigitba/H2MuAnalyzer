#! /usr/bin/env python

########################################################################
##                    Workspace and Datacard Maker                    ##
########################################################################
##  Makes .root file and datacard needed for shape limit setting.     ##
##  Output .root and .txt files to be used with Higgs Combine.        ##
##  Modeled after Viktor Khristenko's code in:                        ##
##  UFDimuAnalysis/python/limit_setting/WorkspaceAndDatacardMaker.py  ##
########################################################################


#============================================
# User-defined settings
#============================================

SOURCE   = 'abrinke1'
INDIR    = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/WH_lep_AWB_2019_01_21_lepMVA_test_v1'
INFILE_1 = 'StackPlots.root'
INFILE_2 = 'HADD/histos_Presel2017_'
CATS     = ['e2mu_looseLepMVA_mt150_noBtag_noZ_mass12', '3mu_looseLepMVA_mt150_noBtag_noZ_mass12',
            'e2mu_tightLepMVA_mt150_noBtag_noZ_mass12', '3mu_tightLepMVA_mt150_noBtag_noZ_mass12']

SOURCE   = 'xzuo'
FILE_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/'
#INDIR    = FILE_DIR + 'WH_ele_loose_ID_loose_iso_loose_mu_iso_v1/plots' 
INDIR    = FILE_DIR + 'VH_selection_2019april/pt10_iso04/WH_ele_high_dimu_pt/plots'
INFILE_1 = 'mass_hists_cut_444.root'
CATS     = ['WH_ele_lepMVA04']


MASS_MIN = 110
MASS_MAX = 160
COUNT_RANGE_MIN = 120
COUNT_RANGE_MAX = 130

signals = ["ttH", "ZH", "WH", "VBF", "ggH"]
bkgs = ["DY",  "WZ", "ttbar", "ZZ", "WW", "ttZ", "tZq", "tW", "triboson", "others"]
BKG_CORRELATED = 'Correlated'

#============================================
# Imports
#============================================

import re
import os
import sys
import string
import argparse
import prettytable

import ROOT as R
import ROOT.RooFit as RF

# sys.path.insert(0, '%s/python' % os.getcwd())
# import BkgModels as pdfB
# import BGSFrun

sys.path.insert(1, '%s/../FitBackground/python' % os.getcwd())
import FitFunctions as FF


#============================================
# Main code
#============================================

class WorkspaceAndDatacardMakerNew:
# Object to make workspace, root files, and datacards needed for
# analytic shape or template limit setting via Higgs Combine.

    in_file_1  = 0
    in_file_2  = 0
    cat        = ''
    sig_hist   = 0
    bkg_hist   = 0
    sig_hists  = {}
    bkg_hists  = {}
    data_hist  = 0
    nuisances = []

    def __init__(self, in_file_name_1, in_file_name_2, category):
        print 'Opening file_1: %s' % in_file_name_1
        self.in_file_1 = R.TFile(in_file_name_1)
        print 'Opening file_2: %s' % in_file_name_2
        self.in_file_2 = R.TFile(in_file_name_2)
        self.cat       = category
        self.sig_hist  = self.getSigHist() .Clone('sig_hist')
        self.bkg_hist  = self.getBkgHist() .Clone('bkg_hist')
        self.data_hist = self.getDataHist().Clone('data_hist')
	for sig in signals:
	    self.sig_hists[sig] = self.getChannelHist(sig).Clone('sig_' + sig)
	for bkg in bkgs:
	    self.bkg_hists[bkg] = self.getChannelHist(bkg).Clone('bkg_' + bkg) 

    def getSigHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Signal')
        if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_zoom_Net_Sig')
	if SOURCE == 'xzuo':	 return self.in_file_1.Get('dimu_mass_Net_Sig')

    def getBkgHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Bkg')
        if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_zoom_Net_Bkg')
	if SOURCE == 'xzuo':     return self.in_file_1.Get('dimu_mass_Net_Bkg')
        # if SOURCE == 'abrinke1': return self.in_file_2.Get('ttW_H_mass_zoom')
        # if SOURCE == 'abrinke1': return self.in_file_2.Get('WZ_3l_H_mass_on')

    def getDataHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Data')
        if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_zoom_Net_Data')
	if SOURCE == 'xzuo':     return self.in_file_1.Get('dimu_mass_Net_Data')
    
    def getChannelHist(self, channel):
	if SOURCE == 'xzuo':     return self.in_file_1.Get('dimu_mass_Net_' + channel)


    #-------------------------------------------------------------------------
    # Make datacard with analytic signal and background fitting functions
    #-------------------------------------------------------------------------
    def makeShapeDatacard(self):

        ## Compute the size of the column
        width = 5 + len(self.cat)

        ## Get the binning of the histograms for normalization, assuming all channels use same binning
	
	iMin = self.sig_hist.FindBin(COUNT_RANGE_MIN+0.1)
	iMax = self.sig_hist.FindBin(COUNT_RANGE_MAX-0.1)

        print '\nIn category %s, sig  integral in bins [%d, %d] = %.3f'   % (self.cat, iMin, iMax, self.sig_hist.Integral(iMin, iMax))
        print   'In category %s, bkg  integral in bins [%d, %d] = %.3f'   % (self.cat, iMin, iMax, self.bkg_hist.Integral(iMin, iMax))
        print   'In category %s, data integral in bins [%d, %d] = %.3f'   % (self.cat, iMin, iMax, self.data_hist.Integral(iMin, iMax))

	print '-----------break down by channels----------'
	for sig in signals:
	    print 'In category %s, signal %s in bins [%d, %d] = %.3f' %(self.cat, sig, iMin, iMax, self.sig_hists[sig].Integral(iMin, iMax))
	for bkg in bkgs:
	    print 'In category %s, bkg %s in bins [%d, %d] = %.3f' %(self.cat, bkg, iMin, iMax, self.bkg_hists[bkg].Integral(iMin, iMax))
	print '\n\n\n'

        card_MC = open('out_files/datacard/'+self.cat+'_Counting_%d_%d_%s.txt'%(COUNT_RANGE_MIN,COUNT_RANGE_MAX,BKG_CORRELATED), 'w') 
	########## there is not much point breaking signal into different channels. Will just use the one inclusive signal #########
	card_MC.write('imax 1\n')
        card_MC.write('jmax %d\n' %len(bkgs))
        card_MC.write('kmax *\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write('bin            '+self.cat+'\n')
        card_MC.write('observation    -1.0\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
	card_MC.write('bin'.ljust(width))
	for channel in ['sig'] + bkgs:
            card_MC.write(self.cat.ljust(width))
	card_MC.write('\nprocess'.ljust(width))
	for channel in ['sig'] + bkgs:
            card_MC.write((channel).ljust(width))
	card_MC.write('\nprocess'.ljust(width))
	for count, channel in enumerate(['sig']+bkgs):
            card_MC.write(('%d'%count).ljust(width))
	card_MC.write('\nrate'.ljust(width))
	card_MC.write(('%f' % self.sig_hist.Integral(iMin,iMax)).ljust(width))
	for bkg in bkgs:
            card_MC.write(('%f' % self.bkg_hists[bkg].Integral(iMin,iMax)).ljust(width))
	card_MC.write('\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
	if BKG_CORRELATED=="Correlated":
	    card_MC.write('prompt_bkg'.ljust(width-5) + ('lnN  ').ljust(2)+('-').ljust(width))
	    for bkg in bkgs:
		if (bkg=='WZ' or bkg=='ZZ' or bkg=='ttZ' or bkg=='tZq' or bkg=='triboson'):
		    card_MC.write('1.2'.ljust(width))
		else:
		    card_MC.write('-'.ljust(width))
	    card_MC.write('\nnonprompt_bkg'.ljust(width-5) + ('lnN  ').ljust(2)+('-').ljust(width))
	    for bkg in bkgs:
		if (bkg=='WZ' or bkg=='ZZ' or bkg=='ttZ' or bkg=='tZq' or bkg=='triboson'):
		    card_MC.write('-'.ljust(width))
		else:
		    card_MC.write('1.4'.ljust(width))
	elif BKG_CORRELATED=="Uncorrelated":
	    for bkg in bkgs:
		card_MC.write((bkg+'yield').ljust(width-5) + ('lnN  ').ljust(2)+('-').ljust(width))
		for channel in bkgs:
		    if (channel != bkg):
			card_MC.write('-'.ljust(width))
		    elif (channel=='WZ' or channel=='ZZ' or bkg=='ttZ' or bkg=='tZq' or bkg=='triboson'):
			card_MC.write('1.2'.ljust(width))
		    else:
			card_MC.write('1.4'.ljust(width))
		card_MC.write('\n')
        #card_MC.write((self.cat).ljust(width-5)+('lnN  ').ljust(2)+('-').ljust(width)+('9.99').ljust(width)+'\n')


print '\n\n*** Inside MakeWorkspaceDatacardNew.py ... and running! ***\n'

## WDM = WorkspaceAndDatacardMakerNew('/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/WH_lep_AWB_2019_01_21_lepMVA_test_v1/plots/e2mu_tightLepMVA_mass12/StackPlots.root', 'e2mu_tightLepMVA_mass12') 

for cat in CATS:
    if SOURCE == 'acarnes':  WDM = WorkspaceAndDatacardMakerNew(INDIR+'/'+INFILE_1, '', cat)
    if SOURCE == 'abrinke1': WDM = WorkspaceAndDatacardMakerNew( INDIR+'/plots/'+cat+'/'+INFILE_1,
                                                                 INDIR+'/files/'+INFILE_2+cat+'.root',
                                                                 'cat_'+cat )

    if SOURCE == 'xzuo' :    WDM = WorkspaceAndDatacardMakerNew( INDIR+'/'+INFILE_1, '', cat)

    WDM.makeShapeDatacard()

