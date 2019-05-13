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

# SOURCE   = 'abrinke1'
# INDIR    = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/ttH_3l_AWB_2019_01_21_lepMVA_test_v1'
# INFILE_1 = 'StackPlots.root'
# INFILE_2 = 'HADD/histos_Presel2017_'
# CATS     = ['e2mu_looseLepMVA_noZ_ge3j_btag_mass12', '3mu_looseLepMVA_noZ_ge3j_btag_mass12',
#             'e2mu_tightLepMVA_noZ_ge3j_btag_mass12', '3mu_tightLepMVA_noZ_ge3j_btag_mass12']

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
BLIND    = [120, 130]  ## Set to "[]" for unblinded, "[120, 130]" to blind signal region


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
    data_hist  = 0
    h_sig_fit  = 0
    h_bkg_fit  = 0
    h_data_fit = 0
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

    def getSigHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Signal')
        if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_zoom_Net_Sig')
	if SOURCE == 'xzuo':	 return self.in_file_1.Get('dimu_mass_Net_Sig')

    def getBkgHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Bkg')
        if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_zoom_Net_Bkg')
	if SOURCE == 'xzuo':     return self.in_file_1.Get('dimu_mass_Net_Bkg')

    def getDataHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Data')
        if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_zoom_Net_Data')
	if SOURCE == 'xzuo':     return self.in_file_1.Get('dimu_mass_Net_Data')
    
    #-------------------------------------------------------------------------
    # Make workspace for signal and background fitting functions
    #-------------------------------------------------------------------------
    def makeShapeWorkspace(self):
        ## Don't print plots to screen while running (faster)
	R.gROOT.SetBatch(R.kTRUE)

	# Save histograms for template fit
        outfile_MC = R.TFile.Open('out_files/workspace/'+self.cat+'_template_MC.root', "RECREATE")
	outfile_MC.cd()
	sig_temp = self.sig_hist.Clone()
	sig_temp.SetName('sig_fit_'+self.cat)
   	sig_temp.Write()
        bkg_temp = self.bkg_hist.Clone()
	bkg_temp.SetName('bkg_fit_'+self.cat)
	bkg_temp.Write()
	data_temp = self.bkg_hist.Clone()
	data_temp.Add(sig_temp)
	data_temp.SetName('data_obs')
	data_temp.Write()
        outfile_MC.Close()

    ## End function: makeShapeWorkspace(self):


    #-------------------------------------------------------------------------
    # Make datacard with analytic signal and background fitting functions
    #-------------------------------------------------------------------------
    def makeShapeDatacard(self):

        ## Compute the size of the column
        width = 20 + len(self.cat)

        ## Get the binning of the histograms for normalization
        for i in range(1, self.sig_hist.GetNbinsX()+2):
            if self.sig_hist.GetBinLowEdge(i) == MASS_MIN: iMin_sig = i
            if self.sig_hist.GetBinLowEdge(i) == MASS_MAX: iMax_sig = i-1
        for i in range(1, self.bkg_hist.GetNbinsX()+2):
            if self.bkg_hist.GetBinLowEdge(i) == MASS_MIN: iMin_bkg = i
            if self.bkg_hist.GetBinLowEdge(i) == MASS_MAX: iMax_bkg = i-1
        for i in range(1, self.data_hist.GetNbinsX()+2):
            if self.data_hist.GetBinLowEdge(i) == MASS_MIN: iMin_data = i
            if self.data_hist.GetBinLowEdge(i) == MASS_MAX: iMax_data = i-1

        print '\nIn category %s, sig  integral in bins [%d, %d] = %.3f'   % (self.cat, iMin_sig, iMax_sig, self.sig_hist.Integral(iMin_sig, iMax_sig))
        print '  *                    Integral from fit         = %.3f'   % (self.h_sig_fit.Integral())
        print   'In category %s, bkg  integral in bins [%d, %d] = %.3f'   % (self.cat, iMin_bkg, iMax_bkg, self.bkg_hist.Integral(iMin_bkg, iMax_bkg))
        print '  *                    Integral from fit         = %.3f'   % (self.h_bkg_fit.Integral())
        print   'In category %s, data integral in bins [%d, %d] = %.3f'   % (self.cat, iMin_data, iMax_data, self.data_hist.Integral(iMin_data, iMax_data))
        print '  *                    Integral from fit         = %.3f\n' % (self.h_data_fit.Integral())


        card_data = open('out_files/datacard/'+self.cat+'_data.txt', 'w') 
        card_data.write('imax *\n')
        card_data.write('jmax *\n')
        card_data.write('kmax *\n')
        card_data.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_data.write('shapes * * out_files/workspace/'+self.cat+'_data.root '+self.cat+'_data:$PROCESS\n')
        card_data.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_data.write('bin            '+self.cat+'_data\n')
        card_data.write('observation    -1.0\n')
        card_data.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_data.write('bin'.ljust(width)+(self.cat+'_data').ljust(width)+(self.cat+'_data').ljust(width)+'\n')
        card_data.write('process'.ljust(width)+('sig_fit_'+self.cat).ljust(width)+('data_fit_'+self.cat).ljust(width)+'\n')
        card_data.write('process'.ljust(width)+'0'.ljust(width)+'1'.ljust(width)+'\n')
        # card_data.write('rate'.ljust(width)+('1.0').ljust(width)+('1.0').ljust(width)+'\n')
        card_data.write('rate'.ljust(width)+('%f' % self.h_sig_fit.Integral()).ljust(width)+('%f' % self.h_data_fit.Integral()).ljust(width)+'\n')
        card_data.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_data.write((self.cat).ljust(width-5)+('lnN  ').ljust(2)+('-').ljust(width)+('9.99').ljust(width)+'\n')
        # card_data.write('alpha'.ljust(width)+'rateParam'.ljust(width)+self.cat.ljust(width)+('bmodel_'+self.cat).ljust(width)+'1'.ljust(width)+'[0.9,1.1]'.ljust(width)+'\n')
        # card_data.write('lumi'.ljust(width)+'lnN'.ljust(width)+'1.05'.ljust(width)+'1.05'.ljust(width)+'\n')
      
        # # get maximum length of the strings to figure out the width of the columns for the systematics section 
        # pwidth = len('param')
        # for n in self.nuisance_params:
        #     if len(n) > pwidth: pwdith = len(n)
        # pwidth+=4
      
        # # write the systematics section
        # for n in self.nuisance_params:
        #     f.write(n.ljust(pwidth)+'param'.ljust(pwidth)+'0.0'.ljust(pwidth)+'0.1'.ljust(pwidth)+'\n')

        card_MC = open('out_files/datacard/'+self.cat+'_MC.txt', 'w') 
        card_MC.write('imax *\n')
        card_MC.write('jmax *\n')
        card_MC.write('kmax *\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write('shapes * * out_files/workspace/'+self.cat+'_MC.root '+self.cat+'_MC:$PROCESS\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write('bin            '+self.cat+'_MC\n')
        card_MC.write('observation    -1.0\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write('bin'.ljust(width)+(self.cat+'_MC').ljust(width)+(self.cat+'_MC').ljust(width)+'\n')
        card_MC.write('process'.ljust(width)+('sig_fit_'+self.cat).ljust(width)+('bkg_fit_'+self.cat).ljust(width)+'\n')
        card_MC.write('process'.ljust(width)+'0'.ljust(width)+'1'.ljust(width)+'\n')
        # card_MC.write('rate'.ljust(width)+('1.0').ljust(width)+('1.0').ljust(width)+'\n')
        card_MC.write('rate'.ljust(width)+('%f' % self.h_sig_fit.Integral()).ljust(width)+('%f' % self.h_bkg_fit.Integral()).ljust(width)+'\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write((self.cat).ljust(width-5)+('lnN  ').ljust(2)+('-').ljust(width)+('9.99').ljust(width)+'\n')

	

	card_MC = open('out_files/datacard/'+self.cat+'_template_MC.txt', 'w')
        card_MC.write('imax *\n')
        card_MC.write('jmax *\n')
        card_MC.write('kmax *\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write('shapes * * out_files/workspace/'+self.cat+'_template_MC.root '+'$PROCESS\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write('bin            '+self.cat+'_MC\n')
        card_MC.write('observation    -1.0\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write('bin'.ljust(width)+(self.cat+'_MC').ljust(width)+(self.cat+'_MC').ljust(width)+'\n')
        card_MC.write('process'.ljust(width)+('sig_fit_'+self.cat).ljust(width)+('bkg_fit_'+self.cat).ljust(width)+'\n')
        card_MC.write('process'.ljust(width)+'0'.ljust(width)+'1'.ljust(width)+'\n')
        # card_MC.write('rate'.ljust(width)+('1.0').ljust(width)+('1.0').ljust(width)+'\n')
        card_MC.write('rate'.ljust(width)+('%f' % self.sig_hist.Integral()).ljust(width)+('%f' % self.bkg_hist.Integral()).ljust(width)+'\n')
        card_MC.write('----------------------------------------------------------------------------------------------------------------------------------\n')
        card_MC.write((self.cat).ljust(width-5)+('lnN  ').ljust(2)+('-').ljust(width)+('9.99').ljust(width)+'\n')



print '\n\n*** Inside MakeWorkspaceDatacardNew.py ... and running! ***\n'

## WDM = WorkspaceAndDatacardMakerNew('/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/WH_lep_AWB_2019_01_21_lepMVA_test_v1/plots/e2mu_tightLepMVA_mass12/StackPlots.root', 'e2mu_tightLepMVA_mass12') 

for cat in CATS:
    if SOURCE == 'acarnes':  WDM = WorkspaceAndDatacardMakerNew(INDIR+'/'+INFILE_1, '', cat)
    if SOURCE == 'abrinke1': WDM = WorkspaceAndDatacardMakerNew( INDIR+'/plots/'+cat+'/'+INFILE_1,
                                                                 INDIR+'/files/'+INFILE_2+cat+'.root',
                                                                 'cat_'+cat )

    if SOURCE == 'xzuo' :    WDM = WorkspaceAndDatacardMakerNew( INDIR+'/'+INFILE_1, '', cat)

    WDM.makeShapeWorkspace()
    WDM.makeShapeDatacard()

