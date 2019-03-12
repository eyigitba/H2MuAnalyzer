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

# SOURCE   = 'acarnes'
# # INDIR    = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/acarnes/h2mumu/rfiles'
# # INFILE_1 = 'validate_UNBLINDED_dimu_mass_Roch_90_200_categories3_tree_bdt_res_mass_n16_mbg25_unc0_np1p0000_sf1_sd0_sb0_PoissonSignificance_36814_dyAMC-J_minpt20_b-4_sig-xlumi0.root'
# # INFILE_1 = 'validate_UNBLINDED_dimu_mass_KaMu_90_200_categories3_tree_categorization_final_36814_dyAMC-J_minpt10_b-4_sig-xlumi1.root'
# INDIR    = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/acarnes/h2mumu/analysis_note_plots/validation'
# INFILE_1 = 'validate_blinded_dimu_mass_KaMu_110_160_categories2_36814_dyAMC-J_minpt20_b1_sig-xlumi1.root'
# SUBDIR   = 'net_histos'  ## Directory containing signal histogram within root file
# # CATS   = ['c0']
# # CATS   = ['c3', 'c4', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c14', 'c15']
# CATS   = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10', 'c11', 'c12', 'c13', 'c14']

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

    def getBkgHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Bkg')
        # if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_on_Net_Bkg')
        # if SOURCE == 'abrinke1': return self.in_file_2.Get('ttW_H_mass_zoom')
        if SOURCE == 'abrinke1': return self.in_file_2.Get('WZ_3l_H_mass_on')

    def getDataHist(self):
        if SOURCE == 'acarnes':  return self.in_file_1.Get(SUBDIR+'/'+self.cat+'_Net_Data')
        if SOURCE == 'abrinke1': return self.in_file_1.Get('h_H_mass_zoom_Net_Data')

    
    #-------------------------------------------------------------------------
    # Make workspace for signal and background fitting functions
    #-------------------------------------------------------------------------
    def makeShapeWorkspace(self):
        ## Don't print plots to screen while running (faster)
	R.gROOT.SetBatch(R.kTRUE)

        sig_fit  = FF.FitFunction('sig_fit_%s'  % self.cat, self.sig_hist,  'Gaus',    3, [MASS_MIN, MASS_MAX], [],    'dimu_mass')
        bkg_fit  = FF.FitFunction('bkg_fit_%s'  % self.cat, self.bkg_hist,  'BWZRed',  1, [MASS_MIN, MASS_MAX], [],    'dimu_mass')
        data_fit = FF.FitFunction('data_fit_%s' % self.cat, self.data_hist, 'BWZRed',  1, [MASS_MIN, MASS_MAX], BLIND, 'dimu_mass')
        # bkg_fit  = FF.FitFunction('bkg_fit_%s'  % self.cat, self.bkg_hist,  'PolyPlusBWZ', 1, [MASS_MIN, MASS_MAX], [],    'dimu_mass')
        # data_fit = FF.FitFunction('data_fit_%s' % self.cat, self.data_hist, 'PolyPlusBWZ', 1, [MASS_MIN, MASS_MAX], BLIND, 'dimu_mass')
        # bkg_fit  = FF.FitFunction('bkg_fit_%s'  % self.cat, self.bkg_hist,  'Bern', 3, [MASS_MIN, MASS_MAX], [],    'dimu_mass')
        # data_fit = FF.FitFunction('data_fit_%s' % self.cat, self.data_hist, 'Bern', 3, [MASS_MIN, MASS_MAX], BLIND, 'dimu_mass')

        FF.DoFit(sig_fit)
        FF.DoFit(bkg_fit)
        FF.DoFit(data_fit)

        self.h_sig_fit  = sig_fit.fit_hist
        self.h_bkg_fit  = bkg_fit.fit_hist
        self.h_data_fit = data_fit.fit_hist

        ## After fitting, we nail the parameters down so that 
        ##   Higgs Combine knows what the SM signal shape is
        for j in range(len(sig_fit.params[0])):
            for i in range(len(sig_fit.params)):
                sig_fit.params[i][j].setConstant(R.kTRUE)

        ## Create a new workspace for this category
        WS_data = R.RooWorkspace(self.cat+'_data')  ## Using data for background
        WS_MC   = R.RooWorkspace(self.cat+'_MC')    ## Using MC for background
        ## "Data" input to workspace must have name "data_obs"
        ## Rename signal and background models for simplicity
        data_fit.dat.SetName('data_obs')
        bkg_fit .dat.SetName('data_obs')
        data_fit.model.SetName('data_fit_'+self.cat)
        bkg_fit .model.SetName('bkg_fit_'+self.cat)
        sig_fit .model.SetName('sig_fit_'+self.cat)
        ## "import" is a keyword so workspace.import() doesn't work in python
        ##   Have to do this instead:
        getattr(WS_data, 'import')(data_fit.dat,   R.RooCmdArg())
        getattr(WS_data, 'import')(data_fit.model, R.RooCmdArg())
        getattr(WS_MC,   'import')(bkg_fit.dat,    R.RooCmdArg())
        getattr(WS_MC,   'import')(bkg_fit.model,  R.RooCmdArg())
        getattr(WS_data, 'import')(sig_fit.model,  R.RooCmdArg())
        getattr(WS_MC,   'import')(sig_fit.model,  R.RooCmdArg())

	# Save workspaces to root file
        WS_data.SaveAs('out_files/workspace/'+self.cat+'_data.root')
        WS_data.Print()
        WS_MC.SaveAs('out_files/workspace/'+self.cat+'_MC.root')
        WS_MC.Print()


	#-------------------------------------------------------------------
        ## Plot data and fits into a frame
        c_sig  = R.TCanvas('c_%s_sig'  % self.cat, 'c_%s_sig'  % self.cat, 800, 600)
        c_bkg  = R.TCanvas('c_%s_bkg'  % self.cat, 'c_%s_bkg'  % self.cat, 800, 600)
        c_data = R.TCanvas('c_%s_data' % self.cat, 'c_%s_data' % self.cat, 800, 600)

        ## Signal
        c_sig.cd()
        fra_sig = sig_fit.var.frame()
        sig_fit.dat.plotOn(fra_sig)
        sig_fit.model.plotOn(fra_sig, RF.LineColor(R.kBlue), RF.Range('FULL'))
        sig_fit.model.plotOn(fra_sig, RF.Components(sig_fit.arg_sets[0]), RF.LineStyle(R.kDashed), RF.LineColor(R.kGreen), RF.Range('FULL'))
        sig_fit.model.plotOn(fra_sig, RF.Components(sig_fit.arg_sets[1]), RF.LineStyle(R.kDashed), RF.LineColor(R.kRed), RF.Range('FULL'))
        sig_fit.model.plotOn(fra_sig, RF.Components(sig_fit.arg_sets[2]), RF.LineStyle(R.kDashed), RF.LineColor(R.kViolet), RF.Range('FULL'))
        sig_fit.model.paramOn(fra_sig, RF.Layout(0.55, 0.90, 0.90))
        fra_sig.Draw()
        sig_chi = R.TLatex(0.7, 0.3, "#chi^{2} = %.2f" % fra_sig.chiSquare())
        sig_chi.SetNDC(R.kTRUE)
        sig_chi.Draw()
        # self.h_sig_fit.Draw('same')
        c_sig.SaveAs('out_files/png/sig_fit_%s.png' % self.cat)

        ## Background
        c_bkg.cd()
        fra_bkg = bkg_fit.var.frame()
        bkg_fit.dat.plotOn(fra_bkg)
        bkg_fit.model.plotOn(fra_bkg, RF.LineColor(R.kBlue), RF.Range('FULL'))
        # bkg_fit.model.plotOn(fra_bkg, RF.Components(bkg_fit.arg_sets[0]), RF.LineStyle(R.kDashed), RF.LineColor(R.kGreen), RF.Range('FULL'))
        # bkg_fit.model.plotOn(fra_bkg, RF.Components(bkg_fit.arg_sets[1]), RF.LineStyle(R.kDashed), RF.LineColor(R.kRed), RF.Range('FULL'))
        bkg_fit.model.paramOn(fra_bkg, RF.Layout(0.55, 0.90, 0.90))
        fra_bkg.Draw()
        bkg_chi = R.TLatex(0.7, 0.5, "#chi^{2} = %.2f" % fra_bkg.chiSquare())
        bkg_chi.SetNDC(R.kTRUE)
        bkg_chi.Draw()
        # self.h_bkg_fit.Draw('same')
        c_bkg.SaveAs('out_files/png/bkg_fit_%s.png' % self.cat)

        ## Data
        c_data.cd()
        fra_data = data_fit.var.frame()
        data_fit.dat.plotOn(fra_data)
        data_fit.model.plotOn(fra_data, RF.LineColor(R.kBlue), RF.Range('FULL'))
        # data_fit.model.plotOn(fra_data, RF.Components(data_fit.arg_sets[0]), RF.LineStyle(R.kDashed), RF.LineColor(R.kGreen), RF.Range('FULL'))
        # data_fit.model.plotOn(fra_data, RF.Components(data_fit.arg_sets[1]), RF.LineStyle(R.kDashed), RF.LineColor(R.kRed), RF.Range('FULL'))
        data_fit.model.paramOn(fra_data, RF.Layout(0.55, 0.90, 0.90))
        fra_data.Draw()
        data_chi = R.TLatex(0.7, 0.5, "#chi^{2} = %.2f" % fra_data.chiSquare())
        data_chi.SetNDC(R.kTRUE)
        data_chi.Draw()
        # self.h_data_fit.Draw('same')
        c_data.SaveAs('out_files/png/data_fit_%s.png' % self.cat)


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



print '\n\n*** Inside MakeWorkspaceDatacardNew.py ... and running! ***\n'

## WDM = WorkspaceAndDatacardMakerNew('/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/WH_lep_AWB_2019_01_21_lepMVA_test_v1/plots/e2mu_tightLepMVA_mass12/StackPlots.root', 'e2mu_tightLepMVA_mass12') 

for cat in CATS:
    if SOURCE == 'acarnes':  WDM = WorkspaceAndDatacardMakerNew(INDIR+'/'+INFILE_1, '', cat)
    if SOURCE == 'abrinke1': WDM = WorkspaceAndDatacardMakerNew( INDIR+'/plots/'+cat+'/'+INFILE_1,
                                                                 INDIR+'/files/'+INFILE_2+cat+'.root',
                                                                 'cat_'+cat )
    WDM.makeShapeWorkspace()
    WDM.makeShapeDatacard()

