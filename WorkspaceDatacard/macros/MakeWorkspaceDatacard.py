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

from shutil import rmtree
import xml.etree.ElementTree as ET

sys.path.insert(1, '%s/python' % os.getcwd())
sys.path.insert(1, '%s/../FitBackground/python' % os.getcwd())
import FitFunctions    as FF  ## From FitBackground/python/FitFunctions.py
import DataLoader      as DL  ## From WorkspaceDatacard/python/DataLoader.py
import PlotHelper      as PH  ## From WorkspaceDatacard/python/PlotHelper.py
import WorkspaceHelper as WH  ## From WorkspaceDatacard/python/WorkspaceHelper.py
import DatacardHelper  as DH  ## From WorkspaceDatacard/python/DatacardHelper.py


#============================================
# User-defined settings
#============================================

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Configure the channels, distributions, and fits to use
## Chosen from XML files in configs/ directory

#====== AWB CONFIG
# CONFIGS = ['WH_lep_AWB_2019_05_15_TMVA_v3']
# CONFIGS = ['WH_lep_AWB_2019_05_20_v1_shape']

#====== XWZ CONGIF
#CONFIGS = ['WH_2018_3lep_AWB_mass_06_19_med']
#CONFIGS = ['WH_2018_3lep_AWB_BDT_06_19_med']
#CONFIGS = ['WH_1718_3lep_AWB_mass_06_25']
#CONFIGS = ['WH_1718_3lep_AWB_BDT_06_25_med']
#CONFIGS = ['WH_1718_3lep_AWB_Bernstein_mass_07_03', 'WH_1718_3lep_AWB_template_mass_07_03']
#CONFIGS = ['ZH_lep_XWZ_mass_05_23_MVAp04_v1']
# CONFIGS = ['WH_lep_XWZ_2019_05_14_TMVA_out_v1']
#CONFIGS = ['ttH_1718_3lep_AWB_mass_07_04']
#CONFIGS = ['ttH_1718_3lep_AWB_mass_med_07_08']
#CONFIGS = ['WH_Run2_3lep_AWB_mass_07_10']
CONFIGS = ['ttH_Run2_3lep_AWB_mass_07_10']

#============================================
# Main code
#============================================

class WorkspaceAndDatacardMaker:
# Class to make workspace, root files, and datacards needed for
# analytic shape or template limit setting via Higgs Combine

    def __init__(self, _source, _in_dir, _in_file, _out_dir, _cat, _cat_loc, _dist, _min_max, _blind, _rebin, _models):

        self.out_dir = _out_dir ## Output directory for datacards, workspaces, plots, etc.
        self.cat     = _cat     ## Category name for this workspace
        self.dist    = _dist    ## Distribution for signal extraction
        self.blind   = _blind   ## Range to blind in data
        self.rebin   = _rebin   ## Rebin input distribution for optimal S^2/
        self.models  = _models  ## Set of analytic or template models to use

        self.in_data = DL.DataLoader(USER, _source, _in_dir, _in_file, _cat_loc, _dist, _min_max, _rebin)

        self.h_sig_fit  = None  ## Histogram of analytic best-fit to signal shape
        self.h_bkg_fit  = None  ## Histogram of analytic best-fit to background shape
        self.h_data_fit = None  ## Histogram of analytic best-fit to data shape
        # self.nuisances = []     ## List of nuisance parameters (not yet implemented - AWB 13.05.19)

    def Clear(self):
        if not self.in_data    is None: self.in_data.Clear()
        if not self.in_data    is None: del self.in_data
        if not self.h_sig_fit  is None: del self.h_sig_fit
        if not self.h_bkg_fit  is None: del self.h_bkg_fit
        if not self.h_data_fit is None: del self.h_data_fit

        self.out_dir    = ''
        self.cat        = ''
        self.dist       = ''
        self.blind      = []
        self.rebin      = []
        self.models     = []
        self.in_data    = None
        self.h_sig_fit  = None
        self.h_bkg_fit  = None
        self.h_data_fit = None
    

    #-------------------------------------------------------------------------
    # Make workspaces for signal and background template-based model
    #-------------------------------------------------------------------------

    def makeTemplateWorkspace(self, model):
        ## Don't print plots to screen while running (faster)
	R.gROOT.SetBatch(R.kTRUE)

	## Save histograms for template fit.  Two separate types of workspace:
        ##  * 'stack' for a single summed stack for signal and background,
        ##  * 'group' for independent signal and background groups by process
        model_str = model.replace('template', 'rebin') if self.rebin != False else model
        WS = R.TFile.Open(self.out_dir+'/workspace/'+self.cat+'_'+self.dist+'_'+model_str+'.root', 'RECREATE')
        WS.cd()

        if 'stack' in model and self.rebin == False:
            self.in_data.data_hist.Write('data_obs')
            self.in_data.sig_hists[0].Write()
            self.in_data.bkg_hists[0].Write()
            for i in range( min(4, len(self.in_data.hists_sys)) ):
                self.in_data.hists_sys[i].Write()
        ## End conditional: if 'stack' in model and not rebin:
        elif 'group' in model and self.rebin == False:
            self.in_data.data_hist.Write('data_obs')
            for i in range(len(self.in_data.sig_hists)):
                self.in_data.sig_hists[i].Write()
            for i in range(len(self.in_data.bkg_hists)):
                self.in_data.bkg_hists[i].Write()
            for i in range(len(self.in_data.hists_sys)):
                self.in_data.hists_sys[i].Write()
        ## End conditional elif 'group' in model and not rebin:
        elif 'stack' in model and len(self.rebin) == 3:
            self.in_data.data_rebin.Write('data_obs')
            self.in_data.sig_rebin[0].Write()
            self.in_data.bkg_rebin[0].Write()
            for i in range( min(4, len(self.in_data.rebin_sys)) ):
                self.in_data.rebin_sys[i].Write()
        ## End conditional: elif 'stack' in model and rebin:
        elif 'group' in model and len(self.rebin) == 3:
            self.in_data.data_rebin.Write('data_obs')
            for i in range(len(self.in_data.sig_rebin)):
                self.in_data.sig_rebin[i].Write()
            for i in range(len(self.in_data.bkg_rebin)):
                self.in_data.bkg_rebin[i].Write()
            for i in range(len(self.in_data.rebin_sys)):
                self.in_data.rebin_sys[i].Write()
        ## End conditional elif 'group' in model and rebin:
        else:
            print 'Model %s not valid!  Exiting.' % model
            sys.exit()

        WS.Print()

    ## End function: def makeTemplateWorkspace(self, model):


    #-------------------------------------------------------------------------
    # Make workspace for signal and background analytic-fit model
    #-------------------------------------------------------------------------

    def makeShapeWorkspaces(self, sig_mod, sig_ord, sig_frz, bkg_mod, bkg_ord, bkg_frz):
        ## Don't print plots to screen while running (faster)
	R.gROOT.SetBatch(R.kTRUE)

        sig_fit  = FF.FitFunction('sig_fit_%s'  % self.cat, self.in_data.sig_hists[0], sig_mod, sig_ord, self.in_data.min_max, [],         'dimu_mass')
        bkg_fit  = FF.FitFunction('bkg_fit_%s'  % self.cat, self.in_data.bkg_hists[0], bkg_mod, bkg_ord, self.in_data.min_max, [],         'dimu_mass')
        data_fit = FF.FitFunction('data_fit_%s' % self.cat, self.in_data.data_hist,    bkg_mod, bkg_ord, self.in_data.min_max, self.blind, 'dimu_mass')

        FF.DoFit(sig_fit)
        FF.DoFit(bkg_fit)
        FF.DoFit(data_fit)

        self.h_sig_fit  = sig_fit.fit_hist
        self.h_bkg_fit  = bkg_fit.fit_hist
        self.h_data_fit = data_fit.fit_hist

        ## Plot data and fits into a frame
        PH.DrawFits(sig_fit, bkg_fit, data_fit, self.cat+'_'+self.dist+'_'+sig_mod+str(sig_ord)+'_'+bkg_mod+str(bkg_ord), self.out_dir)

        ## After fitting, we freeze all parameters in the signal fit
        WH.FreezeParams(sig_fit, sig_frz)
        ## Also freeze all the fit parameters in the background fit (to be revisited - AWB 09.05.2019)
        WH.FreezeParams(bkg_fit, bkg_frz)
        ## And all the fit parameters in the data fit (also to be revisited - AWB 09.05.2019)
        WH.FreezeParams(data_fit, bkg_frz)

        sig_frz_str = '_frz'.join(f for f in sig_frz)
        if len(sig_frz_str) > 0: sig_frz_str = '_frz'+sig_frz_str
        bkg_frz_str = '_frz'.join(f for f in bkg_frz)
        if len(bkg_frz_str) > 0: bkg_frz_str = '_frz'+bkg_frz_str
        shape_str = sig_mod+str(sig_ord)+sig_frz_str+'_'+bkg_mod+str(bkg_ord)+bkg_frz_str+'_shape'

        ## Create a new workspace for this category
        WS_data = R.RooWorkspace(self.cat+'_'+self.dist+'_'+shape_str+'_data')  ## Using data for background
        WS_MC   = R.RooWorkspace(self.cat+'_'+self.dist+'_'+shape_str+'_MC')    ## Using MC for background
        ## "Data" input to workspace must have name "data_obs"
        ## Rename signal and background models for simplicity
        data_fit.dat  .SetName('data_obs')
        bkg_fit .dat  .SetName('data_obs')
        data_fit.model.SetName('data_fit')
        bkg_fit .model.SetName('bkg_fit')
        sig_fit .model.SetName('sig_fit')
        ## "import" is a keyword so workspace.import() doesn't work in python
        ##   Have to do this instead:
        getattr(WS_data, 'import')(data_fit.dat,   R.RooCmdArg())
        getattr(WS_data, 'import')(data_fit.model, R.RooCmdArg())
        getattr(WS_MC,   'import')(bkg_fit.dat,    R.RooCmdArg())
        getattr(WS_MC,   'import')(bkg_fit.model,  R.RooCmdArg())
        getattr(WS_data, 'import')(sig_fit.model,  R.RooCmdArg())
        getattr(WS_MC,   'import')(sig_fit.model,  R.RooCmdArg())

	# Save workspaces to root file
        WS_data.SaveAs(self.out_dir+'/workspace/'+WS_data.GetName()+'.root')
        WS_data.Print()
        WS_MC.SaveAs(self.out_dir+'/workspace/'+WS_MC.GetName()+'.root')
        WS_MC.Print()

    ## End function: def makeShapeWorkspaces(self, sig_mod, sig_ord, sig_frz, bkg_mod, bkg_ord, bkg_frz):


    #---------------------------------------------------------------------------------------
    # Make datacard with template-based and analytic signal and background fitting functions
    #---------------------------------------------------------------------------------------

    def makeShapeDatacards(self, sig_mod, sig_ord, sig_frz, bkg_mod, bkg_ord, bkg_frz):

        ## Compute the size of the column
        width = max(15, 3 + len(self.cat))

        sig_frz_str = '_frz'.join(f for f in sig_frz)
        if len(sig_frz_str) > 0: sig_frz_str = '_frz'+sig_frz_str
        bkg_frz_str = '_frz'.join(f for f in bkg_frz)
        if len(bkg_frz_str) > 0: bkg_frz_str = '_frz'+bkg_frz_str
        shape_str = sig_mod+str(sig_ord)+sig_frz_str+'_'+bkg_mod+str(bkg_ord)+bkg_frz_str+'_shape'

        card_data = open(self.out_dir+'/datacard/'+self.cat+'_'+self.dist+'_'+shape_str+'_data.txt', 'w')
        DH.WriteHeader    (card_data, self.cat, self.out_dir, self.cat+'_'+self.dist+'_'+shape_str+'_data')
        DH.WriteSigBkgBody(card_data, self.cat, self.dist, 'shape_data', width, self.h_sig_fit.Integral(), self.h_data_fit.Integral())
        
        card_MC = open(self.out_dir+'/datacard/'+self.cat+'_'+self.dist+'_'+shape_str+'_MC.txt', 'w') 
        DH.WriteHeader    (card_MC, self.cat, self.out_dir, self.cat+'_'+self.dist+'_'+shape_str+'_MC')
        DH.WriteSigBkgBody(card_MC, self.cat, self.dist, 'shape_MC', width, self.h_sig_fit.Integral(), self.h_bkg_fit.Integral())

    ## End function: def makeShapeDatacards(self, sig_mod, sig_ord, sig_frz, bkg_mod, bkg_ord, bkg_frz):

    def makeTemplateDatacards(self, model):

        ## Compute the size of the column
        width = max(20, 8 + len(self.cat))

        model_str = model.replace('template', 'rebin') if self.rebin else model

        card = open(self.out_dir+'/datacard/'+self.cat+'_'+self.dist+'_'+model_str+'.txt', 'w')
        DH.WriteHeader(card, self.cat, self.out_dir, self.cat+'_'+self.dist+'_'+model_str)

        if model_str == 'template_stack':
            DH.WriteSigBkgBody(card, self.cat, self.dist, 'template_stack', width, self.in_data.sig_hists[0].Integral(), self.in_data.bkg_hists[0].Integral())
        elif model_str == 'template_group':
            DH.WriteGroupBody(card, self.cat, self.dist, 'template_group', width, self.in_data.sig_hists, self.in_data.bkg_hists, len(self.in_data.hists_sys) > 0)
        elif model_str == 'rebin_stack':
            DH.WriteSigBkgBody(card, self.cat, self.dist, 'rebin_stack', width, self.in_data.sig_rebin[0].Integral(), self.in_data.bkg_rebin[0].Integral())
        elif model_str == 'rebin_group':
            DH.WriteGroupBody(card, self.cat, self.dist, 'rebin_group', width, self.in_data.sig_rebin, self.in_data.bkg_rebin, len(self.in_data.rebin_sys) > 0)
        else:
            print 'Invalid model string %s!!!  Exiting.' % model_str
            sys.exit()

    ## End function: def makeTemplateDatacards(self, model):

    def makeCutAndCountDatacards(self):

	## Compute the size of the column
        width = max(15, 3 + len(self.cat))
	MASS_WINDOW = [120, 130]  # currently hard coded here, could use MASS_WINDOW = self.blind

	card = open(self.out_dir+'/datacard/'+self.cat+'_'+self.dist+'_CutAndCount'+'.txt', 'w')
	DH.WriteCutAndCount(card, self.cat, self.out_dir, self.dist, width, MASS_WINDOW, self.in_data.sig_hists, self.in_data.bkg_hists)
    ## End function: makeCutAndCountDatacards(self):

#####################
##  Main function  ##
#####################

def main():

    print '\n\n*** Inside MakeWorkspaceDatacardNew.py ... and running! ***\n'

    for config in CONFIGS:

        ## Create new output sub-directory
        out_dir = 'out_files/'+config
        if os.path.exists(out_dir):
            delete_dir = raw_input('\n*** Directory %s already exists!!! ***\nType "Y" to delete and continue, "N" to exit.\n\n' % out_dir)
            if delete_dir == 'Y':
                rmtree(out_dir)
                print '\nDeleted %s\n' % out_dir
            else:
                print 'You typed %s, not "Y" - exiting\n' % delete_dir.lower()
        os.makedirs(out_dir+'/workspace')
        os.makedirs(out_dir+'/datacard')
        os.makedirs(out_dir+'/plot')


        ## Parse the XML configuration
        XML = ET.parse('configs/'+config+'.xml').getroot()
        source  = (XML.find('source') .text).replace(' ','')
        in_dir  = (XML.find('in_dir') .text).replace(' ','')
        in_file = (XML.find('in_file').text).replace(' ','')
        print '\nParsed XML from configs/'+config+'.xml specifying:'
        print '  * source  = %s' % source
        print '  * in_dir  = %s' % in_dir
        print '  * in_file = %s' % in_file

        ## Loop over categories in XML
        for cat in XML.find('categories'):
            cat_name = cat.get('name')
            cat_loc  = cat.get('loc')
            print '\nInitializing workspace and datacard maker for category %s (from %s)' % (cat_name, cat_loc)

            ## Loop over discriminants
            for discr in XML.find('discriminants'):

                ###############################################################################
                ##  Parse info from XML (should move to separate function - AWB 13.05.2019)  ##
                ###############################################################################
                dist     = (discr.find('dist').text).replace(' ','')
                blind    = (discr.find('blind').text).replace(' ','')
                min_max  = (discr.find('min_max').text).replace(' ','')
                rebin    = (discr.find('rebin').text).replace(' ','')
                models   = (discr.find('models').text).replace(' ','')

                if (not min_max.startswith('[') or not min_max.endswith(']')):
                    print '\n\nInvalid option for min_max = %s' % min_max
                elif min_max == '[]': min_max = []
                else: min_max = [float(min_max.split(',')[0].replace('[','')), float(min_max.split(',')[1].replace(']',''))]

                if (not blind.startswith('[') or not blind.endswith(']')):
                    print '\n\nInvalid option for blind = %s' % blind
                elif blind == '[]': blind = []
                else: blind = [float(blind.split(',')[0].replace('[','')), float(blind.split(',')[1].replace(']',''))]

                if (rebin != 'False' and not (rebin.startswith('[') and rebin.endswith(']'))):
                    print '\n\nInvalid option for rebin = %s' % rebin
                    sys.exit()
                elif rebin == 'False': rebin = False
                else: rebin = [str(rebin.split(',')[0].replace('[','')), float(rebin.split(',')[1]), float(rebin.split(',')[2].replace(']',''))]
                if rebin != False and len(rebin) != 3:
                    print '\n\nInvalid option rebin with !=3 elements: ['+','.join(str(i) for i in rebin)+']'
                    sys.exit()

                if (not models.startswith('[') or not models.endswith(']')):
                    print '\n\nInvalid option for models = %s' % models
                else: models = models.replace('[','').replace(']','')
                model_list = []
                for mod in models.split(','): model_list.append(mod)
                models = model_list

                ## Create new instance of WorkspaceAndDatacardMaker
                print '  * Using discriminant %s, rebin = %s' % (dist, ('False' if rebin == False else ', '.join(str(i) for i in rebin)))
                print '    Bin range = %s, blinding data in %s' % (' to '.join(str(i) for i in min_max), ' to '.join(str(j) for j in blind))
                print '    Models = %s' % ','.join(k for k in models)

                ## Initialize using specified configuration
                WDM = WorkspaceAndDatacardMaker(source, in_dir, in_file, out_dir, cat_name, cat_loc, dist, min_max, blind, rebin, models)

                ## Create MC template-based workspaces with grouped processes
                ##  or inclusive signal and background stacks
                for model in models:
                    if not 'template' in model: continue
                    WDM.makeTemplateWorkspace(model)
                    WDM.makeTemplateDatacards(model)

                ## Loop over signal and background fit settings

                ###############################################################################
                ##  Parse info from XML (should move to separate function - AWB 13.05.2019)  ##
                ###############################################################################
                if not 'shape' in models:
                    WDM.Clear()
                    continue

                for sig_fit in discr.find('sig_fits'):
                    sig_mod  = (sig_fit.find('shape').text).replace(' ','')
                    sig_ord  = int((sig_fit.find('order').text).replace(' ',''))
                    sig_frz = (sig_fit.find('freeze').text).replace(' ','').replace('[','').replace(']','')
                    sig_frz_list = []
                    for frz in sig_frz.split(','): sig_frz_list.append(frz)
                    sig_frz = sig_frz_list
                    
                    for bkg_fit in discr.find('bkg_fits'):
                        bkg_mod  = (bkg_fit.find('shape').text).replace(' ','')
                        bkg_ord  = int((bkg_fit.find('order').text).replace(' ',''))
                        bkg_frz = (bkg_fit.find('freeze').text).replace(' ','').replace('[','').replace(']','')
                        bkg_frz_list = []
                        for frz in bkg_frz.split(','): bkg_frz_list.append(frz)
                        bkg_frz = bkg_frz_list

                        ## Create analytic shape-based workspaces on data and MC
                        for model in models:
                            if not model == 'shape': continue
                            print '\nCreating shape workspace with signal model %s (order %d), background model %s (order %d)' % (sig_mod, sig_ord, bkg_mod, bkg_ord)
                            WDM.makeShapeWorkspaces(sig_mod, sig_ord, sig_frz, bkg_mod, bkg_ord, bkg_frz)
                            WDM.makeShapeDatacards (sig_mod, sig_ord, sig_frz, bkg_mod, bkg_ord, bkg_frz)

                        ## End loop: for model in models:
                    ## End loop: for bkg_fit in discr.find('bkg_fits'):
                ## End loop: for sig_fit in discr.find('sig_fits'):

                WDM.Clear()

            ## End loop: for discr in XML.find('discriminants'):
        ## End loop: for cat in XML.find('categories'):
    ## End loop: for config in CONFIGS:


## End function: def main()

if __name__ == '__main__':
    main()
