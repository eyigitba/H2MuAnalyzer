#
# grabbing systematics. 
# SysConfig contains general systematics from yellow report
# getSys get files 
#
# general recommandations for datacard conventions can be found at
# https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationConventions


import os
import sys

import ROOT as R
import ROOT.RooFit as RF

sys.path.insert(1, '%s/../FitBackground/python' % os.getcwd())



# function to convert names from UF convention to HiggsCombine convention
def SysNameConversion(sys_name): 
    if   sys_name == 'JES': 	  return 'CMS_scale_j'
    elif sys_name == 'PU_wgt':    return 'CMS_pu'
    elif sys_name == 'IsoMu_SF':  return 'CMS_eff_m_trg'
    elif sys_name == 'LepMVA_SF': return 'CMS_eff_m_mva'
    elif sys_name == 'Btag_SF':   return 'CMS_eff_b'
    else:	return 'invalid'
## End function SysNameConversion(sys_name):




#################################################
##  class to store values for all systematics  ##
#################################################
class SystematicsConfig:
    def __init__(self, _exp_sys, _hists_sys):

	self.signals    = [] # complete list of relavant signals,     it is a superset of what is in the xml config files
	self.sys_names  = [] # complete list of relavant systematics, it is a superset of what is in the xml config files
	self.sys_values = {}

	## should add the possibility to run systematics on inclusive signal -- XWZ 2019.08.18

	if _exp_sys[0] == 'Norminal':  self.sys_names = ['Norminal']
	else:  
	  self.LoadSystematics(_exp_sys, _hists_sys)

    def LoadSystematics(self, exp_sys, hists_sys):

	################################
	##  load all signal channels  ##
	################################	
	for year in ['2016', '2017', '2018']:
	  self.signals.append('ggH_'     + year)
	  self.signals.append('VBF_'     + year)
	  self.signals.append('WH_pos_'  + year)
	  self.signals.append('WH_neg_'  + year)
	  self.signals.append('ZH_'      + year)
	  self.signals.append('ttH_'     + year)


	##########################################
	## get values for different systematics ##
	##########################################

	#-----------------------------------
	# luminosity systematics
	# from Lumi POG https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiComb
	#-----------------------------------

	# luminosity sum of the uncorrelated parts for 2016
	self.sys_names.append('lumi_2016_uncor')
	self.sys_values['lumi_2016_uncor'] = {}
	for signal in self.signals:
	    if '2016' in signal: self.sys_values['lumi_2016_uncor'][signal] = '1.022'
	    else:		 self.sys_values['lumi_2016_uncor'][signal] = '-'

	# luminosity sum of the uncorrelated parts for 2017	
	self.sys_names.append('lumi_2017_uncor')
        self.sys_values['lumi_2017_uncor'] = {}
        for signal in self.signals:
            if '2017' in signal: self.sys_values['lumi_2017_uncor'][signal] = '1.020'
	    else:		 self.sys_values['lumi_2017_uncor'][signal] = '-'

	# luminosity sum of the uncorrelated parts for 2018      
        self.sys_names.append('lumi_2018_uncor')
        self.sys_values['lumi_2018_uncor'] = {}
        for signal in self.signals:
            if '2018' in signal: self.sys_values['lumi_2018_uncor'][signal] = '1.015'
	    else:		 self.sys_values['lumi_2018_uncor'][signal] = '-'

	# X-Y factorization      
        self.sys_names.append('lumi_XY')
        self.sys_values['lumi_XY'] = {}
        for signal in self.signals:
            if   '2016' in signal: self.sys_values['lumi_XY'][signal] = '1.009'
	    elif '2017' in signal: self.sys_values['lumi_XY'][signal] = '1.008'
	    elif '2018' in signal: self.sys_values['lumi_XY'][signal] = '1.020'
	    else:		   self.sys_values['lumi_XY'][signal] = '-'

	# length scale 17-18      
        self.sys_names.append('lumi_LengScale')
        self.sys_values['lumi_LengScale'] = {}
        for signal in self.signals:
	    if   '2017' in signal: self.sys_values['lumi_LengScale'][signal] = '1.003'
	    elif '2018' in signal: self.sys_values['lumi_LengScale'][signal] = '1.002'
	    else:		   self.sys_values['lumi_LengScale'][signal] = '-'

	# Beam-beam deflection 15-17     
        self.sys_names.append('lumi_BeamDef')
        self.sys_values['lumi_BeamDef'] = {}
        for signal in self.signals:
	    if   '2016' in signal: self.sys_values['lumi_BeamDef'][signal] = '1.004'
	    elif '2017' in signal: self.sys_values['lumi_BeamDef'][signal] = '1.004'
	    else: 		   self.sys_values['lumi_BeamDef'][signal] = '-'           

	# dynamic beta 15-17      
        self.sys_names.append('lumi_DynaBeta')
        self.sys_values['lumi_DynaBeta'] = {}
        for signal in self.signals:
            if   '2016' in signal: self.sys_values['lumi_DynaBeta'][signal] = '1.005'
	    elif '2017' in signal: self.sys_values['lumi_DynaBeta'][signal] = '1.005'
	    else:		   self.sys_values['lumi_DynaBeta'][signal] = '-'

 	# Beam current calibration 17-18      
        self.sys_names.append('lumi_BeamCal')
        self.sys_values['lumi_BeamCal'] = {}
        for signal in self.signals:
            if   '2017' in signal: self.sys_values['lumi_BeamCal'][signal] = '1.003'
	    elif '2018' in signal: self.sys_values['lumi_BeamCal'][signal] = '1.002'
	    else:		   self.sys_values['lumi_BeamCal'][signal] = '-'

	# Ghosts and satellites 15-17      
        self.sys_names.append('lumi_GhostSat')
        self.sys_values['lumi_GhostSat'] = {}
        for signal in self.signals:
            if   '2016' in signal: self.sys_values['lumi_GhostSat'][signal] = '1.004'
	    elif '2017' in signal: self.sys_values['lumi_GhostSat'][signal] = '1.001'
	    else:		   self.sys_values['lumi_GhostSat'][signal] = '-'



	#---------------------------------
	# theoretical and phenomenological systematics
	# from CERN Yellow Report CERN-2017-002-M
	#---------------------------------

	# Branching Rato
	self.sys_names.append('BR_hmm')
	self.sys_values['BR_hmm'] = {}
	for signal in self.signals:
	    self.sys_values['BR_hmm'][signal] = '1.017'
	
	# cross sections
	self.sys_names.append('QCDscale_ggH')
	self.sys_values['QCDscale_ggH'] = {}
	for signal in self.signals:
	    if 'ggH' in signal: self.sys_values['QCDscale_ggH'][signal] = '1.039'
	    else:		self.sys_values['QCDscale_ggH'][signal] = '-'

	#
	self.sys_names.append('QCDscale_VBF')
        self.sys_values['QCDscale_VBF'] = {}
        for signal in self.signals:
            if 'VBF' in signal: self.sys_values['QCDscale_VBF'][signal] = '0.997/1.004'
            else:               self.sys_values['QCDscale_VBF'][signal] = '-'

	#
	self.sys_names.append('QCDscale_WH')
        self.sys_values['QCDscale_WH'] = {}
        for signal in self.signals:
            if 'WH' in signal: self.sys_values['QCDscale_WH'][signal] = '0.993/1.005'  # correlated for WH_pos and WH_neg
            else:              self.sys_values['QCDscale_WH'][signal] = '-'

	#
	self.sys_names.append('QCDscale_ZH')
        self.sys_values['QCDscale_ZH'] = {}
        for signal in self.signals:
            if 'ZH' in signal: self.sys_values['QCDscale_ZH'][signal] = '0.969/1.038'
            else:              self.sys_values['QCDscale_ZH'][signal] = '-'

	#
	self.sys_names.append('QCDscale_ttH')
        self.sys_values['QCDscale_ttH'] = {}
        for signal in self.signals:
            if 'ttH' in signal: self.sys_values['QCDscale_ttH'][signal] = '0.908/1.058'
            else:               self.sys_values['QCDscale_ttH'][signal] = '-'

	#
	self.sys_names.append('pdf_Higgs_gg')
        self.sys_values['pdf_Higgs_gg'] = {}
        for signal in self.signals:
            if 'ggH' in signal: self.sys_values['pdf_Higgs_gg'][signal] = '1.032'
            else:               self.sys_values['pdf_Higgs_gg'][signal] = '-'

	#
	self.sys_names.append('pdf_Higgs_qqbar')
        self.sys_values['pdf_Higgs_qqbar'] = {}
        for signal in self.signals:
            if  'VBF' in signal: self.sys_values['pdf_Higgs_qqbar'][signal] = '0.979/1.021'
	    elif 'WH' in signal: self.sys_values['pdf_Higgs_qqbar'][signal] = '1.019'
	    elif 'ZH' in signal: self.sys_values['pdf_Higgs_qqbar'][signal] = '1.016'
            else:                self.sys_values['pdf_Higgs_qqbar'][signal] = '-'

	#
	self.sys_names.append('pdf_Higgs_ttH')
        self.sys_values['pdf_Higgs_ttH'] = {}
        for signal in self.signals:
            if 'ttH' in signal: self.sys_values['pdf_Higgs_ttH'][signal] = '1.036'
            else:               self.sys_values['pdf_Higgs_ttH'][signal] = '-'


        #---------------------------------
        # experimental systematics
        # from up/down shifted histograms
        #---------------------------------
	for sys_name in exp_sys:
	    self.getSystematicShifts(sys_name, hists_sys)

    ## End function: LoadSystematics(self):



    ###############################################
    ##  function to calculate systematic values  ##
    ##  from up/down histogram yields            ##
    ###############################################

    def getSystematicShifts(self, sys_name, hists_sys): # hists_sys should be a dictionary, intended to be the same as in DataLoader
	HC_sys_name = SysNameConversion(sys_name)
	if HC_sys_name == 'invalid':
	    print 'Not valid systematic name: %s'  %sys_name
	    return

	hist_keys = hists_sys.keys()
	print hist_keys
	print '\n\n'
	for year in ['2016', '2017', '2018']:
	  self.sys_names.append(HC_sys_name + '_' + year)
	  self.sys_values[HC_sys_name + '_' + year] = {}

	  for signal in self.signals:
	    if year not in signal:	
	    # systematics in this year does not effect signals in other years
		self.sys_values[HC_sys_name + '_' + year][signal] = '-'  

	    elif (signal + '_noSys') not in hist_keys or (signal + '_' + sys_name + '_up') not in hist_keys or (signal + '_' + sys_name + '_down') not in hist_keys:
	    # some hists do not exist in the given collection
		print 'cannot find %s hists\n' %signal
		self.sys_values[HC_sys_name + '_' + year][signal] = '-'
		
	    else:  # all shifted hists found for this signal
	        noSys_intg = hists_sys[signal + '_noSys'].Integral()
		up_intg	   = hists_sys[signal + '_' + sys_name + '_up'].Integral()
		down_intg  = hists_sys[signal + '_' + sys_name + '_down'].Integral()

		if hists_sys[signal + '_noSys'].GetEntries() < 100:   # hist do not have enough events to make a reasonable signal modeling
		    print 'hist %s only have %d entries' %(signal, hists_sys[signal + '_noSys'].GetEntries())
		    self.sys_values[HC_sys_name + '_' + year][signal] = '-'
		else:     # found proper values to set
		    self.sys_values[HC_sys_name + '_' + year][signal] = '%.3f/%.3f' %( down_intg/noSys_intg, up_intg/noSys_intg )

	    # End of if, elif, else
	  # End of for signal in self.signals:
	# End of for year in ['2016', '2017', '2018']:

    ## End function: getSystematicShifts(self, sys_name, hists_sys):











