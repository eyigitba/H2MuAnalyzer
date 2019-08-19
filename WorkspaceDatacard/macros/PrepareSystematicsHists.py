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

import os
import string

import ROOT as R



#============================================
# User-defined settings
#============================================

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

IN_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/'
LABEL  = 'WH_lep_AWB_2019_07_23_signal_sys_v1'
CATE   = 'noZ10_noBtag'

SYS_NAMES = ['noSys', 'JES_up', 'JES_down', 'PU_wgt_up', 'PU_wgt_down', 'IsoMu_SF_up', 'IsoMu_SF_down', 'LepMVA_SF_up', 'LepMVA_SF_down']
YEARS = ['2016', '2017', '2018']



def GetOneFile(file_name, year, sys_name, hists):
    in_file = R.TFile.Open(file_name, 'READ')

    for hist in in_file.GetListOfKeys():
	hist_name = hist.GetName()
	if '_120' in hist_name or '_130' in hist_name: 	continue
	elif '_H_pair_mass_' not in hist_name: 		continue
	elif 'zoomH' not in hist_name: 			continue

	else:
	  hists.append( in_file.Get(hist_name).Clone() )
	  hists[-1].SetName( hist_name.replace('_H_pair_mass_','_%s_%s_H_pair_mass_' %(year, sys_name)).replace('_125', '').replace('H2Mu_gg_', 'ggH_').replace('H2Mu_','') )
	  hists[-1].SetDirectory(0)
	  print hists[-1].GetName()

    in_file.Close()
## End of function GetOneFile(file_name, year, sys_name, hists):



def main():

    hists = []
    print len(hists)
    for year in YEARS:
      for sys_name in SYS_NAMES:
	in_file_name = IN_DIR + year + '/Histograms/' + LABEL + '/files/HADD/' + 'signals_' + CATE + '_' + sys_name + '.root'
	GetOneFile(in_file_name, year, sys_name, hists)

    out_file_name = IN_DIR + 'Run2' + '/Histograms/' + LABEL + '/' + 'Systematics_Run2.root'
    out_file = R.TFile.Open(out_file_name, 'RECREATE')

    out_file.cd()
    print len(hists)
    for hist in hists:
	print hist.GetName()
	hist.Write()



main()

