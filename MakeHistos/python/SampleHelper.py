###################################################################
###     macro for calculating normalization for each sample     ###
###                                                             ###
###  currently envisioned to be used in GenerateBatchScript,py  ###
###  out_file for storing normalization info which will be      ###
###  taken by batch jobs as an arg input                        ###
###                                                             ###
###             Xunwu Zuo   01.09.2018                          ###
###################################################################
 

import ROOT as R
import os
import sys

def GetNormForSample(out_file, samp_name, xsec, lumi, in_dir_name, file_list):
# file_list as taken from GenerateBatchScript.py  (name, size)
    if "SingleMu" in samp_name:
	print "No need to renormalize data"
	return 0

    file_chain = R.TChain("dimuons/metadata","chain_"+samp_name)
    for in_file in file_list:
        # print 'Adding file root://eoscms.cern.ch/%s/%s' % (in_dir_name, in_file[0])
        file_tmp = R.TFile.Open( 'root://eoscms.cern.ch/%s/%s' % (in_dir_name, in_file[0]) )
        if (not file_tmp):
            print 'COULD NOT OPEN FILE!!!'
            print 'root://eoscms.cern.ch/%s/%s' % (in_dir_name, in_file[0])
            sys.exit()
	file_chain.Add( 'root://eoscms.cern.ch/%s/%s' % (in_dir_name, in_file[0]) )

    total_event = 0
    for i in range(file_chain.GetEntries()):
	file_chain.GetEntry(i)
	total_event += file_chain.sumEventWeights
        
#    print "sample = %s " %(samp_name) + "\tnumber of files = %d " %(len(file_list)) + "\tsumEventWeights = %f , xsec = %f , Norm = %f" %(total_event, xsec, xsec * lumi / total_event)

    return xsec * lumi / total_event

# assign an ID to the sample based on its parent particle PDG ID   -- XWZ 20.11.2018
# generic notations: 00 for quark, 10 for lepton, 20 for boson
def GetSampleID(name):
    # 0 for data
    if "SingleMu" in samp_name: 	return 0
    # positive value for signal
    elif "H2Mu_ggH" in samp_name: 	return 25
    elif "H2Mu_VBF" in samp_name:  	return 010225  # represented by one up and one down quark
    elif "H2Mu_WH"  in samp_name:  	return 2425
    elif "H2Mu_ZH"  in samp_name:  	return 2325
    elif "H2Mu_ttH" in samp_name:  	return 060625
    # negative value for bkg
    elif "ZJets" in samp_name:  	return -23; # can add more entries if need to distinguish between AMC/MG, 0j/1j/2j samples
    elif "tt_ll" in samp_name or samp_name == "tt": 	  return -0606;
    elif "tW_pos" in samp_name or "tW_neg" in samp_name:  return -0624;
    elif samp_name == "tZq": 		return -062300; # 00 for quark in general
    elif samp_name == "tZW": 		return -062324;
    elif "ttW" in samp_name: 		return -060624;
    elif samp_name == "ttZ": 		return -060623;
    elif samp_name == "ttH": 		return -060625;
    elif "WW" in samp_name:  		return -2424;
    elif "WZ" in samp_name: 		return -2423;
    elif "ZZ" in samp_name:  		return -2323;
    elif "WWW" in samp_name:  		return -242424;
    elif "WWZ" in samp_name: 		return -242423;
    elif "WZZ" in samp_name:  		return -242323;
    elif "ZZZ" in samp_name: 		return -232323;

    return -999;  #just so that this function can be compiled. do not expect -999 to show in any case


