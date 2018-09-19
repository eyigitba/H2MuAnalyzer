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

#def GetNormForFile():


def GetNormForSample(out_file, samp_name, xsec, in_dir_name, file_list):
# file_list as taken from GenerateBatchScript.py  (name, size)
    if "SingleMu" in samp_name:
	print "no need for renormalization for data"
	return 0

    lumi = 41000.0

    file_chain = R.TChain("dimuons/metadata","chain_"+samp_name)
    for in_file in file_list:
	file_chain.Add(  "root://eoscms.cern.ch/%s/%s" %(in_dir_name, in_file[0]) )

    total_event = 0
    for i in range(file_chain.GetEntries()):
	file_chain.GetEntry(i)
	total_event += file_chain.sumEventWeights
        
#    print "sample = %s " %(samp_name) + "\tnumber of files = %d " %(len(file_list)) + "\tsumEventWeights = %f , xsec = %f , Norm = %f" %(total_event, xsec, xsec * lumi / total_event)

    return xsec * lumi / total_event
