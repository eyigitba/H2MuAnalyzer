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


def GetNormForSample(out_file, samp_name, xsec, file_list):
# file_list as taken from GenerateBatchScript.py  (name, size)
    if "SingleMu" in samp_name:
	print "no need for renormalization for data"
	return

    lumi = 41.0

    file_chain = TChain("dimuons/metadata","chain_"+samp_name)
    for in_file in file_list:
	file_chain.Add(in_file[0])

    total_event = 0
    for i in range(file_chain.GetEntries()):
	total_event += file_chain.sumEventWeights
        
    out_file.write( "\nsample = %s " %(samp_name) + "\tnumber of files = %d " %(len(file_list)) + "\tsumEventWeights = %f , xsec = %f , Norm = %f" %(total_event, xsec, xsec * lumi / total_event) )
    print "wrote: " + "sample = %s " %(samp_name) + "\tnumber of files = %d " %(len(file_list)) + "\tsumEventWeights = %f , xsec = %f , Norm = %f" %(total_event, xsec, xsec * lumi / total_event)


