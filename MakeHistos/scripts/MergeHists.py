#! /usr/bin/env python

######################################################
###    Merge different histograms from the same    ###
###    MC process, weigthed by effective entries   ###
######################################################

## Basic python includes
import os
import sys
import copy

## ROOT includes
import ROOT as R

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Settings for this stack-drawing job
if USER == 'abrinke1':
    # YEARS    = ['2016', '2017', '2018']  ## Dataset years (2016, 2017, 2018)
    YEARS    = ['2018']
    PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/YEAR/Histograms'

    SAMPLES = {}  ## Sets of samples to merge
    SAMPLES['H2Mu_VBF_120'] = ['H2Mu_VBF_120_NLO_1', 'H2Mu_VBF_120_NLO_2']
    SAMPLES['H2Mu_VBF']     = ['H2Mu_VBF_125_NLO_1', 'H2Mu_VBF_125_NLO_2']
    SAMPLES['H2Mu_VBF_130'] = ['H2Mu_VBF_130_NLO_1', 'H2Mu_VBF_130_NLO_2']
    SAMPLES['ZJets']        = ['ZJets_AMC', 'ZJets_MG_1', 'ZJets_MG_2']
    SAMPLES['ZJets_hiM']    = ['ZJets_hiM_AMC', 'ZJets_hiM_MG']
    SAMPLES['tt_ll']        = ['tt_ll_POW', 'tt_ll_MG', 'tt_ll_AMC']
    SAMPLES['ttZ']          = ['ttZ', 'ttZ_1', 'ttZ_2']

    # LABEL      = 'WH_lep_AWB_2019_07_08_signal_v1'  ## Sub-folder within PLOT_DIR containing histograms
    # CATEGORIES = ['3lep_looseLepMVA_noZ5_noBtag',   ## Categories for which to merge histograms
    #               '3lep_medLepMVA_noZ10_noBtag',
    #               '3lep_hiPt_lepW20_medLepMVA_noZ10_noBtag']

    # LABEL      = 'WH_lep_AWB_2019_07_08_sideband_v1'      ## Sub-folder within PLOT_DIR containing histograms
    # CATEGORIES = ['3lep_allMass_medLepMVA_noZ10_noBtag',  ## Categories for which to merge histograms 
    #               '3lep_allMass_medLepMVA_onZ10_noBtag',
    #               '3mu_allMass_medLepMVA_noZ10_noBtag',
    #               '3mu_allMass_medLepMVA_onZ10_noBtag',
    #               'e2mu_allMass_medLepMVA_noZ10_noBtag',
    #               'e2mu_allMass_medLepMVA_onZ10_noBtag']

    LABEL      = 'ttH_3l_AWB_2019_07_12_signal_v1'  ## Sub-folder within PLOT_DIR containing histograms
    # CATEGORIES = ['3lep_looseLepMVA_ge2j_btag',     ## Categories for which to merge histograms
    #               '3lep_medLepMVA_noZ10_ge2j_btag',
    #               '3lep_hiPt_lepW20_medLepMVA_noZ10_ge2j_btag',
    #               '3lep_hiPt_lep20_tightLepMVA_noZ10_btag']
    CATEGORIES = ['3lep_looseLepMVA_ge2j_btag']

    # LABEL      = 'ttH_3l_AWB_2019_07_08_sideband_v1'         ## Sub-folder within PLOT_DIR containing histograms
    # CATEGORIES = ['3lep_allMass_medLepMVA_noZ10_ge2j_btag',  ## Categories for which to merge histograms 
    #               '3lep_allMass_medLepMVA_onZ10_ge2j_btag',
    #               '3mu_allMass_medLepMVA_noZ10_ge2j_btag',
    #               '3mu_allMass_medLepMVA_onZ10_ge2j_btag',
    #               'e2mu_allMass_medLepMVA_noZ10_ge2j_btag',
    #               'e2mu_allMass_medLepMVA_onZ10_ge2j_btag']


elif USER == 'xzuo':
    PLOT_DIR = 'NONE'
elif USER == 'bortigno':
    PLOT_DIR = 'NONE'
else: print 'Invalid USER = %s' % USER


def MergeOneCat(category, year):

    print '\n\n********************************'
    print '    Inside MergeOneCat(%s, %s)' % (category, year)
    print '********************************'

    #######################################
    ###  Set up input and output files  ###
    #######################################

    if USER == 'abrinke1': in_file_dir = PLOT_DIR.replace('YEAR', year)+'/'+LABEL+'/files/HADD'
    else:                  in_file_dir = PLOT_DIR.replace('YEAR', year)+'/'+LABEL+'/files/sum'

    in_file_name  = in_file_dir+'/histos_PreselRun2_%s.root' % category  ## File with input histograms
    out_file_name = in_file_name.replace('.root', '_merged.root')
    
    print 'Opening input file %s' % in_file_name
    in_file = R.TFile.Open(in_file_name, 'READ')

    print 'Creating output file %s' % out_file_name
    out_file = R.TFile.Open(out_file_name, 'RECREATE')
    out_file.cd()

    global_dists = []
    
    for new_samp in SAMPLES.keys():

        old_samps     = copy.deepcopy(SAMPLES[new_samp])
        old_samps_mod = copy.deepcopy(SAMPLES[new_samp])

        print '\n\n*** Attempting to merge '+', '.join(old_samps)+' into '+new_samp+' ***\n'

        dists      = []  ## All distributions contained in the file
        hists_excl = []  ## Excluded (non-existent) histograms
        entries    = {}  ## Maximum number of effective entries corresponding to each old sample
        entries['TOTAL'] = 0

        ## Pick the sample with the longest name as the unique "key" sample
        ## At least one histogram from this sample must exist in the file
        in_file.cd()
        key_samp = ''
        for samp in old_samps:
            samp_exists = False
            for hist in R.gDirectory.GetListOfKeys():
                ## See if histogram name starts with the "key" sample name
                hist_name = hist.GetName()
                if hist_name.startswith(samp+'_'):
                    samp_exists = True
                    break
            if not samp_exists:
                print 'Sample %s does not exist in this file: skipping' % samp
                old_samps_mod.remove(samp)
            elif len(samp) > len(key_samp):
                key_samp = samp

        if len(key_samp) > 0:
            print '\nSet "key" sample to %s' % key_samp
        else:
            print '\nNone of the input samples intended for merging exist in this file.  Moving on.'
            continue

        old_samps = copy.deepcopy(old_samps_mod)
        for samp in old_samps:
            entries[samp] = 0

        ## Loop over all the histograms in the file
        ## Histogram name composed of sample+'_'+distribution
        for hist in R.gDirectory.GetListOfKeys():
            ## See if histogram name starts with some sample name
            hist_name = hist.GetName()
            found_dist = False
            for samp in old_samps:
                ## If histogram starts with "key" sample name, then it exists
                if hist_name.startswith(key_samp+'_'):
                    dists.append(hist_name.replace(key_samp+'_', ''))
                    found_dist = True
                    break
                ## If histogram starts with another sample name, but not any of the *other*
                ##    sample names, then the histogram name = samp_dist
                elif hist_name.startswith(samp+'_'):
                    found_other_dist = False
                    for other_samp in old_samps:
                        if hist_name.startswith(other_samp+'_') and not other_samp is samp:
                            found_other_dist = True
                    if not found_other_dist:
                        dists.append(hist_name.replace(samp+'_', ''))
                        found_dist = True
                        break

            if not found_dist: continue
            ## Check if this histogram has the maximum number of effective entries
            for samp in old_samps:
                try:
                    entries[samp] = max(entries[samp], in_file.Get(samp+'_'+dists[-1]).GetEffectiveEntries())
                except:
                    print '\nIt appears that histogram %s does not exist: skipping' % (samp+'_'+dists[-1])
                    hists_excl.append(samp+'_'+dists[-1])

        ## End loop: for hist in R.gDirectory.GetListOfKeys():

        for samp in old_samps:
            entries['TOTAL'] += entries[samp]

        print '\nFound %d distributions, where the maximum number of effective entries was:' % len(dists)
        for samp in old_samps:
            print '  * %s : %.2f, weight by %.4f' % (samp, entries[samp], entries[samp] / entries['TOTAL'])
            ## Reset the "key" sample to the one with the maximum number of effective entries
            if entries[samp] > entries[key_samp]:
                key_samp = samp
                print '    - Reset "key" sample to %s' % key_samp

        print '\nLooping over distributions and adding together component histograms'
        out_file.cd()
        for dist in dists:
            if not dist in global_dists: global_dists.append(dist)
            # print '  * Looking at distribution %s' % dist
            found_one = False
            for samp in old_samps:
                if samp+'_'+dist in hists_excl:
                    continue
                if not found_one:
                    out_hist = in_file.Get(samp+'_'+dist).Clone(new_samp+'_'+dist)
                    out_hist.Scale(entries[samp] / entries['TOTAL'])
                    out_hist.SetTitle(new_samp+'_'+dist)
                    found_one = True
                else:
                    out_hist.Add( in_file.Get(samp+'_'+dist), entries[samp] / entries['TOTAL'] )

            out_file.cd()
            out_hist.Write()

    ## End loop: for new_samp in SAMPLES.keys():

    print '\n\nEnd loop over samples to merge: copying other histograms over unchanged'
    in_file.cd()
    for hist in R.gDirectory.GetListOfKeys():
        ## See if histogram name starts with one of the sample names
        hist_name = hist.GetName()
        merged_sample = False
        for new_samp in SAMPLES.keys():
            old_samps = SAMPLES[new_samp]
            for samp in old_samps:
                if hist_name.startswith(samp+'_') and hist_name.replace(samp+'_', '') in global_dists:
                    merged_sample = True
        if not merged_sample:
            out_file.cd()
            copy_hist = in_file.Get(hist_name).Clone()
            copy_hist.Write()

    ## Double-check you didn't loose anything you wanted to keep
    in_file.cd()
    in_file_hists = []
    for in_hist in R.gDirectory.GetListOfKeys():
        in_file_hists.append(in_hist.GetName())
    out_file.cd()
    out_file_hists = []
    for out_hist in R.gDirectory.GetListOfKeys():
        out_file_hists.append(out_hist.GetName())

    print '\nChecking for input file histograms missing from output file'
    for in_hist in in_file_hists:
        not_found_at_all = True
        if in_hist in out_file_hists:
            not_found_at_all = False
        for new_samp in SAMPLES.keys():
            old_samps = SAMPLES[new_samp]
            for samp in old_samps:
                if in_hist.startswith(samp+'_') and in_hist.replace(samp, new_samp) in out_file_hists:
                    not_found_at_all = False
        if not_found_at_all:
            print in_hist

    # print '\nChecking for output file histograms missing from input file'
    # for out_hist in out_file_hists:
    #     if not out_hist in in_file_hists:
    #         print out_hist

    ## Write and close output file
    print '\nWriting output file %s' % out_file.GetName()
    out_file.Write()
    out_file.Close()

    print '\nFinished writing all distributions.  Done!'

## End function: def MergeOneCat(category)


## Main function: loop over categories and merge histograms
def main():

    for category in CATEGORIES:
        for year in YEARS:
            MergeOneCat(category, year)

main()  ## End of main function

