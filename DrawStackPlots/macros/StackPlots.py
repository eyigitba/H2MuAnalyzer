#! /usr/bin/env python

####################################################
###    Make stack and ratio from added histos    ###
####################################################

## N.B. For some reason frequently segfaults after plotting ~20 stacks - AWB 21.01.2019
## Can run as ./macros/StackPlots.py i j to plots stacks with indices [i, j]
## Script prints index of each stack so you can see the last one which succeded

## Basic python includes for manipulating files
import os
import sys

## ROOT includes
import ROOT as R
R.gROOT.SetBatch(True)  ## Don't display canvases to screen while running
R.gStyle.SetOptStat(0)

## Info about data and MC NTuples from H2MuAnalyzer/MakeHistos/python/SampleDatabase.py
sys.path.insert(0, '%s/../MakeHistos/python' % os.getcwd())
from SampleDatabase import GetSamples
## Configuration for stack plots from H2MuAnalyzer/DrawStackPlots/python/StackPlotConfig.py
sys.path.insert(0, '%s/python' % os.getcwd())
from StackPlotConfig import ConfigStackPlot

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Settings for this stack-drawing job
if USER == 'abrinke1':
    PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms'
    CONFIG   = 'ttH_3l'   ## Pre-defined stack configuration from python/StackPlotConfig.py
    YEAR     = '2017'     ## Dataset year (2016 or 2017)
    # LABEL    = 'WH_lep_AWB_2019_01_19_lepMVA_test_v1'  ## Sub-folder within PLOT_DIR containing histograms
    LABEL    = 'ttH_3l_AWB_2019_01_19_lepMVA_test_v1'  ## Sub-folder within PLOT_DIR containing histograms
    # CATEGORY = '3mu_0b_mt150_mass12_noZ'  ## Category for which to draw plots
    CATEGORY = 'e2mu_noZ_btag_ge3j_mass12'  ## Category for which to draw plots
    # LABEL    = 'WH_lep_CERN_hiM_11_10_2018_v1'  ## Sub-folder within PLOT_DIR containing histograms
    # CATEGORY = 'e2mu_NONE'  ## Category for which to draw plots
    # CATEGORY = '3mu_tight_0b_mt150_mass12_noZ'  ## Category for which to draw plots
    IN_FILE  = 'histos_Presel2017_%s.root' % CATEGORY  ## File with input histograms
    SCALE     = 'lin' ## 'log' or 'lin' scaling of y-axis
    RATIO_MIN = 0.0   ## Minimum value in ratio plot
    RATIO_MAX = 2.0   ## Maximum value in ratio plot
elif USER == 'xzuo':
    PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'
elif USER == 'bortigno':
    PLOT_DIR = 'NONE'
else: print 'Invalid USER = %s' % USER



## Function to draw each stack plot and ratio plot on the same canvas
def DrawOneStack( dist, sig_stack, all_stack, h_data, legend, out_file_name ):   # Do not use TRatioPlot! It is a devil! - XWZ 19.09.2018

    ## Create a new TCanvas
    canv = R.TCanvas('can_'+dist, 'can_'+dist, 1)
    canv.Clear()
    canv.cd()
    
    ## Draw the upper pad, with the stack histogram
    upper_pad = R.TPad('upperPad_'+dist, 'upperPad_'+dist, 0, 0.3, 1, 1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.SetRightMargin(0.20)
    upper_pad.Draw()
    upper_pad.cd()

    ## Draw signal + background stack
    if SCALE == 'log':
        upper_pad.SetLogy()
        all_stack.SetMinimum( 0.01 )
        if h_data: all_stack.SetMaximum( 10.0*max(h_data.GetMaximum(), all_stack.GetMaximum()) )
        else:      all_stack.SetMaximum( 10.0*all_stack.GetMaximum() )
    elif SCALE == 'lin':
        all_stack.SetMinimum(0)
        if h_data: all_stack.SetMaximum( 1.2*max(h_data.GetMaximum(), all_stack.GetMaximum()) )
        else:      all_stack.SetMaximum( 1.2*all_stack.GetMaximum() )
    else: print 'Invalid SCALE = %s. Exiting.' % SCALE, sys.exit()
        
    all_stack.Draw('HIST')

    ## Draw the data histogram
    if h_data:
        ## Blind 120 - 130 GeV range of mass plots
        if 'mass' in h_data.GetName():
            for i in range(1, h_data.GetNbinsX()+1):
                if h_data.GetXaxis().GetBinCenter(i) > 120 and h_data.GetXaxis().GetBinCenter(i) < 130:
                    h_data.SetBinContent(i, 0)
        h_data.SetMarkerStyle(20)
        h_data.Draw('SAME')

    ## Draw the signal histogram, normalized to total MC area
    h_sig = sig_stack.GetStack().Last().Clone('tmp')
    if h_sig.Integral() > 0: h_sig.Scale( all_stack.GetStack().Last().Integral() / h_sig.Integral() )
    h_sig.Draw('HISTSAME')

    ## Draw the legend
    legend.Draw()

    ## Draw the lower pad, with the ratio histogram
    canv.cd()

    lower_pad = R.TPad('lowerPad_'+dist, 'lowerPad_'+dist, 0, 0.05, 1, 0.3)
    lower_pad.SetTopMargin(0.05)
    lower_pad.SetRightMargin(0.20)
    lower_pad.Draw()
    lower_pad.cd()

    ## Create the ratio histogram of data / MC
    if h_data:
        ratio_hist = h_data.Clone('ratioHist')
        ratio_hist.Divide(all_stack.GetStack().Last())
        ratio_hist.SetTitle('')
        ratio_hist.SetMinimum(RATIO_MIN)
        ratio_hist.SetMaximum(RATIO_MAX)
        # ratio_hist.GetYaxis().SetNdivisions(502)  ## For some reason causes segfault after a dozen or so plots - AWB 09.10.2018
        ratio_hist.SetMarkerStyle(20)
        ratio_hist.Draw()

    canv.Update()
    canv.SaveAs(PLOT_DIR+'/'+LABEL+'/plots/'+CATEGORY+'/'+dist+'.png')

    ## Open output root file and save canvas
    out_file_loc = R.TFile.Open(out_file_name, 'UPDATE')
    out_file_loc.cd()
    canv.Write()
    out_file_loc.Close()

    ## Delete objects created in DrawOneStack()
    del canv, upper_pad, lower_pad, h_sig, out_file_loc
    if h_data: del ratio_hist

## End function: DrawOneStack()



def main():

    print '\nInside StackPlots.py, main()'

    #######################################
    ###  Set up input and output files  ###
    #######################################

    if USER == 'abrinke1': in_file_dir = PLOT_DIR+'/'+LABEL+'/files/HADD'
    else:                  in_file_dir = PLOT_DIR+'/'+LABEL+'/files/sum'

    in_file_name  = in_file_dir+'/'+IN_FILE
    out_dir       = PLOT_DIR+'/'+LABEL+'/plots/'+CATEGORY
    out_file_name = out_dir+'/StackPlots.root'
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    
    ## If this is the first time running StackPlots.py, create new output file
    if ( len(sys.argv) == 1 or int(sys.argv[1]) <= 1 ):
        print 'Creating output file %s' % out_file_name
        out_file = R.TFile.Open(out_file_name, 'RECREATE')
    else:  ## If StackPlots.py crashed and you are finishing the plots, update existing file
        print 'Re-opening output file %s' % out_file_name
        out_file = R.TFile.Open(out_file_name, 'UPDATE')
    out_file.Close()  ## Re-open later when saving stack
    del out_file
    print 'Opening input file %s' % in_file_name
    in_file = R.TFile.Open(in_file_name, 'READ')


    #########################################
    ###  Configure stack plot properties  ###
    #########################################

    cfg = ConfigStackPlot(CONFIG, YEAR)

    mass       = cfg.sig_mass
    loc        = cfg.ntuple_loc
    groups     = cfg.groups
    excl_samps = cfg.excl_samps
    colors     = cfg.colors


    ###################################################
    ###  Set up groups of samples to plot together  ###
    ###################################################

    ## List of expected samples from SampleDatabase.py
    DB_samps = GetSamples(loc, int(YEAR))

    ## Exclude samples with incorrect signal mass
    for samp in [S for S in DB_samps if S.evt_type == 'Sig']:
        if    mass == '125' and ('120' in samp.name or '130' in samp.name): excl_samps.append(samp.name)
        elif (mass == '120' or mass == '130') and not mass in samp.name:    excl_samps.append(samp.name)
        elif mass != '120' and mass != '125' and mass != '130': print 'Invalid signal mass choice %s - quitting' % mass, sys.exit()


    ##################################################################
    ###  Find all the samples and distributions in the input file  ###
    ##################################################################

    dists       = []  ## All distributions contained in the file
    found_samps = []  ## All samples found in the file
    ## Loop over all the histograms in the file
    ## Histogram name composed of sample+'_'+distribution
    for key in R.gDirectory.GetListOfKeys():
        ## See if histogram name starts with a known sample name
        for samp in [S for S in DB_samps if key.GetName().startswith(S.name+'_')]:
            ## Make sure sample name is not a sub-string of another sample
            unique_samp = True
            for samp2 in [S for S in DB_samps if key.GetName().startswith(S.name+'_')]:
                if samp2.name != samp.name and samp.name in samp2.name: unique_samp = False
            if not unique_samp: continue

            dist = key.GetName().replace(samp.name+'_', '')  ## Remove sample name to get distribution name
            if not dist in dists: dists.append(dist)         ## Append the distribution to the list of distributions

    print '\nList of found distributions:'
    print ', '.join(dists)

    ## Clean out excluded samples from the Database list
    DB_samps = [S for S in DB_samps if not S.name in excl_samps]

    ## Check to see if there are any samples not known to SampleDatabase.py
    for key in R.gDirectory.GetListOfKeys():
        for dist in [D for D in dists if key.GetName().endswith('_'+D)]:
            ## Make sure distribution name is not a sub-string of another distribution
            unique_dist = True
            for dist2 in dists:
                if dist2 != dist and dist in dist2: unique_dist = False
            if not unique_dist: continue

            samp = key.GetName().replace('_'+dist, '')            ## Remove distribution name to get sample name
            if samp in excl_samps: continue                       ## Don't consider explicitly excluded samples
            if not samp in found_samps: found_samps.append(samp)  ## Append the sample to the list of found samples
            if not any(samp == S.name for S in DB_samps):
                print 'Found histogram %s with distribution %s and sample %s, not found in SampleDatabase.py or excl_samps:'
                print excl_samps
                print 'Check that %s is in SampleDatabase.py or excl_samps, and that GetSamples(%s, %s) is configured properly' % (samp, loc, YEAR)
                sys.exit()

    print '\nList of found samples:'
    print ', '.join(found_samps)

    ## Check to see if any samples in SampleDatabase.py do not appear in the found list
    not_found_samps = []
    for samp in [S for S in DB_samps if not S.name in found_samps]: not_found_samps.append(samp.name)
    if len(not_found_samps) > 0:
        print '\nSURPRISING SITUATION!!!  The following samples not found in the input file:'
        print ', '.join(not_found_samps)

    ## Clean out not found samples from the Database list
    DB_samps = [S for S in DB_samps if not S.name in not_found_samps]
                

    ##################################
    ###  Assign samples to groups  ###
    ##################################

    ## If not otherwise assigned, put sample in 'Data', 'Signal', or 'Other' category
    for samp in DB_samps:
        if not any(samp.name in samps for group, samps in groups[samp.evt_type].items()):
            if samp.evt_type  == 'Data': groups['Data']['Data']  .append(samp.name)
            elif samp.evt_type == 'Sig': groups['Sig'] ['Signal'].append(samp.name)
            elif samp.evt_type == 'Bkg': groups['Bkg'] ['Other'] .append(samp.name)

    ## Sanity check that found_samps and DB_samps contain the same list of samples
    for samp in found_samps:
        if not any(samp == S.name for S in DB_samps): print '\nSample %s in found_samps not in SampleDatabase.py. Exiting', sys.exit()
    for samp in DB_samps:
        if not samp.name in found_samps: print '\nSample %s in SampleDatabase.py not in found_samps. Exiting', sys.exit()

    ## Check to see if we're missing any samples we were expecting
    for evt_type in groups.keys():
        for group, samps in groups[evt_type].items():
            for samp in [S for S in samps if not S in found_samps and not S in not_found_samps]:
                print '\nSURPRISING SITUATION!!!  Sample %s from %s group %s not found in final list of samples:' % (samp, evt_type, group)
                print ', '.join(found_samps)
            ## Drop samples from group if they do not appear in found_samps
            groups[evt_type][group] = [S for S in samps if S in found_samps]

    print '\nFinal list of data samples:'
    for group in groups['Data'].keys(): print group+': '+', '.join(groups['Data'][group])
    print '\nFinal list of signal samples:'
    for group in groups['Sig'].keys(): print group+': '+', '.join(groups['Sig'][group])
    print '\nFinal list of background samples:'
    for group in groups['Bkg'].keys(): print group+': '+', '.join(groups['Bkg'][group])


    ##########################################
    ###  Histograms, stacks, and canvases  ###
    ##########################################

    print '\nAbout to start creating and filling the plots'
    print 'Looping over samples and distributions to get histograms'
    iDist = 0
    for dist in dists:
        iDist += 1
        ## Only plot stacks in certain index range, if speficied by the user
        if ( len(sys.argv) > 1 and (iDist < int(sys.argv[1]) or iDist > int(sys.argv[2])) ): continue
        print '  * Looking at distribution %s (#%d)' % (dist, iDist)

	group_hist = {}  ## Summed sample histograms by group
	stack_all  = R.THStack('all_stack_'+dist, dist+' signal + background')
	stack_sig  = R.THStack('sig_stack_'+dist, dist+' signal')
	stack_dat  = R.THStack('dat_stack_'+dist, dist+' data')

        ## Construct histogram for each group from samples in group
        for evt_type in groups.keys():
            for group, samps in groups[evt_type].items():
                for samp in samps:
                    try:
                        hist = in_file.Get(samp+'_'+dist).Clone('tmp')

                        ## Include multiple samples scaled by available statistics
                        if YEAR == '2016':
                            if 'ZJets_AMC' in samp: hist.Scale(0.6)
                            if 'ZJets_MG'  in samp: hist.Scale(0.4)
                            # if 'tt_ll_MG'  in samp: hist.Scale(0.4)  ## Experimental SF from 3LooseMu_ttbar_3l_val_mu category
                        if YEAR == '2017':
                            if 'ZJets_AMC' in samp: hist.Scale(0.5)
                            if 'ZJets_MG'  in samp: hist.Scale(0.25)  ## Using both MG_1 and MG_2
                            if 'tt_ll_POW' in samp: hist.Scale(0.7)
                            if 'tt_ll_MG'  in samp: hist.Scale(0.3)

                        if not group in group_hist.keys():
                            group_hist[group] = hist.Clone('hist_'+dist+'_'+group)
                        else:
                            group_hist[group].Add(hist)
                        del hist
                    except:
                        print '  - Could not find histogram '+samp+'_'+dist

        ## Fill the signal stack
	for group in groups['Sig'].keys():
            group_hist[group].SetLineColor(colors[group])
            group_hist[group].SetLineWidth(2)
	    stack_sig.Add(group_hist[group])
	group_hist['Sig'] = stack_sig.GetStack().Last()
        ## Fill the signal + background stack
        for group in (groups['Bkg'].keys() + groups['Sig'].keys()):
            group_hist[group].SetFillColor(colors[group])
            stack_all.Add(group_hist[group])
	group_hist['MC'] = stack_all.GetStack().Last()
        ## Fill the data stack
        MC_only = True  ## If we don't find any data histograms, fill 'data' as sum of MC
        for group in groups['Data'].keys():
            if group in group_hist.keys():
                stack_dat.Add(group_hist[group])
                MC_only = False
        if not MC_only:
            group_hist['Dat'] = stack_dat.GetStack().Last()
        else:
            print '\nFor distribution %s, could not find any data.  Filling with sum of MC.\n' % dist
            group_hist['Dat'] = 0


        ## Create TLegend
        legend = R.TLegend(0.82, 0.3, 0.98, 0.9)
        for evt_type in groups.keys():
            for group in groups[evt_type].keys():
                if MC_only and group == 'Data': continue
                legend.AddEntry(group_hist[group], group)


        ## Draw stack plot
	DrawOneStack( dist, stack_sig, stack_all, group_hist['Dat'], legend, out_file_name )
        ## Delete objects created in loop over dists
        del group_hist, stack_all, stack_sig, stack_dat, legend
        
    ## End loop: for dist in dists

    print '\nFinished looping over all distributions.  Done!'

main()  ## End of main function

