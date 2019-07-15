
import sys

## Python includes
from collections import OrderedDict

## ROOT includes
import ROOT as R


class StackPlotConfig:

    def __init__(self):
        self.year       = 'NONE',  ## '2016' or '2017'
        self.ntuple_loc = 'NONE',  ## 'CERN', 'CERN_hiM', 'CERN_3l', or 'UF'
        self.sig_mass   = 'NONE',  ## '120', '125', or '130'
        self.groups     = {}       ## 'Data', 'Sig', and 'Bkg' defined below
        self.weight     = {}       ##  By individual sample, defined below
        self.colors     = {}       ##  Colors for stacked data in each group

## End class StackPlotConfig


## Configure stack plot based on pre-defined settings
def ConfigStackPlot(known_config, year):

    print '\nConfiguring stack plots for %s with year %s' % (known_config, year)

    cfg = StackPlotConfig()

    cfg.groups = {}
    cfg.groups['Data'] = {}
    cfg.groups['Sig']  = {}
    cfg.groups['Bkg']  = {}

    ##################################################################################
    ###  Configuration for validation of WH backgrounds (Z+jets, ttbar, ttW, ttZ)  ###
    ##################################################################################

    if (known_config == 'WH_3l_bkg_val'):
        
        cfg.year       = year
        cfg.ntuple_loc = 'CERN_hiM'
        cfg.sig_mass   = '125'

        if (year == 'Leg2016'):

            cfg.groups['Data']['Data']  = []  ## By default, all data samples go into 'Data'
            cfg.groups['Sig']['Signal'] = []  ## By default, all signal samples go into 'Signal'

            cfg.groups['Bkg']['ZJets'] = ['ZJets_AMC_0j_A', 'ZJets_AMC_1j_A','ZJets_AMC_2j_A', 'ZJets_MG',
                                          'ZJets_MG_HT_70_100','ZJets_MG_HT_100_200_A', 'ZJets_MG_HT_100_200_B',
                                          'ZJets_MG_HT_200_400_A', 'ZJets_MG_HT_200_400_B', 'ZJets_MG_HT_400_600_A',
                                          'ZJets_MG_HT_400_600_B', 'ZJets_MG_HT_600_800', 'ZJets_MG_HT_800_1200',
                                          'ZJets_MG_HT_1200_2500', 'ZJets_MG_HT_2500_inf']
            ## After GEN_wgt (+/-1) factor, AMC_0j sample has ~60% of the total stats,
            ##   vs. ~40% for ZJets_MG.  Weight by 0.6 and 0.4, respectively.
            for samp in cfg.groups['Bkg']['ZJets']:
                if 'ZJets_AMC' in samp: cfg.weight[samp] = 0.6
                if 'ZJets_MG'  in samp: cfg.weight[samp] = 0.4

            cfg.groups['Bkg']['WZ']         = ['WZ_3l_AMC']
            cfg.groups['Bkg']['ttbar']      = ['tt_ll_MG_1', 'tt_ll_MG_2']
            cfg.groups['Bkg']['single top'] = ['tW_pos_1', 'tW_pos_2', 'tW_neg_1', 'tW_neg_2']
            cfg.groups['Bkg']['ttX']        = ['ttW_1', 'ttW_2', 'ttZ', 'ttH', 'tZq', 'tZW']
            cfg.groups['Bkg']['Other']      = []  ## Any background sample not of a specified type goes into 'Other'
            
            ## Preserve ordering of groups
            cfg.groups['Data'] = OrderedDict(cfg.groups['Data'])
            cfg.groups['Sig']  = OrderedDict(cfg.groups['Sig'])
            cfg.groups['Bkg']  = OrderedDict(cfg.groups['Bkg'])

            ## Samples to exclude from consideration
            cfg.excl_samps = ['tt_ll_AMC',  ## Has ~80% of the tt_ll_MG stats, but ~1/3 are negative weight, so effectively ~1/4 the stats in the end
                              'ZJets_AMC',  ## Inclusive sample unnecessary with ZJets_AMC_Nj samples (N = 0/1/2)
                              'ZJets_hiM']  ## Doesn't cover 3-muon phase space, where Higgs candidate pair may include a non-prompt muon not from Z
            
        else: print 'Year %s not valid for config %s. Exiting.' % year, sys.exit()

        ## Colors for each group
        cfg.colors = {}
        cfg.colors['Signal']     = R.kRed
        cfg.colors['Data']       = R.kBlack
        cfg.colors['ZJets']      = R.kAzure
        cfg.colors['WZ']         = R.kViolet
        cfg.colors['ttbar']      = R.kOrange+2
        cfg.colors['single top'] = R.kOrange
        cfg.colors['ttX']        = R.kSpring
        cfg.colors['Other']      = R.kCyan

    ## End conditional: if (known_config == 'WH_3l_bkg_val')


    ##################################################################################
    ###  Configuration for validation of WH backgrounds (Z+jets, ttbar, ttW, ttZ)  ###
    ##################################################################################

    if (known_config == 'WH_lep' or known_config == 'WH_lep_allMass' or
        known_config == 'ttH_3l' or known_config == 'ttH_3l_allMass'):

        cfg.year       = year
        cfg.ntuple_loc = 'CERN_3l'
        cfg.sig_mass   = '125'

        if (year == '2016' or year == '2017' or year == '2018' or year == 'Run2'):

            cfg.groups['Data']['Data']  = []  ## By default, all data samples go into 'Data'
            cfg.groups['Sig']['Signal'] = []  ## By default, all signal samples go into 'Signal'
            cfg.excl_samps              = []  ## Samples to exclude from consideration

            ## Start by excluding component histograms merged in MakeHistos/scripts/MergeHists.py
            cfg.excl_samps += ['H2Mu_VBF_120_NLO_1', 'H2Mu_VBF_120_NLO_2']  ## Merged into H2Mu_VBF_120
            cfg.excl_samps += ['H2Mu_VBF_125_NLO_1', 'H2Mu_VBF_125_NLO_2']  ## Merged into H2Mu_VBF
            cfg.excl_samps += ['H2Mu_VBF_130_NLO_1', 'H2Mu_VBF_130_NLO_2']  ## Merged into H2Mu_VBF_130
            cfg.excl_samps += ['ZJets_AMC', 'ZJets_MG_1', 'ZJets_MG_2']     ## Merged into ZJets
            cfg.excl_samps += ['ZJets_hiM_AMC', 'ZJets_hiM_MG']             ## Merged into ZJets_hiM
            cfg.excl_samps += ['tt_ll_POW', 'tt_ll_MG', 'tt_ll_AMC']        ## Merged into tt_ll
            cfg.excl_samps += ['ttZ_1', 'ttZ_2']                            ## Merged into ttZ

            ## Additional samples to exclude from consideration (redundant)
            cfg.excl_samps += ['WWW_lep', 'ZZ_4l_amc']


            cfg.groups['Bkg']['WZ'] = ['WZ_3l']

            cfg.groups['Bkg']['ZZ'] = ['ZZ_4l']
            if (year == '2017' or year == 'Run2'):
                cfg.groups['Bkg']['ZZ'] += ['ZZ_4l_gg_4mu', 'ZZ_4l_gg_2e2mu']
            if (year == '2018' or year == 'Run2'):
                cfg.groups['Bkg']['ZZ'] += ['ggZZ_4mu', 'ggZZ_4tau', 'ggZZ_2e2mu', 'ggZZ_2mu2tau']

            cfg.groups['Bkg']['ttbar'] = ['tt_ll', 'tW_pos', 'tW_neg']

            if ('allMass' in known_config):
                ## Inclusive (mass > 50 GeV) Z+jets samples, to use for Z mass background validation plots
                cfg.groups['Bkg']['ZJets'] = ['ZJets']
                cfg.excl_samps            += ['ZJets_hiM']  ## Using inclusive samples, don't need 'hiM'

            elif not ('allMass' in known_config):
                ## High-mass (105 - 160 GeV) Z+jets samples, to use for Higgs mass signal region plots
                cfg.groups['Bkg']['ZJets']  = ['ZJets_hiM']
                cfg.excl_samps             += ['ZJets']  ## Using 'hiM' samples, don't need inclusive

            else: print 'Config %s is not valid for year %s. Exiting.' % (known_config, year), sys.exit()

            cfg.groups['Bkg']['Other'] = []  ## Any background sample not of a specified type goes into 'Other'


            ## Special modifications for WH or ttH plots
            if known_config == 'WH_lep' or known_config == 'WH_lep_allMass':

                cfg.groups['Bkg']['ttX'] = ['ttZ', 'ttW', 'ttH', 'tZq', 'tHq', 'tHW', 'ttZ_lowM', 'ttWW']

            elif known_config == 'ttH_3l' or known_config == 'ttH_3l_allMass':
                cfg.groups['Bkg']['top+X'] = ['tZq', 'tHq', 'tHW', 'ttWW'] ## Missing tZW
                cfg.groups['Bkg']['ttH']   = ['ttH']
                cfg.groups['Bkg']['ttW']   = ['ttW']
                cfg.groups['Bkg']['ttZ']   = ['ttZ', 'ttZ_lowM']

            else: print 'Config %s is not valid for year %s. Exiting.' % (known_config, year), sys.exit()

            ## Preserve ordering of groups
            cfg.groups['Data'] = OrderedDict(cfg.groups['Data'])
            cfg.groups['Sig']  = OrderedDict(cfg.groups['Sig'])
            cfg.groups['Bkg']  = OrderedDict(cfg.groups['Bkg'])

        else: print 'Year %s not valid for config %s. Exiting.' % (year, known_config), sys.exit()


    if (year == '2016' or year == '2017' or year == '2018' or year == 'Run2'):

        ## Colors for each group
        cfg.colors = {}
        cfg.colors['Signal']     = R.kRed
        cfg.colors['Data']       = R.kBlack
        cfg.colors['WZ']         = R.kViolet
        cfg.colors['ZZ']         = R.kPink+7
        cfg.colors['ZJets']      = R.kAzure
        cfg.colors['ttbar']      = R.kOrange+2
        cfg.colors['Other']      = R.kCyan

        if known_config == 'WH_lep' or known_config == 'WH_lep_allMass':
            cfg.colors['ttX']        = R.kSpring
        else:
            cfg.colors['top+X']      = R.kOrange
            cfg.colors['ttW']        = R.kSpring
            cfg.colors['ttZ']        = R.kGreen+3
            cfg.colors['ttH']        = R.kYellow-7


    ## End conditional: if (known_config == 'ttH_3l')


    ##############################
    ###  Other configurations  ###
    ##############################

    # elif (known_config == 'other'):
    #     cfg.properties = 'other'

    else: print 'Inside StackPlotConfig.py, unknown configuration %s. Exiting.' % known_config, sys.exit()


    return cfg

## End function: ConfigStackPlot(known_config)
