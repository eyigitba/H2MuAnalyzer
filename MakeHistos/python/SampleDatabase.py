
import sys

#################################
##      SampleDatabase.py      ##
##  Return location of NTuple  ##
##  files from lxplus ('CERN') ##
##  or UF HPC/IHPEA ('UF') for ##
##  Leg2016, 2016, 2017, 2018  ##
#################################

class SampleInfo:

    def __init__(self, name, DAS_name, xsec, year, in_dir, evt_type):
        self.name     = name       ## Name of sample
        self.DAS_name = DAS_name   ## Name of sample from DAS, designates location of input NTuple files
        self.xsec     = xsec       ## Cross section in pb
        self.year     = year       ## Leg2016, 2016, 2017, or 2018
        self.in_dir   = in_dir     ## Input directory in lxplus or UF HPC
        self.evt_type = evt_type   ## 'Data', 'Sig', or 'Bkg'


def GetSamples(location = 'CERN', year = '2017'):

    print '\nGetting samples: location = %s, year = %s\n' % (location, year)

    ## Vector of sample infos to return
    samples = []
    ## Initialize in_dir to null
    in_dir  = ''

    #####################################################
    ###  Set location of samples in lxplus or UF HPC  ###
    #####################################################
    if   (location == 'UF'):
        if (year == 'Leg2016'):
            in_dir = '/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/simplified/'
    elif (location == 'CERN'):
        if   (year == 'Leg2016'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13'
        elif (year == '2016'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2016/94X_v3/prod-v16.0.7'
        elif (year == '2017'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17'
        elif (year == '2018'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/prod-v18-pre-tag'
        else:
            print 'Invalid location (%s) and/or year (%s)!!!  Exiting.' % (location, year)
            sys.exit
    elif (location == 'CERN_3l'):
        if   (year == '2016'):
            in_dir = '/eos/cms/store/user/bortigno/h2mm/ntuples/2016/94X_v3/STR'
        elif (year == '2017'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1'
        elif (year == '2018'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/prod-v18.1.6.skim3l'
        else:
            print 'Invalid location (%s) and/or year (%s)!!!  Exiting.' % (location, year)
            sys.exit
    elif (location == 'CERN_hiM'):
        if (year == 'Leg2016'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM'
    elif (location == 'CERN_lepMVA_test_v1'):
        if (year == '2017'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2018_12_13_LepMVA_2l_test_v1'
    elif (location == 'CERN_lepMVA_test_v2'):
        if (year == '2017'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_14_LepMVA_2l_test_v2'
    elif (location == 'CERN_lepMVA_hiM_test_v2'):
        if (year == '2017'):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_14_LepMVA_2l_hiM_test_v2'



    if ( in_dir == '' ):
        print 'Invalid location (%s) and/or year (%s)!!!  Exiting.' % (location, year)
        sys.exit


    ##########################
    ###  Get data samples  ###
    ##########################
    if (year == 'Leg2016' or year == '2016'):
        samples.append( SampleInfo('SingleMu_2016B',   'SingleMuon', 5800, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016C',   'SingleMuon', 2600, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016D',   'SingleMuon', 4300, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016E',   'SingleMuon', 4100, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016G',   'SingleMuon', 7800, year, in_dir, 'Data') )
    if (year == 'Leg2016'):
        samples.append( SampleInfo('SingleMu_2016F_1', 'SingleMuon', 1600, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016F_2', 'SingleMuon', 1600, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016H_1', 'SingleMuon', 4507, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016H_2', 'SingleMuon', 4507, year, in_dir, 'Data') )
    if (year == '2016'):
        samples.append( SampleInfo('SingleMu_2016F',   'SingleMuon', 3200, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016H',   'SingleMuon', 9014, year, in_dir, 'Data') )

    if (year == '2017'):
        samples.append( SampleInfo('SingleMu_2017B', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017C', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017D', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017E', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017F', 'SingleMuon', 1, year, in_dir, 'Data') )

    if (year == '2018'):
        samples.append( SampleInfo('SingleMu_2018A', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2018B', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2018C', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2018D', 'SingleMuon', 1, year, in_dir, 'Data') )

    ########################
    ###  Get MC samples  ###
    ########################

    ################
    ###  Signal  ###
    ################
    # if (year == '2017'): mH = ''

    if (year == '2017'): mH = '_125'
    else:                mH = ''

    if (year == 'Leg2016'): mVH = ''
    else:                   mVH = '_125'

    if (year == '2017'): ORD = '_NLO'
    else:                ORD = ''

    if   (year == 'Leg2016'): ggH_gen = '_13TeV_powheg_pythia8'
    elif (year == '2017'):    ggH_gen = '_13TeV_amcatnloFXFX_pythia8'
    else:                     ggH_gen = '_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8'

    if (year == 'Leg2016'): VBF_gen = '_13TeV_powheg_pythia8'
    else:                   VBF_gen = '_13TeV_amcatnlo_pythia8'

    if (year == 'Leg2016'): ggH_str = 'GluGlu_HToMuMu_'
    else:                   ggH_str = 'GluGluHToMuMu_'

    if (year == 'Leg2016'): VBF_str = 'VBF_'
    else:                   VBF_str = 'VBF'

    if (year == 'Leg2016'): WH_pos_str = 'WPlusH_HToMuMu_'
    else:                   WH_pos_str = 'WplusH_HToMuMu_WToAll_'

    if (year == 'Leg2016'): WH_neg_str = 'WMinusH_HToMuMu_'
    else:                   WH_neg_str = 'WminusH_HToMuMu_WToAll_'

    if (year == 'Leg2016'): ZH_str = 'ZH_HToMuMu_'
    else:                   ZH_str = 'ZH_HToMuMu_ZToAll_'

    if (year == '2016' or year == '2018'): VH_gen = '_TuneCP5_PSweights_13TeV_powheg_pythia8'
    else:                                  VH_gen = '_13TeV_powheg_pythia8'

    if (year == '2017'): ttH_gen = '_TuneCP5_13TeV-powheg-pythia8'
    else:                ttH_gen = '_TuneCP5_PSweights_13TeV-powheg-pythia8'


    ## H2Mu_gg
    if (year == '2017'):
      samples.append( SampleInfo('H2Mu_gg', ggH_str+'M125'+ggH_gen, 0.009618, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_gg'+mH +ORD, ggH_str+'M125'+ggH_gen, 0.009618, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_gg_120'+ORD, ggH_str+'M120'+ggH_gen, 0.009618, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_gg_130'+ORD, ggH_str+'M130'+ggH_gen, 0.009618, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018

    ## H2Mu_VBF
    if (year == '2017'):
        samples.append( SampleInfo('H2Mu_VBF_125_NLO_1', VBF_str+'HToMuMu_M125'+VBF_gen, 0.0008208, year, in_dir, 'Sig') )
        samples.append( SampleInfo('H2Mu_VBF_125_NLO_2', VBF_str+'HToMuMu_M125'+VBF_gen, 0.0008208, year, in_dir, 'Sig') )
        samples.append( SampleInfo('H2Mu_VBF_120_NLO_2', VBF_str+'HToMuMu_M120'+VBF_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
        samples.append( SampleInfo('H2Mu_VBF_130_NLO_2', VBF_str+'HToMuMu_M130'+VBF_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
        samples.append( SampleInfo('H2Mu_VBF',           VBF_str+'HToMuMu_M125'+VBF_gen, 0.0008208, year, in_dir, 'Sig') )
    else:
        samples.append( SampleInfo('H2Mu_VBF'+mH,        VBF_str+'HToMuMu_M125'+VBF_gen, 0.0008208, year, in_dir, 'Sig') )
        samples.append( SampleInfo('H2Mu_VBF_120',       VBF_str+'HToMuMu_M120'+VBF_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
        samples.append( SampleInfo('H2Mu_VBF_130',       VBF_str+'HToMuMu_M130'+VBF_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    ## H2Mu_VH
    samples.append( SampleInfo('H2Mu_ZH'+mVH, ZH_str+'M125'+VH_gen, 0.0001923, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_ZH_120', ZH_str+'M120'+VH_gen, 0.0001923, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_ZH_130', ZH_str+'M130'+VH_gen, 0.0001923, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    samples.append( SampleInfo('H2Mu_WH_pos'+mVH, WH_pos_str+'M125'+VH_gen, 0.0001858, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_WH_pos_120', WH_pos_str+'M120'+VH_gen, 0.0001858, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_WH_pos_130', WH_pos_str+'M130'+VH_gen, 0.0001858, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    samples.append( SampleInfo('H2Mu_WH_neg'+mVH, WH_neg_str+'M125'+VH_gen, 0.0001164, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_WH_neg_120', WH_neg_str+'M120'+VH_gen, 0.0001164, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_WH_neg_130', WH_neg_str+'M130'+VH_gen, 0.0001164, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    ## H2Mu_ttH
    if (year != 'Leg2016'):
        samples.append( SampleInfo('H2Mu_ttH_120', 'ttHToMuMu_M120'+ttH_gen, 0.0001103, year, in_dir, 'Sig') )  ## 0.5071 x 0.0002176 from YR4 (125 GeV)
        samples.append( SampleInfo('H2Mu_ttH_125', 'ttHToMuMu_M125'+ttH_gen, 0.0001103, year, in_dir, 'Sig') )  ## 0.5071 x 0.0002176 from YR4 (125 GeV)
        samples.append( SampleInfo('H2Mu_ttH_130', 'ttHToMuMu_M130'+ttH_gen, 0.0001103, year, in_dir, 'Sig') )  ## 0.5071 x 0.0002176 from YR4 (125 GeV)
        samples.append( SampleInfo('H2Mu_ttH', 'ttHToMuMu_M125'+ttH_gen, 0.0001103, year, in_dir, 'Sig') )  ## 0.5071 x 0.0002176 from YR4 (125 GeV)

    
    ####################
    ###  Background  ###
    ####################

    if ('2016' in year): py_tune = '_TuneCUETP8M1_'
    else:                py_tune = '_TuneCP5_'
 
    ## DYJetsToLL
    samples.append(     SampleInfo('ZJets',          'Combination_of_ZJets_AMC_ZJets_MG_1_ZJets_MG_2',                             5765.4,             year, in_dir, 'Bkg') )
    samples.append(     SampleInfo('ZJets_hiM',      'Combination_of_ZJets_hiM_AMC_ZJets_hiM_MG',                                    46.948,           year, in_dir, 'Bkg') )
    samples.append(     SampleInfo('ZJets_AMC',      'DYJetsToLL_M-50'+py_tune+'13TeV-amcatnloFXFX-pythia8',                       5765.4,             year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (6020.85) - AWB 28.09.2018
    samples.append(     SampleInfo('ZJets_hiM_MG',   'DYJetsToLL_M-105To160_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',            46.948,           year, in_dir, 'Bkg') ) ## xsec taken from 2017
    if (not '2016' in year):
        samples.append( SampleInfo('ZJets_MG_1',     'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',                          5765.4,             year, in_dir, 'Bkg') )
    if (year == '2016'):
        samples.append( SampleInfo('ZJets_MG_1',     'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',                     5765.4,             year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_2',     'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',                     5765.4,             year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_hiM_AMC',  'DYJetsToLL_M-105To160_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',                46.948,           year, in_dir, 'Bkg') ) ## xsec taken from 2017
    if (not '2016' in year):
        samples.append( SampleInfo('ZJets_hiM_AMC',  'DYJetsToLL_M-105To160_TuneCP5_PSweights_13TeV-amcatnloFXFX-pythia8',           46.948,           year, in_dir, 'Bkg') ) ## xsec taken from 2017
    if (year == '2017'):
        samples.append( SampleInfo('ZJets_MG_2',      'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',                         5765.4,             year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_m_10_50',   'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8',                    18610.0,             year, in_dir, 'Bkg') )

    if (year == 'Leg2016'):
        samples.append( SampleInfo('ZJets_AMC_0j_A', 'DYToLL_0J_13TeV-amcatnloFXFX-pythia8',                               4754   * 0.96,              year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_AMC_1j_A', 'DYToLL_1J_13TeV-amcatnloFXFX-pythia8',                                888.9 * 0.86*0.985*0.995,  year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_AMC_2j_A', 'DYToLL_2J_13TeV-amcatnloFXFX-pythia8',                                348.8 * 0.88*0.975*0.992,  year, in_dir, 'Bkg') )
        ## samples.append( SampleInfo('ZJets_AMC_0j_B', 'DYToLL_0J_13TeV-amcatnloFXFX-pythia8',                              4754   * 0.96,               year, in_dir, 'Bkg') ) ## Not used? - AWB 16.09.2018
        ## samples.append( SampleInfo('ZJets_AMC_1j_B', 'DYToLL_1J_13TeV-amcatnloFXFX-pythia8',                                888.9 * 0.86*0.985*0.995,  year, in_dir, 'Bkg') ) ## Not used? - AWB 16.09.2018
        ## samples.append( SampleInfo('ZJets_AMC_2j_B', 'DYToLL_2J_13TeV-amcatnloFXFX-pythia8',                                348.8 * 0.88*0.975*0.992,  year, in_dir, 'Bkg') ) ## Not used? - AWB 16.09.2018
        samples.append( SampleInfo('ZJets_MG',              'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',              5765.4,             year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_70_100',    'DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    178.952    * 0.98, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_100_200_A', 'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   181.302    * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_100_200_B', 'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   181.302    * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_200_400_A', 'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    50.4177   * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_200_400_B', 'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    50.4177   * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_400_600_B', 'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',     6.98394  * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_400_600_A', 'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',     6.98394  * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_600_800',   'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',     1.68141  * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_800_1200',  'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    0.775392 * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_1200_2500', 'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   0.186222 * 0.96, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_HT_2500_inf',  'DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    0.004385 * 0.96, year, in_dir, 'Bkg') )

    ## TTJets
    samples.append(     SampleInfo('tt_ll',     'Combination_of_tt_ll_MG_tt_ll_POW_tt_ll_AMC',             85.656,       year, in_dir, 'Bkg') )
    if (year != 'Leg2016'):
        samples.append( SampleInfo('tt_ll_MG',  'TTJets_DiLept'+py_tune+'13TeV-madgraphMLM-pythia8',       85.656,       year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018
    if (year == '2016' or year == '2017'):
        samples.append( SampleInfo('tt_ll_POW', 'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',        85.656,       year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018
    if (year == '2018'):
        samples.append( SampleInfo('tt_ll_AMC',  'TT_DiLept_TuneCP5_13TeV-amcatnlo-pythia8',               85.656,       year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tt_ll_POW',  'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',		   85.656,       year, in_dir, 'Bkg') )
    if (year == 'Leg2016'):
        samples.append( SampleInfo('tt_ll_AMC', 'TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8', 85.656 * 0.9, year, in_dir, 'Bkg') ) ## Why the factor of 0.9?!? - AWB 28.09.2018
        samples.append( SampleInfo('tt_ll_MG_1', 'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   85.656,       year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018
        samples.append( SampleInfo('tt_ll_MG_2', 'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   85.656,       year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018

    ## Single top (+X)
    if (year == '2016'):
        samples.append( SampleInfo('tW_pos',   'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',     66.02,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_neg',   'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8', 66.02,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tZq',      'tZq_ll_4f_13TeV-amcatnlo-pythia8',                                         0.09418, year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('tHq',      'THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1',                            0.07096, year, in_dir, 'Bkg') )  ## AN-2016/378
        samples.append( SampleInfo('tHW',      'THW_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1',                            0.01561, year, in_dir, 'Bkg') )  ## AN-2016/378
    if (year == '2017' or year == '2018'):
        samples.append( SampleInfo('tW_pos',   'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',         35.85,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_neg',   'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',     35.85,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tHq',      'THQ_4f_Hincl_13TeV_madgraph_pythia8',                                      0.07096, year, in_dir, 'Bkg') )  ## AN-2016/378
        samples.append( SampleInfo('tHW',      'THW_5f_Hincl_13TeV_madgraph_pythia8',                                      0.01561, year, in_dir, 'Bkg') )  ## AN-2016/378
    if (year == '2017'):
        samples.append( SampleInfo('tZq',      'tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8',               0.09418, year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('tZW',      'ST_tWll_5f_LO_TuneCP5_PSweights_13TeV-madgraph-pythia8',                   0.01123, year, in_dir, 'Bkg') )  ## AN-2018/025
    if (year == '2018'):
        samples.append( SampleInfo('tZq',      'tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-madgraph-pythia8',                         0.09418, year, in_dir, 'Bkg') )  ## AN-2018/025

    if (year == 'Leg2016'):
        samples.append( SampleInfo('tW_pos_1', 'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',            35.85,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_pos_2', 'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',            35.85,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_neg_1', 'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',        35.85,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_neg_2', 'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',        35.85,    year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tZq',      'tZq_ll_4f_13TeV-amcatnlo-pythia8',                                         0.09418, year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('tZW',      'ST_tWll_5f_LO_13TeV-MadGraph-pythia8',                                     0.01123, year, in_dir, 'Bkg') )  ## AN-2018/025

    ## TTX
    samples.append(     SampleInfo('ttZ',      'TTZToLLNuNu_M-10'+py_tune+'13TeV-amcatnlo-pythia8',                 0.2728,  year, in_dir, 'Bkg') )  ## AN-2018/025
    if (year == '2016'):
        samples.append( SampleInfo('ttZ_1',    'TTZToLLNuNu_M-10'+py_tune+'13TeV-amcatnlo-pythia8',                 0.2728,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttZ_2',    'TTZToLLNuNu_M-10'+py_tune+'13TeV-amcatnlo-pythia8',                 0.2728,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttWW',     'TTWW_TuneCUETP8M2T4_13TeV-madgraph-pythia8',                        0.00783, year, in_dir, 'Bkg') )  ## AN-2018/025 (arXiv 1405.0301 has 0.0099)
    if (not '2016' in year):
        samples.append( SampleInfo('ttWW',     'TTWW_TuneCP5_13TeV-madgraph-pythia8',                               0.00783, year, in_dir, 'Bkg') )  ## AN-2018/025 (arXiv 1405.0301 has 0.0099)
    if (year == '2016' or year == '2018'):
        samples.append( SampleInfo('ttW',      'TTWJetsToLNu'+py_tune+'13TeV-amcatnloFXFX-madspin-pythia8',         0.2043,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttH',      'ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8',                      0.2151,  year, in_dir, 'Bkg') )
    if (year == '2017' or year == '2018'):
        samples.append( SampleInfo('ttZ_lowM', 'TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8',                    0.0493,  year, in_dir, 'Bkg') )  ## AN-2018/025
    if (year == '2017'):
        samples.append( SampleInfo('ttW',      'TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8', 0.2043,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttH',      'ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8',     0.2151,  year, in_dir, 'Bkg') )
    if (year == 2018):
	samples.append( SampleInfo('ttW',      'TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8', 0.2043,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttWW',     'TTWW_TuneCP5_13TeV-madgraph-pythia8',                               0.00783, year, in_dir, 'Bkg') )  ## AN-2018/025 (arXiv 1405.0301 has 0.0099)
        samples.append( SampleInfo('ttZ',      'TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8',                   0.2728,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttZ_lowM', 'TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8',                    0.0493,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttH',      'ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8',     		    0.2151,  year, in_dir, 'Bkg') )

    if (year == 'Leg2016'):
        samples.append( SampleInfo('ttW_1',    'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',      0.2043,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttW_2',    'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',      0.2043,  year, in_dir, 'Bkg') )  ## AN-2018/025
        samples.append( SampleInfo('ttH',      'ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8',     0.2151,  year, in_dir, 'Bkg') )

    ## Diboson
    samples.append(     SampleInfo('WZ_2l',          'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',       4.409,   year, in_dir, 'Bkg') ) ## Consistent with current WZ_3l, see below - AWB 28.09.2018
    samples.append(     SampleInfo('ZZ_2l_2q',       'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',       3.22,    year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
    if (year != 'Leg2016'):
        samples.append( SampleInfo('WZ_3l',          'WZTo3LNu'+py_tune+'13TeV-amcatnloFXFX-pythia8',     4.666,   year, in_dir, 'Bkg') ) ## We used 2.113 in 2016 - AWB 09.10.2018    4.430(2016) from TOP-18-008, 4.666(2016) from AN-2017/277, 5.063(2017) from TOP-18-008. use 4.666 to align with Hamburg - XZW 13.05.2019
    if (year == '2016' or year == '2018'):
        samples.append( SampleInfo('ggZZ_4mu',       'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8',     0.00159, year, in_dir, 'Bkg') ) ## From AN-2018/340
        samples.append( SampleInfo('ggZZ_4tau',      'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8',    0.00159, year, in_dir, 'Bkg') ) ## From AN-2018/340
        samples.append( SampleInfo('ggZZ_2e2mu',     'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8',   0.00319, year, in_dir, 'Bkg') ) ## From AN-2018/340
        samples.append( SampleInfo('ggZZ_2e2tau',    'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8',  0.00319, year, in_dir, 'Bkg') ) ## From AN-2018/340
        samples.append( SampleInfo('ggZZ_2mu2tau',   'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8', 0.00319, year, in_dir, 'Bkg') ) ## From AN-2018/340
    if (year == '2016'):
        samples.append( SampleInfo('WW_2l_1',        'WWTo2L2Nu_13TeV-powheg',                           12.46,    year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
        samples.append( SampleInfo('ZZ_4l',          'ZZTo4L_13TeV_powheg_pythia8_ext1',                  1.212,   year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
        samples.append( SampleInfo('ZZ_4l_amc',      'ZZTo4L_13TeV-amcatnloFXFX-pythia8',                 1.212,   year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
    if (year == '2017'):
        samples.append( SampleInfo('WW',             'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8',   12.46,    year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
        samples.append( SampleInfo('ZZ_4l',          'ZZTo4L_13TeV_powheg_pythia8',                       1.212,   year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
        samples.append( SampleInfo('ZZ_4l_gg_4mu',   'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8',     0.00159, year, in_dir, 'Bkg') ) ## From AN-2018/340
        samples.append( SampleInfo('ZZ_4l_gg_2e2mu', 'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8',   0.00319, year, in_dir, 'Bkg') ) ## From AN-2018/340
    if (year == '2018'):
        samples.append( SampleInfo('WW_2l_1',        'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8',   12.46,    year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
        samples.append( SampleInfo('ZZ_4l',          'ZZTo4L_TuneCP5_13TeV_powheg_pythia8',               1.212,   year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns

    if (year == 'Leg2016'):
        samples.append( SampleInfo('WW',             'WWTo2L2Nu_13TeV-powheg',                           12.46,    year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
        samples.append( SampleInfo('WZ_3l_AMC',      'WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8',  4.666,   year, in_dir, 'Bkg') ) ## We used 2.113 in 2016 - AWB 09.10.2018    4.430(2016) from TOP-18-008, 4.666(2016) from AN-2017/277, 5.063(2017) from TOP-18-008. use 4.666 to align with Hamburg - XZW 13.05.2019
        samples.append( SampleInfo('ZZ_2l_2v',       'ZZTo2L2Nu_13TeV_powheg_pythia8',                    0.564,   year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns
        samples.append( SampleInfo('ZZ_4l_AMC',      'ZZTo4L_13TeV-amcatnloFXFX-pythia8',                 1.212,   year, in_dir, 'Bkg') ) ## https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns


    ## Triboson
    if (year != '2017'): WWZ_str = 'WWZ'
    if (year == '2017'): WWZ_str = 'WWZ_4F'

    samples.append( SampleInfo('WWW', 'WWW_4F'+py_tune+'13TeV-amcatnlo-pythia8', 0.2086,  year, in_dir, 'Bkg') ) ## From AN-2017/277
    samples.append( SampleInfo('WWZ', WWZ_str+py_tune+'13TeV-amcatnlo-pythia8',  0.1651,  year, in_dir, 'Bkg') ) ## From AN-2017/277
    samples.append( SampleInfo('WZZ', 'WZZ'+py_tune+'13TeV-amcatnlo-pythia8',    0.05565, year, in_dir, 'Bkg') ) ## From AN-2017/277
    samples.append( SampleInfo('ZZZ', 'ZZZ'+py_tune+'13TeV-amcatnlo-pythia8',    0.01398, year, in_dir, 'Bkg') ) ## From AN-2017/277
    if (year == '2016'):
        samples.append( SampleInfo('WWW_lep', 'WWW_4F_DiLeptonFilter_TuneCUETP8M1_13TeV-amcatnlo-pythia8', 0.0515,  year, in_dir, 'Bkg') ) ## Scale WWW cross section by 0.324^3 + 3*0.676*0.324^2 = 0.247

    ## VH with Higgs to WW, tau-tau, and ZZ
    if (year == '2017'):
        samples.append( SampleInfo('H2W_ZH_125',       'GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5',                     0.1888,  year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2W_WH_neg_125',   'HWminusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5',                     0.1138,  year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2W_WH_pos_125',   'HWplusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5',                      0.1796,  year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Tau_WH_neg_125', 'WminusHToTauTau_M125_13TeV_powheg_pythia8',                            0.03339, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Tau_WH_pos_125', 'WplusHToTauTau_M125_13TeV_powheg_pythia8',                             0.05268, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Tau_ZH_125',     'ZHToTauTau_M125_13TeV_powheg_pythia8',                                 0.05541, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Z_WH_neg_125',   'WminusH_HToZZTo2L2X_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8', 0.01395, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Z_WH_pos_125',   'WplusH_HToZZTo2L2X_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8',  0.02200, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Z_ZH_125',       'ZH_HToZZ_2LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8',   0.02314, year, in_dir, 'Bkg') )


    ####################################################
    ###  Shorter lists of samples for special tests  ###
    ####################################################

    # CERN_lepMVA_test_v1_list = [ 'SingleMu_2017B', 'SingleMu_2017C', 'SingleMu_2017D', 'SingleMu_2017E', 'SingleMu_2017F',
    #                              'H2Mu_gg', 'H2Mu_ZH_120', 'H2Mu_ZH_125', 'H2Mu_ZH_130', 'H2Mu_WH_pos_120',
    #                              'H2Mu_WH_pos_125', 'H2Mu_WH_pos_130' 'H2Mu_WH_neg_125', 'H2Mu_WH_neg_130', 'H2Mu_ttH_125',
    #                              'ZJets_AMC', 'tt_ll_MG', 'tt_ll_POW', 'tt_ljj_POW_1', 'tt_ljj_POW_2', 'tW_neg', 'tW_pos',
    #                              'WZ_3l', 'ZZ_4l' ]
    
    # CERN_lepMVA_test_v2_list = [ 'SingleMu_2017B', 'SingleMu_2017C', 'SingleMu_2017D', 'SingleMu_2017E', 'SingleMu_2017F',
    #                              'ZJets_AMC', 'tt_ll_MG', 'tt_ll_POW', 'tt_ljj_POW_1', 'tt_ljj_POW_2', 'tW_neg', 'tW_pos' ]

    # new_samples = []
    # if (location == 'CERN_lepMVA_test_v1' and year == '2017'):
    #     for i in range( len(CERN_lepMVA_test_v1_list) ):
    #         if samples.at(i).name in CERN_lepMVA_test_v1_list:
    #             new_samples.append( samples.at(i) )
    #     samples = new_samples
    # if (location == 'CERN_lepMVA_test_v2' and year == '2017'):
    #     for i in range( len(CERN_lepMVA_test_v2_list) ):
    #         if samples.at(i).name in CERN_lepMVA_test_v2_list:
    #             new_samples.append( samples.at(i) )
    #     samples = new_samples


    # ##################################################################################
    # ###  Double-Higgs signal and background samples: not needed for H2Mu analysis  ###
    # ##################################################################################

    # ## '/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1'

    # ## Double Higgs - different "nodes" correspond to non-SM HH couplings
    # if (year == '2017'):
    #     samples.append( SampleInfo('HH_n12_incl_2', 'GluGluToHHTo2B2Tau_node_12_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
    #     samples.append( SampleInfo('HH_n2_incl_2',  'GluGluToHHTo2B2Tau_node_2_13TeV-madgraph',  0.00243, year, in_dir, 'Bkg') )
    #     samples.append( SampleInfo('HH_n7_incl_2',  'GluGluToHHTo2B2Tau_node_7_13TeV-madgraph',  0.00243, year, in_dir, 'Bkg') )
    #     samples.append( SampleInfo('HH_SM_incl_2',  'GluGluToHHTo2B2Tau_node_SM_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
    #     samples.append( SampleInfo('HH_n12_incl_1', 'GluGluToHHTo4B_node_12_13TeV-madgraph',     0.00243, year, in_dir, 'Bkg') )
    #     samples.append( SampleInfo('HH_n2_incl_1',  'GluGluToHHTo4B_node_2_13TeV-madgraph',      0.00243, year, in_dir, 'Bkg') )
    #     samples.append( SampleInfo('HH_n7_incl_1',  'GluGluToHHTo4B_node_7_13TeV-madgraph',      0.00243, year, in_dir, 'Bkg') )
    #     samples.append( SampleInfo('HH_SM_incl_1',  'GluGluToHHTo4B_node_SM_13TeV-madgraph',     0.00243, year, in_dir, 'Bkg') )


    ## /////// ***** END OF ALL SAMPLE INITIALIZATION ***** ///////

    return samples

## End function: GetSamples(location = 'CERN', year = '2017')
