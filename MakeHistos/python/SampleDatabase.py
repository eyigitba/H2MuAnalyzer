

#################################
##      SampleDatabase.py      ##
##  Return location of NTuple  ##
##  files from lxplus ('CERN') ##
##  or UF HPC/IHPEA ('UF') for ##
##        2016 or 2017         ##
#################################

class SampleInfo:

    def __init__(self, name, DAS_name, xsec, year, in_dir, evt_type):
        self.name     = name       ## Name of sample
        self.DAS_name = DAS_name   ## Name of sample from DAS, designates location of input NTuple files
        self.xsec     = xsec       ## Cross section in pb
        self.year     = year       ## 2016 or 2017
        self.in_dir   = in_dir     ## Input directory in lxplus or UF HPC
        self.evt_type = evt_type   ## 'Data', 'Sig', or 'Bkg'


def GetSamples(location = 'CERN', year = '2017'):

    print '\nGetting samples: location = %s, year = %d\n' % (location, year)

    ## Vector of sample infos to return
    samples = []

    #####################################################
    ###  Set location of samples in lxplus or UF HPC  ###
    #####################################################
    if (location == 'UF'):
        if (year == 2016):
            in_dir = '/cms/data/store/user/t2/users/acarnes/h2mumu/awb_samples/simplified/'
        else:
            print 'Invalid location (%s) and/or year (%d)!!!  Exiting.' % (location, year)
            sys.exit
    elif (location == 'CERN'):
        if   (year == 2016):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13'
        elif (year == 2017):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17'
        elif (year == 2018):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2018/102X/'
        else:
            print 'Invalid location (%s) and/or year (%d)!!!  Exiting.' % (location, year)
            sys.exit
    elif (location == 'CERN_hiM'):
        if (year == 2016):
            in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13_hiM'
        else:
            print 'Invalid location (%s) and/or year (%d)!!!  Exiting.' % (location, year)
            sys.exit
    elif (location == 'CERN_lepMVA_test_v1'):
	if (year == 2017):
	    in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2018_12_13_LepMVA_2l_test_v1'
	else:
	    print 'Invalid location (%s) and/or year (%d) !!! Exiting.' % (location, year)
	    sys.exit
    elif (location == 'CERN_lepMVA_test_v2'):
	if (year == 2017):
	    in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_14_LepMVA_2l_test_v2'
	else:
	    print 'Invalid location (%s) and/or year (%d) !!! Exiting.' % (location, year)
    elif (location == 'CERN_lepMVA_3l_test_v1_with_extra_samples'): #Be cautious! For now only works if in_dir is overwritten in GenerateBatch.py. Will change later. 
	if (year == 2017):
	    in_dir = '/eos/cms/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17'
	else: 
	    print 'Invalid location (%s) and/or year (%d) !!! Exiting.' % (location, year)
    else:
        print 'Invalid location (%s) and/or year (%d)!!!  Exiting.' % (location, year)
        sys.exit


    ##########################
    ###  Get data samples  ###
    ##########################
    if (year == 2016):
        samples.append( SampleInfo('SingleMu_2016B',   'SingleMuon', 5800, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016C',   'SingleMuon', 2600, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016D',   'SingleMuon', 4300, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016E',   'SingleMuon', 4100, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016F_1', 'SingleMuon', 1600, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016F_2', 'SingleMuon', 1600, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016G',   'SingleMuon', 7800, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016H_1', 'SingleMuon', 4507, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2016H_2', 'SingleMuon', 4507, year, in_dir, 'Data') )

    if (year == 2017):
        samples.append( SampleInfo('SingleMu_2017B', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017C', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017D', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017E', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2017F', 'SingleMuon', 1, year, in_dir, 'Data') )

    if (year == 2018):
        samples.append( SampleInfo('SingleMu_2018A', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2018B', 'SingleMuon', 1, year, in_dir, 'Data') )
        samples.append( SampleInfo('SingleMu_2018C', 'SingleMuon', 1, year, in_dir, 'Data') )


    ########################
    ###  Get MC samples  ###
    ########################

    ################
    ###  Signal  ###
    ################

    if (year == 2016): mH = ''
    if (year == 2017 or year == 2018): mH = '_125'

    if (year == 2016): ORD = ''
    if (year == 2017 or year == 2018): ORD = '_NLO'

    if (year == 2016): sig_gen = '_powheg_pythia8'
    if (year == 2017 or year == 2018): sig_gen = '_amcatnloFXFX_pythia8'
    
    if (year == 2016): ggH_str = 'GluGlu_HToMuMu_'
    if (year == 2017 or year == 2018): ggH_str = 'GluGluHToMuMu_'

    if (year == 2016): VBF_str = 'VBF_'
    if (year == 2017 or year == 2018): VBF_str = 'VBFH_'

    if (year == 2016): WH_pos_str = 'WPlusH_HToMuMu_'
    if (year == 2017 or year == 2018): WH_pos_str = 'WplusH_HToMuMu_WToAll_'

    if (year == 2016): WH_neg_str = 'WMinusH_HToMuMu_'
    if (year == 2017 or year == 2018): WH_neg_str = 'WminusH_HToMuMu_WToAll_'

    if (year == 2016): ZH_str = 'ZH_HToMuMu_'
    if (year == 2017 or year == 2018): ZH_str = 'ZH_HToMuMu_ZToAll_'

    ## H2Mu_gg
    samples.append( SampleInfo('H2Mu_gg'+mH+ORD,  ggH_str+'M125_13TeV'+sig_gen, 0.009618, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_gg_120'+ORD, ggH_str+'M120_13TeV'+sig_gen, 0.009618, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_gg_130'+ORD, ggH_str+'M130_13TeV'+sig_gen, 0.009618, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    ## H2Mu_VBF
    if (year == 2016):
        samples.append( SampleInfo('H2Mu_VBF'+mH,  VBF_str+'HToMuMu_M125_13TeV'+sig_gen, 0.0008208, year, in_dir, 'Sig') )
        samples.append( SampleInfo('H2Mu_VBF_120', VBF_str+'HToMuMu_M120_13TeV'+sig_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
        samples.append( SampleInfo('H2Mu_VBF_130', VBF_str+'HToMuMu_M130_13TeV'+sig_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('H2Mu_VBF_125_NLO_1', 'VBFHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnlo_pythia8', 0.0008208, year, in_dir, 'Sig') )
        samples.append( SampleInfo('H2Mu_VBF_125_NLO_2', 'VBFHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnlo_pythia8', 0.0008208, year, in_dir, 'Sig') )
        samples.append( SampleInfo('H2Mu_VBF_120_NLO_2', 'VBFHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnlo_pythia8', 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
        samples.append( SampleInfo('H2Mu_VBF_130_NLO_2', 'VBFHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnlo_pythia8', 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    ## H2Mu_VH
    samples.append( SampleInfo('H2Mu_ZH'+mH,  ZH_str+'M125_13TeV_powheg_pythia8', 0.0001923, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_ZH_120', ZH_str+'M120_13TeV_powheg_pythia8', 0.0001923, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_ZH_130', ZH_str+'M130_13TeV_powheg_pythia8', 0.0001923, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    samples.append( SampleInfo('H2Mu_WH_pos'+mH,  WH_pos_str+'M125_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_WH_pos_120', WH_pos_str+'M120_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_WH_pos_130', WH_pos_str+'M130_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    samples.append( SampleInfo('H2Mu_WH_neg'+mH,  WH_neg_str+'M125_13TeV_powheg_pythia8', 0.0001164, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_WH_neg_120', WH_neg_str+'M120_13TeV_powheg_pythia8', 0.0001164, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    samples.append( SampleInfo('H2Mu_WH_neg_130', WH_neg_str+'M130_13TeV_powheg_pythia8', 0.0001164, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.09.2018
    
    ## H2Mu_ttH
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('H2Mu_ttH'+mH, 'ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8', 0.0001103, year, in_dir, 'Sig') )  ## 0.5071 x 0.0002176 from YR4 (125 GeV)

    
    ####################
    ###  Background  ###
    ####################

    if (year == 2016): py_tune = '_TuneCUETP8M1_'
    if (year == 2017 or year ==2018): py_tune = '_TuneCP5_'
    
    ## DYJetsToLL
    samples.append( SampleInfo('ZJets_AMC',      'DYJetsToLL_M-50'+py_tune+'13TeV-amcatnloFXFX-pythia8', 5765.4,        year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (6020.85) - AWB 28.09.2018
    ## if (year == 2017):
        ## samples.append( SampleInfo('ZJets_AMC_2','DYJetsToLL_M-50'+py_tune+'13TeV-amcatnloFXFX-pythia8', 5765.4,        year, in_dir, 'Bkg') )
    if (year == 2016):
        samples.append( SampleInfo('ZJets_AMC_0j_A', 'DYToLL_0J_13TeV-amcatnloFXFX-pythia8',                 4754   * 0.96, year, in_dir, 'Bkg') )
        ## samples.append( SampleInfo('ZJets_AMC_0j_B', 'DYToLL_0J_13TeV-amcatnloFXFX-pythia8',                4754   * 0.96, year, in_dir, 'Bkg') ) ## Not used? - AWB 16.09.2018
        samples.append( SampleInfo('ZJets_AMC_1j_A', 'DYToLL_1J_13TeV-amcatnloFXFX-pythia8',                  888.9 * 0.86*0.985*0.995, year, in_dir, 'Bkg') )
        ## samples.append( SampleInfo('ZJets_AMC_1j_B', 'DYToLL_1J_13TeV-amcatnloFXFX-pythia8',                 888.9 * 0.86*0.985*0.995, year, in_dir, 'Bkg') ) ## Not used? - AWB 16.09.2018
        samples.append( SampleInfo('ZJets_AMC_2j_A', 'DYToLL_2J_13TeV-amcatnloFXFX-pythia8',                  348.8 * 0.88*0.975*0.992, year, in_dir, 'Bkg') )
        ## samples.append( SampleInfo('ZJets_AMC_2j_B', 'DYToLL_2J_13TeV-amcatnloFXFX-pythia8',                 348.8 * 0.88*0.975*0.992, year, in_dir, 'Bkg') ) ## Not used? - AWB 16.09.2018

    if (year == 2016):
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
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('ZJets_MG_1',            'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',                   5765.4,             year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_MG_2',            'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',                   5765.4,             year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZJets_hiM_AMC',         'DYJetsToLL_M-105To160_TuneCP5_PSweights_13TeV-amcatnloFXFX-pythia8',   46.948,            year, in_dir, 'Bkg') ) ## A. Marini
        samples.append( SampleInfo('ZJets_hiM_MG',          'DYJetsToLL_M-105To160_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',    46.948,            year, in_dir, 'Bkg') ) ## A. Marini
        ## samples.append( SampleInfo('ZJets_m_10_50',         'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8',              18610.0,             year, in_dir, 'Bkg') )
    
    ## TTJets
    if (year == 2016):
        samples.append( SampleInfo('tt_ll_AMC', 'TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8', 85.656 * 0.9, year, in_dir, 'Bkg') ) ## Why the factor of 0.9?!? - AWB 28.09.2018
        samples.append( SampleInfo('tt_ll_MG_1', 'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   85.656,       year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018
        samples.append( SampleInfo('tt_ll_MG_2', 'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   85.656,       year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('tt_ll_MG',  'TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8',       85.656, year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018
        samples.append( SampleInfo('tt_ll_POW', 'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',      85.656, year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.09.2018

    ## SingleTop
    if (year == 2016): tZq_DAS = 'tZq_ll_4f_13TeV-amcatnlo-pythia8'
    if (year == 2017 or year == 2018): tZq_DAS = 'tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8'

    if (year == 2016):
        samples.append( SampleInfo('tW_pos_1', 'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',        35.85,     year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_pos_2', 'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',        35.85,     year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_neg_1', 'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',    35.85,     year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_neg_2', 'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',    35.85,     year, in_dir, 'Bkg') )
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('tW_pos',   'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',     35.85,     year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tW_neg',   'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8', 35.85,     year, in_dir, 'Bkg') )
    samples.append(     SampleInfo('tZq',       tZq_DAS,                                                               0.0942,   year, in_dir, 'Bkg') )
    samples.append(     SampleInfo('tZW',      'ST_tWll_5f_LO_13TeV-MadGraph-pythia8',                                 0.01103,  year, in_dir, 'Bkg') )

    ## TTX
    if (year == 2016): ttH_tune = '_TuneCUETP8M2_ttHtranche3_'
    if (year == 2017 or year == 2018): ttH_tune = '_TuneCP5_'

    if (year == 2016):
        samples.append( SampleInfo('ttW_1', 'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',      0.2043, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttW_2', 'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8',      0.2043, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttZ',   'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8',              0.2529, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttH',   'ttHToNonbb_M125'+ttH_tune+'13TeV-powheg-pythia8',                   0.2151, year, in_dir, 'Bkg') )
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('ttW',   'TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8', 0.2043, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttZ',   'TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8',                   0.2529, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttH',   'ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8',     0.2151, year, in_dir, 'Bkg') )

    ## Diboson
    if (year == 2016): WW_DAS = 'WWTo2L2Nu_13TeV-powheg'
    if (year == 2017 or year == 2018): WW_DAS = 'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8'

    samples.append(     SampleInfo('WW',              WW_DAS,                                          12.46,    year, in_dir, 'Bkg') )
    samples.append(     SampleInfo('WZ_2l',          'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',     4.409,   year, in_dir, 'Bkg') ) ## Consistent with current WZ_3l, see below - AWB 28.09.2018
    if (year == 2016):
        samples.append( SampleInfo('WZ_3l_AMC',      'WZTo3LNu'+py_tune+'13TeV-amcatnloFXFX-pythia8',   4.430,   year, in_dir, 'Bkg') ) ## From TOP-18-008. We used 2.113 in 2016 - AWB 09.10.2018
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('WZ_3l',          'WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8',     4.430,   year, in_dir, 'Bkg') ) ## From TOP-18-008. We used 2.113 in 2016 - AWB 09.10.2018
    samples.append(     SampleInfo('ZZ_2l_2q',       'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',     3.22,    year, in_dir, 'Bkg') )
    if (year == 2016):
        samples.append( SampleInfo('ZZ_2l_2v',       'ZZTo2L2Nu_13TeV_powheg_pythia8',                  0.564,   year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ZZ_4l_AMC',      'ZZTo4L_13TeV-amcatnloFXFX-pythia8',               1.212,   year, in_dir, 'Bkg') ) ## From AN-2018/340
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('ZZ_4l',          'ZZTo4L_13TeV_powheg_pythia8',                     1.256,   year, in_dir, 'Bkg') ) ## From AN-2018/340
    if (year == 2017 or year == 2018):
        samples.append( SampleInfo('ZZ_4l_gg_2e2mu', 'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8', 0.00319, year, in_dir, 'Bkg') ) ## From AN-2018/340
        samples.append( SampleInfo('ZZ_4l_gg_4mu',   'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8',   0.00159, year, in_dir, 'Bkg') ) ## From AN-2018/340

    ## Triboson
    if (year == 2016): WWZ_str = 'WWZ'
    if (year == 2017 or year == 2018): WWZ_str = 'WWZ_4F'

    samples.append( SampleInfo('WWW', 'WWW_4F'+py_tune+'13TeV-amcatnlo-pythia8', 0.2086,  year, in_dir, 'Bkg') )
    samples.append( SampleInfo('WWZ', WWZ_str+py_tune+'13TeV-amcatnlo-pythia8',  0.1651,  year, in_dir, 'Bkg') )
    samples.append( SampleInfo('WZZ', 'WZZ'+py_tune+'13TeV-amcatnlo-pythia8',    0.05565, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZZZ', 'ZZZ'+py_tune+'13TeV-amcatnlo-pythia8',    0.01398, year, in_dir, 'Bkg') )

    # ## To do - find cross sections and include later
    # if (year == 2017 or year == 2018):
    #     samples.append( SampleInfo('tHq') )
    #     samples.append( SampleInfo('tHW') )
    #     samples.append( SampleInfo('ttWW') )
    #     samples.append( SampleInfo('ttZ_lowM') )

    ## add lepMVA samples for tests. first empty samples vector to keep the integraty of the original list -- XWZ 27.12.2018
    if (location == 'CERN_lepMVA_test_v1' and year == 2017):
	samples = []
	samples.append( SampleInfo('SingleMu_2017B', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017C', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017D', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017E', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017F', 'SingleMuon', 1, year, in_dir, 'Data') )

	samples.append( SampleInfo('H2Mu_gg', 'GluGluHToMuMu_M125_13TeV_amcatnloFXFX_pythia8', 0.009618, year, in_dir, 'Bkg' ) )
        samples.append( SampleInfo('H2Mu_ZH_120', 'ZH_HToMuMu_ZToAll_M120_13TeV_powheg_pythia8', 0.0002136, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Mu_ZH_125', 'ZH_HToMuMu_ZToAll_M125_13TeV_powheg_pythia8', 0.0002136, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Mu_ZH_130', 'ZH_HToMuMu_ZToAll_M130_13TeV_powheg_pythia8', 0.0002136, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Mu_WH_pos_120', 'WplusH_HToMuMu_WToAll_M120_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Mu_WH_pos_125', 'WplusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Mu_WH_pos_130', 'WplusH_HToMuMu_WToAll_M130_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('H2Mu_WH_neg_125', 'WminusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8',  0.0001164, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Mu_WH_neg_130', 'WminusH_HToMuMu_WToAll_M130_13TeV_powheg_pythia8',  0.0001164, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Mu_ttH_125', 'ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8', 0.0001103, year, in_dir, 'Bkg' ) )

	samples.append( SampleInfo('ZJets_AMC', 'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8', 5765.4, year, in_dir, 'Bkg') )

 	samples.append( SampleInfo('tt_ll_MG',  'TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8',    85.656, year ,in_dir, 'Bkg') )
	samples.append( SampleInfo('tt_ll_POW', 'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',   85.656, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('tt_ljj_POW_1', 'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8', 358.46, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('tt_ljj_POW_2', 'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8', 358.46, year, in_dir, 'Bkg') )

	samples.append( SampleInfo('tW_neg', 'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8', 35.85, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('tW_pos', 'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',     35.85, year, in_dir, 'Bkg') )

	samples.append( SampleInfo('WZ_3l',  'WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8', 4.430, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('ZZ_4l',  'ZZTo4L_13TeV_powheg_pythia8', 1.212, year, in_dir, 'Bkg') )
	#samples.append( SampleInfo() )


    if (location == 'CERN_lepMVA_test_v2' and year == 2017):
	samples = []
	samples.append( SampleInfo('SingleMu_2017B', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017C', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017D', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017E', 'SingleMuon', 1, year, in_dir, 'Data') )
	samples.append( SampleInfo('SingleMu_2017F', 'SingleMuon', 1, year, in_dir, 'Data') )

	samples.append( SampleInfo('ZJets_AMC', 'DY', 							5765.4, year, in_dir, 'Bkg') )

        samples.append( SampleInfo('tt_ll_MG', 	'TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8', 	85.656, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tt_ll_POW', 'TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8', 	85.656, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tt_ljj_POW_1', 'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8', 358.46, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('tt_ljj_POW_2', 'TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8', 358.46, year, in_dir, 'Bkg') )

	samples.append( SampleInfo('tW_neg', 'ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8', 35.85, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('tW_pos', 'ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8',     35.85, year, in_dir, 'Bkg') )


    if (location == 'CERN_lepMVA_3l_test_v1_with_extra_samples' and year == 2017):
	in_dir = '/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/2017/94X_v2/2019_01_15_LepMVA_3l_test_v1' #sample from loc = 'CERN' not cleared, only extra samples are listed here
	# for whose xsec is not found yet, use 1 as a placeholder
	# diHiggs
	samples.append( SampleInfo('HH_n12_incl_2', 'GluGluToHHTo2B2Tau_node_12_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') ) # HH To BBTauTau, not sure meaning of node
	samples.append( SampleInfo('HH_n2_incl_2', 'GluGluToHHTo2B2Tau_node_2_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('HH_n7_incl_2', 'GluGluToHHTo2B2Tau_node_7_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('HH_SM_incl_2', 'GluGluToHHTo2B2Tau_node_SM_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('HH_n12_incl_1', 'GluGluToHHTo4B_node_12_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('HH_n2_incl_1', 'GluGluToHHTo4B_node_2_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('HH_n7_incl_1', 'GluGluToHHTo4B_node_7_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('HH_SM_incl_1', 'GluGluToHHTo4B_node_SM_13TeV-madgraph', 0.00243, year, in_dir, 'Bkg') )
	# H2W
	samples.append( SampleInfo('H2W_ZH_125', 'GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5', 0.1888, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2W_WH_neg_125', 'HWminusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5', 0.1138, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2W_WH_pos_125', 'HWplusJ_HToWW_M125_13TeV_powheg_pythia8_TuneCP5', 0.1796, year, in_dir, 'Bkg') )
	# THX
	samples.append( SampleInfo('tHq', 'THQ_4f_Hincl_13TeV_madgraph_pythia8', 1, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('tHW', 'THW_5f_Hincl_13TeV_madgraph_pythia8', 1, year, in_dir, 'Bkg') )
	# ttX
	samples.append( SampleInfo('ttWW', 'TTWW_TuneCP5_13TeV-madgraph-pythia8', 1, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('ttZ_lowM', 'TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8', 1, year, in_dir, 'Bkg') )
	# HtoOthers
	samples.append( SampleInfo('H2Tau_WH_neg_125', 'WminusHToTauTau_M125_13TeV_powheg_pythia8', 0.03339, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Tau_WH_pos_125', 'WplusHToTauTau_M125_13TeV_powheg_pythia8', 0.05268, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Z_WH_neg_125', 'WminusH_HToZZTo2L2X_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8', 0.01395, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Z_WH_pos_125', 'WplusH_HToZZTo2L2X_M125_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8', 0.02200, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Tau_ZH_125', 'ZHToTauTau_M125_13TeV_powheg_pythia8', 0.05541, year, in_dir, 'Bkg') )
	samples.append( SampleInfo('H2Z_ZH_125', 'ZH_HToZZ_2LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8', 0.02314, year, in_dir, 'Bkg') )

    return samples

## End function: GetSamples(location = 'CERN', year = '2017')
