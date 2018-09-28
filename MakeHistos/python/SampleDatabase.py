

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
            in_dir = '/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/Moriond17/Mar13'
        elif (year == 2017):
            in_dir = '/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17'
        else:
            print 'Invalid location (%s) and/or year (%d)!!!  Exiting.' % (location, year)
            sys.exit
    else:
        print 'Invalid location (%s) and/or year (%d)!!!  Exiting.' % (location, year)


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


    ########################
    ###  Get MC samples  ###
    ########################

    ################
    ###  Signal  ###
    ################

    if (year == 2016): sig_gen = '_powheg_pythia8'
    if (year == 2017): sig_gen = '_amcatnloFXFX_pythia8'

    if (year == 2016): VBF_str = 'VBF'
    if (year == 2017): VBF_str = 'VBFH_'

    ## H2Mu_gg
    samples.append( SampleInfo('H2Mu_gg',     'GluGlu_HToMuMu_M125_13TeV'+sig_gen, 0.009618, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_gg_120', 'GluGlu_HToMuMu_M120_13TeV'+sig_gen, 0.009618, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    samples.append( SampleInfo('H2Mu_gg_130', 'GluGlu_HToMuMu_M130_13TeV'+sig_gen, 0.009618, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    
    ## H2Mu_VBF
    samples.append( SampleInfo('H2Mu_VBF',     VBF_str+'HToMuMu_M125_13TeV'+sig_gen, 0.0008208, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_VBF_120', VBF_str+'HToMuMu_M120_13TeV'+sig_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    samples.append( SampleInfo('H2Mu_VBF_130', VBF_str+'HToMuMu_M130_13TeV'+sig_gen, 0.0008208, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    
    ## H2Mu_VH
    samples.append( SampleInfo('H2Mu_ZH',     'ZH_HToMuMu_M125_13TeV_powheg_pythia8', 0.0002136, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_ZH_120', 'ZH_HToMuMu_M120_13TeV_powheg_pythia8', 0.0002136, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    samples.append( SampleInfo('H2Mu_ZH_130', 'ZH_HToMuMu_M130_13TeV_powheg_pythia8', 0.0002136, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    
    samples.append( SampleInfo('H2Mu_WH_pos',     'WPlusH_HToMuMu_M125_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_WH_pos_120', 'WPlusH_HToMuMu_M120_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    samples.append( SampleInfo('H2Mu_WH_pos_130', 'WPlusH_HToMuMu_M130_13TeV_powheg_pythia8', 0.0001858, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    
    samples.append( SampleInfo('H2Mu_WH_neg',     'WMinusH_HToMuMu_M125_13TeV_powheg_pythia8', 0.0001164, year, in_dir, 'Sig') )
    samples.append( SampleInfo('H2Mu_WH_neg_120', 'WMinusH_HToMuMu_M120_13TeV_powheg_pythia8', 0.0001164, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    samples.append( SampleInfo('H2Mu_WH_neg_130', 'WMinusH_HToMuMu_M130_13TeV_powheg_pythia8', 0.0001164, year, in_dir, 'Sig') ) ## Presumably incorrect? - AWB 16.08.2018
    
    ## H2Mu_ttH
    if (year == 2017):
        samples.append( SampleInfo('H2Mu_ttH', 'ttH_HToMuMu_M125_13TeV_powheg_pythia8', 0.0001103, year, in_dir, 'Sig') )  ## 0.5071 x 0.0002176 from YR4 (125 GeV)

    
    ####################
    ###  Background  ###
    ####################

    if (year == 2016): py_tune = '_TuneCUETP8M1_'
    if (year == 2017): py_tune = '_TuneCP5_'
    
    ## DYJetsToLL
    samples.append( SampleInfo('ZJets_AMC',      'DYJetsToLL_M-50'+py_tune+'13TeV-amcatnloFXFX-pythia8', 5765.4,        year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (6020.85) - AWB 28.08.2018
    if (year == 2017):
        samples.append( SampleInfo('ZJets_AMC_2','DYJetsToLL_M-50'+py_tune+'13TeV-amcatnloFXFX-pythia8', 5765.4,        year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_AMC_0j_A', 'DYToLL_0J_13TeV-amcatnloFXFX-pythia8',                 4754   * 0.96, year, in_dir, 'Bkg') )
    ## samples.append( SampleInfo('ZJets_AMC_0j_B', 'DYToLL_0J_13TeV-amcatnloFXFX-pythia8',                4754   * 0.96, year, in_dir, 'Bkg') ) ## Not used? - AWB 16.08.2018
    samples.append( SampleInfo('ZJets_AMC_1j_A', 'DYToLL_1J_13TeV-amcatnloFXFX-pythia8',                  888.9 * 0.86*0.985*0.995, year, in_dir, 'Bkg') )
    ## samples.append( SampleInfo('ZJets_AMC_1j_B', 'DYToLL_1J_13TeV-amcatnloFXFX-pythia8',                 888.9 * 0.86*0.985*0.995, year, in_dir, 'Bkg') ) ## Not used? - AWB 16.08.2018
    samples.append( SampleInfo('ZJets_AMC_2j_A', 'DYToLL_2J_13TeV-amcatnloFXFX-pythia8',                  348.8 * 0.88*0.975*0.992, year, in_dir, 'Bkg') )
    ## samples.append( SampleInfo('ZJets_AMC_2j_B', 'DYToLL_2J_13TeV-amcatnloFXFX-pythia8',                 348.8 * 0.88*0.975*0.992, year, in_dir, 'Bkg') ) ## Not used? - AWB 16.08.2018
    
    samples.append( SampleInfo('ZJets_MG',              'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',              5765.4,             year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_70_100',    'DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    178.952    * 0.98, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_100_200',   'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   181.302    * 0.96, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_200_400',   'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    50.4177   * 0.96, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_400_600',   'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',     6.98394  * 0.96, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_600_800',   'DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',     1.68141  * 0.96, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_800_1200',  'DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    0.775392 * 0.96, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_1200_2500', 'DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',   0.186222 * 0.96, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_MG_HT_2500_inf',  'DYJetsToLL_M-50_HT-2500toinf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    0.004385 * 0.96, year, in_dir, 'Bkg') )
    
    samples.append( SampleInfo('ZJets_hiM',          'DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8', 5765.4 * 1.235, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZJets_hiM_SpringPU', 'DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8', 5765.4 * 1.235, year, in_dir, 'Bkg') )
    if (year == 2017):
        samples.append( SampleInfo('ZJets_m_10_50', 'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8', 18610.0, year, in_dir, 'Bkg') )
    
    ## TTJets
    samples.append( SampleInfo('tt_ll_AMC', 'TTJets_Dilept_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8', 85.656 * 0.9, year, in_dir, 'Bkg') ) ## Why the factor of 0.9?!? - AWB 28.08.2018
    samples.append( SampleInfo('tt_ll_MG',  'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',    85.656,       year, in_dir, 'Bkg') ) ## A bit lower than TOP-18-008 (87.315) - AWB 28.08.2018
    if (year == 2017):
        samples.append( SampleInfo('tt', 'TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8', 815.96, year, in_dir, 'Bkg') )

    ## SingleTop
    if (year == 2016): tZq_DAS = 'tZq_ll_4f_13TeV-amcatnlo-pythia8'
    if (year == 2017): tZq_DAS = 'tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8'

    samples.append( SampleInfo('tW_pos', 'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1',     35.85,    year, in_dir, 'Bkg') )
    samples.append( SampleInfo('tW_neg', 'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1', 35.85,    year, in_dir, 'Bkg') )
    samples.append( SampleInfo('tZq',    tZq_DAS,                                                             0.0942,  year, in_dir, 'Bkg') )
    samples.append( SampleInfo('tZW',    'ST_tWll_5f_LO_13TeV-MadGraph-pythia8',                              0.01103, year, in_dir, 'Bkg') )

    ## TTX
    if (year == 2016): ttH_tune = '_TuneCUETP8M2_ttHtranche3_'
    if (year == 2017): ttH_tune = '_TuneCP5_'

    if (year == 2016):
        samples.append( SampleInfo('ttW_1', 'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8', 0.2043, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttW_2', 'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8', 0.2043, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttZ',   'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8',         0.2529, year, in_dir, 'Bkg') )
    if (year == 2017):
        samples.append( SampleInfo('ttW',   'ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8',                    0.6244, year, in_dir, 'Bkg') )
        samples.append( SampleInfo('ttZ',   'ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8',                    0.8462, year, in_dir, 'Bkg') )
    samples.append(     SampleInfo('ttH',   'ttHToNonbb_M125'+ttH_tune+'13TeV-powheg-pythia8',              0.2151, year, in_dir, 'Bkg') )

    ## Diboson
    if (year == 2016): WW_DAS = 'WWTo2L2Nu_13TeV-powheg'
    if (year == 2017): WW_DAS = 'WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8'

    if (year == 2016): ZZ_4l_DAS = 'ZZTo4L_13TeV-amcatnloFXFX-pythia8'
    if (year == 2017): ZZ_4l_DAS = 'ZZTo4L_13TeV_powheg_pythia8'

    samples.append( SampleInfo('WW',        WW_DAS,                                         12.46,  year, in_dir, 'Bkg') )
    samples.append( SampleInfo('WZ_2l',     'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',   4.409, year, in_dir, 'Bkg') ) ## Consistent with current WZ_3l, see below - AWB 28.08.2018
    samples.append( SampleInfo('WZ_3l_AMC', 'WZTo3LNu'+py_tune+'13TeV-amcatnloFXFX-pythia8', 2.113, year, in_dir, 'Bkg') ) ## Much lower than TOP-18-008 (4.4297) - AWB 28.08.2018
    samples.append( SampleInfo('ZZ_2l_2v',  'ZZTo2L2Nu_13TeV_powheg_pythia8',                0.564, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZZ_2l_2q',  'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',   3.22,  year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZZ_4l',     ZZ_4l_DAS,                                       1.212, year, in_dir, 'Bkg') )

    ## Triboson
    if (year == 2016): WWZ_str = 'WWZ'
    if (year == 2017): WWZ_str = 'WWZ_4F'

    samples.append( SampleInfo('WWW', 'WWW_4F'+py_tune+'13TeV-amcatnlo-pythia8', 0.2086,  year, in_dir, 'Bkg') )
    samples.append( SampleInfo('WWZ', WWZ_str+py_tune+'13TeV-amcatnlo-pythia8',  0.1651,  year, in_dir, 'Bkg') )
    samples.append( SampleInfo('WZZ', 'WZZ'+py_tune+'13TeV-amcatnlo-pythia8',    0.05565, year, in_dir, 'Bkg') )
    samples.append( SampleInfo('ZZZ', 'ZZZ'+py_tune+'13TeV-amcatnlo-pythia8',    0.01398, year, in_dir, 'Bkg') )

    return samples

## End function: GetSamples(location = 'CERN', year = '2017')
