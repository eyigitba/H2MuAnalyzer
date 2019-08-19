
## Current directory
HOMEDIR=$PWD
## CMSSW release containing Higgs Combine
COMBINEDIR="/afs/cern.ch/work/x/xzuo/combine/CMSSW_8_1_0/src"
## Directory containing workspace directory
WSDIR="$HOMEDIR"
## Datacard name
declare -a DCARDS=(

############################################
#    "naive_case_test.txt"
#    "test_card.txt"

#    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/test_combine6_MC.txt"

############################################

    "out_files/ZH_Run2_lep_BDT_08_14/datacard/ZH_4l_BDT_and_mass_template_stack.txt"
    "out_files/ZH_Run2_lep_BDT_08_14/datacard/ZH_4l_BDT_mass_min_template_stack.txt"
    "out_files/ZH_Run2_lep_BDT_08_14/datacard/ZH_4l_BDT_mass_more_template_stack.txt"


    "out_files/ZH_Run2_lep_08_14/datacard/ZH_2cats_BWZ_MC_0814.txt"
    "out_files/ZH_Run2_lep_08_14/datacard/ZH_4l_dimu_mass_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
    "out_files/ZH_Run2_lep_08_14/datacard/ZH_4l_dimu_mass_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
    "out_files/ZH_Run2_lep_08_14/datacard/ZH_4l_dimu_mass_Gaus3_frzCoefs_frzParams_Bern3_shape_MC.txt"
    "out_files/ZH_Run2_lep_08_14/datacard/ZH_4l_dimu_mass_BDT_n10_p05_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
    "out_files/ZH_Run2_lep_08_14/datacard/ZH_4l_dimu_mass_BDT_p05_p10_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"

#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZRed2_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern3_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_expo2_shape_MC.txt"

#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_BDT_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_BDT_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_BDT_p02_p06_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_BDT_p06_p10_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_BDT_p060_p068_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_BDT_p068_p076_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_3lep_H_pair_mass_BDT_p076_p10_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"

#    "out_files/WH_Run2_3lep_07_31/datacard/WH_4cats_BWZ_MC_0731.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_4cats_Bern2_MC_0731.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_4cats_Bern3_MC_0731.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_4cats_expo2_MC_0731.txt"
#    "out_files/WH_Run2_3lep_07_31/datacard/WH_6cats_BWZ_MC_0731.txt"


#    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern3_shape_MC.txt"
#    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_expo2_shape_MC.txt"
##    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_BDT_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
##    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_BDT_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
##    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_BDT_p02_p06_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
##    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_BDT_p06_p10_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
##    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_BDT_p060_p068_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
##    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_BDT_p068_p076_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_3lep_H_pair_mass_BDT_p076_p10_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_4cats_BWZ_MC.txt"
#    "out_files/WH_Run2_3lep_AWB_mass_07_10/datacard/WH_6cats_BWZ_MC.txt"


#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_expo2_shape_MC.txt"
##    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
##    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
##    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p02_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
##    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p02_p05_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
##    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p05_p06_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
##    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p06_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_3cats_Bern2_MC.txt"
#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_4cats_Bern2_MC.txt"
#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/ttH_6cats_Bern2_MC.txt"

#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/WH_ttH_BWZ_MC.txt"
#    "out_files/ttH_Run2_3lep_AWB_mass_07_10/datacard/WH_ttH_subcats_MC.txt"



#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern3_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_med_3cats_Bern1.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_med_3cats_Bern2.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_med_4cats_Bern1.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_med_4cats_Bern2.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_p02_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_p02_p05_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_p05_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_p02_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_p02_p05_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_med_07_08/datacard/ttH_3lep_med_H_pair_mass_BDT_v1_med_p05_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"




#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern3_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_tight_slim_3cats_Bern1.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_tight_slim_3cats_Bern2.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_tight_slim_4cats_Bern1.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_tight_slim_4cats_Bern2.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p02_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p02_p05_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p05_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p02_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p02_p05_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_tight_07_08/datacard/ttH_3lep_tight_H_pair_mass_BDT_v1_tight_noBDT_p05_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"


#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern3_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3cats_Bern1.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3cats_Bern2.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_4cats_Bern1.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_4cats_Bern2.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_p02_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_p02_p05_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_p05_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern1_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_n10_n02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_n02_p02_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_p02_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_p02_p05_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"
#    "out_files/ttH_1718_3lep_AWB_mass_07_04/datacard/ttH_3lep_H_pair_mass_BDT_v1_med_p05_p10_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"




#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_AWB_v3_combo_rebin_group.txt"
#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_AWB_v3_combo_rebin_stack.txt"
#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_AWB_v3_combo_template_group.txt"
#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_AWB_v3_combo_template_stack.txt"
   
#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_XWZ_withMass_rebin_group.txt"
#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_XWZ_withMass_rebin_stack.txt"
#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_XWZ_withMass_template_group.txt"
#    "out_files/WH_1718_3lep_AWB_BDT_06_25_med/datacard/WH_3lep_BDT_XWZ_withMass_template_stack.txt"



#    "out_files/WH_1718_3lep_AWB_template_mass_07_03/datacard/WH_3lep_H_pair_mass_zoomH_template_group.txt"

#    "out_files/WH_1718_3lep_AWB_Bernstein_mass_07_03/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern3_frzCoefs_frzParams_shape_MC.txt"
#    "out_files/WH_1718_3lep_AWB_Bernstein_mass_07_03/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_Bern2_shape_MC.txt"


###    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/combine6cats_data.txt "
#    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/combine2cats_MC.txt"
#    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/combine6cats_MC.txt"

#    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/WH_3lep_H_pair_mass_BDT_p06_p10_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/WH_3lep_H_pair_mass_BDT_p076_p10_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_1718_3lep_AWB_mass_06_25/datacard/WH_3lep_H_pair_mass_BDT_p068_p076_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"



#    "out_files/WH_combine_2017_2018/datacard/combine_6cats_2017_2018_BWZ1_MC.txt"
#    "out_files/WH_combine_2017_2018/datacard/combine_BDT_2017_2018_rebin_stack.txt"

#    "out_files/WH_2018_3lep_AWB_mass_06_19_med/datacard/WH_3lep_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
#    "out_files/WH_2018_3lep_AWB_mass_06_19_med/datacard/combine2_2018.txt"
#    "out_files/WH_2018_3lep_AWB_mass_06_19_med/datacard/combine6_2018.txt"

#    "out_files/WH_2018_3lep_AWB_BDT_06_19_med/datacard/WH_3lep_BDT_AWB_v3_combo_rebin_group.txt"
#    "out_files/WH_2018_3lep_AWB_BDT_06_19_med/datacard/WH_3lep_BDT_AWB_v3_combo_template_group.txt"
#    "out_files/WH_2018_3lep_AWB_BDT_06_19_med/datacard/WH_3lep_BDT_AWB_v3_combo_rebin_stack.txt"
#    "out_files/WH_2018_3lep_AWB_BDT_06_19_med/datacard/WH_3lep_BDT_AWB_v3_combo_template_stack.txt"
#    "out_files/WH_2018_3lep_AWB_BDT_06_19_med/datacard/WH_3lep_BDT_XWZ_withMass_rebin_group.txt"
#    "out_files/WH_2018_3lep_AWB_BDT_06_19_med/datacard/WH_3lep_BDT_XWZ_withMass_rebin_stack.txt"


    # "out_files/WH_lep_AWB_2019_05_20_v1_shape/datacard/combine_2cats.txt"
    #"out_files/WH_lep_AWB_2019_05_20_v1_shape/datacard/combine_6cats.txt "
    # "out_files/WH_lep_AWB_2019_05_20_v1_shape/datacard/lep3_H_pair_mass_BDT_p06_p10_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"
    # "out_files/WH_lep_AWB_2019_05_20_v1_shape/datacard/lep3_H_pair_mass_zoomH_Gaus3_frzCoefs_frzParams_BWZ1_shape_MC.txt"



#    "out_files/ZH_lep_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_4l_BDT_and_mass_rebin_group.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_4l_BDT_mass_more_rebin_group.txt"
#    "out_files/ZH_lep_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_4l_dimu_mass_template_group.txt"
#    "out_files/ZH_lep_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_4l_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_shape_MC.txt"

#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_and_mass_rebin_group.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_and_mass_rebin_stack.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_and_mass_template_group.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_and_mass_template_stack.txt"
#
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_more_rebin_group.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_more_rebin_stack.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_more_template_group.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_more_template_stack.txt"
#
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_min_rebin_group.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_min_rebin_stack.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_min_template_group.txt"
#    "out_files/ZH_lep_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4l_BDT_mass_min_template_stack.txt"


#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_and_mass_rebin_group.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_and_mass_rebin_stack.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_and_mass_template_group.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_and_mass_template_stack.txt"
#
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_more_rebin_group.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_more_rebin_stack.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_more_template_group.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_more_template_stack.txt"
#
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_min_rebin_group.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_min_rebin_stack.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_min_template_group.txt"
#    "out_files/ZH_mu_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_4mu_BDT_mass_min_template_stack.txt"

#    "out_files/ZH_mu_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_4mu_dimu_mass_template_group.txt"
#    "out_files/ZH_mu_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_4mu_dimu_mass_template_stack.txt"
#    "out_files/ZH_mu_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_4mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_shape_MC.txt"
#    "out_files/ZH_mu_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_4mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_frzCoefs_frzParams_shape_MC.txt"


#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_and_mass_rebin_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_and_mass_rebin_stack.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_and_mass_template_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_and_mass_template_stack.txt"
#
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_more_rebin_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_more_rebin_stack.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_more_template_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_more_template_stack.txt"
#
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_min_rebin_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_min_rebin_stack.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_min_template_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_noMVA_v1/datacard/ZH_2e2mu_BDT_mass_min_template_stack.txt"

#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_and_mass_rebin_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_and_mass_rebin_stack.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_and_mass_template_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_and_mass_template_stack.txt"

#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_more_rebin_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_more_rebin_stack.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_more_template_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_more_template_stack.txt"

#   "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_min_rebin_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_min_rebin_stack.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_min_template_group.txt"
#    "out_files/ZH_ele_XWZ_BDTs_05_23_MVAp04_v1/datacard/ZH_2e2mu_BDT_mass_min_template_stack.txt"


#    "out_files/ZH_ele_XWZ_mass_05_23_noMVA_v1/datacard/ZH_2e2mu_dimu_mass_template_group.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_noMVA_v1/datacard/ZH_2e2mu_dimu_mass_template_stack.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_noMVA_v1/datacard/ZH_2e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_shape_MC.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_noMVA_v1/datacard/ZH_2e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_frzCoefs_frzParams_shape_MC.txt"
#
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAn04_v1/datacard/ZH_2e2mu_dimu_mass_template_group.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAn04_v1/datacard/ZH_2e2mu_dimu_mass_template_stack.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAn04_v1/datacard/ZH_2e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_shape_MC.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAn04_v1/datacard/ZH_2e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_frzCoefs_frzParams_shape_MC.txt"
#
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_2e2mu_dimu_mass_template_group.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_2e2mu_dimu_mass_template_stack.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_2e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_shape_MC.txt"
#    "out_files/ZH_ele_XWZ_mass_05_23_MVAp04_v1/datacard/ZH_2e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_frzCoefs_frzParams_shape_MC.txt"


#    "out_files/WH_mu_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_template_stack.txt"
#    "out_files/WH_mu_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_template_group.txt"
#    "out_files/WH_mu_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_rebin_stack.txt"
#    "out_files/WH_mu_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_rebin_group.txt"

#    "out_files/WH_ele_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_template_stack.txt"
#    "out_files/WH_ele_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_template_group.txt"
#    "out_files/WH_ele_XWZ_2019_05_15_v1_merge3/datacard/e2mu_BDTG_UF_v2_template_stack.txt"
#    "out_files/WH_ele_XWZ_2019_05_15_v1_merge3/datacard/e2mu_BDTG_UF_v2_template_group.txt"

#    "out_files/WH_ele_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_rebin_stack.txt"
#    "out_files/WH_ele_XWZ_2019_05_15_v1/datacard/e2mu_BDTG_UF_v2_rebin_group.txt"

#    "out_files/WH_ele_XWZ_mass_05_16_v1/datacard/e2mu_dimu_mass_template_group.txt"
#    "out_files/WH_ele_XWZ_mass_05_16_v1/datacard/e2mu_dimu_mass_template_stack.txt"
#    "out_files/WH_ele_XWZ_mass_05_16_v1/datacard/e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_frzCoefs_frzParams_shape_MC.txt"
#    "out_files/WH_ele_XWZ_mass_05_16_v1/datacard/e2mu_dimu_mass_Gaus3_frzCoefs_frzParams_BWZRed3_shape_MC.txt"

#     "out_files/datacard/WH_ele_lepMVA04_BDT_neg_template_MC.txt"
#     "out_files/datacard/WH_ele_lepMVA04_BDT_pos_template_MC.txt"
#     "out_files/datacard/WH_ele_lepMVA04_BDT_template_no_mass.txt"
#     "out_files/datacard/WH_ele_lepMVA04_BDT_template_no_mass_wrong.txt"
  #   "out_files/datacard/WH_ele_lepMVA04_BDT_template.txt"
#     "out_files/datacard/WH_ele_lepMVA04_BDT_template_wrong.txt"
#     "out_files/datacard/WH_ele_lepMVA04_BDT_template_50.txt"
#     "out_files/datacard/WH_ele_lepMVA04_BDT_template_wrong_50.txt"
  #   "out_files/datacard/WH_ele_lepMVA04_template_MC.txt"
  #   "out_files/datacard/WH_ele_lepMVA04_MC.txt"
  #   "out_files/datacard/WH_ele_lepMVA04_data.txt"
#     "out_files/datacard/WH_ele_lepMVA04_Counting_120_130_Correlated.txt"
#     "out_files/datacard/WH_ele_lepMVA04_Counting_120_130_Uncorrelated.txt"     


#    "out_files/datacard/WH_mu_444_template_MC.txt"
#    "out_files/datacard/WH_mu_884_template_MC.txt"
#    "out_files/datacard/WH_mu_888_template_MC.txt"
#    "out_files/datacard/WH_mu_no_template_MC.txt"

#    "out_files/datacard/WH_ele_n4n4n4_template_MC.txt"
#    "out_files/datacard/WH_ele_n4n44_template_MC.txt"
#    "out_files/datacard/WH_ele_4n44_template_MC.txt"
#    "out_files/datacard/WH_ele_444_template_MC.txt"
#    "out_files/datacard/WH_ele_no_template_MC.txt"

#    "out_files/datacard/cut_cat_sum20_template_MC.txt"
#    "out_files/datacard/cut_cat_sum_template_MC.txt"
#    "out_files/datacard/cut_cat_template_MC.txt"
#    "out_files/datacard/cut_cat_MC.txt"
#    "out_files/datacard/cut_cat_data.txt"
)


## Unique label name
#LABEL="-n _T400_syst_0p3_fit"
#LABEL="-n _T10_syst_0p3_fit"
#LABEL="-n _test_XWZ_BDT_temp_v1"
#LABEL="-n _test_T2500_bias_dev"
#LABEL="-n _test_6cats_T2500"
#LABEL="-n _Bern3_inc_fixShapeMC_T2500"
#LABEL="-n _WH_mass_template_T2500"
#LABEL="-n _BWZ_against_Bern3"
LABEL="-n _ttH_3lep"
## Expected signal strength for limits and significance
MULIM="--expectSignal=0"
MUSIG="--expectSignal=1"
## Indicate whether to include systematics
SYST="-S 1"  ## "-S 0" or "-S 1"
## Indicate range of allowed r-values
RANGE="--rMin -100 --rMax 100"
## Specify number of toys (-t N), use a random seed (-s -1), and save toys
TOYS="-t 2500 --saveToys"
TEMPLATE_TOYS="--toysFile higgsCombine_WH_mass_template_T2500.FitDiagnostics.mH120.123456.root" 
FIXBERN_TOYS="--toysFile higgsCombine_Bern3_inc_fixShapeMC_T2500.FitDiagnostics.mH120.123456.root"
FREEBERN_TOYS="--toysFile higgsCombine_Bern_inc_freeShapeMC_T800.FitDiagnostics.mH120.123456.root"
## Generate toys with unconstrained systematics (""), no systematics,
##   or systematics constrained from fit to data ("toysFrequentist")
FREQ="--toysNoSystematics"  ## "", "--toysNoSystematics", or "--toysFrequentist"
## Save output
##   * "--saveShapes" and "--plots" only save info from the last toy
##   * "--saveShapes" and "--saveWithUncertainties" can't be run simultaneously
SAVE="--trackParameters rgx{.*} --saveNormalizations --saveWithUncertainties --savePredictionsPerToy"




cd $COMBINEDIR
echo "Running 'cmsenv' in CMSSW release with Higgs Combine: $COMBINEDIR"
eval `scramv1 runtime -sh`
echo "Navigating to workspace directory: $WSDIR"
cd $WSDIR
for DCARD in "${DCARDS[@]}"
do
    echo "---------------------------------------------------------------------------------------------------"
    echo "-------------------------------- $DCARD --------------------------------"
    echo "---------------------------------------------------------------------------------------------------"

    echo "Running combine -M AsymptoticLimits $DCARD --run blind $MUSIG $SYST $FREQ $LABEL"
    combine -M AsymptoticLimits $DCARD --run blind $MUSIG $SYST $FREQ $LABEL
    echo "---------------------------------------------"

    echo "---------------------------------------------"
    echo "Running combine -M Significance --uncapped=true $DCARD -t -1  $MUSIG $SYST $FREQ $LABEL"
    combine -M Significance --uncapped=true $DCARD -t -1  $MUSIG $SYST $FREQ $LABEL
    echo "---------------------------------------------"


#    echo "---------------------------------------------"
#    echo "combine -M GenerateOnly $TOYS $DCARD $MUSIG $SYST $FREQ $LABEL $RANGE"
#    combine -M GenerateOnly $TOYS $DCARD $MUSIG $SYST $FREQ $LABEL $RANGE
#    echo "---------------------------------------------"


#    echo "---------------------------------------------"
#    echo "combine -M FitDiagnostics $TOYS $DCARD $MUSIG $SYST $FREQ $LABEL $RANGE"
#    combine -M FitDiagnostics $TOYS $DCARD $MUSIG $SYST $FREQ $LABEL $RANGE
#    echo "---------------------------------------------"


##    echo "---------------------------------------------"
##    echo "combine -M FitDiagnostics $DCARD $FIXBERN_TOYS $TOYS $SYST $LABEL $RANGE"
##    combine -M FitDiagnostics $DCARD $FIXBERN_TOYS $TOYS $SYST $LABEL $RANGE
##    echo "---------------------------------------------"





#    echo "---------------------------------------------"
#    echo "Running combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=10"
#    combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=10
#    echo "---------------------------------------------"

#    echo "---------------------------------------------"
#    echo "Running combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=50"
#    combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=50
#    echo "---------------------------------------------"

#    echo "---------------------------------------------"
#    echo "Running combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=1 --toysFreq"
#    combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=1 --toysFreq
#    echo "---------------------------------------------"
#    echo "Running combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=0"
#    combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=0
#    echo "---------------------------------------------"
#    echo "combine -M MaxLikelihoodFit --uncapped=true $DCARD"
#    combine -M MaxLikelihoodFit --uncapped=true $DCARD
#    echo "---------------------------------------------"
#    echo "combine -M FitDiagnostics $DCARD --rMin -100 --rMax 100"
#    combine -M FitDiagnostics $DCARD --rMin -100 --rMax 100
#    echo "---------------------------------------------"
#    echo "combine -M FitDiagnostics $DCARD -t 1000 --rMin -100 --rMax 100"
#    combine -M FitDiagnostics $DCARD -t 1000 --rMin -100 --rMax 100
#    echo "---------------------------------------------"
#    echo "Running combine -M Significance --signif $DCARD -t 100 --expectSignal=1"
#    combine -M Significance --signif $DCARD -t 100 --expectSignal=1

#    echo "Running combine -M AsymptoticLimits $DCARD --run blind -S 0 --expectSignal 0"
#    combine -M AsymptoticLimits $DCARD --run blind -S 0 --expectSignal 0
#    echo "Running combine -M AsymptoticLimits $DCARD --run blind --expectSignal 0"
#    combine -M AsymptoticLimits $DCARD --run blind --expectSignal 0
    # echo "Running combine -M MultiDimFit $DCARD -S 0 --preFitValue 0 --do95 1"
    # combine -M MultiDimFit $DCARD -S 0 --preFitValue 0 --do95 1
    # echo "Running combine -M FitDiagnostics $DCARD -S 0 --preFitValue 0"
    # combine -M FitDiagnostics $DCARD -S 0 --preFitValue 0
done
echo "Returning to home directory and running 'cmsenv': $HOMEDIR"
cd $HOMEDIR
eval `scramv1 runtime -sh`
echo "Finished running run_test_limits.sh!"
