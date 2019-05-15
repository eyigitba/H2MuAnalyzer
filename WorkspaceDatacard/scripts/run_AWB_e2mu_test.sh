
## Current directory
HOMEDIR=$PWD
## CMSSW release containing Higgs Combine
COMBINEDIR="/afs/cern.ch/user/a/abrinke1/HiggsToMuMu/Combine/CMSSW_8_1_0/src"
## Directory containing workspace directory
WSDIR="$HOMEDIR"
## Datacard name
declare -a DCARDS=(

    "out_files/WH_lep_AWB_2019_05_14_TMVA_retrain_v1/datacard/e2mu_BDTG_UF_v1_rebin_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_14_TMVA_retrain_v1/datacard/e2mu_BDTG_UF_v1_template_stack.txt"

    "out_files/WH_lep_XWZ_2019_05_14_TMVA_out_v1/datacard/e2mu_BDTG_UF_v2_rebin_stack.txt"
    # "out_files/WH_lep_XWZ_2019_05_14_TMVA_out_v1/datacard/e2mu_BDTG_UF_v2_template_stack.txt"

    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_BDT_mass_template_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v3/datacard/mu3_BDT_mass_rebin_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v3/datacard/mu3_BDT_mass_template_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_on_template_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_on_template_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_zoom_Gaus3_frzCoefs_frzParams_BWZRed1_frzCoefs_frzParams_shape_MC.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_zoom_Gaus3_frzCoefs_frzParams_BWZRed1_shape_MC.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_zoom_rebin_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_zoom_rebin_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_zoom_template_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/mu3_H_mass_zoom_template_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_BDT_mass_template_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_BDT_mass_template_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_on_template_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_on_template_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_zoom_Gaus3_frzCoefs_frzParams_BWZRed1_frzCoefs_frzParams_shape_MC.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_zoom_Gaus3_frzCoefs_frzParams_BWZRed1_shape_MC.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_zoom_rebin_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_zoom_rebin_stack.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_zoom_template_group.txt"
    # "out_files/WH_lep_AWB_2019_05_01_v1/datacard/e2mu_H_mass_zoom_template_stack.txt"

    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_H_mass_zoom_template_group.txt"
    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_H_mass_zoom_template_stack.txt"
    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_H_mass_zoom_shape_MC.txt"
    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_H_mass_zoom_shape_data.txt"
    
    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_BDT_mass_template_stack.txt"
    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_BDT_mass_template_group.txt"
    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_BDT_mass_rebin_stack.txt"
    # "out_files/datacard/e2mu_medLepMVA_noZ_noBtag_mass12_BDT_mass_rebin_group.txt"

    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_MC.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_h_H_mass_on_template_MC_autoMCStats-10.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_h_H_mass_on_template_MC.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_ge2j_btag_mass12_h_H_mass_on_template_MC_autoMCStats0.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_template_MC_autoMCStats500.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_template_MC_autoMCStats0.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_template_MC.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_ge2j_btag_mass12_MC.txt"
    # "out_files/datacard/cat_e2mu_medLepMVA_noZ_ge2j_btag_mass12_template_MC.txt"
)

## Unique label name
#LABEL="-n _T400_syst_0p3_fit"
#LABEL="-n _T10_syst_0p3_fit"
LABEL="-n _test_XWZ_new_script_v1"
## Expected signal strength for limits and significance
MULIM="--expectSignal=0"
MUSIG="--expectSignal=1"
## Indicate whether to include systematics
SYST="-S 1"  ## "-S 0" or "-S 1"
## Indicate range of allowed r-values
RANGE="--rMin -100 --rMax 100"
## Specify number of toys (-t N), use a random seed (-s -1), and save toys
TOYS="-t 400 --saveToys"
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
    echo "------------------ $DCARD ------------------"
    echo "---------------------------------------------------------------------------------------------------"
    # echo "combine -M FitDiagnostics $DCARD $RANGE $MULIM $SYST $TOYS $FREQ $SAVE $LABEL"
    #       combine -M FitDiagnostics $DCARD $RANGE $MULIM $SYST $TOYS $FREQ $SAVE $LABEL
    # echo "---------------------------------------------------------------------------------------------------"
    # echo "combine -M AsymptoticLimits $DCARD --run blind $MULIM $SYST $FREQ $LABEL"
    #       combine -M AsymptoticLimits $DCARD --run blind $MULIM $SYST $FREQ $LABEL
    echo "---------------------------------------------------------------------------------------------------"
    echo "combine -M Significance $DCARD --uncapped=true -t -1 $MUSIG $SYST $FREQ $LABEL"
          combine -M Significance $DCARD --uncapped=true -t -1 $MUSIG $SYST $FREQ $LABEL
    echo "---------------------------------------------------------------------------------------------------"
    # ## deltaNLL = 1.0 <--> 68% CL, 3.85 <--> 95% CL
    # echo "text2workspace.py $DCARD"
    #       text2workspace.py $DCARD
    # echo "combine out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_template_MC.root -M MultiDimFit -t 100 --algo grid --points 601 --setParameterRanges r=-20.05,40.05"
    #       combine out_files/datacard/cat_e2mu_medLepMVA_noZ_noBtag_mass12_template_MC.root -M MultiDimFit -t 100 --algo grid --points 601 --setParameterRanges r=-20.05,40.05
    # echo "combine higgsCombine_noSyst.FitDiagnostics.mH120.123456.root -M MultiDimFit --algo grid --points 601 --setParameterRanges r=-20.05,20.05"
    # 	  combine higgsCombine_noSyst.FitDiagnostics.mH120.123456.root -M MultiDimFit --algo grid --points 601 --setParameterRanges r=-20.05,40.05
    echo "---------------------------------------------------------------------------------------------------"
done
echo "Returning to home directory and running 'cmsenv': $HOMEDIR"
cd $HOMEDIR
eval `scramv1 runtime -sh`
echo "Finished running run_test_limits.sh!"
