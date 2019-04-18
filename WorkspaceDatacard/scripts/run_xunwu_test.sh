
## Current directory
HOMEDIR=$PWD
## CMSSW release containing Higgs Combine
COMBINEDIR="/afs/cern.ch/work/x/xzuo/combine/CMSSW_8_1_0/src"
## Directory containing workspace directory
WSDIR="$HOMEDIR"
## Datacard name
declare -a DCARDS=(
    "out_files/datacard/WH_mu_444_template_MC.txt"
    "out_files/datacard/WH_mu_884_template_MC.txt"
    "out_files/datacard/WH_mu_888_template_MC.txt"
    "out_files/datacard/WH_mu_no_template_MC.txt"

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
    echo "Running combine -M AsymptoticLimits $DCARD --run blind"
    combine -M AsymptoticLimits $DCARD --run blind
    echo "---------------------------------------------"
    echo "Running combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=1"
    combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=1
    echo "---------------------------------------------"
    echo "Running combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=1 --toysFreq"
    combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=1 --toysFreq
    echo "---------------------------------------------"
    echo "Running combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=0"
    combine -M Significance --uncapped=true $DCARD -t -1 --expectSignal=0
#    echo "---------------------------------------------"
#    echo "combine -M MaxLikelihoodFit --uncapped=true $DCARD"
#    combine -M MaxLikelihoodFit --uncapped=true $DCARD
    echo "---------------------------------------------"
    echo "combine -M FitDiagnostics $DCARD"
    combine -M FitDiagnostics $DCARD
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
