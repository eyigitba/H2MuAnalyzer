
## Current directory
HOMEDIR=$PWD
## CMSSW release containing Higgs Combine
COMBINEDIR="/afs/cern.ch/user/a/abrinke1/HiggsToMuMu/Combine/CMSSW_8_1_0/src"
## Directory containing workspace directory
WSDIR="$HOMEDIR"
## Datacard name
# declare -a DCARDS=("out_files/datacard/c14_data.txt")
# declare -a DCARDS=(
#     "out_files/datacard/c0_data.txt"
#     "out_files/datacard/c1_data.txt"
#     "out_files/datacard/c2_data.txt"
#     "out_files/datacard/c3_data.txt"
#     "out_files/datacard/c4_data.txt"
#     "out_files/datacard/c5_data.txt"
#     "out_files/datacard/c6_data.txt"
#     "out_files/datacard/c7_data.txt"
#     "out_files/datacard/c8_data.txt"
#     "out_files/datacard/c9_data.txt"
#     "out_files/datacard/c10_data.txt"
#     "out_files/datacard/c11_data.txt"
#     "out_files/datacard/c12_data.txt"
#     "out_files/datacard/c13_data.txt"
#     "out_files/datacard/c14_data.txt"
# )
# declare -a DCARDS=(
#     # "out_files/datacard/cat_e2mu_looseLepMVA_noZ_ge3j_btag_mass12_data.txt"
#     # "out_files/datacard/cat_3mu_looseLepMVA_noZ_ge3j_btag_mass12_data.txt"
#     "out_files/datacard/cat_e2mu_looseLepMVA_noZ_ge3j_btag_mass12_MC.txt"
#     "out_files/datacard/cat_3mu_looseLepMVA_noZ_ge3j_btag_mass12_MC.txt"
# )
declare -a DCARDS=(
    # "out_files/datacard/cat_e2mu_looseLepMVA_mt150_noBtag_noZ_mass12_data.txt"
    # "out_files/datacard/cat_3mu_looseLepMVA_mt150_noBtag_noZ_mass12_data.txt"
    "out_files/datacard/cat_e2mu_looseLepMVA_mt150_noBtag_noZ_mass12_MC.txt"
    "out_files/datacard/cat_3mu_looseLepMVA_mt150_noBtag_noZ_mass12_MC.txt"
)

cd $COMBINEDIR
echo "Running 'cmsenv' in CMSSW release with Higgs Combine: $COMBINEDIR"
eval `scramv1 runtime -sh`
echo "Navigating to workspace directory: $WSDIR"
cd $WSDIR
for DCARD in "${DCARDS[@]}"
do
    echo "Running combine -M AsymptoticLimits $DCARD --run blind -S 0 --expectSignal 0"
    combine -M AsymptoticLimits $DCARD --run blind -S 0 --expectSignal 0
    echo "Running combine -M AsymptoticLimits $DCARD --run blind --expectSignal 0"
    combine -M AsymptoticLimits $DCARD --run blind --expectSignal 0
    # echo "Running combine -M MultiDimFit $DCARD -S 0 --preFitValue 0 --do95 1"
    # combine -M MultiDimFit $DCARD -S 0 --preFitValue 0 --do95 1
    # echo "Running combine -M FitDiagnostics $DCARD -S 0 --preFitValue 0"
    # combine -M FitDiagnostics $DCARD -S 0 --preFitValue 0
done
echo "Returning to home directory and running 'cmsenv': $HOMEDIR"
cd $HOMEDIR
eval `scramv1 runtime -sh`
echo "Finished running run_test_limits.sh!"
