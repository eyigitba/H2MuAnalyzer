
## Current directory
HOMEDIR=$PWD
## CMSSW release containing Higgs Combine
COMBINEDIR="/afs/cern.ch/work/x/xzuo/combine/CMSSW_8_1_0/src"
## Directory containing workspace directory
WSDIR="$HOMEDIR/hig-17-019/paper/13TeV"
## Datacard name
# DCARD="cms_datacard_hmumu_cat14_fixed_v2.txt"
DCARD="cms_datacard_hmumu_cat14_fixed_v2_cat12_only.txt"

cd $COMBINEDIR
echo "Running 'cmsenv' in CMSSW release with Higgs Combine: $COMBINEDIR"
eval `scramv1 runtime -sh`
echo "Navigating to workspace directory: $WSDIR"
cd $WSDIR
echo "Running combine -M AsymptoticLimits $DCARD --run blind"
combine -M AsymptoticLimits $DCARD --run blind
echo "Returning to home directory and running 'cmsenv': $HOMEDIR"
cd $HOMEDIR
eval `scramv1 runtime -sh`
echo "Finished running run_2016_limits.sh!"

