
pwd_cmd="/bin/pwd"
run_dir=`${pwd_cmd}`
in_dir_17="/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/WH_lep_AWB_2019_06_24_v1/files/HADD/"
in_dir_18="/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms/WH_lep_AWB_2019_06_24_v1/files/HADD/"
out_dir="/afs/cern.ch/work/a/abrinke1/public/H2Mu/Run2/Histograms/WH_lep_AWB_2019_06_24_v1/files/HADD/"
hadd_cmd="/cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw-patch/CMSSW_10_2_11_patch1/external/slc6_amd64_gcc700/bin/hadd -f"

echo $in_dir_17
echo $in_dir_18
echo $out_dir
echo $hadd_cmd

## Loop over files from 2017
for file17 in `ls $in_dir_17`; do
    ## Make sure this is a ROOT file
    if test "${file17#*root}" != "$file17"; then
        ## File from 2018 should have the same name, modulo a change from Presel2017 --> PreselRun2
	file18="${file17/Presel2017/PreselRun2}"
	## For some reason, some names got corrupted and have "Presel2017e2mu" instead of "Presel2017_e2mu"
	file18="${file18/PreselRun2e2mu/PreselRun2_e2mu}"
	hadd_str="$hadd_cmd $out_dir$file18 $in_dir_17$file17 $in_dir_18$file18"
	echo ""
	echo $hadd_str
	`$hadd_str`
	echo ""
    fi
done
