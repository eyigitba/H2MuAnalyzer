
pwd_cmd="/bin/pwd"
run_dir=`${pwd_cmd}`
in_dir_16="/afs/cern.ch/work/a/abrinke1/public/H2Mu/2016/Histograms/ttH_3l_AWB_2019_07_12_signal_v1/files/HADD/"
in_dir_17="/afs/cern.ch/work/a/abrinke1/public/H2Mu/2017/Histograms/ttH_3l_AWB_2019_07_12_signal_v1/files/HADD/"
in_dir_18="/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms/ttH_3l_AWB_2019_07_12_signal_v1/files/HADD/"
out_dir="/afs/cern.ch/work/a/abrinke1/public/H2Mu/Run2/Histograms/ttH_3l_AWB_2019_07_12_signal_v1/files/HADD/"
hadd_cmd="/cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw-patch/CMSSW_10_2_11_patch1/external/slc6_amd64_gcc700/bin/hadd -f"

echo $in_dir_16
echo $in_dir_17
echo $in_dir_18
echo $out_dir
echo $hadd_cmd

## Loop over files from 2018
for file18 in `ls $in_dir_18`; do
    ## Make sure this is a ROOT file
    if test "${file18#*_merged.root}" != "$file18"; then
	## File from 2016 should have the same name
	file16="${file18}"
	## File from 2017 should have the same name
	file17="${file18}"

	hadd_str="$hadd_cmd $out_dir$file18 $in_dir_16$file16 $in_dir_17$file17 $in_dir_18$file18"
	echo ""
	echo $hadd_str
	`$hadd_str`
	echo ""
    fi
done
