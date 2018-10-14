
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("SingleMu_2017E", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/SingleMuon/SingleMu_2017E/180802_164036", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/WH_cat_2017_v4_v4/files", {"0000/tuple_8.root", "0000/tuple_80.root"}, "252", -1, 1000, 0.000000)'