
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("SingleMu_2017C", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/SingleMuon/SingleMu_2017C/180802_163916", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/VH_toy_2017_v4_v2/files", {"0000/tuple_144.root", "0000/tuple_145.root", "0000/tuple_146.root", "0000/tuple_147.root", "0000/tuple_148.root"}, "30", -1, 1000, 0.000000)'