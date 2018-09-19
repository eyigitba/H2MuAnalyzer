
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("SingleMu_2017C", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/SingleMuon/SingleMu_2017C/180802_163916", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/VH_toy_2017_v4_v2/files", {"0000/tuple_187.root", "0000/tuple_188.root", "0000/tuple_189.root", "0000/tuple_19.root", "0000/tuple_190.root", "0000/tuple_191.root"}, "38", -1, 1000, 0.000000)'