
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("tt", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/tt/180802_165355", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/MC_data_comparison_2017_v4_v2/files", {"0000/tuple_6.root", "0000/tuple_60.root", "0000/tuple_61.root", "0000/tuple_62.root", "0000/tuple_63.root", "0000/tuple_64.root", "0000/tuple_65.root", "0000/tuple_66.root", "0000/tuple_67.root", "0000/tuple_68.root", "0000/tuple_69.root", "0000/tuple_7.root", "0000/tuple_70.root", "0000/tuple_71.root", "0000/tuple_72.root", "0000/tuple_73.root", "0000/tuple_74.root", "0000/tuple_75.root", "0000/tuple_76.root", "0000/tuple_77.root", "0000/tuple_78.root", "0000/tuple_79.root", "0000/tuple_8.root"}, "229", -1, 1000, -0.000714)'