
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("tt", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/tt/180802_165355", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/MC_data_comparison_2017_v4_v2/files", {"0000/tuple_161.root", "0000/tuple_162.root", "0000/tuple_163.root", "0000/tuple_164.root", "0000/tuple_165.root", "0000/tuple_166.root", "0000/tuple_167.root", "0000/tuple_168.root", "0000/tuple_169.root", "0000/tuple_17.root", "0000/tuple_170.root", "0000/tuple_171.root", "0000/tuple_172.root", "0000/tuple_173.root", "0000/tuple_174.root", "0000/tuple_175.root", "0000/tuple_176.root", "0000/tuple_177.root", "0000/tuple_178.root", "0000/tuple_179.root", "0000/tuple_18.root", "0000/tuple_180.root", "0000/tuple_181.root", "0000/tuple_182.root", "0000/tuple_183.root", "0000/tuple_184.root"}, "210", -1, 1000, -0.000714)'