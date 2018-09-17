
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("tt", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/tt/180802_165355", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/MC_data_comparison_2017_v4_v2/files", {"0000/tuple_34.root", "0000/tuple_340.root", "0000/tuple_341.root", "0000/tuple_342.root", "0000/tuple_343.root", "0000/tuple_344.root", "0000/tuple_345.root", "0000/tuple_346.root", "0000/tuple_347.root", "0000/tuple_348.root", "0000/tuple_349.root", "0000/tuple_35.root", "0000/tuple_350.root", "0000/tuple_351.root", "0000/tuple_352.root", "0000/tuple_353.root", "0000/tuple_354.root", "0000/tuple_355.root", "0000/tuple_356.root", "0000/tuple_357.root", "0000/tuple_358.root", "0000/tuple_359.root", "0000/tuple_36.root"}, "218", -1, 1000, -0.000714)'