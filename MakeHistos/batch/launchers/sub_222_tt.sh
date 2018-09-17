
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("tt", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/tt/180802_165355", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/MC_data_comparison_2017_v4_v2/files", {"0000/tuple_425.root", "0000/tuple_426.root", "0000/tuple_427.root", "0000/tuple_428.root", "0000/tuple_429.root", "0000/tuple_43.root", "0000/tuple_430.root", "0000/tuple_431.root", "0000/tuple_432.root", "0000/tuple_433.root", "0000/tuple_434.root", "0000/tuple_435.root", "0000/tuple_436.root", "0000/tuple_437.root", "0000/tuple_438.root", "0000/tuple_439.root", "0000/tuple_44.root", "0000/tuple_440.root", "0000/tuple_441.root", "0000/tuple_442.root", "0000/tuple_443.root", "0000/tuple_444.root", "0000/tuple_445.root"}, "222", -1, 1000, -0.000714)'