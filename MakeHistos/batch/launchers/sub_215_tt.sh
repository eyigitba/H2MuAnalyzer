
run_dir="/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos"
cd ${run_dir}
eval `scramv1 runtime -sh`
root -b -l -q '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/MakeHistos/macros/MC_data_comparison.C("tt", "/store/group/phys_higgs/HiggsExo/H2Mu/UF/ntuples/data_2017_and_mc_fall17/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/tt/180802_165355", "/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/MC_data_comparison_2017_v4_v2/files", {"0000/tuple_276.root", "0000/tuple_277.root", "0000/tuple_278.root", "0000/tuple_279.root", "0000/tuple_28.root", "0000/tuple_280.root", "0000/tuple_281.root", "0000/tuple_282.root", "0000/tuple_283.root", "0000/tuple_284.root", "0000/tuple_285.root", "0000/tuple_286.root", "0000/tuple_287.root", "0000/tuple_288.root", "0000/tuple_289.root", "0000/tuple_29.root", "0000/tuple_290.root", "0000/tuple_291.root", "0000/tuple_292.root", "0000/tuple_293.root", "0000/tuple_294.root", "0000/tuple_295.root", "0000/tuple_296.root", "0000/tuple_297.root"}, "215", -1, 1000, -0.000714)'