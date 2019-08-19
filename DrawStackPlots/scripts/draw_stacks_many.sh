
MAX_PLOT=221  ## Maximum number of stack plots to create
NUM_PLOT=10   ## Number of stack plots to create per job

# for category in "3lep_looseLepMVA_noZ5_noBtag" "3lep_medLepMVA_noZ10_noBtag" "3lep_hiPt_lepW20_medLepMVA_noZ10_noBtag"; do
#     for year in "Run2" "2016" "2017" "2018"; do
# 	for i in $(eval echo "{0..$MAX_PLOT..$NUM_PLOT}"); do
# 	    echo "Running macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category WH_lep_AWB_2019_07_08_signal_v1 WH_lep"
# 	    `macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category WH_lep_AWB_2019_07_08_signal_v1 WH_lep`
# 	done
#     done
# done

for category in "3lep_looseLepMVA_ge2j_btag" "3lep_medLepMVA_noZ10_ge2j_btag" "3lep_hiPt_lepW20_medLepMVA_noZ10_ge2j_btag" "3lep_hiPt_lep20_tightLepMVA_noZ10_btag"; do
    for year in "Run2" "2016" "2017" "2018"; do
	for i in $(eval echo "{0..$MAX_PLOT..$NUM_PLOT}"); do
	    echo "Running macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category ttH_3l_AWB_2019_07_12_signal_v1 ttH_3l"
	    `macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category ttH_3l_AWB_2019_07_12_signal_v1 ttH_3l`
	done
    done
done

# for category in "3lep_allMass_medLepMVA_onZ10_noBtag" "3mu_allMass_medLepMVA_onZ10_noBtag" "e2mu_allMass_medLepMVA_onZ10_noBtag" "3lep_allMass_medLepMVA_noZ10_noBtag" "3mu_allMass_medLepMVA_noZ10_noBtag" "e2mu_allMass_medLepMVA_noZ10_noBtag"; do
# # for category in "3lep_allMass_medLepMVA_onZ10_noBtag" "3lep_allMass_medLepMVA_noZ10_noBtag"; do
#     for year in "Run2" "2016" "2017" "2018"; do
# 	for i in $(eval echo "{0..$MAX_PLOT..$NUM_PLOT}"); do
# 	    echo "Running macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category WH_lep_AWB_2019_07_08_sideband_v1 WH_lep_allMass"
# 	    `macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category WH_lep_AWB_2019_07_08_sideband_v1 WH_lep_allMass`
# 	done
#     done
# done

# for category in "3lep_allMass_medLepMVA_onZ10_ge2j_btag" "3mu_allMass_medLepMVA_onZ10_ge2j_btag" "e2mu_allMass_medLepMVA_onZ10_ge2j_btag" "3lep_allMass_medLepMVA_noZ10_ge2j_btag" "3mu_allMass_medLepMVA_noZ10_ge2j_btag" "e2mu_allMass_medLepMVA_noZ10_ge2j_btag"; do
# # for category in "3lep_allMass_medLepMVA_onZ10_ge2j_btag" "3lep_allMass_medLepMVA_noZ10_ge2j_btag"; do
#     for year in "Run2" "2016" "2017" "2018"; do
# 	for i in $(eval echo "{0..$MAX_PLOT..$NUM_PLOT}"); do
# 	    echo "Running macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category ttH_3l_AWB_2019_07_08_sideband_v1 ttH_3l_allMass"
# 	    `macros/StackPlots.py $i $((i+NUM_PLOT-1)) $year $category ttH_3l_AWB_2019_07_08_sideband_v1 ttH_3l_allMass`
# 	done
#     done
# done

