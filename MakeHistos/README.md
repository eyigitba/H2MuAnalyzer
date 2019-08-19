### updates on systematic uncertainties
#XWZ 2019.08.19

general instructions:
To prepare the histograms in a channel with systematics, one needs to run the same macro twice.
The first time run the macro on all samples with SYS_SHIFTS = [] in the batch, to get the StackPlot for signal, bkg, data
The second time run the macro on signals only with SYS_SHIFTS specified in the batch script, to get the up/down shifted signal histograms.

A detailed description of the changes are listed below

- histo making part 
JES, PU_wgt, IsoMu_SF, LepMVA_SF are added to the histo making library.

JES is achieved in src/LoadNTupleBranches.cc by writing up/down shifted objects/values to the jet related branches.
PU and IsoMu are achieved in src/EventWeight.cc by using up/down shifted weights instead of non-shifted ones.
lepMVA is achieved in src/ObjectHelper.cc by loading the up/down shifted SF instead of the non-shifted ones.

MuID and MuIso systematics are also available, in case lepMVA is not in use.

All these changes are compatible with older scripts. The only change needed for old macros is to add 'std::string SYS = ""' in the macro arguments after 'TString hist_tree = ""'


- batch part
Added hard-coded configuration SYS_SHIFTS, which is to be taken as an argument in the macros.
Modified names of sub_files and hadd_files to run different systematics setups on the same sample


- bug fix
1. for 2017 and 2018, use Roch as default pt, instead of KaMu
2. In some macros, if lepMVA in use, set evt_wgt.muon_ID  = false  and  evt_wgt.muon_Iso = false
3. In some macros, when getting lepMVA SF for muons, use the absolute value of muon eta. (reference muon SF has eta range 0~2.4)
