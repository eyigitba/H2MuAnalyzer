####################################################
###    make stack and ratio from miniNTuple      ###
####################################################

## Basic python includes for manipulating files
import os
import sys

sys.path.insert(0, '%s/lib' % os.getcwd() )
from ROOT import *
from MNT_Helper import LoadColors, LinearStack, RatioPlot, FillHistTerm, GetSF
#R.gROOT.SetBatch(True)

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/h2mm_944/src/H2MuAnalyzer/TrainMVA/output'

BDT_WITH_MASS    = '2017_WH_ele_high_pt_dimu_lepMVA04_train_with_mass.root'
BDT_WITHOUT_MASS = '2017_WH_ele_high_pt_dimu_lepMVA04_train_without_mass.root'
DO_BDT_CATS = True

if USER == 'xzuo':     
    OUT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'
    OUT_DIR = OUT_DIR + '/VH_selection_2019april/pt10_iso04/WH_ele_high_dimu_pt/plots'
    MASS_FILE = 'mass_hists_cut_444.root'

def BDT_Template_Limit():
    out_name = "BDT_lepMVA04_with_mass.root"
    out_file = TFile( OUT_DIR + "/" + out_name , "RECREATE")

    in_file_BDT  = TFile.Open( PLOT_DIR + "/" + BDT_WITH_MASS, "READ")
    in_file_norm = TFile.Open( OUT_DIR + "/" + MASS_FILE, "READ")

    bdt_name = "2017_WH_ele_with_mass_all_sig_all_bkg_ge0j"
    method_name = "BDTG_UF_v2"
    sig_hist_name = bdt_name + "/" + "Method_" + method_name + "/" + method_name + "/" + "MVA_" + method_name + "_S"
    bkg_hist_name = bdt_name + "/" + "Method_" + method_name + "/" + method_name + "/" + "MVA_" + method_name + "_B"
    BDT_sig = in_file_BDT.Get(sig_hist_name).Clone()
    BDT_sig.SetNameTitle("BDT_Net_Sig", "BDT_Net_Sig")
    BDT_bkg = in_file_BDT.Get(bkg_hist_name).Clone()
    BDT_bkg.SetNameTitle("BDT_Net_Bkg", "BDT_Net_Bkg")

    BDT_data = TH1D("BDT_Net_Data", "BDT_Net_Data", 100,-1,1) #placeholder empty data plot

    sig_norm = in_file_norm.Get("dimu_mass_Net_Sig").Integral()
    bkg_norm = in_file_norm.Get("dimu_mass_Net_Bkg").Integral()
    print "signal norm = %f" %sig_norm
    print "bkg    norm = %f" %bkg_norm

    print "signal BDT Integral = %f" %BDT_sig.Integral()
    print "bkg    BDT Integral = %f" %BDT_bkg.Integral()
    BDT_sig.Scale( 1/BDT_sig.Integral() )
    BDT_sig.Scale( sig_norm )
    BDT_bkg.Scale( 1/BDT_bkg.Integral() )
    BDT_bkg.Scale( bkg_norm )


    out_file.cd()
    BDT_sig.Write()
    BDT_bkg.Write()
    BDT_data.Write()
    out_file.Close()


def BDT_Cut_Mass_Template():
    out_name_neg = "mass_BDT_neg_lepMVA04.root"
    out_name_pos = "mass_BDT_pos_lepMVA04.root"
    out_file_neg = TFile( OUT_DIR + "/" + out_name_neg , "RECREATE")
    out_file_pos = TFile( OUT_DIR + "/" + out_name_pos , "RECREATE")

    in_file_BDT  = TFile.Open( PLOT_DIR + "/" + BDT_WITHOUT_MASS, "READ")

    TestTree = in_file_BDT.Get("2017_WH_ele_without_mass_all_sig_all_bkg_ge0j/TestTree")
    TestTree.Draw("dimu_mass>>bkg_neg(12,100,160)", "2 * event_wgt * xsec_norm * ( BDTG_UF_v2 < 0 && Sample_ID < 0)")
    bkg_neg = gDirectory.Get("bkg_neg")
    bkg_neg.SetName("dimu_mass_Net_Bkg")
    bkg_pos = TestTree.Draw("dimu_mass>>bkg_pos(12,100,160)", "2 * event_wgt * xsec_norm * ( BDTG_UF_v2 >= 0 && Sample_ID < 0)")
    bkg_pos = gDirectory.Get("bkg_pos")
    bkg_pos.SetName("dimu_mass_Net_Bkg")
    sig_neg = TestTree.Draw("dimu_mass>>sig_neg(12,100,160)", "2 * event_wgt * xsec_norm * ( BDTG_UF_v2 < 0 && Sample_ID > 0)")
    sig_neg = gDirectory.Get("sig_neg")
    sig_neg.SetName("dimu_mass_Net_Sig")
    sig_pos = TestTree.Draw("dimu_mass>>sig_pos(12,100,160)", "2 * event_wgt * xsec_norm * ( BDTG_UF_v2 >= 0 && Sample_ID > 0)")
    sig_pos = gDirectory.Get("sig_pos")
    sig_pos.SetName("dimu_mass_Net_Sig")

    data_neg = bkg_neg.Clone()   # place holder
    data_neg.SetName("dimu_mass_Net_Data")
    data_pos = bkg_pos.Clone()
    data_pos.SetName("dimu_mass_Net_Data")

    out_file_neg.cd()
    bkg_neg.Write()
    sig_neg.Write()
    data_neg.Write()
    out_file_neg.Close()

    out_file_pos.cd()
    bkg_pos.Write()
    sig_pos.Write()
    data_pos.Write()
    out_file_pos.Close()

def main():
    BDT_Template_Limit()
    if DO_BDT_CATS: BDT_Cut_Mass_Template()
main()











	
