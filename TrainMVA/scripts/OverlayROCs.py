#! /usr/bin/env python 

#============================================
# Imports
#============================================

import array
import ROOT as R

XWZ_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2017/Histograms/VH_selection_2019april/pt10_iso04/WH_ele_with_BDT/plots/'
AWB_DIR = '/afs/cern.ch/user/a/abrinke1/HiggsToMuMu/2017/CMSSW_9_4_10/src/H2MuAnalyzer/TrainMVA/output/'

XWZ_1 = '2017_WH_ele_against_inclu_trimvar_with_mass_all_sig_all_bkg_ge0j'
AWB_1 = 'f_Opt_withMassBDT_lepMVA_v1_all_sig_all_bkg_ge0j'
AWB_2 = 'f_Opt_noMassBDT_mass_v1_all_sig_all_bkg_ge0j'

AWB_2 = 'f_Opt_AWB_noMass_v3_all_sig_all_bkg_'
AWB_3 = 'f_Opt_AWB_noMass_v3_resWgt_all_sig_all_bkg_resWgt'


def OverlayROCs():

    print '\nInside OverlayROCs.py'

    ## Don't print plots to screen while running (faster)
    R.gROOT.SetBatch(R.kTRUE)

    inputs = {}
    # inputs['XWZ_e2mu_BDT_withMass'] = [XWZ_DIR+'2017_WH_ele_against_inclu_lepMVA04.root',
    #                                    XWZ_1, 'BDTG_UF_v2', R.kBlack]
    # inputs['AWB_e2mu_LD_withMassBDT_lepMVA'] = [AWB_DIR+'TMVA_retrain_WH_lep_e2mu_2019_05_15.root',
    #                                             AWB_1, 'LD', R.kViolet]
    # inputs['AWB_e2mu_BDT_withMassBDT_lepMVA'] = [AWB_DIR+'TMVA_retrain_WH_lep_e2mu_2019_05_15.root',
    #                                              AWB_1, 'BDTG_UF_v1', R.kBlue]
    # inputs['AWB_e2mu_LD_noMassBDT_mass'] = [AWB_DIR+'TMVA_retrain_WH_lep_e2mu_2019_05_15.root',
    #                                         AWB_2, 'LD', R.kGreen]
    # inputs['AWB_e2mu_BDT_noMassBDT_mass'] = [AWB_DIR+'TMVA_retrain_WH_lep_e2mu_2019_05_15.root',
    #                                          AWB_2, 'BDTG_UF_v1', R.kRed]

    inputs['AWB_lep_noMassBDT_UF_v1'] = [AWB_DIR+'TMVA_BDT_2017_WH_lep_all_vs_all_v2.root',
                                         AWB_2, 'BDTG_UF_v1', R.kBlack]
    inputs['AWB_lep_noMassBDT_resWgt_UF_v1'] = [AWB_DIR+'TMVA_BDT_2017_WH_lep_all_vs_all_v2.root',
                                                AWB_2, 'BDTG_UF_v1', R.kBlue]
    inputs['AWB_lep_noMassBDT_resWgt_UF_v2'] = [AWB_DIR+'TMVA_BDT_2017_WH_lep_all_vs_all_v2.root',
                                                AWB_2, 'BDTG_UF_v2', R.kViolet]
    inputs['AWB_lep_noMassBDT_resWgt_UF_v3'] = [AWB_DIR+'TMVA_BDT_2017_WH_lep_all_vs_all_v2.root',
                                                AWB_2, 'BDTG_UF_v3', R.kRed]

    
    dists = ['trainingRejBvsS',
             'trainingEffBvsS',
             'B', 'Train_B']
             # 'S', 'Train_S']

    ## Output file
    out_file = R.TFile('plots/OverlayROCs.root', 'recreate')
    out_file.cd()
    cans  = {}
    hists = {}

    ## Loop over distributions
    for dist in dists:
        out_file.cd()
        cans[dist] = R.TCanvas('c_%s' % dist, 'c_%s' % dist, 800, 600)
        hists[dist] = {}

        ## Loop over input methods
        for key in inputs.keys():
            hist_dir = '/Method_{0}/{0}/MVA_{0}_'.format(inputs[key][2])
            print '\n  * Opening file %s' % inputs[key][0]
            print   '  * Getting histogram %s' % (inputs[key][1]+hist_dir+dist)
            in_file = R.TFile(inputs[key][0])

            out_file.cd()
            hists[dist][key] = in_file.Get(inputs[key][1]+hist_dir+dist).Clone(dist+'_'+key)
            hists[dist][key].SetLineWidth(2)
            hists[dist][key].SetLineColor(inputs[key][3])

            ## If input is BDT histogram, replace with cumulative integral version
            if dist == 'B' or dist == 'S' or dist == 'Train_B' or dist == 'Train_S':
                ## hists[dist][key] = hists[dist][key].GetCumulative(False)  ## Doesn't handle error bars properly
                nBins = hists[dist][key].GetNbinsX()
                for i in range(1, nBins+1):
                    err = R.Double(-99)
                    cum = hists[dist][key].IntegralAndError(i, nBins, err)
                    hists[dist][key].SetBinContent(i, cum)
                    hists[dist][key].SetBinError(i, err)

            ## Draw BDT output from training as ratio to testing
            if dist == 'Train_B' or dist == 'Train_S':
                hists[dist][key].Divide(hists[dist[-1]][key])
                can = R.TCanvas('c_%s_%s' % (dist, key), 'c_%s_%s' % (dist, key), 800, 600)
                hists[dist][key].GetYaxis().SetRangeUser(0, 2)
                hists[dist][key].Draw('histe')
                leg = R.TLegend(0.11, 0.11, 0.45, 0.35)
                leg.AddEntry(hists[dist][key], key)
                leg.Draw('same')
                can.Write()
                can.SaveAs( 'plots/%s.png' % can.GetName() )

        ## End loop: for key in inputs.keys():
            
        ## Don't draw BDT output histogram overlays
        if dist == 'B' or dist == 'S' or dist == 'Train_B' or dist == 'Train_S': continue

        cans[dist].cd()
        if (dist == 'trainingRejBvsS'): leg = R.TLegend(0.11, 0.11, 0.45, 0.55)
        else:                           leg = R.TLegend(0.11, 0.55, 0.45, 0.89)
        for key in hists[dist].keys():
            hists[dist][key].Draw('same')
            leg.AddEntry(hists[dist][key], key)
        leg.Draw('same')
        cans[dist].Write()
        cans[dist].SaveAs( 'plots/%s.png' % cans[dist].GetName() )
        cans[dist].SetLogy()
        cans[dist].SaveAs( 'plots/%s_log.png' % cans[dist].GetName() )

    ## End loop: for dist in dists:

    out_file.Close()

    print '\nDone!'

## End function: OverlayROCs():


def main():

    OverlayROCs()

main()  ## End of main function


