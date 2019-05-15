
## Functions to print analytic fits to data and MC

import os
import sys

import ROOT as R
import ROOT.RooFit as RF

sys.path.insert(1, '%s/../FitBackground/python' % os.getcwd())
import FitFunctions as FF  ## From FitBackground/python/FitFunctions.py


def DrawFits(sig_fit, bkg_fit, data_fit, cat, out_dir):

    #-------------------------------------------------------------------
    ## Plot data and fit into a frame

    c_sig  = R.TCanvas('c_%s_sig'  % cat, 'c_%s_sig'  % cat, 800, 600)
    c_bkg  = R.TCanvas('c_%s_bkg'  % cat, 'c_%s_bkg'  % cat, 800, 600)
    c_data = R.TCanvas('c_%s_data' % cat, 'c_%s_data' % cat, 800, 600)

    ## Signal
    c_sig.cd()
    fra_sig = sig_fit.var.frame()
    sig_fit.dat  .plotOn (fra_sig)
    sig_fit.model.plotOn (fra_sig, RF.LineColor(R.kBlue), RF.Range('FULL'))
    ## Plot sub-components of signal model
    if len(sig_fit.arg_sets) > 1:
        sig_fit.model.plotOn (fra_sig, RF.Components(sig_fit.arg_sets[0]), RF.LineStyle(R.kDashed), RF.LineColor(R.kGreen),  RF.Range('FULL'))
        sig_fit.model.plotOn (fra_sig, RF.Components(sig_fit.arg_sets[1]), RF.LineStyle(R.kDashed), RF.LineColor(R.kRed),    RF.Range('FULL'))
    if len(sig_fit.arg_sets) > 2:
        sig_fit.model.plotOn (fra_sig, RF.Components(sig_fit.arg_sets[2]), RF.LineStyle(R.kDashed), RF.LineColor(R.kViolet), RF.Range('FULL'))
    sig_fit.model.paramOn(fra_sig, RF.Layout(0.55, 0.90, 0.90))
    fra_sig.Draw()
    sig_chi = R.TLatex(0.7, 0.3, "#chi^{2} = %.2f" % fra_sig.chiSquare())
    sig_chi.SetNDC(R.kTRUE)
    sig_chi.Draw()
    c_sig.SaveAs(out_dir+'/plot/sig_fit_%s_%s_%d.png' % (cat, sig_fit.fit_type, sig_fit.order))

    ## Background
    c_bkg.cd()
    fra_bkg = bkg_fit.var.frame()
    bkg_fit.dat  .plotOn(fra_bkg)
    bkg_fit.model.plotOn(fra_bkg, RF.LineColor(R.kBlue), RF.Range('FULL'))
    ## Plot sub-components of background model
    if len(bkg_fit.arg_sets) > 1:
        bkg_fit.model.plotOn(fra_bkg, RF.Components(bkg_fit.arg_sets[0]), RF.LineStyle(R.kDashed), RF.LineColor(R.kGreen),  RF.Range('FULL'))
        bkg_fit.model.plotOn(fra_bkg, RF.Components(bkg_fit.arg_sets[1]), RF.LineStyle(R.kDashed), RF.LineColor(R.kRed),    RF.Range('FULL'))
    if len(bkg_fit.arg_sets) > 2:
        bkg_fit.model.plotOn(fra_bkg, RF.Components(bkg_fit.arg_sets[2]), RF.LineStyle(R.kDashed), RF.LineColor(R.kViolet), RF.Range('FULL'))
    bkg_fit.model.paramOn(fra_bkg, RF.Layout(0.55, 0.90, 0.90))
    fra_bkg.Draw()
    bkg_chi = R.TLatex(0.7, 0.5, "#chi^{2} = %.2f" % fra_bkg.chiSquare())
    bkg_chi.SetNDC(R.kTRUE)
    bkg_chi.Draw()
    c_bkg.SaveAs(out_dir+'/plot/bkg_fit_%s_%s_%d.png' % (cat, bkg_fit.fit_type, bkg_fit.order))

    ## Data
    c_data.cd()
    fra_data = data_fit.var.frame()
    data_fit.dat  .plotOn(fra_data)
    data_fit.model.plotOn(fra_data, RF.LineColor(R.kBlue), RF.Range('FULL'))
    ## Plot sub-components of data model
    if len(bkg_fit.arg_sets) > 1:
        data_fit.model.plotOn(fra_data, RF.Components(data_fit.arg_sets[0]), RF.LineStyle(R.kDashed), RF.LineColor(R.kGreen),  RF.Range('FULL'))
        data_fit.model.plotOn(fra_data, RF.Components(data_fit.arg_sets[1]), RF.LineStyle(R.kDashed), RF.LineColor(R.kRed),    RF.Range('FULL'))
    if len(data_fit.arg_sets) > 2:
        data_fit.model.plotOn(fra_data, RF.Components(data_fit.arg_sets[2]), RF.LineStyle(R.kDashed), RF.LineColor(R.kViolet), RF.Range('FULL'))
    data_fit.model.paramOn(fra_data, RF.Layout(0.55, 0.90, 0.90))
    fra_data.Draw()
    data_chi = R.TLatex(0.7, 0.5, "#chi^{2} = %.2f" % fra_data.chiSquare())
    data_chi.SetNDC(R.kTRUE)
    data_chi.Draw()
    c_data.SaveAs(out_dir+'/plot/data_fit_%s_%s_%d.png' % (cat, data_fit.fit_type, data_fit.order))

## End function: def DrawFits(sig_fit, bkg_fit, data_fit, cat):

