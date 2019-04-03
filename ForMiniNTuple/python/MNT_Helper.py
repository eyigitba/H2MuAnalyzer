##################################################
#####  functions for making plots from MNT   #####
##################################################

from ROOT import *

def LinearStack( name, all_stack, scaled_signal, h_data, legend, plot_dir):
    canv = TCanvas("Stack_" + name, "Stack_" + name, 600,600)
    canv.Clear()

    canv.cd()
    all_stack.Draw("HIST")
    scaled_signal.Draw("HISTSAME")
    h_data.Draw("SAMEPE")
    legend.Draw()
    canv.Update()
    canv.Write()
#    canv.SaveAs( plot_dir + "/linear_stack_" + name + ".png")




def RatioPlot( name, all_stack, scaled_signal, h_data, ratio_graph, legend, plot_dir, log_plot = False ):   # Do not use TRatioPlot! It is a devil!    -XWZ 19.09.2018
    canv = TCanvas("ratio_"+name, "ratio_"+name, 600,600)
    canv.Clear()

    upper_pad = TPad("upperpad_"+name, "upperpad_"+name, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.Draw()
    upper_pad.cd()
    if log_plot:
        upper_pad.SetLogy()
#       all_stack.SetMinimum(1e-5)
#       all_stack.SetMaximum(1e6)
    	all_stack.SetMinimum(1e-3)
    	all_stack.SetMaximum(1e8)
    all_stack.Draw("HIST")
    scaled_signal.Draw("HISTSAME")
    h_data.SetMarkerStyle(20)
    h_data.Draw("SAMEPE")
    legend.Draw()

    canv.cd()
    lower_pad = TPad("lowerpad_"+name, "lowerpad_"+name, 0,0.05, 1,0.2)
    lower_pad.SetTopMargin(0.05)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
    ratio_graph.SetMinimum(0.5)
    ratio_graph.SetMaximum(1.5)
#    ratio_graph.SetMinimum(0.5)
#    ratio_graph.SetMaximum(2.5)

    ratio_graph.GetXaxis().SetRangeUser( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_graph.SetMarkerStyle(20)
#    ratio_graph.GetYaxis().SetNdivisions(510)
    ratio_graph.GetYaxis().SetNdivisions(505)
    ratio_graph.GetXaxis().SetLabelSize(0.15)
    ratio_graph.GetYaxis().SetLabelSize(0.13)
    ratio_graph.GetYaxis().SetTitle("data/MC")
    ratio_graph.GetYaxis().SetTitleSize(0.20)
    ratio_graph.GetYaxis().SetTitleOffset(0.2)
    ratio_graph.Draw()

    canv.Update()
    canv.Write()
#    canv.SaveAs( plot_dir + "/linear_stack_" + name + ".png")


def PassLepMVA(lepMVA_cut, lep_value):
    Pass = False
    cut_value = float(lepMVA_cut)
    if lep_value > cut_value:
        Pass = True
    return Pass


def FillHistTerm(histos, term, signals, bkgs, value, Sample_ID, event_wgt):
    #blind
    if Sample_ID == 0 and term == "dimu_mass" and value > 120 and value < 130:
        return
    #data
    if Sample_ID == 0:
        histos[term]["data"].Fill(value, 1)
    #signal
    elif Sample_ID == 25:
        histos[term]["ggH"].Fill(value, event_wgt )
    elif Sample_ID == 010225:
        histos[term]["VBF"].Fill(value, event_wgt )
    elif Sample_ID == 2425:
        histos[term]["WH"].Fill(value, event_wgt )
    elif Sample_ID == 2325:
        histos[term]["ZH"].Fill(value, event_wgt )
    elif Sample_ID == 060625:
        histos[term]["ttH"].Fill(value, event_wgt )

    #background
    elif Sample_ID == -23:
        histos[term]["DY"].Fill(value, event_wgt / 3)  # used two DY samples
    elif Sample_ID == -0606:
        histos[term]["ttbar"].Fill(value, event_wgt / 2)  # used two ttbar samples
    elif Sample_ID == -2423:
        histos[term]["WZ"].Fill(value, event_wgt)
    elif Sample_ID == -2323:  #ZZ_4l and ZZ_2l
        histos[term]["ZZ"].Fill(value, event_wgt)
    elif Sample_ID == -2424: # WW
	histos[term]["WW"].Fill(value, event_wgt)
    elif Sample_ID/10000 < -23:   # triboson
        histos[term]["triboson"].Fill(value, event_wgt)
    elif Sample_ID == -062300:
        histos[term]["tZq"].Fill(value, event_wgt)
    elif Sample_ID == -060623:
        histos[term]["ttZ"].Fill(value, event_wgt)
    else:
        histos[term]["tX+ttX"].Fill(value, event_wgt)
