####################################################
###    make stack and ratio from added histos    ###
####################################################

## Basic python includes for manipulating files
import os

from ROOT import *

## Configure the script user
if 'abrinke1' in os.getcwd(): USER = 'abrinke1'
if 'bortigno' in os.getcwd(): USER = 'bortigno'
if 'xzuo'     in os.getcwd(): USER = 'xzuo'

## Directory for input histograms and output plots
if USER == 'abrinke1': PLOT_DIR = '/afs/cern.ch/work/a/abrinke1/public/H2Mu/2018/Histograms'
if USER == 'xzuo':     PLOT_DIR = '/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms'

LABEL = 'MC_data_comparison_2017_v4_v7'  ## Sub-folder within PLOT_DIR containing histograms


def ratioplot( term, sig_stack, all_stack, h_data, ratio_graph, legend ):   # Do not use TRatioPlot! It is a devil!    -XWZ 19.09.2018
    canv = TCanvas("ratio_"+term, "ratio_"+term, 1)
    canv.Clear()
    
    upper_pad = TPad("upperpad_"+term, "upperpad_"+term, 0,0.3, 1,1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.Draw()
    upper_pad.cd()
    upper_pad.SetLogy()
    all_stack.SetMinimum(1e-4)
    all_stack.Draw("HIST")
    h_data.SetMarkerStyle(20)
    h_data.Draw("SAME")
    for histo in sig_stack:
	histo.Draw("SAMEHIST")
    legend.Draw()

    canv.cd()
    lower_pad = TPad("lowerpad_"+term, "lowerpad_"+term, 0,0.05, 1,0.3)
    lower_pad.SetTopMargin(0.05)
    lower_pad.Draw()
    lower_pad.cd()
    ratio_graph.SetMinimum(0.5)
    ratio_graph.SetMaximum(1.3)
    ratio_graph.GetXaxis().SetRangeUser( h_data.GetXaxis().GetXmin(), h_data.GetXaxis().GetXmax() )
    ratio_graph.SetMarkerStyle(20)
    ratio_graph.Draw()

    canv.Update()
    canv.Write()
    canv.SaveAs(PLOT_DIR+"/"+LABEL+"/files/sum" + "/plots" + "/ratio_" + term + ".png")


def main():

    file_dir=PLOT_DIR+"/"+LABEL+"/files/sum"
    out_file = TFile( file_dir + "/stack" + ".root", "RECREATE")
    in_file = TFile.Open( file_dir + "/all.root", "READ")

    terms = ["dimuon_mass", "dimuon_pt","dimuon_eta", "dimuon_delta_eta", "dimuon_delta_phi",
	     "leading_muon_pt", "leading_muon_eta", "subleading_muon_pt", "subleading_muon_eta",
	     "dijet_mass", "dijet_pt", "dijet_eta", "dijet_delta_eta", "dijet_delta_phi",
	     "leading_jet_pt", "leading_jet_eta", "subleading_jet_pt", "subleading_jet_eta",
	     "MET", "nJets", "nBjets", "nVertices", "nElectrons",  "nMuons"]
    signals = ["H2Mu_WH_neg", "H2Mu_WH_pos", "H2Mu_ZH", "H2Mu_VBF", "H2Mu_gg"]
    bkgs = ["tt", "ZJets_AMC"]
    data = ["SingleMu_2017B", "SingleMu_2017C", "SingleMu_2017D", "SingleMu_2017E", "SingleMu_2017F"]
    samples = signals + bkgs + data

    color = {}
    color["H2Mu_gg"] = kRed
    color["H2Mu_VBF"] = kViolet - 6
    color["H2Mu_ZH"] = kGreen+2
    color["H2Mu_WH_pos"] = kMagenta
    color["H2Mu_WH_neg"] = kCyan+1

    color["ZJets_AMC"] = kOrange - 3
    color["tt"] = kAzure + 8
#    color[] =

    directories = {}
    canvas = {}
    histos = {}
    stack_all = {}
    stack_sig = {}
    stack_data = {}
    ratios = {}
    legend = {}
    for term in terms:
	directories[term] = out_file.mkdir(term)
	histos[term] = {}
	stack_all[term] = THStack("h_stack_"+term, term)
	stack_sig[term] = THStack("sig_stack_"+term, "sig_"+term)
	stack_data[term] = THStack("data_stack"+term, "data_"+term)

    for sample in samples:
	for term in terms:
	    if in_file.Get(sample + "_" + term) is None:
		histos[term][sample] = TH1F()
	    else:
	        histos[term][sample] = in_file.Get(sample+ "_" +term).Clone()

    out_file.cd()
    for term in terms:
	for sample in signals:
	    histos[term][sample].SetLineColor(color[sample])
	    histos[term][sample].SetLineWidth(2)
	    stack_sig[term].Add(histos[term][sample])

	histos[term]["signal"] = stack_sig[term].GetStack().Last()
	stack_all[term].Add(histos[term]["signal"])
        histos[term]["signal"].SetFillColor( kGray+2 )
	histos[term]["signal"].SetLineWidth(0)
	for sample in bkgs:
	    histos[term][sample].SetFillColor(color[sample])
	    stack_all[term].Add(histos[term][sample])
	histos[term]["MC"] = stack_all[term].GetStack().Last()
	for sample in data:
	    stack_data[term].Add(histos[term][sample])
	histos[term]["data"]= stack_data[term].GetStack().Last()

	ratios[term] = TGraphAsymmErrors()
        ratios[term].Divide(histos[term]["data"], histos[term]["MC"], "pois")
        ratios[term].SetName("ratiograph_"+term)

        legend[term] = TLegend(0.9,0.5,1,1)
        legend[term].AddEntry(histos[term]["data"], "data")
        legend[term].AddEntry(histos[term]["signal"], "signal sum")
        for sample in bkgs + signals:
            legend[term].AddEntry(histos[term][sample], sample )
	ratioplot( term, stack_sig[term], stack_all[term], histos[term]["data"], ratios[term], legend[term])


#	canvas["ratio_"+term].SaveAs(PLOT_DIR+"/"+LABEL+"/files/sum/plots/ratio"+term+".png")




    out_file.Close()
main()

