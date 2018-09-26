####################################################
###    make stack and ratio from added histos    ###
####################################################

from ROOT import *

def ratioplot( direc, term, sig_stack, all_stack, h_data, ratio_graph, legend ):   # Do not use TRatioPlot! It is a devil!    -XWZ 19.09.2018
#    direc.cd()
    canv = TCanvas("ratio_"+term, "ratio_"+term, 600, 600)
    canv.Clear()
    
    upper_pad = TPad("upperpad_"+term, "upperpad_"+term, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.Draw()
    upper_pad.cd()
    upper_pad.SetLogy()
    all_stack.SetMinimum(1e-5)
    all_stack.SetMaximum(1e6)
#    all_stack.SetMinimum(1e-3)
#    all_stack.SetMaximum(1e10)
    all_stack.Draw("HIST")
    h_data.SetMarkerStyle(20)
    h_data.Draw("SAME")
    for histo in sig_stack:
	histo.Draw("SAMEHIST")
    legend.Draw()

    canv.cd()
    lower_pad = TPad("lowerpad_"+term, "lowerpad_"+term, 0,0.05, 1,0.2)
    lower_pad.SetTopMargin(0.05)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
#    ratio_graph.SetMinimum(0.6)
#    ratio_graph.SetMaximum(1.7)
    ratio_graph.SetMinimum(0.5)
    ratio_graph.SetMaximum(2.5)
    
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
    canv.SaveAs("/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/VH_toy_2017_v4_v7/files/sum" + "/plots_WHlep" + "/ratio_" + term + ".png")


def main():

    file_dir="/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/VH_toy_2017_v4_v7/files/sum"
    out_file = TFile( file_dir + "/WHlep_stack" + ".root", "RECREATE")
    in_file = TFile.Open( file_dir + "/all_WHlep.root", "READ")

    terms = ["dimuon_mass", "dimuon_pt","dimuon_eta", "dimuon_delta_eta", "dimuon_delta_phi",
	     "leading_muon_pt", "leading_muon_eta", "subleading_muon_pt", "subleading_muon_eta",
	     "dijet_mass", "dijet_pt", "dijet_eta", "dijet_delta_eta", "dijet_delta_phi",
	     "leading_jet_pt", "leading_jet_eta", "subleading_jet_pt", "subleading_jet_eta",
	     "MET", "nJets", "nBJets", "nVertices", "nElectrons",  "nMuons"]
    signals = ["H2Mu_ttH", "H2Mu_WH_neg", "H2Mu_WH_pos", "H2Mu_ZH", "H2Mu_VBF", "H2Mu_gg"]
    bkgs = ["diboson", "tt", "ZJets_AMC"]
    diboson = ['WW', 'WZ_3l_AMC', 'ZZ_2l_2v', 'ZZ_4l']
    data = ["SingleMu_2017B", "SingleMu_2017C", "SingleMu_2017D", "SingleMu_2017E", "SingleMu_2017F"]
    samples = signals + bkgs + data

    color = {}
    color["H2Mu_gg"] = kRed
    color["H2Mu_VBF"] =  kBlue + 1
    color["H2Mu_ZH"] =  kOrange + 7
    color["H2Mu_WH_pos"] = kGreen + 2
    color["H2Mu_WH_neg"] = kViolet + 1 
    color["H2Mu_ttH"]  = kPink + 6 

    color["ZJets_AMC"] =  kAzure + 7
    color["tt"] =       kGreen - 9
    color["diboson"] =  kCyan - 7
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

    for term in terms:
	histos[term]["diboson"] = None
	for sample in diboson:
	    if histos[term]["diboson"] is None:
		print "copying first diboson  " + term
		if not in_file.GetListOfKeys().Contains(sample + "_" + term) :
			print sample + "_" + term + "\n"
		
	    	histos[term]["diboson"] = in_file.Get(sample+ "_" +term).Clone()
		histos[term]["diboson"].SetName(sample+ "_" + "diboson")
	    else:
		histos[term]["diboson"].Add( in_file.Get(sample+ "_" +term) )

    for sample in samples:
	for term in terms:
	    if not in_file.GetListOfKeys().Contains(sample + "_" + term) :
		print sample 
		print term + "\n\n"
#		histos[term][sample] = TH1F()
#		histos[term][sample].SetName(sample + "_" + term)
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
        histos[term]["signal"].SetFillColor( kGray )
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

        legend[term] = TLegend(0.7,0.7,1,1)
        legend[term].AddEntry(histos[term]["data"], "data")
        legend[term].AddEntry(histos[term]["signal"], "signal sum")
        for sample in bkgs + signals:
            legend[term].AddEntry(histos[term][sample], sample )
	ratioplot( directories[term], term, stack_sig[term], stack_all[term], histos[term]["data"], ratios[term], legend[term])


#	canvas["ratio_"+term].SaveAs("/afs/cern.ch/work/x/xzuo/public/H2Mu/2018/Histograms/MC_data_comparison_2017_v4_v4/files/sum/plots/ratio"+term+".png")




    out_file.Close()
main()

