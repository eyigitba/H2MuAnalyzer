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

LABEL = 'lepMVA_variables_nonprompt_v1'  ## Sub-folder within PLOT_DIR containing histograms


def ratioplot( term, MC_stack, h_data, ratio_graph, legend ):   # Do not use TRatioPlot! It is a devil!    -XWZ 19.09.2018
    canv = TCanvas("ratio_"+term, "ratio_"+term, 600,600)
    canv.Clear()
    
    upper_pad = TPad("upperpad_"+term, "upperpad_"+term, 0,0.2, 1,1)
    upper_pad.SetBottomMargin(0.05);
    upper_pad.Draw()
    upper_pad.cd()
    upper_pad.SetLogy()
#    MC_stack.SetMinimum(1e-5)
#    MC_stack.SetMaximum(1e6)
    MC_stack.SetMinimum(1e-3)
    MC_stack.SetMaximum(1e8)
    MC_stack.Draw("HIST")
    h_data.SetMarkerStyle(20)
    h_data.Draw("SAME")
    legend.Draw()

    canv.cd()
    lower_pad = TPad("lowerpad_"+term, "lowerpad_"+term, 0,0.05, 1,0.2)
    lower_pad.SetTopMargin(0.05)
    lower_pad.SetGridy()
    lower_pad.Draw()
    lower_pad.cd()
#    ratio_graph.SetMinimum(0.5)
#    ratio_graph.SetMaximum(1.5)
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
#    canv.SaveAs(PLOT_DIR+"/"+LABEL+"/files/sum" + "/plots_none" + "/ratio_" + term + ".png")

def count_event(histos, term, samples):
    print "counting events in plot " + term
    sum_num = 0
    for sample in samples:
	print sample + "\t %6.3f" %histos[term][sample].Integral()
        sum_num += histos[term][sample].Integral()
    print "sum" + "\t %6.3f \n" %sum_num
 

def grab_hists(in_file, histos, hist_name, samples, tW, tt, tt_ljj):
    for sample in samples:
	if sample is "tW":
	    histos[hist_name]["tW"] = None
	    for thing in tW:
		if histos[hist_name]["tW"] is None:
		    if not in_file.GetListOfKeys().Contains(thing + "_" + hist_name) :
		        print thing + "_" + hist_name + "  not in in_file\n"
		        histos[hist_name]["tW"] = TH1D()
		    else:
		    	histos[hist_name]["tW"] = in_file.Get(thing + "_" + hist_name).Clone()
		    histos[hist_name]["tW"].SetName("tW" + "_" + hist_name)
		else:
		    histos[hist_name]["tW"].Add(in_file.Get(thing + "_" + hist_name))
	elif sample is "tt":
	    histos[hist_name]["tt"] = None
            for thing in tt:
                if histos[hist_name]["tt"] is None:
                    if not in_file.GetListOfKeys().Contains(thing + "_" + hist_name) :
                        print thing + "_" + hist_name + "  not in in_file\n"
                        histos[hist_name]["tt"] = TH1D()
                    else:
                        histos[hist_name]["tt"] = in_file.Get(thing + "_" + hist_name).Clone()
                    histos[hist_name]["tt"].SetName("tt" + "_" + hist_name)
                else:
                    histos[hist_name]["tt"].Add(in_file.Get(thing + "_" + hist_name))
	    histos[hist_name]["tt"].Scale(0.5)
	elif sample is "tt_ljj":
            histos[hist_name]["tt_ljj"] = None
            for thing in tt_ljj:
                if histos[hist_name]["tt_ljj"] is None:
                    if not in_file.GetListOfKeys().Contains(thing + "_" + hist_name) :
                        print thing + "_" + hist_name + "  not in in_file\n"
                        histos[hist_name]["tt_ljj"] = TH1D()
                    else:
                        histos[hist_name]["tt_ljj"] = in_file.Get(thing + "_" + hist_name).Clone()
                    histos[hist_name]["tt_ljj"].SetName("tt_ljj" + "_" + hist_name)
                else:
                    histos[hist_name]["tt_ljj"].Add(in_file.Get(thing + "_" + hist_name))
            histos[hist_name]["tt_ljj"].Scale(0.5)

	else:
	    if not in_file.GetListOfKeys().Contains(sample + "_" + hist_name):
		print sample + "_" + hist_name + "   not in in_file\n"
		histos[hist_name][sample] = TH1D()
	    else:
		histos[hist_name][sample] = in_file.Get(sample + "_" + hist_name).Clone()


def main():

    file_dir=PLOT_DIR+"/"+LABEL
    out_file = TFile( file_dir + "/stack_all_5" + ".root", "RECREATE") # use update to handle seg breaks. "recreate" recommanded for normal cases
    in_file = TFile.Open( file_dir + "/all_samples.root", "READ")

    terms = {}
    terms["extra"] = []#["leading_mu_pt", "leading_mu_abs_eta", "trailing_mu_pt", "trailing_mu_abs_eta"]
    terms["mu"] =  [#"pt_ratio", "pt_rel", 
		#"jet_deepCSV", "trk_mult",
		    #"relIso_03", "miniIso", "miniIsoCharged", 
		    #"SIP_3D", "abs_dxy", "abs_dz", 
		"segCompat", "lepMVA"]
    terms["ele"] = [#"pt", "abs_eta", "pt_ratio", 
		#"pt_rel", "jet_deepCSV", "trk_mult", 
		    #"relIso_03", "miniIso", "miniIsoCharged",
		    #"SIP_3D", "abs_dxy", "abs_dz", 
		"muSegComp", "mvaID", "lepMVA"]

#    signals = ["H2Mu_ttH", "H2Mu_WH_neg", "H2Mu_WH_pos", "H2Mu_ZH", "H2Mu_VBF", "H2Mu_gg"]
    bkgs = ["ZZ_4l", "WZ_3l", "tW", "tt_ljj", "tt", "ZJets_AMC"]
    bkgs = ["tW", "tt", "tt_ljj", "ZJets_AMC"]  # CERN_lepMVA_v2
    tW = ["tW_neg", "tW_pos"]
    tt = ["tt_ll_MG", "tt_ll_POW"] 
    tt_ljj = ["tt_ljj_POW_1", "tt_ljj_POW_2"]
    data = ["SingleMu_2017B", "SingleMu_2017C", "SingleMu_2017D", "SingleMu_2017E", "SingleMu_2017F"]
    samples =  bkgs + data

    color = {}
    color["H2Mu_gg"] = kRed
    color["H2Mu_VBF"] =  kBlue + 1
    color["H2Mu_ZH"] =  kOrange + 7
    color["H2Mu_WH_pos"] = kGreen + 2
    color["H2Mu_WH_neg"] = kViolet + 1 
    color["H2Mu_ttH"]  = kPink + 6 

    color["ZJets_AMC"] =  kAzure + 7
    color["tt"] =       kGreen - 9
    color["tW"] =  kCyan - 7
    color["tt_ljj"] = kGray
    color["ZZ_4l"] = kYellow - 9
    color["WZ_3l"] = kViolet - 9
#    color[] =

    directories = {}
    directories["mu"] = out_file.mkdir("mu")
    directories["ele"] = out_file.mkdir("ele")

    histos = {}
    stack_MC = {}
    stack_data = {}
    ratios = {}
    legend = {}

    for term in terms["extra"]:
	histos[term] = {}
	stack_MC[term] = THStack("MC_stack_"+term, term)
	stack_data[term] = THStack("data_stack"+term, "data_"+term)
	grab_hists(in_file, histos, term, samples, tW, tt, tt_ljj)

    for lepton in ["mu", "ele"]:
    	for term in terms[lepton]:
	    for mva in ["tight","neg", "all"]:
		hist_name = lepton + "_" + term + "_" + mva
		histos[hist_name] = {}
		stack_MC[hist_name] =   THStack("MC_stack_"   + hist_name, "MC_"   + hist_name)
		stack_data[hist_name] = THStack("data_stack_" + hist_name, "data_" + hist_name)
		grab_hists(in_file, histos, hist_name, samples, tW, tt, tt_ljj)

#    if "none" not in in_file.GetName():
#      for term in terms:
#	for sample in samples:
#	    if histos[term][sample].GetNbinsX() == 50:   histos[term][sample].Rebin(2)


#    count_event(histos, "dimuon_mass", signals)
#    count_event(histos, "dimuon_mass", bkgs)


    out_file.cd()
    for term in histos.keys():
#	for sample in signals:
#	    histos[term][sample].SetLineColor(color[sample])
#	    histos[term][sample].SetLineWidth(2)

	for sample in bkgs:
	    histos[term][sample].SetFillColor(color[sample])
	    stack_MC[term].Add(histos[term][sample])
	histos[term]["MC"] = stack_MC[term].GetStack().Last()
	for sample in data:
	    stack_data[term].Add(histos[term][sample])
	histos[term]["data"]= stack_data[term].GetStack().Last()

	ratios[term] = TGraphAsymmErrors()
        ratios[term].Divide(histos[term]["data"], histos[term]["MC"], "pois")
        ratios[term].SetName("ratiograph_"+term)

        legend[term] = TLegend(0.7,0.7,1,1)
        legend[term].AddEntry(histos[term]["data"], "data")
#        legend[term].AddEntry(histos[term]["signal"], "signal sum")
        for sample in bkgs:
            legend[term].AddEntry(histos[term][sample], sample )
	if "ele" in term:
	    directories["ele"].cd()
	else:
	    directories["mu"].cd()
	ratioplot( term, stack_MC[term], histos[term]["data"], ratios[term], legend[term])


    out_file.Close()
main()

