
## Class and functions to load input data for making workspaces and datacards

import os
import sys
import array

import ROOT as R

sys.path.insert(1, '%s/../AutoCat/python' % os.getcwd())
import AutoCatTools as ACT  ## From AutoCat/python/AutoCatTools.py


class DataLoader:

    ## ============================================ ##
    ##  Initialize based on user and configuration  ##
    ## ============================================ ##

    def __init__(self, _user, _source, _in_dir, _in_file, _cat_loc, _dist, _min_max, _rebin):

        print '\nInstantiating DataLoader with user = %s and source = %s' % (_user, _source)

        self.user       = _user     ## Who is running the script: abrinke1, xzuo, pbortigno
        self.source     = _source   ## Data source: abrinke1, xzuo, pbortigno, acarnes
        self.in_dir     = _in_dir   ## Data input directory
        self.in_file    = _in_file  ## Data input file (first initialized to the file name)
        self.cat_loc    = _cat_loc  ## Full category name from input file to run over
        self.dist       = _dist     ## Distribution to be used for signal extraction
        self.min_max    = _min_max  ## Minimum and maximum variable values considered
        self.sig_hists  = []        ## Signal histograms, including total stack in [0]
        self.bkg_hists  = []        ## Background histograms, including total stack in [0]
        self.data_hist  = 0         ## Rebinned data histogram
        self.rebin      = _rebin    ## Whether to rebin maximizing S^2/B
        self.sig_rebin  = []        ## Rebinned signal histograms, including total stack in [0]
        self.bkg_rebin  = []        ## Rebinned background histograms, including total stack in [0]
        self.data_rebin = 0         ## Rebinned data histogram

        self.LoadHistograms()

    ## End function: def __init__(self, _user, _source, _in_dir, _in_file, _cat_loc, _dist, _min_max, _rebin):


    ## ============================== ##
    ##  Clear values and free memory
    ## ============================== ##
        
    def Clear(self):

        if not self.in_file is None: self.in_file.Close()
        for hist in self.sig_hists: del hist
        for hist in self.bkg_hists: del hist
        if not self.data_hist is None: del self.data_hist
        for hist in self.sig_rebin: del hist
        for hist in self.bkg_rebin: del hist
        if not self.data_rebin is None: del self.data_rebin
        
        self.user       = ''
        self.source     = ''
        self.in_dir     = ''
        self.in_file    = None
        self.cat_loc    = ''
        self.dist       = ''
        self.min_max    = []
        self.sig_hists  = []
        self.bkg_hists  = []
        self.data_hist  = 0
        self.rebin      = []
        self.sig_rebin  = []
        self.bkg_rebin  = []
        self.data_rebin = 0
        

    ## ================================= ##
    ##  Load histograms from input file
    ## ================================= ##

    def LoadHistograms(self):

        if   (self.source == 'abrinke1'):
            in_file_name = self.in_dir+'/plots/'+self.cat_loc+'/'+self.in_file
        elif (self.source == 'abrinke1_TMVA'):
            in_file_name = self.in_dir+'/'+self.in_file
        elif (self.source == 'xzuo_TMVA'):
            in_file_name = self.in_dir+'/'+self.in_file
        elif (self.source == 'XXX'):
            in_file_name = 'YYY'
        else:
            print '\n\nNo input file location for source %s!!!\nExiting!' % self.source
            sys.exit()

        print '\nOpening input file %s' % in_file_name
        self.in_file = R.TFile(in_file_name)


        ## Configuration from Andrew's StackPlots.root files
        if (self.source == 'abrinke1'):
            print '  * Loading %s_Net_Data, %s_Net_Sig, and %s_Net_Bkg' % (self.dist, self.dist, self.dist)
            self.data_hist =      self.in_file.Get(self.dist+'_Net_Data').Clone('Net_Data')
            self.sig_hists.append(self.in_file.Get(self.dist+'_Net_Sig').Clone('Net_Sig'))
            self.bkg_hists.append(self.in_file.Get(self.dist+'_Net_Bkg').Clone('Net_Bkg'))
            self.in_file.cd('groups')
            ## Loop over group histograms that went into the stack
            for hist in R.gDirectory.GetListOfKeys():
                hstr = hist.GetName()
                if hstr.startswith(self.dist):
                    print '    - Loading group %s' % hstr
                    if ('H2Mu' in hstr or 'Signal' in hstr):
                        self.sig_hists.append(self.in_file.Get('groups/'+hstr).Clone())
                        self.sig_hists[-1].SetName(hstr.replace(self.dist+'_', '').replace(' ', '_'))
                    else:
                        self.bkg_hists.append(self.in_file.Get('groups/'+hstr).Clone())
                        self.bkg_hists[-1].SetName(hstr.replace(self.dist+'_', '').replace(' ', '_'))
        ## End conditional: if (self.source == 'abrinke1'):

        ## Configuration from Andrew's TMVA training output file
        elif (self.source == 'abrinke1_TMVA'):
            h_str = '%s/Method_%s/%s/MVA_%s' % (self.cat_loc, self.dist, self.dist, self.dist)
            print '  * Loading '+h_str+'_S and _B'
            self.data_hist =      self.in_file.Get(h_str+'_B_high').Clone('Net_Data')
            self.sig_hists.append(self.in_file.Get(h_str+'_S_high').Clone('Net_Sig'))
            self.bkg_hists.append(self.in_file.Get(h_str+'_B_high').Clone('Net_Bkg'))

            ## Scale signal and background to expected yields
            self.data_hist   .Scale(254.7 / self.data_hist.Integral())
            self.sig_hists[0].Scale(1.368 / self.sig_hists[0].Integral())
            self.bkg_hists[0].Scale(254.7 / self.bkg_hists[0].Integral())

        ## End conditional: if (self.source == 'xzuo_TMVA'):

        ## Configuration from Xunwu's TMVA training output file
        elif (self.source == 'xzuo_TMVA'):
            h_str = '%s/Method_%s/%s/MVA_%s' % (self.cat_loc, self.dist, self.dist, self.dist)
            print '  * Loading '+h_str+'_S and _B'
            self.data_hist =      self.in_file.Get(h_str+'_B').Clone('Net_Data')
            self.sig_hists.append(self.in_file.Get(h_str+'_S').Clone('Net_Sig'))
            self.bkg_hists.append(self.in_file.Get(h_str+'_B').Clone('Net_Bkg'))

            ## Scale signal and background to expected yields
            self.data_hist   .Scale(97.71 / self.data_hist.Integral())
            self.sig_hists[0].Scale(0.643 / self.sig_hists[0].Integral())
            self.bkg_hists[0].Scale(97.71 / self.bkg_hists[0].Integral())

        ## End conditional: if (self.source == 'xzuo_TMVA'):

        ## Configuration for input files from XXX
        elif (self.source == 'XXX'):
            in_file_name = 'YYY'
        else:
            print '\n\nNo input file configuration for source %s!!!\nExiting!' % self.source
            sys.exit()


        ## Set content for bins outside specified range to 0
        if len(self.min_max) == 2:
            for iBin in range(1, self.data_hist.GetNbinsX()+1):
                if ( self.data_hist.GetBinLowEdge(iBin)   < self.min_max[0] or
                     self.data_hist.GetBinLowEdge(iBin+1) > self.min_max[1] ):

                    print '      ** Blinding bin %d (%f to %f)' % ( iBin, self.data_hist.GetBinLowEdge(iBin),
                                                                    self.data_hist.GetBinLowEdge(iBin+1) )
                    self.data_hist.SetBinContent(iBin, 0)
                    self.data_hist.SetBinError  (iBin, 0)
                    for hist in self.sig_hists:
                        hist.SetBinContent(iBin, 0)
                        hist.SetBinError  (iBin, 0)
                    for hist in self.bkg_hists:
                        hist.SetBinContent(iBin, 0)
                        hist.SetBinError  (iBin, 0)
        
        ## Rebin histograms for maximum S/sqrt(B) significance with a conservative estimate of uncertainties
        if self.rebin != False:
            ## rebin[0] sets the significance calculation strategy: 'noSyst', 'nominal', or 'conserv'
            ## rebin[1] sets the minimum allowed background in a single bin
            ## rebin[2] sets the maximum allowed loss in total significance from rebinning
            new_bins = array.array('d', ACT.MergeBins(self.sig_hists[0], self.bkg_hists[0],  self.rebin[0], self.rebin[1], self.rebin[2], False) )
            self.data_rebin = self.data_hist.Rebin( len(new_bins)-1, self.data_hist.GetName(), new_bins )
            for hist in self.sig_hists:
                self.sig_rebin.append( hist.Rebin( len(new_bins)-1, hist.GetName(), new_bins ) )
            for hist in self.bkg_hists:
                self.bkg_rebin.append( hist.Rebin( len(new_bins)-1, hist.GetName(), new_bins ) )

    ## End function: def LoadHistograms(self):
