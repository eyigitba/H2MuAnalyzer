

##############################
##      FitFunctions.py     ##
##  Standard functions for  ##
##  signal and background   ##
##############################

import ROOT as R
import ROOT.RooFit as RF


class FitFunction:

    def __init__(self, name, hist, fit_type, order, x_range, x_blind):
        self.name      = name            ## Name of fit
        self.hist      = hist.Clone()    ## Input data histogram
        self.fit_type  = fit_type        ## Type of fit (exponential, Breit-Wigner, etc.)
        self.order     = order           ## Order of the fit
        self.x_range   = x_range         ## Range of data to fit
        self.x_blind   = x_blind         ## Range of data to blind
        self.var       = 0               ## RooRealVar from the hist
        self.dat       = 0               ## RooDataHist from the hist
        self.params    = [[]]            ## RooRealVar parameters to the fit
        self.funcs     = []              ## RooAbsPDF component functions of the fit
        self.arg_sets  = []              ## RooArgSets, each containing one RooAbsPdf from funcs
        self.amp_vars  = []              ## RooRealVars with the amplitude of each function
        self.coef_list = R.RooArgList()  ## Coefficients of each RooAbsPDF component from funcs
        self.func_list = R.RooArgList()  ## Container of RooArgSets from arg_sets
        self.model     = 0

        InitHist(self)
        InitVarDat(self)
        if fit_type == 'expo':
            InitExpo(self)


## Modify the histogram range and blind the signal region
def InitHist(FF):
    
    ## Zoom in on a given mass region
    if len(FF.x_range) == 2:
        FF.hist.GetXaxis().SetRangeUser(FF.x_range[0], FF.x_range[1])
    elif len(FF.x_range) != 0:
        print '\nIn InitHist(), trying to set the following invalid range:'
        print FF.x_range
        print 'Will continue without changing the original histogram range'
    
    ## Blind signal region, typically  120 - 130 GeV
    nBins = FF.hist.GetNbinsX()
    if len(FF.x_blind) == 2:
        for i in range(1, nBins):
            if FF.hist.GetBinLowEdge(i) >= FF.x_blind[0] and FF.hist.GetBinLowEdge(i) < FF.x_blind[1]:
                FF.hist.SetBinContent(i, 0)
    elif len(FF.x_blind) != 0:
        print '\nIn InitHist(), trying to blind the following invalid range:'
        print FF.x_blind
        print 'Will continue without blinding anything'

## End function: InitHist()


## Create a RooRealVar and RooDataHist from the input data histogram
def InitVarDat(FF):
    
    if len(FF.x_range) == 2:
        x_min = FF.x_range[0]
        x_max = FF.x_range[1]
    else:
        x_min = FF.hist.GetBinLowEdge(1)
        x_max = FF.hist.GetBinLowEdge( FF.hist.GetNbinsX() + 1 )
    FF.var = R.RooRealVar('var_'+FF.name, 'var_'+FF.name, x_min, x_max)
    if len(FF.x_blind) == 2:
        FF.var.setRange('loM', x_min, FF.x_blind[0])
        FF.var.setRange('hiM', FF.x_blind[1], x_max)

    FF.dat = R.RooDataHist('dat_'+FF.name, 'dat_'+FF.name, R.RooArgList(R.RooArgSet(FF.var)), FF.hist)

## End function InitVarDat()


## Create an exponential model
def InitExpo(FF):

    for i in range(FF.order):
        FF.params.append([])
        FF.params[i].append( R.RooRealVar('slope%d' % i, 'Slope of exponential #%d'  % i, -0.01, -1.0, -0.0001) )
        FF.funcs    .append( R.RooExponential('expo%d'  % i, 'e^(x*slope)', FF.var, FF.params[i][0]) )
        FF.arg_sets .append( R.RooArgSet(FF.funcs[i]) )
        FF.amp_vars .append( R.RooRealVar('amp%d' % i, 'Amplitude of exponential #%d' % i, 1.0/FF.order, 0.0, 1.0) )

        FF.func_list.add(FF.funcs[i])
        if (i < FF.order - 1):
            FF.coef_list.add(FF.amp_vars[i])

    ## Use a "recursive" fit (kTRUE) to ensure all-positive coefficients
    FF.model = R.RooAddPdf('mod_'+FF.name, 'Sum of %d exponentials' % FF.order, FF.func_list, FF.coef_list, R.kTRUE)

## End function InitExpo()


## Fit the model to the data
def DoFit(FF):
    
    if len(FF.x_blind) == 2:
        FF.model.fitTo(FF.dat, RF.Save(), RF.Range('loM,hiM'))
    else:
        FF.model.fitTo(FF.dat, RF.Save())

## End function DoFit()
