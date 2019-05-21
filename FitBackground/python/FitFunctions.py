

##############################
##      FitFunctions.py     ##
##  Standard functions for  ##
##  signal and background   ##
##############################

import sys

import ROOT as R
import ROOT.RooFit as RF
import ROOT.TMath as TM


class FitFunction:

    def __init__(self, name, hist, fit_type, order, x_range, x_blind = [], var_name = None):
        self.name      = name            ## Name of fit
        self.hist_orig = hist.Clone()    ## Input data histogram - original, unblinded
        self.hist      = hist.Clone()    ## Input data histogram
        self.fit_type  = fit_type        ## Type of fit (exponential, Breit-Wigner, etc.)
        self.order     = order           ## Order of the fit
        self.x_range   = x_range         ## Range of data to fit
        self.x_blind   = x_blind         ## Range of data to blind
        self.x_bins    = [0,0,0]         ## Binning of data in fit range: # of bins, low bin, high bin
        self.var       = 0               ## RooRealVar from the hist
        self.var_name  = var_name        ## RooRealVar from the hist
        self.dat       = 0               ## RooDataHist from the hist
        self.params    = [[]]            ## RooRealVar parameters to the fit
        self.funcs     = []              ## RooAbsPDF component functions of the fit
        self.arg_sets  = []              ## RooArgSets, each containing one RooAbsPdf from funcs
        self.amp_vars  = []              ## RooRealVars with the amplitude of each function
        self.coef_list = R.RooArgList()  ## Coefficients of each RooAbsPDF component from funcs
        self.func_list = R.RooArgList()  ## Container of RooArgSets from arg_sets
        self.model     = 0               ## RooAddHist model from func_list and coef_list
        self.fit_hist  = 0               ## TH1 version of fit

        if var_name is None:
            var_name = 'var_'+name
        InitHist(self)
        InitVarDat(self)

        ## Different functions from 2016 with constant ranges found here:
        ## https://github.com/uiowahep/Analysis/blob/master/Configuration/higgs/UF_AWB_settings.py#L261
        if fit_type == 'poly':
            InitPoly(self)
        elif fit_type == 'expo':
            InitExpo(self)
        elif fit_type == 'Bern':
            InitBern(self)
        elif fit_type == 'BWZ':
            InitBWZ(self)
        elif fit_type == 'BWZRed':
            InitBWZRed(self)
        elif fit_type == 'PolyPlusBWZ':
            InitPolyPlusBWZ(self)
        elif fit_type == 'Gaus':
            InitGaus(self)
        else:
            print 'Fit type %s does not match any valid option!!! Exiting.' % fit_type
            sys.exit()

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
        for i in range(1, nBins+1):
            if FF.hist.GetBinLowEdge(i) >= FF.x_blind[0] and FF.hist.GetBinLowEdge(i) < FF.x_blind[1]:
                FF.hist.SetBinContent(i, 0)
    elif len(FF.x_blind) != 0:
        print '\nIn InitHist(), trying to blind the following invalid range:'
        print FF.x_blind
        print 'Will continue without blinding anything'

    ## Get the bins corresponding to the fit range (nBins, low bin, high bin)
    for i in range(1, nBins+1):
        if (FF.hist.GetBinLowEdge(i) == FF.x_range[0]):
            FF.x_bins[1] = i
        if (FF.hist.GetBinLowEdge(i+1) == FF.x_range[1]):
            FF.x_bins[2] = i
    FF.x_bins[0] = FF.x_bins[2] - FF.x_bins[1] + 1

## End function: InitHist()


## Create a RooRealVar and RooDataHist from the input data histogram
def InitVarDat(FF):

    ## Set x axis range if it was not done manually
    if len(FF.x_range) != 2:
        FF.x_range = [FF.hist.GetBinLowEdge(1), FF.hist.GetNbinsX() + 1]
    ## Get minimum and maximum values in x range
    x_min = FF.x_range[0]
    x_max = FF.x_range[1]

    FF.var = R.RooRealVar(FF.var_name, FF.var_name, x_min, x_max)
    if len(FF.x_blind) == 2:
        FF.var.setRange('loM', x_min, FF.x_blind[0])
        FF.var.setRange('hiM', FF.x_blind[1], x_max)

    FF.dat = R.RooDataHist('dat_'+FF.name, 'dat_'+FF.name, R.RooArgList(R.RooArgSet(FF.var)), FF.hist)

## End function InitVarDat()


## Fit the model to the data
def DoFit(FF):
    
    print '\n\nAbout to fit model %s to data!' % FF.model.GetName()

    ## RooFit or Minuit strategy, error level, etc?
    ## https://root.cern.ch/doc/master/classRooAbsPdf.html#a8f802a3a93467d5b7b089e3ccaec0fa8
    ## https://root.cern.ch/doc/master/classRooMinuit.html#a605d27ee6cfbd36d5a61e8085bed0539
    if len(FF.x_blind) == 2:
        FF.model.fitTo(FF.dat, RF.Save(), RF.Range('loM,hiM'))
    else:
        FF.model.fitTo(FF.dat, RF.Save())

    ## Creates a histogram of the model, with arbitrary normalization
    FF.fit_hist = FF.model.createHistogram(FF.var_name, FF.x_bins[0])
    ## Find the normalization of the fitted range, and scale model histogram
    f_norm = FF.fit_hist.Integral()
    h_norm = FF.hist_orig.Integral(FF.x_bins[1], FF.x_bins[2])
    ## If a range of the histogram had been pre-blinded, need to use only fitted range
    if len(FF.x_blind) == 2:
        f_norm = 0
        h_norm = 0
        nBins = FF.hist.GetNbinsX()
        for i in range(FF.x_bins[1], FF.x_bins[2]+1):
            if FF.hist.GetBinLowEdge(i+1) <= FF.x_blind[0] or FF.hist.GetBinLowEdge(i) >= FF.x_blind[1]:
                f_norm += FF.fit_hist.GetBinContent(i)
                h_norm += FF.hist.GetBinContent(i)

    ## Scale fit histogram to data normalization in the fitted range
    FF.fit_hist.Scale( h_norm / f_norm )
    ## Set line color red, remove error bars
    FF.fit_hist.SetLineColor( R.kRed )
    for i in range(FF.fit_hist.GetNbinsX()):
        FF.fit_hist.SetBinError(i+1, 0)

    print 'Post-fit amplitude values are:'
    for i in range(len(FF.coef_list)):
        print '  * %s: %f' % (FF.coef_list[i].GetName(), FF.coef_list[i].getValV())
    print 'Post-fit parameter values are:'
    for j in range(len(FF.params[0])):
        for i in range(len(FF.params)):
            print '  * %s: %f' % (FF.params[i][j].GetName(), FF.params[i][j].getValV())

## End function DoFit()



############################################################
#####    *** Initialization of different models ***    #####
############################################################


###************************************###
###  Create a simple polynomial model  ###
###************************************###
## Re-scale x-axis so x = 0 lies to the right of the data being fit (i.e. all data has negative x-values)
## This means all polynomial terms will be "falling"
## The only variables are the normalization and x = 0 offset of each polynomial term
def InitPoly(FF):
    print '\nInitializing a polynomial of order %s' % FF.order

    x_min = FF.x_range[0]
    x_max = FF.x_range[1]
    x_wid = x_max - x_min

    x_var = '(%s - %.2f)/%.2f' % (FF.var_name, x_min, x_wid)

    for i in range(FF.order):
        if i != 0: FF.params.append([])
        ## "x = 0" is always >= x_max.  Use sqrt(off*off) to allow smoth crossing across off = 0 boundary.
        FF.params[i].append( R.RooRealVar('x%d' % (i+1), 'Offset term #%d' % (i+1), 1.0, -999, 999) )
        ## Start with linear term: constant term can be picked up with linear term + offset
        FF.funcs    .append( R.RooGenericPdf( 'poly%d' % (i+1), 'poly%d' % (i+1),
                                              'pow(sqrt(x%d*x%d) + 1 - %s, %d)' % ((i+1), (i+1), x_var, (i+1)),
                                              R.RooArgList(R.RooArgSet(FF.var, FF.params[i][0])) ) )
        FF.arg_sets .append( R.RooArgSet(FF.funcs[i]) )
        ## Each term has an amplitude somewhere between 0 and 1.0
        FF.amp_vars .append( R.RooRealVar('poly_amp%d' % (i+1), 'Amplitude of term #%d' % (i+1), 1.0/FF.order, 0.0, 1.0) )
        
        FF.func_list.add(FF.funcs[i])
        if (i < FF.order - 1):
            FF.coef_list.add(FF.amp_vars[i])

    ## Model composed of the sum of the polynomial terms
    ## Use a "recursive" fit (kTRUE) to ensure all-positive coefficients
    FF.model = R.RooAddPdf('mod_'+FF.name, 'Polynomial of order %d' % FF.order, FF.func_list, FF.coef_list, R.kTRUE)

## End function InitPoly()


###**************************************###
###  Create a sum of exponentials model  ###
###**************************************###
def InitExpo(FF):

    for i in range(FF.order):
        if i != 0: FF.params.append([])
        FF.params[i].append( R.RooRealVar('slope%d' % i, 'Slope of exponential #%d' % i, -0.01, -1.0, -0.0001) )
        FF.funcs    .append( R.RooExponential('expo%d'  % i, 'e^(x*slope)', FF.var, FF.params[i][0]) )
        FF.arg_sets .append( R.RooArgSet(FF.funcs[i]) )
        FF.amp_vars .append( R.RooRealVar('expo_amp%d' % i, 'Amplitude of exponential #%d' % i, 1.0/FF.order, 0.0, 1.0) )

        FF.func_list.add(FF.funcs[i])
        if (i < FF.order - 1):
            FF.coef_list.add(FF.amp_vars[i])

    ## Model composed of the sum of the exponential terms
    ## Use a "recursive" fit (kTRUE) to ensure all-positive coefficients
    FF.model = R.RooAddPdf('mod_'+FF.name, 'Sum of %d exponentials' % FF.order, FF.func_list, FF.coef_list, R.kTRUE)

## End function InitExpo()


## Create a sum of Bernstein polynomials model
def InitBern(FF):

    x_str = FF.var_name
    x_var = '(%s - %.2f)/%.2f' % (x_str, FF.x_range[0], FF.x_range[1] - FF.x_range[0])
    params_arg_list = R.RooArgList( R.RooArgSet(FF.var) )

    for i in range(FF.order+1):
        binom = int(TM.Binomial(FF.order, i))
        # print 'Appending Bernstein #%d : %d*pow(%s, %d)*pow(1 - %s, %d)' % (i, binom, x_var, i, x_var, FF.order-i)
        FF.funcs   .append( R.RooGenericPdf('Bern%d' % i, 'Bern%d' % i, '%d*pow(%s, %d)*pow(1 - %s, %d)' % (binom, x_var, i, x_var, FF.order-i), params_arg_list) )
        FF.arg_sets.append( R.RooArgSet(FF.funcs[i]) )
        FF.amp_vars.append( R.RooRealVar('Bern_amp%d' % i, 'Amplitude of Bernstein #%d' % i, 1.0/(FF.order+1), 0.0, 1.0) )

        FF.func_list.add(FF.funcs[i])
        if (i < FF.order):
            FF.coef_list.add(FF.amp_vars[i])

    ## Use a "recursive" fit (kTRUE) to ensure all-positive coefficients
    FF.model = R.RooAddPdf('mod_'+FF.name, '%d degree Bernstein polynomial' % FF.order, FF.func_list, FF.coef_list, R.kTRUE)

## End function InitBern()


## Create a Breit-Wigner model
def InitBWZ(FF):

    FF.order = 1
    FF.params[0].append( R.RooRealVar('a1', 'a1', 7.5, -99.9, 99.9) )
    params_arg_list = R.RooArgList( R.RooArgSet(FF.var, FF.params[0][0]) )
    x_str = FF.var_name
    FF.funcs    .append( R.RooGenericPdf('BWZ', 'BWZ', '2.5*exp(pow(a1/100.0, 2)*'+x_str+')/(pow('+x_str+'-91.2, 2) + pow(2.5/2, 2))', params_arg_list) )
    FF.arg_sets .append( R.RooArgSet(FF.funcs[0]) )
    FF.amp_vars .append( R.RooRealVar('BWZ_amp', 'Amplitude of BWZ', 1.0, 0.0, 1.0) )

    FF.func_list.add(FF.funcs[0])
    FF.coef_list.add(FF.amp_vars[0])

    # FF.model = R.RooAddPdf('mod_'+FF.name, 'BWZ', FF.func_list, FF.coef_list, R.kFALSE)
    FF.model = FF.funcs[0]

## End function InitBWZ()


## Create a reduced Breit-Wigner model
def InitBWZRed(FF):

    if (FF.order < 2 or FF.order > 3):
        print '\n\nMAJOR ERROR!!! Cannot initialize BWZRed with order %d' % FF.order
        sys.exit()

    FF.params[0].append( R.RooRealVar('a1', 'a1', 1.4, -99.9, 99.9) )
    FF.params[0].append( R.RooRealVar('a2', 'a2', 7.5, -99.9, 99.9) )
    if (FF.order >= 3):
        FF.params[0].append( R.RooRealVar('a3', 'a3', 0.1, -99.9, 99.9) )

    params_arg_list = R.RooArgList( R.RooArgSet(FF.var, FF.params[0][0], FF.params[0][1], FF.params[0][2]) )
    x_str = FF.var_name
    # FF.funcs.append( R.RooGenericPdf('BWZRed', 'BWZRed', 'exp(a2*'+x_str+' + a3*'+x_str+'*'+x_str+')/(pow('+x_str+'-91.1876,a1) + pow(2.5/2,a1))', params_arg_list) )
    if (FF.order == 3):
        FF.funcs.append( R.RooGenericPdf('BWZRed', 'BWZRed', 'exp(pow(a2/100.0, 2)*'+x_str+' + pow(a3/100.0, 2)*pow('+x_str+', 2)) / (pow('+x_str+'-91.2, a1*a1) + pow(2.5/2, a1*a1))', params_arg_list) )
    if (FF.order == 2):
        FF.funcs.append( R.RooGenericPdf('BWZRed', 'BWZRed', 'exp(pow(a2/100.0, 2)*'+x_str+' + pow(0.0/100.0, 2)*pow('+x_str+', 2)) / (pow('+x_str+'-91.2, a1*a1) + pow(2.5/2, a1*a1))', params_arg_list) )
    FF.arg_sets .append( R.RooArgSet(FF.funcs[0]) )
    FF.amp_vars .append( R.RooRealVar('BWZRed_amp', 'Amplitude of BWZRed', 1.0, 0.0, 1.0) )

    FF.func_list.add(FF.funcs[0])

    # FF.model = R.RooAddPdf('mod_'+FF.name, 'BWZRed', FF.func_list)
    FF.model = FF.funcs[0]

## End function InitBWZRed()


## Create a Breit-Wigner + linear
def InitPolyPlusBWZ(FF):

    ## BWZ piece
    FF.params[0].append( R.RooRealVar('a1', 'a1', 7.5, -99.9, 99.9) )
    FF.funcs    .append( R.RooGenericPdf('BWZ', 'BWZ', '2.5*exp(pow(a1/100.0, 2)*'+FF.var_name+')/(pow('+FF.var_name+'-91.2, 2) + pow(2.5/2, 2))', R.RooArgList( R.RooArgSet(FF.var, FF.params[0][0]) )) )
    FF.arg_sets .append( R.RooArgSet(FF.funcs[0]) )
    FF.amp_vars .append( R.RooRealVar('BWZ_amp', 'Amplitude of BWZ', 0.5, 0.0, 1.0) )
    FF.func_list.add(FF.funcs[0])
    FF.coef_list.add(FF.amp_vars[0])

    ## Polynomial piece
    x_min = FF.x_range[0]
    x_max = FF.x_range[1]
    x_wid = x_max - x_min

    x_var = '(%s - %.2f)/%.2f' % (FF.var_name, x_min, x_wid)

    for i in range(1, FF.order+1):
        FF.params.append([])
        ## "x = 0" is always >= x_max.  Use sqrt(off*off) to allow smoth crossing across off = 0 boundary.
        FF.params[i].append( R.RooRealVar('x%d' % i, 'Offset term #%d' % i, 1.0, -999, 999) )
        ## Start with linear term: constant term can be picked up with linear term + offset
        FF.funcs    .append( R.RooGenericPdf( 'poly%d' % i, 'poly%d' % i,
                                              'pow(sqrt(x%d*x%d) + 1 - %s, %d)' % (i, i, x_var, i),
                                              R.RooArgList(R.RooArgSet(FF.var, FF.params[i][0])) ) )
        FF.arg_sets .append( R.RooArgSet(FF.funcs[i]) )
        ## Each term has an amplitude somewhere between 0 and 1.0
        FF.amp_vars .append( R.RooRealVar('poly_amp%d' % i, 'Amplitude of polynomial term #%d' % i, 0.5/FF.order, 0.0, 1.0) )
        
        FF.func_list.add(FF.funcs[i])
        if (i < FF.order):
            FF.coef_list.add(FF.amp_vars[i])

    ## Model composed of the Breit-Wigner plus the sum of the polynomial terms
    ## Use a "recursive" fit (kTRUE) to ensure all-positive coefficients
    FF.model = R.RooAddPdf('mod_'+FF.name, 'Breit-Wigner + polynomial of order %d' % FF.order, FF.func_list, FF.coef_list, R.kTRUE)

## End function InitPolyPlusBWZ()


###***********************************###
###  Create a sum of Gaussians model  ###
###***********************************###
def InitGaus(FF):

    ## Fit starting points and ranges, for mean + width + amplitude
    coef = {}
    for i in range(FF.order):
        coef[i] = [125, 115, 135, 3.0, 0.9, 99.9, 1.0/FF.order, 0.0, 1.0]

    ## For triple-gaussian H2Mu signal, order narrow to wide
    if FF.order == 3:
        mean = FF.hist.GetMean()
        rms  = FF.hist.GetRMS()
        coef[0] = [mean, mean - 0.2*rms, mean + 0.2*rms, 0.6*rms, 0.2*rms, 1.2*rms, 0.7, 0.5, 1.0]
        coef[1] = [mean, mean - 0.4*rms, mean + 0.4*rms, 1.0*rms, 0.4*rms, 1.8*rms, 0.2, 0.0, 1.0]
        coef[2] = [mean, mean - 2.4*rms, mean + 0.8*rms, 2.4*rms, 1.2*rms, 4.8*rms, 0.1, 0.0, 1.0]

    for i in range(FF.order):
        if i != 0: FF.params.append([])
        FF.params[i].append( R.RooRealVar('mean%d'  % (i+1), 'Mean of Gaussian #%d'  % (i+1), coef[i][0], coef[i][1], coef[i][2]) )
        FF.params[i].append( R.RooRealVar('width%d' % (i+1), 'Width of Gaussian #%d' % (i+1), coef[i][3], coef[i][4], coef[i][5]) )
        FF.funcs    .append( R.RooGaussian('gaus%d' % (i+1), 'gaus(x - mean, width) #%d' % (i+1), FF.var, FF.params[i][0], FF.params[i][1]) )
        FF.arg_sets .append( R.RooArgSet(FF.funcs[i]) )
        FF.amp_vars .append( R.RooRealVar('Gaus_amp%d' % (i+1), 'Amplitude of Gaussian #%d' % (i+1), coef[i][6], coef[i][7], coef[i][8]) )
        
        FF.func_list.add(FF.funcs[i])
        if (i < FF.order - 1):
            FF.coef_list.add(FF.amp_vars[i])
            
    ## Model composed of the sum of the Gaussian terms
    ## Use a "recursive" fit (kTRUE) to ensure all-positive coefficients
    FF.model = R.RooAddPdf('mod_'+FF.name, 'Sum of %d Gaussians' % FF.order, FF.func_list, FF.coef_list, R.kTRUE)

## End function InitGaus()
