
## Helper functions to modify existing fits

import os
import sys

import ROOT as R
import ROOT.RooFit as RF

sys.path.insert(1, '%s/../FitBackground/python' % os.getcwd())
import FitFunctions as FF  ## From FitBackground/python/FitFunctions.py


def FreezeParams(fit, to_freeze):

    ## Freeze fit coefficients
    if 'Coefs' in to_freeze:
        # print '\nFreezing coefficients:'
        for i in range(len(fit.coef_list)):
            # print fit.coef_list[i].GetName()
            fit.coef_list[i].setConstant(R.kTRUE)

    ## Freeze other shape parameters
    if 'Params' in to_freeze:
        # print '\nFreezing parameters:'
        for j in range(len(fit.params[0])):
            for i in range(len(fit.params)):
                # print fit.params[i][j].GetName()
                fit.params[i][j].setConstant(R.kTRUE)

## End function: def FreezeParams(fit, mask):


