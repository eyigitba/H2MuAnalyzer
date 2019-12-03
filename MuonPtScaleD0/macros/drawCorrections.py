#! /usr/bin/env python

###############################################
###          drawCorrections.py             ###
###                                         ###
###    Draw percent corrections on muon pT  ###
###     vs d0_BS and pT to show the effect  ###
###                                         ###
###              Efe Yigitbasi              ###
###               2.11.2019                 ###
###############################################

## Basic python includes for manipulating files
import sys
import os

## More python tools
import math
import numpy

## ROOT includes
import ROOT as R
import ROOT.RooFit as RF
R.gROOT.SetBatch(True)  ## Don't display histograms or canvases when drawn

OUTDIR = 'plots_corrections_d0_BS'

## Main function executed by ./macros/FindPeaks.py
def main():

    print 'Drawing 2D correction functions in d0_BS and pT'
    canv = {}

    print '\nCreating output directory %s' % OUTDIR
    if not os.path.exists(OUTDIR):
      os.makedirs(OUTDIR)

    for eta in ['eta_0_0p9', 'eta_0p9_1p7', 'eta_1p7_inf']: # |eta| binned
        for year in ['2016', '2017', '2018']:
          for sample in ['DY', 'ttH']:
            p_str = 'pT_corrections' + '_' + sample + '_' + eta + '_' + year
            canv[p_str] = R.TCanvas(p_str, p_str, 1600, 1200)
            if sample == 'DY':
              if eta == 'eta_0_0p9':
                if year == "2016": 
                  f =  R.TF2(p_str,"411.343/10000 * abs(x) * y",-0.01,0.01,20,300);
                if year == "2017": 
                  f =  R.TF2(p_str,"650.839/10000 * abs(x) * y",-0.01,0.01,20,300);
                if year == "2018": 
                  f =  R.TF2(p_str,"582.32/10000 * abs(x) * y",-0.01,0.01,20,300);
              if eta == 'eta_0p9_1p7':
                if year == "2016": 
                  f =  R.TF2(p_str,"673.398/10000 * abs(x) * y",-0.01,0.01,20,300);
                if year == "2017": 
                  f =  R.TF2(p_str,"988.369/10000 * abs(x) * y",-0.01,0.01,20,300);
                if year == "2018": 
                  f =  R.TF2(p_str,"974.047/10000 * abs(x) * y",-0.01,0.01,20,300);
              if eta == 'eta_1p7_inf':
                if year == "2016": 
                  f =  R.TF2(p_str,"1098.984/10000 * abs(x) * y",-0.01,0.01,20,300);
                if year == "2017": 
                  f =  R.TF2(p_str,"1484.616/10000 * abs(x) * y",-0.01,0.01,20,300);
                if year == "2018": 
                  f =  R.TF2(p_str,"1263.388/10000 * abs(x) * y",-0.01,0.01,20,300);
            else:
              continue
            f.Draw('colz')
            f.SetMaximum(1)
            f.SetMinimum(0.0001)            # contours = {}
            # contours[0] = 0.05
            # contours[1] = 0.1
            # contours[2] = 0.2
            # f.SetContour(4)
            # f.SetContourLevel(1,0.05)
            # f.SetContourLevel(2,0.1)
            # f.SetContourLevel(3,0.2)
            # f.DrawCopy('cont3 same')
            # f.SetMaximum(1)
            # f.SetMinimum(0.0001)
            R.gPad.SetLogz()
            canv[p_str].SaveAs('%s/%s.png' % (OUTDIR, p_str))

            
    print '... and, done!'
        
## End function: main()
    

if __name__ == '__main__':
    main()
