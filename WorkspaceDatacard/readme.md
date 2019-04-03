### Limit Setting
1. MakeWorkspaceDatacard.py creates workspace containing data, signal and background model and
   datacard for Higgs combine tool. You need to give a file path for your histogram root file and
   category at the end of the script. Then run it by doing 'python macros/MakeWorkspaceDatacard.py'
2. Run combine by doing 'combine resulting_datacard_filename.txt + your-options'.

* To read your data/MC/signal rootfile into MakeWorkspaceDatacard.py you need to follow some
  convention for naming your histograms.
    - For background: cateogry + '_Net_Bkg'
    - For signal: category + '_Net_signal'
    - Or you can change the script so that the script can find your histograms.
        ** Look at "def setNetBackgroundHist(self):" and below in MakeWorkspaceDatacard.py
           and modify to match to your histogram names.

* BGSFrun.py performs fitting taking functions from PDFDatabase.py
* PDFDatabase.py has various defintions of functions. You can add your functions.

Let me know if anything not working.
