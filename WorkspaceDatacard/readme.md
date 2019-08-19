## updates with systematic uncertainties (XWZ 2019.08.19)

**general instructions:**


With the current update, datacards are written with several signal channels (by year and production mode), and one net bkg.

This update is made only on mass-shape workspace and datacard.

To run the limits:

First run macros/PrepareSystematicsHists.py to prepare the histograms for systematics

Then write a config file 

Then run the macros/MakeWorkspaceDatacard.py with that config file (just like before)

**Note that:** in DataLoader.py, the __ source == 'Sys_test' __ option is used so the other options are intact for Andrew's possible usage. Later it shoule be merged into the others.

**more details are below**
1. Added options for signals and sys_names in the config file. 
It will be used in the class WorkspaceAndDatacardMaker and DataLoader
2. Added class SystematicsConfig
It contains all systematic uncertainties to be used in the datacard making, including theoretical, phenominological, and experimetal uncertainties
3. changes in the DatacardHelper to write uncertainty values



## Limit Setting (initial push)
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
