// Sample.h

#ifndef SAMPLE
#define SAMPLE

#include <iostream>

#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"

/* #include "TMath.h" */
/* #include "TTree.h" */
/* #include "TEntryList.h" */

/* #include "LumiReweightingStandAlone.h" */
/* #include "VarSet.h" */
/* #include "BranchSet.h" */

class Sample
{
    public:
        Sample();
	Sample(TString _name, TString _sampleType = "NONE", std::vector<TString> _filenames = {});
        ~Sample();

        TString name;
        TString filename;
	std::vector<TString> filenames;
        TString treename;
	TChain* chain;
        /* reweight::LumiReWeighting* lumiWeights;  // Information for pileup reweighting  */

        /* TString dir;           // DAS directory */
        /* TString pileupfile;    // used for pile up reweighting through lumiWeights */
        TString sampleType;    // "data", "signal", "background"

        /* int plotColor;         // the color used when plotting the sample */
        int nOriginal;         // the number of events run over to get this sample
        int nOriginalWeighted; // the number of original events run over to get this sample accounting for genWeights
        /* int N;                 // the number of events in the sample */

        float xsec;            // xsec in pb
        float lumi;            // the luminosity of the data or effective luminosity of MC

        /* VarSet vars;           // all of the variables from the ttree */
        /* BranchSet branches;    // all of the branches from the ttree, which we link to the vars */

        /* int getEntry(int i);                    // load the ith event from the ttree into vars */
        /* int getEntry(int i, TEntryList* list);  // load ith event from the list into vars */
        /*                                         // the ith event in the list maps to the jth tree entry */

        void calculateNoriginal(const bool verbose = false);                      // calculate nOriginal and nOriginalWeighted
        /* void setBranchAddresses(TString options = "");  // link the values in the tree to vars */
        /* double getWeight();                             // get the weight for the histogram based upon the pileup weight and the MC gen weight */

        // scale by xsec*lumi/N_weighted for MC and 1.0 for data
        float getLumiScaleFactor(const float luminosity, const bool verbose = false);

    protected:
	std::vector<TFile*> files; // files with the ttree

};



#endif
