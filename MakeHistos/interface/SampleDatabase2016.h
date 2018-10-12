#ifndef SAMPLE_DATABASE_2016_H
#define SAMPLE_DATABASE_2016_H

#include "H2MuAnalyzer/MakeHistos/interface/Sample.h"

std::map<TString, Sample*> GetSamples2016(std::map<TString, Sample*>& samples, TString location, TString select = "ALL", bool info_only = false);

#endif
