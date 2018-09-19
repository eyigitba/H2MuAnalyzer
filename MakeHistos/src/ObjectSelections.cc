#include "H2MuAnalyzer/MakeHistos/interface/ObjectSelections.h"


bool MuonPass(MuonInfo & muon, int pt_cut , std::string id , bool verbose )
{
	bool muon_pass = false;

	if( id == "medium" and muon.pt > pt_cut and abs(muon.eta) < 2.4 and muon.relIso < 0.25 and muon.isMediumID == 1) muon_pass = true; 
	else if(id == "tight") muon_pass = false; // for now
	
	return muon_pass;
}


bool JetPass(JetInfo & jet, int pt_cut , std::string PU_ID , bool verbose )
{
	bool jet_pass = false;
	
	if( jet.pt <= pt_cut or abs(jet.eta) >= 4.7 ) return jet_pass = false;
	else
	{
	    if(PU_ID == "loose")
	    {
		if(abs(jet.eta) < 2.5)
		{
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.89)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.97)  return false;
		}
                else if(abs(jet.eta) < 2.75)
                {
		        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.52)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.68)  return false;
                }
		else if(abs(jet.eta) < 3.0)
		{
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.38)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.53)  return false;
		}
                else if(abs(jet.eta) < 5)
		{
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.30)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.47)  return false;
		}
	    }
	    else if(PU_ID == "medium")
	    {
		if(abs(jet.eta) < 2.5)
		{
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID < 0.61)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID < 0.18)  return false;
		}
                else if(abs(jet.eta) < 2.75)
                {
		        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.35)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.55)  return false;
		}
                else if(abs(jet.eta) < 3.0)
                {
		        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.23)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.42)  return false;
		}
                else if(abs(jet.eta) < 5)
                {
		        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.17)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.36)  return false;
		}
	    }
	    else if(PU_ID == "tight")
	    {
		if(abs(jet.eta) < 2.5)
                {
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID < 0.86)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID < 0.69)  return false;
                }
                else if(abs(jet.eta) < 2.75)
                {
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.10)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.35)  return false;
                }
                else if(abs(jet.eta) < 3.0)
                {
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.23)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.42)  return false;
                }
                else if(abs(jet.eta) < 5)
                {
                        if (jet.pt >=30 and jet.pt< 50 and jet.puID <-0.17)  return false;
                        if (jet.pt >=10 and jet.pt< 30 and jet.puID <-0.36)  return false;
                }
	    }
	}
	return jet_pass = true;
}

