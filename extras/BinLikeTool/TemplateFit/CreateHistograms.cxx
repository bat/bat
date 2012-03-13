#include <vector> 

#include <TROOT.h> 
#include <TFile.h>
#include <TH1D.h> 
#include <TMath.h> 
#include <TRandom3.h> 

double expectation(double m, double x); 

int CreateHistograms()
{
	int nbins = 20; 
	int nev = 1000; 

	gRandom = new TRandom3(1000); 

	TFile * file = new TFile("histograms.root", "recreate"); 

	int nm = 6;
	double m[6] = {150, 160, 170, 180, 190, 200}; 
	double nevents[6] = {0.1*nev, 0.5*nev, nev, 0.5*nev, 0.3*nev, 0.1*nev}; 

	std::vector<TH1D*> HistCont; 

	for (int i = 0; i < nm; ++i) {
		TH1D * hist = new TH1D(Form("temp_%i", i), "", nbins, 100.0, 250.0); 
		
		for (int j = 1; j <= 100; ++j) {
			hist -> SetBinContent(j, 
														expectation(m[i], hist->GetBinCenter(j))); 
		}
		HistCont.push_back(hist); 
	}

	for (int i = 0; i < nm; ++i) {
		TH1D * hist = new TH1D(Form("hist_%i", i), "", nbins, 100.0, 250.0); 

		hist->FillRandom(HistCont.at(i), nevents[i]); 

		hist->Write(); 
	}

	file->Close(); 

	return 0; 
}


double expectation(double m, double x)
{
	double mean = 150.0 + 0.2 * m; 
	double rms  =   0.0 + 0.1 * m; 

	return TMath::Gaus(x, mean, rms, true);
}
