#include <iostream>

#include "BinLikeModel.h"

#include <TROOT.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TLatex.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
BinLikeModel::BinLikeModel() : BCModel()
{
	// default MCMC settings for this model
	MCMCSetNLag(10);
	//	MCMCSetNIterationsRun(1000000);
}

// ---------------------------------------------------------
BinLikeModel::~BinLikeModel()
{
	int nhist = GetNHistograms(); 
	for (int i = 0; i < nhist; ++i)
		delete fErrHistContainer.at(i); 
}

// ---------------------------------------------------------
double BinLikeModel::LogLikelihood(std::vector<double> parameters)
{
	double logprob = 0.;

	// loop over histograms
	int nhist = GetNHistograms(); 
	for (int ihist = 0; ihist < nhist; ++ihist){

		// get histogram
		TH1D* hist = fHistogramContainer.at(ihist); 
		
		// calculate integral 
		double integral = hist->Integral(); 

		// get bin width (of first bin)
		double binwidth = hist->GetBinWidth(1); 

		// get leading parameter value 
		double parval = GetLeadParValue(ihist); 

		// loop over bins 
		int nbins = GetNBins(ihist); 
		for (int ibin = 1; ibin <= nbins; ++ibin) {

			// calculate expectation 
			double expectation = integral*Expectation(parameters, parval, hist->GetBinCenter(ibin))*binwidth; 

			// get observed number of events 
			double observed = hist->GetBinContent(ibin); 

			// calculate Poisson probability for this bin 
			logprob += BCMath::LogPoisson(observed, expectation); 
		}
	}
	
	return logprob;
}

// ---------------------------------------------------------
double BinLikeModel::LogAPrioriProbability(std::vector<double> parameters)
{
	double logprob = 0.;

	return logprob;
}

// ---------------------------------------------------------
int BinLikeModel::AddHistogram(TH1D* hist, double parvalue, double weight)
{
	// set histogram properties
	hist->SetStats(kFALSE); 

	// add histogram to container
	fHistogramContainer.push_back(hist); 

	// add leading parameter value to container 
	fLeadParContainer.push_back(parvalue); 

	// add histogram weight to container 
	fHistWeightContainer.push_back(weight); 
	
	// create histograms for uncertainty determination
	double maxi = 1.5 * hist->GetMaximum();

	TH2D* temphist = new TH2D(Form("UncertaintyExp_%i", BCLog::GetHIndex()), "",
														hist->GetNbinsX(),
														hist->GetXaxis()->GetXmin(),
														hist->GetXaxis()->GetXmax(),
														100,
														0.0,
														maxi);

	fErrHistContainer.push_back(temphist); 
	
	// no errors
	return 1; 
}

// ---------------------------------------------------------
int BinLikeModel::GetNBins(int index)
{
	if (index < 0 || index >= GetNHistograms())
		return 0; 

	else 
		return fHistogramContainer.at(index)->GetNbinsX(); 
}

// ---------------------------------------------------------
double BinLikeModel::GetLeadParValue(int index)
{
	if (index < 0 || index >= GetNHistograms())
		return 0; 

	else 
		return fLeadParContainer.at(index); 
}

// ---------------------------------------------------------
void BinLikeModel::PrintHistograms(const char* path)
{
	// create canvas
	TCanvas * c1 = new TCanvas("c1"); 
	c1->cd(); 

	// get number of histograms 
	int nhist = GetNHistograms(); 

	// loop over histograms 
	for (int ihist = 0; ihist < nhist; ++ihist){ 

		// get histogram 
		TH1D* hist = fHistogramContainer.at(ihist); 

		// get number of bins
		int nbins = hist->GetNbinsX(); 

		// get error histograms
		TH2D* hist_errexp = fErrHistContainer.at(ihist);

		// fill graph
		int npoints = 1000; 
		TGraph gr(npoints); 
		gr.SetLineWidth(1);
		gr.SetLineColor(kRed);

		double xmin = hist->GetXaxis()->GetXmin(); 
		double xmax = hist->GetXaxis()->GetXmax(); 

		for (int i = 0; i <npoints; ++i) {
			double x = xmin + double(i) * (xmax-xmin)/double(npoints-1); 
			double expectation = hist->Integral()*Expectation(GetBestFitParameters(), GetLeadParValue(ihist), x) * hist->GetBinWidth(1); 
			gr.SetPoint(i, x, expectation); 
		}

		// define error graphs
		TGraphAsymmErrors * graph_error_exp = new TGraphAsymmErrors(nbins);
		graph_error_exp->SetLineWidth(0);
		graph_error_exp->SetFillStyle(1001);
		graph_error_exp->SetFillColor(kYellow);
		
		TGraphAsymmErrors * graph_error_obs = new TGraphAsymmErrors(nbins);
		graph_error_obs->SetMarkerStyle(0);
		
		bool flag_error2 = true;

		if (flag_error2)
			for (int i = 1; i <= nbins; ++i)
				{
					TH1D * proj = hist_errexp->ProjectionY("_py", i, i);
					if (proj->Integral() > 0)
						proj->Scale(1.0 / proj->Integral());
					double quantiles[3];
					double sums[3] = {0.16, 0.5, 0.84};
					proj->GetQuantiles(3, quantiles, sums);
					graph_error_exp->SetPoint(i-1, hist->GetBinCenter(i), quantiles[1]);
					graph_error_exp->SetPointError(i-1, 0.0, 0.0, quantiles[1] - quantiles[0], quantiles[2]-quantiles[1]);
					delete proj;
				}
		
		// draw 
		c1->cd(); 
		hist->Draw("PE"); 		
		graph_error_exp->Draw("SAME3");
		gr.Draw("SAMEC"); 
		hist->Draw("SAMEPE"); 		

		// print to file 
		c1 -> Print(Form("%shistogram_%.0f.eps", path, GetLeadParValue(ihist))); 
	}

	// clean up
	delete c1; 
}

// ---------------------------------------------------------
void BinLikeModel::PrintChi2Summary()
{
	// create canvas
	TCanvas * c1 = new TCanvas("c1"); 
	c1->cd(); 

	int nhist = GetNHistograms();

	// define histogram 
	double xmin = fLeadParContainer.at(0)-0.5*(fLeadParContainer.at(1) - fLeadParContainer.at(0));
	double xmax = fLeadParContainer.at(nhist-1)+0.5*(fLeadParContainer.at(nhist-1) - fLeadParContainer.at(nhist-2));
	int nbinsmax = 0; 
	double ymin = 0; 
	double ymax = 0; 

	for (int ihist = 0; ihist < nhist; ++ihist){
		if (fHistogramContainer.at(ihist)->GetNbinsX()>nbinsmax) {
			nbinsmax = fHistogramContainer.at(ihist)->GetNbinsX(); 
			ymax = 2.0 * nbinsmax - GetNParameters();
		}
		if (CalculateChi2Hist(ihist) > ymax) {
			ymax = 1.1*CalculateChi2Hist(ihist); 
		}
	}


	TH2D * hist_axes = new TH2D("hist_axes", ";lead. par. value;#chi^{2};", 
															1, xmin, xmax, 
															1, ymin, ymax);
	hist_axes->SetStats(kFALSE); 

	// define graphs with ndf 
	TGraphErrors * graph_ndf = new TGraphErrors(nhist); 
	graph_ndf->SetMarkerSize(0); 
	graph_ndf->SetFillStyle(1001);
	graph_ndf->SetFillColor(kYellow); 

	// fill estimated ndf
	for (int ihist = 0; ihist < nhist; ++ihist){
		graph_ndf->SetPoint(ihist, fLeadParContainer.at(ihist), fHistogramContainer.at(ihist)->GetNbinsX()-0.5*GetNParameters());
		graph_ndf->SetPointError(ihist, 0.0, 0.5*GetNParameters());
	}

	// define graph with number of events 
	TGraph * graph_nev = new TGraph(nhist);
	graph_nev->SetFillColor(21);
	graph_nev->SetFillStyle(1001); 

	// fill number of events 
	double sum = 0; 
	for (int ihist = 0; ihist < nhist; ++ihist)
		sum+= fHistogramContainer.at(ihist)->GetEntries(); 
	for (int ihist = 0; ihist < nhist; ++ihist){
		graph_nev->SetPoint(ihist, fLeadParContainer.at(ihist), fHistogramContainer.at(ihist)->GetEntries() / sum * 0.3 * (ymax-ymin) + ymin);
	}

	// define graph with chi2
	TGraphErrors * graph_chi2 = new TGraphErrors(nhist); 
	graph_chi2->SetMarkerStyle(21);
	graph_chi2->SetMarkerColor(kBlack);
	graph_chi2->SetMarkerSize(1);

	// fill chi2
	for (int ihist = 0; ihist < nhist; ++ihist){
		graph_chi2->SetPoint(ihist, fLeadParContainer.at(ihist), CalculateChi2Hist(ihist));
		graph_chi2->SetPointError(ihist, 0.0, sqrt(2.0*fHistogramContainer.at(ihist)->GetNbinsX()));
	}

	// define legend
	TLegend * legend = new TLegend(0.51, 0.69, 0.86, 0.86); 
	legend->SetBorderSize(0); 
	legend->SetFillColor(kWhite);
	legend->AddEntry(graph_ndf, "estimated ndf", "F");
	legend->AddEntry(graph_chi2, "#chi^{2}", "P"); 
	legend->AddEntry(graph_nev, "statistics indicator","F"); 

	// define text 
	TLatex * latex = new TLatex(); 
	latex->SetTextSize(0.03); 

	// draw 
	hist_axes->Draw("");
	graph_nev->Draw("SAMEB");
	graph_ndf->Draw("SAME3");
	graph_chi2->Draw("SAMEPE");
	legend->Draw("SAME");
	latex->DrawLatex(xmin+0.0 *(xmax-xmin), ymin+1.1*(ymax-ymin), TString::Format("#chi^{2}/ndf (#chi^{2}-prob.) = %4.3g/%3.3g (%3.2g)", CalculateChi2(), double(GetNDF()), CalculateChi2Prob()));

	// print
	c1->Print("summary_fit.eps");
}

// ---------------------------------------------------------
void BinLikeModel::MCMCUserIterationInterface()
{
	// loop over histograms
	int nhist = GetNHistograms(); 
	for (int ihist = 0; ihist < nhist; ++ihist){

		// get histogram
		TH1D* hist = fHistogramContainer.at(ihist); 
		
		// get bin width (of first bin)
		double binwidth = hist->GetBinWidth(1); 

		// get leading parameter value 
		double parval = GetLeadParValue(ihist); 

		// loop over bins 
		int nbins = hist->GetNbinsX(); 
		for (int ibin = 1; ibin <= nbins; ++ibin) {

			// get x value
			double x = hist->GetBinCenter(ibin); 

			// calculate expectation 
			double expectation = hist->Integral()*Expectation(fMCMCx, parval, x) * binwidth; 

			// set bin content
			fErrHistContainer.at(ihist)->Fill(x, expectation);
		}
	}
}

// ---------------------------------------------------------
int BinLikeModel::GetNDF()
{
	int ndf = 0; 

	// loop over histograms
	int nhist = GetNHistograms(); 
	for (int ihist = 0; ihist < nhist; ++ihist){
		ndf += fHistogramContainer.at(ihist) -> GetNbinsX(); 
	}

	return ndf - GetNParameters(); 
}

// ---------------------------------------------------------
double BinLikeModel::CalculateChi2()
{
	// initialize chi2
	double chi2 = 0; 
	
	// get best fit parameters
	std::vector<double> parameters = GetBestFitParameters();

	// loop over histograms
	int nhist = GetNHistograms(); 
	for (int ihist = 0; ihist < nhist; ++ihist){

		// calculate chi2 of histogram
		chi2 += CalculateChi2Hist(ihist); 
	}
	
	// return chi2
	return chi2; 
}

// ---------------------------------------------------------
double BinLikeModel::CalculateChi2Hist(int index)
{
	// initialize chi2
	double chi2 = 0; 
	
	// get best fit parameters
	std::vector<double> parameters = GetBestFitParameters();

	// get histogram
	TH1D* hist = fHistogramContainer.at(index); 
		
	// get bin width (of first bin)
	double binwidth = hist->GetBinWidth(1); 

	// get leading parameter value 
	double parval = GetLeadParValue(index); 

	// loop over bins 
	int nbins = hist->GetNbinsX(); 
	for (int ibin = 1; ibin <= nbins; ++ibin) {
		
		// get x value
		double x = hist->GetBinCenter(ibin); 

		// calculate expectation 
		double expected = hist->Integral()*Expectation(parameters, parval, x) * binwidth; 

		// get observation 
		double observed = hist->GetBinContent(ibin); 

		// increase chi2
		chi2 += (observed - expected)*(observed - expected)/expected; 
	}
	
	// return chi2
	return chi2; 
}

// ---------------------------------------------------------
double BinLikeModel::CalculateChi2Prob()
{
	double chi2 = CalculateChi2();
	int ndf = GetNDF();

	// return chi2 probability
	return TMath::Prob(chi2, ndf);
}

// ---------------------------------------------------------
