#include <iostream>

#include "StackModel.h"

#include <TROOT.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
StackModel::StackModel() : BCModel()
{
	fUncertaintyHistogramExp = 0;
	fUncertaintyHistogramObsPosterior = 0;
}

// ---------------------------------------------------------
StackModel::StackModel(const char * name) : BCModel(name)
{
	fUncertaintyHistogramExp = 0;
	fUncertaintyHistogramObsPosterior = 0;
}

// ---------------------------------------------------------
StackModel::~StackModel()
{
	if (fUncertaintyHistogramExp)
		delete fUncertaintyHistogramExp;

	if (fUncertaintyHistogramObsPosterior)
		delete fUncertaintyHistogramObsPosterior;
}

// ---------------------------------------------------------
double StackModel::LogLikelihood(std::vector <double> parameters)
{
	double logprob = 0.;

	int nbins      = fDataHistogram.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double nexp = 0;
		double ndata = fDataHistogram.GetBinContent(ibin);

		for (int itemp = 0; itemp < ntemplates; ++itemp)
			nexp += parameters.at(itemp) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

		logprob += BCMath::LogPoisson(ndata, nexp);
	}

	return logprob;
}

// ---------------------------------------------------------
double StackModel::LogAPrioriProbability(std::vector <double> parameters)
{
	double logprob = 0.;

	for (int i = 0; i < int(GetNParameters()); ++i)
		logprob -= log(GetParameter(i)->GetUpperLimit() - GetParameter(i)->GetLowerLimit());

	return logprob;
}

// ---------------------------------------------------------
int StackModel::SetDataHistogram(TH1D hist)
{
	// set histogram
	fDataHistogram = hist;

	// set histogram style
	fDataHistogram.SetMarkerStyle(20);
	fDataHistogram.SetMarkerSize(1.1);
	fDataHistogram.SetStats(kFALSE);

	// create histograms for uncertainty determination
	double maxi = 1.5 * fDataHistogram.GetMaximum();

	fUncertaintyHistogramExp = new TH2D(Form("UncertaintyExp_%i", BCLog::GetHIndex()), "",
																			fDataHistogram.GetNbinsX(),
																			fDataHistogram.GetXaxis()->GetXmin(),
																			fDataHistogram.GetXaxis()->GetXmax(),
																			100,
																			0.0,
																			maxi);

	fUncertaintyHistogramObsPosterior = new TH2D(Form("UncertaintyObsPosterior_%i", BCLog::GetHIndex()), "",
																							 fDataHistogram.GetNbinsX(),
																							 fDataHistogram.GetXaxis()->GetXmin(),
																							 fDataHistogram.GetXaxis()->GetXmax(),
																							 int(maxi) + 1,
																							 -0.5,
																							 double(int(maxi))+0.5);

	// no errors
	return 1;
}

// ---------------------------------------------------------
int StackModel::AddTemplateHistogram(TH1D hist, const char * name, double Nmin, double Nmax)
{
	// check if histogram if filled
	if (hist.Integral() <= 0.)
		return 0;

	// compare histogram properties
	if (CompareHistogramProperties(fDataHistogram, hist) != 1)
		return 0;

	// set histogram color
	hist.SetFillColor( 2 + int(fTemplateHistogramContainer.size()) );
	hist.SetFillStyle(1001);

	// add histogram if it is the first
	if ( int(fTemplateHistogramContainer.size()) == 0)
	{
		hist.Scale(1.0 / hist.Integral());
		fTemplateHistogramContainer.push_back(hist);
		fTemplateNameContainer.push_back(name);
	}

	else
	{
		// add histogram if it has the same properties than those in the
		// container
		if (CompareHistogramProperties(fTemplateHistogramContainer.at(0), hist))
		{
			hist.Scale(1.0 / hist.Integral());
			fTemplateHistogramContainer.push_back(hist);
			fTemplateNameContainer.push_back(name);
		}
		else
			return 0;
	}

	// add a parameter for the fraction
	int index = int(fTemplateHistogramContainer.size()) - 1;

	if (Nmax == 0)
	{
		double sum = fDataHistogram.Integral();
		if (sum > 0. && sum < 10.)
			Nmax = 20;
		else if (sum > 0. && sum >= 10.)
			Nmax = sum + 5.0 * sqrt(sum);
		else
			Nmax = 5.;
	}

	AddParameter(Form("N_%i", index), Nmin, Nmax);

	// successfully added histogram to container
	return 1;
}

// ---------------------------------------------------------
int StackModel::CompareHistogramProperties(TH1D hist1, TH1D hist2)
{
	// compare number of bins
	if (hist1.GetNbinsX() != hist2.GetNbinsX())
		return 0;

	// compare minimum x-values
	if (hist1.GetXaxis()->GetXmin() != hist2.GetXaxis()->GetXmin())
		return 0;

	// compare maximum x-values
	if (hist1.GetXaxis()->GetXmax() != hist2.GetXaxis()->GetXmax())
		return 0;

	// conclusion: they have the same properties
	return 1;
}

// ---------------------------------------------------------
void StackModel::PrintStack(const char * filename, const char * options)
{
	int nbins = fDataHistogram.GetNbinsX();
	int ntemplates = int( fTemplateHistogramContainer.size() );

	// check options
	bool flag_legend = false;
	bool flag_error0 = false; // symm. poisson error for data
	bool flag_error1 = false; // symm. poisson error for exp.
	bool flag_error2 = false; // asymm. poisson error of expectation value
	bool flag_error3 = false; // asymm. poisson error of expected no. of events

	if (std::string(options).find("L") < std::string(options).size())
		flag_legend = true;

	if (std::string(options).find("E0") < std::string(options).size())
		flag_error0 = true;

	if (std::string(options).find("E1") < std::string(options).size())
		flag_error1 = true;

	if (std::string(options).find("E2") < std::string(options).size())
		flag_error2 = true;

	if (std::string(options).find("E3") < std::string(options).size())
		flag_error3 = true;

	// create canvas
	TCanvas c1("c1");
	c1.cd();

	// draw data
	double ymin = 0.0;
	double ymax = 1.1 * fDataHistogram.GetMaximum();
	fDataHistogram.GetYaxis()->SetRangeUser(ymin, ymax);
	fDataHistogram.Draw("P");

	// create stack
	THStack stack("histostack","");

	// create legend
	TLegend legend(0.7, 0.7, 0.95, 0.95);
	legend.SetBorderSize(0);
	legend.SetFillColor(kWhite);

	// create a histogram with the sum of all contributions
 	TH1D * histsum = (TH1D*) fDataHistogram.Clone("temp");

	for (int i = 0; i < ntemplates; ++i)
	{
		TH1D histtemp = fTemplateHistogramContainer.at(i);
		fTemplateHistogramContainer.at(i).Scale(GetBestFitParameter(i) / histtemp.Integral());
		stack.Add(&(fTemplateHistogramContainer.at(i)));
		legend.AddEntry(&(fTemplateHistogramContainer.at(i)), fTemplateNameContainer.at(i).data(), "F");
	}

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double bincontent = 0;

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp)
			bincontent +=fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

		// set bin content
		histsum->SetBinContent(ibin, bincontent);
	}

	// define error graph
	TGraphAsymmErrors * graph_error_exp = new TGraphAsymmErrors(nbins);
	graph_error_exp->SetLineWidth(2);
	graph_error_exp->SetMarkerStyle(0);

	TGraphAsymmErrors * graph_error_obs = new TGraphAsymmErrors(nbins);
	graph_error_obs->SetMarkerStyle(0);

	// calculate uncertainty
	if (flag_error1)
		for (int i = 1; i <= nbins; ++i)
		{
			double nexp = histsum->GetBinContent(i);
			histsum->SetBinError(i, sqrt(nexp));
			histsum->SetMarkerStyle(0);
		}

	if (flag_error2)
		for (int i = 1; i <= nbins; ++i)
		{
			TH1D * proj = fUncertaintyHistogramExp->ProjectionY("_py", i, i);
			if (proj->Integral() > 0)
				proj->Scale(1.0 / proj->Integral());
			double quantiles[3];
			double sums[3] = {0.16, 0.5, 0.84};
			proj->GetQuantiles(3, quantiles, sums);
			graph_error_exp->SetPoint(i-1, histsum->GetBinCenter(i), quantiles[1]);
			graph_error_exp->SetPointError(i-1, 0.0, 0.0, quantiles[1] - quantiles[0], quantiles[2]-quantiles[1]);
			delete proj;
		}

	if (flag_error3)
		for (int i = 1; i <= nbins; ++i)
		{
			TH1D * proj = fUncertaintyHistogramObsPosterior->ProjectionY("_py", i, i);
			if (proj->Integral() > 0)
				proj->Scale(1.0 / proj->Integral());
			double quantiles[3];
			double sums[3] = {0.16, 0.5, 0.84};
			proj->GetQuantiles(3, quantiles, sums);
			graph_error_obs->SetPoint(i-1, histsum->GetBinCenter(i), quantiles[1]);
			graph_error_obs->SetPointError(i-1, 0.0, 0.0, quantiles[1] - TMath::Floor(quantiles[0]), TMath::Ceil(quantiles[2])-quantiles[1]);
			delete proj;
		}

	// draw histograms
	stack.Draw("SAMEA");
	stack.GetHistogram() -> SetXTitle("");
	stack.GetHistogram() -> SetYTitle("");
	stack.GetHistogram() -> GetXaxis() -> SetLabelSize(0);
	stack.GetHistogram() -> GetYaxis() -> SetLabelSize(0);
	stack.Draw("SAME");
	fDataHistogram.Draw("SAMEP");

	if (flag_error0)
		fDataHistogram.Draw("SAMEPE");

	if (flag_error1)
		histsum->Draw("SAMEE");

	if (flag_error3)
		graph_error_obs->Draw("SAMEZ");

	if (flag_error2)
		graph_error_exp->Draw("SAME||");

	if (flag_legend)
		legend.Draw();

	c1.Print(filename);

	// rescale
	for (int i = 0; i < ntemplates; ++i)
		fTemplateHistogramContainer.at(i).Scale(1.0 / fTemplateHistogramContainer.at(i).Integral());

	// delete temporary histograms
 	delete histsum;
	delete graph_error_exp;
	delete graph_error_obs;
}

// ---------------------------------------------------------
double StackModel::CalculateChi2()
{
	int nbins = fDataHistogram.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	std::vector <double> parameters = GetBestFitParameters();

	double chi2 = 0;

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double nexp = 0;
		double ndata = fDataHistogram.GetBinContent(ibin);

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp)
			nexp += parameters.at(itemp) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

		// add to chi2
		chi2 += (nexp - ndata) * (nexp - ndata) / nexp;
	}

	// return chi2
	return chi2;
}

// ---------------------------------------------------------
double StackModel::CalculateChi2Prob()
{
	double chi2 = CalculateChi2();
	int ndf = GetNDF();

	// return chi2 probability
	return TMath::Prob(chi2, ndf);
}

// ---------------------------------------------------------
double StackModel::CalculateMaxLike()
{
	// return maximum likelihood
	return Eval( GetBestFitParameters() );
}

// ---------------------------------------------------------
double StackModel::CalculateKSProb()
{
	// create a histogram with the sum of all contributions
	TH1 * histsum = (TH1D*)(fTemplateHistogramContainer.at(0)).Clone("temp");

	int nbins = fDataHistogram.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	std::vector <double> parameters = GetBestFitParameters();

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double bincontent = 0;

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp)
			bincontent += parameters.at(itemp) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

		// set bin content
		histsum->SetBinContent(ibin, bincontent);
	}

	// perform KS test
	double ksprob = histsum->KolmogorovTest(&fDataHistogram);

	// delete histogram
	delete histsum;

	return ksprob;
}

// ---------------------------------------------------------
double StackModel::CalculatePValue()
{
	// get best fit parameters
	std::vector<double> par = GetBestFitParameters();

	// check size of parameter vector
	if (par.size() != GetNParameters())
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCStackModel::CalculatePValueFast() : Number of parameters is inconsistent.");
		return -1;
	}

	// define temporary variables
	int nbins = fDataHistogram.GetNbinsX();
	int ntemplates = int( fTemplateHistogramContainer.size() );

	std::vector <int> histogram;
	std::vector <double> expectation;
	histogram.assign(nbins, 0);
	expectation.assign(nbins, 0);

	double logp = 0;
	double logp_start = 0;
	int counter_pvalue = 0;

	// define starting distribution
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double nexp = 0;

		for (int itemp = 0; itemp < ntemplates; ++itemp)
			nexp += par.at(itemp) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

		histogram[ibin-1]   = int(nexp);
		expectation[ibin-1] = nexp;

		// calculate p;
		logp += BCMath::LogPoisson(double(int(nexp)), nexp);
		logp_start += BCMath::LogPoisson(fDataHistogram.GetBinContent(ibin), nexp);
	}

	int niter = 100000;

	// loop over iterations
	for (int iiter = 0; iiter < niter; ++iiter)
	{
		// loop over bins
		for (int ibin = 0; ibin < nbins; ++ibin)
		{
			// random step up or down in statistics for this bin
			double ptest = fRandom->Rndm() - 0.5;

			// increase statistics by 1
			if (ptest > 0)
			{
				// calculate factor of probability
				double r = expectation.at(ibin) / double(histogram.at(ibin) + 1);

				// walk, or don't (this is the Metropolis part)
				if (fRandom->Rndm() < r)
				{
					histogram[ibin] = histogram.at(ibin) + 1;
					logp += log(r);
				}
			}

			// decrease statistics by 1
			else if (ptest <= 0 && histogram[ibin] > 0)
			{
				// calculate factor of probability
				double r = double(histogram.at(ibin)) / expectation.at(ibin);

				// walk, or don't (this is the Metropolis part)
				if (fRandom->Rndm() < r)
				{
					histogram[ibin] = histogram.at(ibin) - 1;
					logp += log(r);
				}
			}
		} // end of looping over bins

		// increase counter
		if (logp < logp_start)
			counter_pvalue++;

	} // end of looping over iterations

	// calculate p-value
	return double(counter_pvalue) / double(niter);
}

// ---------------------------------------------------------
void StackModel::MCMCUserIterationInterface()
{
	int nbins = fDataHistogram.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double bincontent = 0;

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp)
			bincontent += fMCMCx.at(itemp) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

		// set bin content
		fUncertaintyHistogramExp->Fill(fDataHistogram.GetBinCenter(ibin), bincontent);

		// loop over bins in the other direction
		int nbinsy = fUncertaintyHistogramObsPosterior->GetNbinsY();
		for (int jbin = 1; jbin <= nbinsy; ++jbin)
		{
			int n = jbin - 1;
			if (fabs(n - bincontent) < 2*sqrt(bincontent))
				fUncertaintyHistogramObsPosterior->Fill(fDataHistogram.GetBinCenter(ibin), n, TMath::Poisson(bincontent, n));
		}
	}
}

// ---------------------------------------------------------
void StackModel::PrintTemp()
{
	TCanvas * c1 = new TCanvas("c1");

	c1->cd();
	fUncertaintyHistogramExp->Draw("COL");
	c1->Print("uncertainty_exp.eps");

	c1->cd();
	fUncertaintyHistogramObsPosterior->Draw("COL");
	c1->Print("uncertainty_obs.eps");

	delete c1;
}

// ---------------------------------------------------------

