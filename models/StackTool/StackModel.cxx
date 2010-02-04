#include <iostream>

#include "StackModel.h"

#include <TROOT.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TPostScript.h>
#include <TPad.h> 
#include <TLine.h> 

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>

#include <iostream>

// ---------------------------------------------------------
StackModel::StackModel() : BCModel()
{
	fHistNorm = 0; 
	fFlagFixNorm = 0;
	fFlagPhysicalLimits = true;
	fNorm = 1; 
	fUncertaintyHistogramExp = 0;
	fUncertaintyHistogramObsPosterior = 0;
}

// ---------------------------------------------------------
StackModel::StackModel(const char * name) : BCModel(name)
{
	fHistNorm = 0; 
	fFlagFixNorm = 0; 
	fFlagPhysicalLimits = true;
	fNorm = 1; 
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

	if (fHistNorm)
		delete fHistNorm; 
}

// ---------------------------------------------------------
double StackModel::LogLikelihood(std::vector <double> parameters)
{
	double logprob = 0.;

	int nbins      = fDataHistogram.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	// define norm and efficiencies 
	std::vector <double> norm; 
	std::vector <double> eff; 

	// set norm and efficiency 
	for (int i = 0; i < ntemplates; ++i) {
		norm.push_back(parameters.at(2*i)); 
		eff.push_back(parameters.at(2*i+1)); 
	}

	// calculate norm
	double sum = 0; 
	for (int i = 0; i < ntemplates; ++i)
		sum += norm.at(i); 

	// fix norm
	if (fFlagFixNorm) {
		for (int i = 0; i < ntemplates; ++i) {
			norm[i] = norm.at(i) / sum; 
		}
	}
	else {
		fNorm = sum; 
	}		

	// loop over bins
	for (int ibin = 1; ibin <= nbins; ++ibin) {
		double nexp = 0;
		double ndata = fDataHistogram.GetBinContent(ibin);
		
		for (int itemp = 0; itemp < ntemplates; ++itemp) {
			nexp += norm.at(itemp) 
				* eff.at(itemp) 
				* fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);
		}		
		
		// check that expectation is larger or equal to zero
		if (nexp < 0)
			nexp = 0; 

		logprob += BCMath::LogPoisson(ndata, nexp);
	}

	return logprob;
}

// ---------------------------------------------------------
double StackModel::LogAPrioriProbability(std::vector <double> parameters)
{
	double logprob = 0.;

	int ntemplates = int(fTemplateHistogramContainer.size());

	// define norm and efficiencies 
	std::vector <double> norm; 
	std::vector <double> eff; 

	// set norm and efficiency 
	for (int i = 0; i < ntemplates; ++i) {
		norm.push_back(parameters.at(2*i)); 
		eff.push_back(parameters.at(2*i+1)); 
	}

	// loop over templates
	for (int i = 0; i < ntemplates; ++i) {
		// efficiency
		if (fTemplateEffErr.at(i)> 0)
			logprob += BCMath::LogGaus(eff.at(i), fTemplateEff.at(i), fTemplateEffErr.at(i)); 

		// norm
		if (fTemplatePriorSigma.at(i) > 0)
			logprob += BCMath::LogGaus(norm.at(i), fTemplatePriorMean.at(i), fTemplatePriorSigma.at(i)); 
	}

	// constraints
	int nconstraints = int(fConstraintSumIndices.size()); 
	if (nconstraints > 0) {

		// loop over constraints
		for (int i = 0; i < nconstraints; ++i) {

			// initialize sum
			double sum = 0; 

			// get number of summands
			int nsummands = int( (fConstraintSumIndices.at(i)).size() ); 

			// loop over summands and add to sum
			for (int j = 0; j < nsummands; ++j) {
				sum += norm.at( (fConstraintSumIndices.at(i)).at(j) ); 
			}

			// add to prior
			logprob += BCMath::LogGaus(sum, fConstraintSumMean.at(i), fConstraintSumRMS.at(i)); 
		}
	}

	return logprob;
}

// ---------------------------------------------------------
int StackModel::SetDataHistogram(TH1D hist)
{
	// create histogram
	int nbins = hist.GetNbinsX();
	double xmin = (hist.GetXaxis())->GetXmin(); 
	double xmax = (hist.GetXaxis())->GetXmax(); 
	fDataHistogram = TH1D("", "", nbins, xmin, xmax); 

	for (int i = 1; i <= nbins; ++i) 
		fDataHistogram.SetBinContent(i, hist.GetBinContent(i)); 

	// set histogram style
	fDataHistogram.SetXTitle((hist.GetXaxis())->GetTitle()); 
	fDataHistogram.SetYTitle((hist.GetYaxis())->GetTitle()); 
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

	// calculate norm
	fNorm = hist.Integral();

	// create histogram containing the normalization
	xmin = 0; 
	xmax = 0; 
	int ntemplates = (fTemplateHistogramContainer.size());

	// calculate the limits on the norm from the sum of all parameter limits 
	for (int i = 0; i < ntemplates; ++i) {
		BCParameter * par = this->GetParameter(2*i); 
		xmin += par->GetLowerLimit(); 
		xmax += par->GetUpperLimit(); 
	}

	// create new histogram for norm
	fHistNorm = new TH1D("", ";N_{norm};dN/dN_{norm}", 100, xmin, xmax); 

	// no errors
	return 1;
}

// ---------------------------------------------------------
int StackModel::AddTemplateHistogram(TH1D hist, const char * name, double Nmin, double Nmax)
{
	// check if histogram if filled
	if (hist.Integral() <= 0.)
		return 0;

	// compare template properties with data
	if (CompareHistogramProperties(fDataHistogram, hist) != 1)
		return 0;

	// check if prior makes sense 
	if (fFlagPhysicalLimits && Nmin < 0)
		Nmin = 0; 

	if (Nmin > Nmax)
		return 0; 

	// get number of templates 
	int ntemplates = int(fTemplateHistogramContainer.size());
	
	// set histogram color and style 
	hist.SetFillColor( 2 + ntemplates );
	hist.SetFillStyle(1001);
	
	// scale histogram
	hist.Scale(1.0 / hist.Integral());

	// check if template is consistent with other templates 
	if (ntemplates > 0)
		if (!CompareHistogramProperties(fTemplateHistogramContainer.at(0), hist))
			return 0; 

	// add template
	fTemplateHistogramContainer.push_back(hist);
	fTemplateNameContainer.push_back(name);
	fTemplateEff.push_back(1.0); 
	fTemplateEffErr.push_back(0.0); 
	fTemplatePriorMean.push_back(0.0);
	fTemplatePriorSigma.push_back(-1.0);

	// add a parameter for the expectation value of the process
	int index = ntemplates;
	
	if (Nmax == 0) {
		double sum = fDataHistogram.Integral();
		if (sum > 0. && sum < 10.)
			Nmax = 20;
		else if (sum > 0. && sum >= 10.)
			Nmax = sum + 5.0 * sqrt(sum);
		else
			Nmax = 5.;
		if (!fFlagPhysicalLimits)
			Nmin = - Nmax;
		else
			Nmin = 0; 
	}

	// a a parameter for the expectation value 
	AddParameter(Form("N_%i", index), Nmin, Nmax);

	// add a parameter for the efficiency
	AddParameter(Form("eff_%i", index), 0.0, 1.0);

	// add prior histogram
	TH1D hist_prior(Form("prior_%i", index), "", 1, Nmin, Nmax); 
	hist_prior.SetBinContent(1, 1); 
	fHistPrior.push_back(hist_prior); 

	// successfully added histogram to container
	return 1;
}


// ---------------------------------------------------------
int StackModel::CalculateRatio(int index, std::vector<int> indices)
{
	// get number of templates
	int ntemplates = int( fTemplateHistogramContainer.size() );
	
	// check index
	if (index < 0 || index >= ntemplates) {
		return 0; 
	} 

	// check indices
	for (int i = 0; i < int(indices.size()); ++i) {
		if (indices.at(i) < 0 || indices.at(i) >= ntemplates) {
			return 0; 
		}
	}

	// create temporary vector
	std::vector<int> tempvec; 
	tempvec.push_back(index); 
	for (int i = 0; i < int(indices.size()); ++i)
		tempvec.push_back(indices.at(i)); 

	// add ratio
	fIndicesRatios1D.push_back(tempvec); 
	
	// get number of ratios
	int nratios = int(fHistRatios1D.size());

	// create histogram
	double fmin = 0.0;
	double fmax = 1.0;
	if (!fFlagPhysicalLimits) {
		fmin = -1.0;
		fmax =  2.0;
	}

	TH1D hist_ratio1d(Form("ratio %i", nratios), ";;", 100, fmin, fmax); 
	hist_ratio1d.SetXTitle(Form("r", nratios)); 
	hist_ratio1d.SetYTitle(Form("p(r|data)", nratios)); 
	fHistRatios1D.push_back(hist_ratio1d); 

	// no error
	return 1;
}

// ---------------------------------------------------------
int StackModel::AddTemplateHistogram(TH1D hist, const char * name, TH1D prior)
{
	// check if prior histogram exists
	if (prior.Integral() <= 0) {
		return 0; 
	}

	// get boundaries
	double nmin = prior.GetXaxis()->GetXmin();
	double nmax = prior.GetXaxis()->GetXmax(); 

	// add template 
	int err = this->AddTemplateHistogram(hist, name, nmin, nmax); 
	
	if (!err) {
		return err; 
	}

	// scale prior histogram
	prior.Scale(1.0/prior.Integral()); 

	// set new prior histogram	
	fHistPrior[int(fHistPrior.size())-1] = prior; 

	// no error
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
	bool flag_diff   = false; // plot difference between data and expectation below stack plot

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

	if (std::string(options).find("D") < std::string(options).size())
		flag_diff = true;

	// create canvas
	TCanvas* c1 = new TCanvas("c1", "", 700, 700);
	c1->cd();
	TPad * pad1; 
	TPad * pad2; 

	double fraction_pads = 0.3; 

	if(!flag_diff)
		fraction_pads=0.0;

	if (flag_diff) {
		pad1 = new TPad("pad1", "", 0.0, fraction_pads, 1.0, 1.0); 
		pad1->SetTopMargin   (0.13/(1.0-fraction_pads));
		pad1->SetBottomMargin(0.0);
    pad1->SetLeftMargin  (0.15);
		pad1->SetRightMargin (0.13);
		pad1->SetFillColor(kWhite); 
		pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, fraction_pads); 
		pad2->SetTopMargin   (0.0);
		pad2->SetBottomMargin(0.15 / fraction_pads);
    pad2->SetLeftMargin  (0.15);
		pad2->SetRightMargin (0.13);
		pad2->SetFillColor(kWhite); 
		pad1->Draw();
		pad2->Draw();
	}
	else {
		pad1 = new TPad("pad1", "",0.0, 0.0, 1.0, 1.0); 
		pad1->SetFillColor(kWhite); 
		pad2 = new TPad(); 
		pad1->Draw(); 
	}

	pad1->cd(); 

	// set style and draw data
	double ymin = 0.01;
	double ymax = 1.1 * (fDataHistogram.GetMaximum() + sqrt(fDataHistogram.GetMaximum()));
	fDataHistogram.GetYaxis()->SetRangeUser(ymin, ymax);
	fDataHistogram.GetXaxis()->SetNdivisions(505); 
	if (flag_diff) {
		fDataHistogram.GetXaxis()->SetLabelSize(fDataHistogram.GetXaxis()->GetLabelSize()/(1.0-fraction_pads)); 
		fDataHistogram.GetXaxis()->SetLabelOffset(fDataHistogram.GetXaxis()->GetLabelOffset()*(1.0-fraction_pads)); 
		fDataHistogram.GetXaxis()->SetTitleSize(fDataHistogram.GetXaxis()->GetTitleSize()/(1.0-fraction_pads)); 
		fDataHistogram.GetXaxis()->SetTitleOffset(fDataHistogram.GetXaxis()->GetTitleOffset()*(1.0-fraction_pads)); 
		fDataHistogram.GetYaxis()->SetLabelSize(fDataHistogram.GetYaxis()->GetLabelSize()/(1.0-fraction_pads)); 
		fDataHistogram.GetYaxis()->SetLabelOffset(fDataHistogram.GetYaxis()->GetLabelOffset()/(fraction_pads)); 
		fDataHistogram.GetYaxis()->SetTitleSize(fDataHistogram.GetYaxis()->GetTitleSize()/(1.0-fraction_pads)); 
		fDataHistogram.GetYaxis()->SetTitleOffset(fDataHistogram.GetYaxis()->GetTitleOffset()*(1.0-fraction_pads)); 
	}
	fDataHistogram.Draw("P");	

	// create a histogram with the sum of all contributions
 	TH1D * histsum = (TH1D*) fDataHistogram.Clone("temp");
	
	// create stack
	THStack stack("histostack","");

	// create legends
	TLegend* legend1; 
	TLegend* legend2;

	if (flag_diff)
		legend1 = new TLegend(0.15, (0.88-fraction_pads)/(1-fraction_pads), 0.50, 0.99);
	else
		legend1 = new TLegend(0.15, 0.88, 0.50, 0.99);
	legend1->SetBorderSize(0);
	legend1->SetFillColor(kWhite);
	legend1->AddEntry(&fDataHistogram, "Data", "LEP");
	legend1->AddEntry(&fDataHistogram, "Total expected uncertainty", "LE");

	double y = 0.99; 
	if (ntemplates > 2 && ntemplates <7)
		y -= 0.11 / 4. * double(ntemplates - 2); 
	legend2 = new TLegend(0.50,(y-fraction_pads)/(1-fraction_pads) , 0.85, 0.99);
	legend2->SetBorderSize(0);
	legend2->SetFillColor(kWhite);

	// scale histograms and add to stack and legend
	for (int i = 0; i < ntemplates; ++i)
	{
		TH1D histtemp = fTemplateHistogramContainer.at(i);
		fTemplateHistogramContainer.at(i).Scale(GetBestFitParameter(2*i) * GetBestFitParameter(2*i+1) / histtemp.Integral());
		stack.Add(&(fTemplateHistogramContainer.at(i)));
		if (i < 2)
			legend1->AddEntry(&(fTemplateHistogramContainer.at(i)), fTemplateNameContainer.at(i).data(), "F");
		else if (i < 6)
			legend2->AddEntry(&(fTemplateHistogramContainer.at(i)), fTemplateNameContainer.at(i).data(), "F");
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

	// create difference histogram
	TH1D* hist_diff = 0; 

	TGraphAsymmErrors * graph_diff_exp = 0; 
	
	if (flag_diff) {
		ymin = 0;
		ymax = 0; 
		hist_diff = new TH1D("hist_diff", "", nbins, histsum->GetXaxis()->GetXmin(), histsum->GetXaxis()->GetXmax() ); 
		hist_diff->GetXaxis()->SetTitle(fDataHistogram.GetXaxis()->GetTitle()); 
		hist_diff->GetYaxis()->SetTitle("#Delta N"); 
		hist_diff->GetXaxis()->SetNdivisions(505); 
		hist_diff->GetXaxis()->SetLabelSize(hist_diff->GetXaxis()->GetLabelSize()/(fraction_pads)); 
		hist_diff->GetXaxis()->SetLabelOffset(hist_diff->GetXaxis()->GetLabelOffset()/fraction_pads*2.); 
		hist_diff->GetXaxis()->SetTitleSize(hist_diff->GetXaxis()->GetTitleSize()/(fraction_pads)); 
		hist_diff->GetXaxis()->SetTitleOffset((hist_diff->GetXaxis()->GetTitleOffset()-(1.0-fraction_pads))/(fraction_pads)); 
		hist_diff->GetYaxis()->SetNdivisions(503); 
		hist_diff->GetYaxis()->SetLabelSize(hist_diff->GetYaxis()->GetLabelSize()/(fraction_pads)); 
		hist_diff->GetYaxis()->SetLabelOffset(hist_diff->GetYaxis()->GetLabelOffset()/(fraction_pads)); 
		hist_diff->GetYaxis()->SetTitleSize(hist_diff->GetYaxis()->GetTitleSize()/(fraction_pads)); 
		hist_diff->GetYaxis()->SetTitleOffset(hist_diff->GetYaxis()->GetTitleOffset()*(fraction_pads)); 
		hist_diff->SetStats(kFALSE); 

		graph_diff_exp = new TGraphAsymmErrors(nbins);
		graph_diff_exp->SetLineWidth(2);
		graph_diff_exp->SetMarkerStyle(0);
		graph_diff_exp->SetFillColor(kYellow);
		for (int i = 0; i < nbins; ++i) {
			hist_diff->SetBinContent(i+1, fDataHistogram.GetBinContent(i+1)-histsum->GetBinContent(i+1)); 
			hist_diff->SetBinError(i+1, fDataHistogram.GetBinError(i+1)); 
			graph_diff_exp->SetPoint(i, (graph_error_exp->GetX())[i], 0.0);
			graph_diff_exp->SetPointEXlow(i, 0.0); 
			graph_diff_exp->SetPointEXhigh(i, 0.0); 
			graph_diff_exp->SetPointEYlow(i, (graph_error_exp->GetEYlow())[i]); 
			graph_diff_exp->SetPointEYhigh(i, (graph_error_exp->GetEYhigh())[i]); 

			if (-(graph_error_exp->GetEYlow())[i] < ymin)
				ymin = -(graph_error_exp->GetEYlow())[i]; 
			if ((graph_error_exp->GetEYhigh())[i] > ymax)
				ymax = (graph_error_exp->GetEYhigh())[i]; 
		}
		if (ymax < (hist_diff->GetMaximum() + hist_diff->GetBinError(hist_diff->GetMaximumBin())))
			ymax = 1.1 * (hist_diff->GetMaximum() + hist_diff->GetBinError(hist_diff->GetMaximumBin())); 
		if (ymin>(hist_diff->GetMinimum() - hist_diff->GetBinError(hist_diff->GetMaximumBin())))
			ymin = 1.1 * (hist_diff->GetMinimum() - hist_diff->GetBinError(hist_diff->GetMaximumBin()));
		(hist_diff->GetYaxis())->SetRangeUser(-1.1*TMath::Max(-ymin, ymax), 1.1*TMath::Max(-ymin, ymax));
		
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

	if (flag_legend) {
		legend1->Draw();
		if (ntemplates>2)
		legend2->Draw();
	}

	TLine * line = 0; 
	if (flag_diff) {
		pad2->cd(); 
		hist_diff->Draw("P"); 
		graph_diff_exp->Draw("SAME4");
		line = new TLine((hist_diff->GetXaxis())->GetXmin(), 0.0, (hist_diff->GetXaxis())->GetXmax(), 0.0); 
		line->SetLineWidth(2); 
		line->SetLineColor(kBlack); 
		line->Draw("SAME"); 
		hist_diff->Draw("SAMEP"); 
	}

	c1->Print(filename);

	// rescale
	for (int i = 0; i < ntemplates; ++i)
		fTemplateHistogramContainer.at(i).Scale(1.0 / fTemplateHistogramContainer.at(i).Integral());

	// delete temporary histograms
	delete pad1;
	delete pad2;
 	delete c1;
	delete legend1;
 	delete legend2;
 	delete graph_error_exp;
 	delete graph_error_obs;
 	delete histsum;
	if (flag_diff) {
 		delete hist_diff;
 		delete graph_diff_exp;
 		delete line; 
	}
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
			nexp += parameters.at(2*itemp) * parameters.at(2*itemp+1) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

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
			bincontent += parameters.at(2*itemp) * parameters.at(2*itemp+1) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

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
			nexp += par.at(2*itemp) * par.at(2*itemp+1) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

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
	int nbins      = fDataHistogram.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double bincontent = 0;

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp)
			bincontent += fMCMCx.at(2 * itemp) * fMCMCx.at(2 * itemp + 1) * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);

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

	// fill normalization
	fHistNorm->Fill(fNorm);

	// fill ratios
	int nratios = int( fIndicesRatios1D.size() ); 

	// loop over fractions to fill 
	for (int i = 0; i < nratios; ++i) {
		int nsum = int( (fIndicesRatios1D.at(i)).size() ) - 1; 
		double sum = 0; 
		for (int j = 1; j <= nsum; ++j) {
			int indexsum = fIndicesRatios1D.at(i).at(j); 
			sum += fMCMCx.at(2*indexsum); 
		}
		fHistRatios1D.at(i).Fill(fMCMCx.at(2*fIndicesRatios1D.at(i).at(0))/sum); 
	}

}

// ---------------------------------------------------------
void StackModel::PrintRatios(const char * filename)
{
	int nratios = int(fHistRatios1D.size());

	TCanvas* c1 = new TCanvas("c1");

	TPostScript * ps = new TPostScript(filename, 112);
	ps->NewPage(); 

	c1->cd(); 
	BCH1D* h1temp = new BCH1D(fHistNorm); 
	h1temp->Draw(); 
	c1->Update(); 
	ps->NewPage(); 
	for (int i = 0; i < nratios; ++i) {
		c1->Update(); 
		ps->NewPage();
		c1->cd(); 
		BCH1D* h1temp = new BCH1D(&fHistRatios1D.at(i));
		h1temp->Draw(); 
	}
	c1->Update();
	ps->Close(); 

	delete c1; 
	delete ps; 
}

// ---------------------------------------------------------
int StackModel::SetTemplateEfficiency(int index, double eff, double err, bool adjust)
{
	// get number of templates
	int ntemplates = int(fTemplateHistogramContainer.size());

	if (index < 0 || index >= ntemplates)
		return 0; 

	// set efficiency and error
	fTemplateEff[index] = eff; 
	fTemplateEffErr[index] = err;

	if (err < 0) {
		int parindex = index*2+1; 
		double effmin = eff;
		double effmax = eff;
		this -> GetParameter(parindex) -> SetLowerLimit(effmin); 
		this -> GetParameter(parindex) -> SetUpperLimit(effmax); 
		fMCMCBoundaryMin[parindex] = eff; 
		fMCMCBoundaryMax[parindex] = eff; 
	}		

	else if (adjust) {
		int parindex = index*2+1; 
		double effmin = TMath::Max(0.0, eff - 5.0*err); 
		double effmax = TMath::Min(1.0, eff + 5.0*err); 
		this -> GetParameter(parindex) -> SetLowerLimit(effmin); 
		this -> GetParameter(parindex) -> SetUpperLimit(effmax); 
		fMCMCBoundaryMin[parindex] = effmin; 
		fMCMCBoundaryMax[parindex] = effmax; 
	}

	// no error 
	return 1;
}

// ---------------------------------------------------------
int StackModel::SetTemplatePrior(int index, double mean, double sigma, bool adjust)
{
	// get number of templates
	int ntemplates = int(fTemplateHistogramContainer.size());

	if (index < 0 || index >= ntemplates)
		return 0; 

	// set efficiency and error
	fTemplatePriorMean[index] = mean; 
	fTemplatePriorSigma[index] = sigma;

	if (adjust) {
		int parindex = index*2; 
		double parmin;
		double parmax; 
		if (fFlagPhysicalLimits) 
			parmin = TMath::Max(0.0, mean - 5.0*sigma); 
		else
			parmin = mean - 5.0 * sigma; 
			parmax = mean + 5.0*sigma; 
		this -> GetParameter(parindex) -> SetLowerLimit(parmin); 
		this -> GetParameter(parindex) -> SetUpperLimit(parmax); 
		fMCMCBoundaryMin[parindex] = parmin; 
		fMCMCBoundaryMax[parindex] = parmax; 
	}

	// no error 
	return 1;
}

// ---------------------------------------------------------
int StackModel::ConstrainSum(std::vector <int> indices, double mean, double rms)
{
	// add contraint to container(s)
	fConstraintSumIndices.push_back(indices); 
	fConstraintSumMean.push_back(mean); 
	fConstraintSumRMS.push_back(rms); 

	// no error
	return 1;
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

