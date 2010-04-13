#include <iostream>

#include "TemplateModel.h"

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
TemplateModel::TemplateModel() : BCModel()
															 , fFlagFixNorm(false)
															 , fFlagPhysicalLimits(true)
															 , fUncertaintyHistogramExp(0)
															 , fUncertaintyHistogramObsPosterior(0)
															 , fNorm(-1)
															 , fNBins(-1)
															 , fXmin(1.)
															 , fXmax(0.)
															 , fPriorNBins(1000)
{
}

// ---------------------------------------------------------
TemplateModel::TemplateModel(const char * name) : BCModel(name)																									
															 , fFlagFixNorm(false)
															 , fFlagPhysicalLimits(true)
															 , fUncertaintyHistogramExp(0)
															 , fUncertaintyHistogramObsPosterior(0)
															 , fNorm(-1)
															 , fNBins(-1)
															 , fXmin(1.)
															 , fXmax(0.)
															 , fPriorNBins(1000)
{
}

// ---------------------------------------------------------
TemplateModel::~TemplateModel()
{
	if (fUncertaintyHistogramExp)
		delete fUncertaintyHistogramExp;

	if (fUncertaintyHistogramObsPosterior)
		delete fUncertaintyHistogramObsPosterior;
}

// ---------------------------------------------------------
double TemplateModel::LogLikelihood(std::vector <double> parameters)
{
	double logprob = 0.;

	// get number pf templates
	int ntemplates = GetNTemplates();

	// get number of sources of systematic uncertainties
	int nsyst = GetNSystErrors(); 
	
	// loop over bins
	for (int ibin = 1; ibin <= fNBins; ++ibin) {
		double nexp = 0;
		double ndata = fHistData.GetBinContent(ibin);

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp) {
			int templateindex = fTemplateParIndexContainer.at(itemp); 
			int effindex = fEffParIndexContainer.at(itemp);

			// get efficiency for the bin
			double efficiency = fEffHistogramContainer.at(itemp).GetBinContent(ibin); 

			// modify efficiency by uncertainty
			double efferr = fEffErrHistogramContainer.at(itemp).GetBinContent(ibin); 

			// check efficiency error 
			if (efferr > 0) 
				efficiency += parameters.at(effindex) * efferr; 

			// loop over sources of systematic uncertainties
			for (int isyst = 0; isyst < nsyst; ++isyst) {
				// get parameter index
				int systindex = fSystErrorParIndexContainer.at(isyst);

				// add efficiency
				double deff = fSystErrorHistogramContainer.at(isyst).at(itemp).GetBinContent(ibin); 

				if (deff > 0)
					efficiency += deff * parameters.at(systindex); 
			}

			// make sure efficiency is positive
			if (efficiency < 0.)
				efficiency = 0.; 

			// calculate expectation nvalue 
			nexp += parameters.at(templateindex) 
				* efficiency
				* fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);
		}		
		
		// check that expectation is larger or equal to zero
		if (nexp < 0) {
			BCLog::OutWarning("TemplateModel::LogLikelihood : Expectation value smaller than 0. Force it to be 0.");
			nexp = 0; 
		}

		// add Poisson term
		logprob += BCMath::LogPoisson(ndata, nexp);
	}

	// return log likelihood
	return logprob;
}

// ---------------------------------------------------------
double TemplateModel::LogAPrioriProbability(std::vector <double> parameters)
{
 	double logprob = 0.;

	// get number of templates
	int ntemplates = GetNTemplates();

	// loop over templates
	for (int i = 0; i < ntemplates; ++i) {

		// prior on process contributions
		double par = parameters.at(fTemplateParIndexContainer.at(i)); 
		int bin = fPriorContainer.at(i).FindBin(par);
		logprob += log( fPriorContainer.at(i).GetBinContent(bin) );

		// prior on efficiences
		par = parameters.at(fEffParIndexContainer.at(i)); 
		logprob += BCMath::LogGaus(par, 0.0, 1.0);
	}
	
	// get number of sources of systematic uncertainties
	int nsyst = GetNSystErrors(); 

	// loop over sources of systematic uncertainties
	for (int i = 0; i < nsyst; ++i) {
		double par = parameters.at(fSystErrorParIndexContainer.at(i)); 
		if (fSystErrorTypeContainer.at(i) == "gauss")
			logprob += BCMath::LogGaus(par, 0.0, 1.0);
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
				int index = fConstraintSumIndices.at(i).at(j);
				double par = parameters.at(fTemplateParIndexContainer.at(index)); 
				sum += par; 
 			}

 			// add to prior
 			logprob += BCMath::LogGaus(sum, fConstraintSumMean.at(i), fConstraintSumRMS.at(i)); 
 		}
 	}
		
	return logprob;
}

// ---------------------------------------------------------
int TemplateModel::SetData(const TH1D& hist)
{
	// create histogram
	fNBins = hist.GetNbinsX();
	fXmin = (hist.GetXaxis())->GetXmin(); 
	fXmax = (hist.GetXaxis())->GetXmax(); 
	fHistData = TH1D("", "", fNBins, fXmin, fXmax); 

	// copy histogram content
	for (int i = 1; i <= fNBins; ++i) 
		fHistData.SetBinContent(i, hist.GetBinContent(i)); 

	// set histogram style
	fHistData.SetXTitle((hist.GetXaxis())->GetTitle()); 
	fHistData.SetYTitle((hist.GetYaxis())->GetTitle()); 
	fHistData.SetMarkerStyle(20);
	fHistData.SetMarkerSize(1.1);
	fHistData.SetStats(kFALSE);

	// calculate norm
	fNorm = hist.Integral();

	// no errors
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::AddTemplate(TH1D hist, const char * name, double Nmin, double Nmax)
{
	// check if histogram if filled
	if (hist.Integral() <= 0.) {
		BCLog::OutError("TemplateModel::AddTemplate : Normalization is zero or less than that."); 
		return 0;
	}

	// compare template properties with data
	if (CompareHistogramProperties(fHistData, hist) != 1) {
		BCLog::OutError("TemplateModel::AddTemplate : Data and template histogram properties are incompatible."); 
		return 0;
	}

	// check if prior makes sense 
	if (fFlagPhysicalLimits && Nmin < 0)
		Nmin = 0; 

	if (Nmin > Nmax) {
		BCLog::OutError("TemplateModel::AddTemplate : Lower limit exceeds upper limit."); 
		return 0; 
	}

	// get number of templates 
	int ntemplates = int(fTemplateHistogramContainer.size());
	
	// set histogram color and style 
	hist.SetFillColor(2 + ntemplates);
	hist.SetFillStyle(1001);
	hist.SetStats(kFALSE);

	// scale histogram
	hist.Scale(1.0 / hist.Integral());

	// check if template is consistent with other templates 
	if (ntemplates > 0)
		if (!CompareHistogramProperties(fTemplateHistogramContainer.at(0), hist)) {
			BCLog::OutError("TemplateModel::AddTemplate : Properties of template histogram is not compatible with older template histograms.");
			return 0; 
		}

	// histograms
	TH1D histprior = TH1D("", "", fPriorNBins, Nmin, Nmax);
	TH1D histsysterror = TH1D("", "", fNBins, fXmin, fXmax);

	// set style
	histprior.SetXTitle(name);

	// fill histograms
	for (int i = 1; i <= fPriorNBins; ++i) 
		histprior.SetBinContent(i, 1.);
	for (int i = 1; i <= fNBins; ++i) 
		histsysterror.SetBinContent(i, 0.);

	// get parameter index
	int parindex = GetNParameters(); 
	int partemplateindex = GetNTemplates(); 

	// add a parameter for the expectation value 
	AddParameter(Form("N_%i", partemplateindex), Nmin, Nmax);

	// get efficiency parameter index
	int effindex = GetNParameters();

	// add a parameter for the efficiency
	AddParameter(Form("eff_%i", partemplateindex), -5.0, 5.0);

	// add histogram, name and index to containers
	fTemplateHistogramContainer.push_back(hist);
	fPriorContainer.push_back(histprior);
	fTemplateNameContainer.push_back(name);
	fTemplateParIndexContainer.push_back(parindex); 
	fEffParIndexContainer.push_back(effindex); 
	
	// set efficiency histograms to one without uncertainty
	fEffHistogramContainer.push_back(TH1D());
	fEffErrHistogramContainer.push_back(TH1D());
	SetTemplateEfficiency(name, 1., 0.);

	// add systematic uncertainties
	for (int i = 0; i < GetNSystErrors(); ++i) {
		std::vector <TH1D> histvector = fSystErrorHistogramContainer.at(i); 
		histvector.push_back(histsysterror);
	}

	// successfully added histogram to container
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::SetTemplatePrior(const char* name, TH1D prior)
{
	// get index
	int parindex = GetIndexTemplate(name); 

	// check parameter index
	if (parindex < 0) {
		BCLog::OutError("TemplateModel::SetTemplatePrior : Did not find parameter."); 
		return 0; 
	}

	// check prior
	if (prior.Integral() <= 0) {
		BCLog::OutError("TemplateModel::SetTemplatePrior : Integral of prior is equal to zero or less than that.");
		return 0; 
	}

	// normalize prior to unity
	prior.Scale(1.0/prior.Integral());

	// replace histogram
	fPriorContainer[parindex] = prior; 

	// no error 
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::SetTemplatePrior(const char* name, double mean, double sigma)
{
	// get index
	int parindex = GetParIndexTemplate(name); 

	// check parameter index
	if (parindex < 0) {
		BCLog::OutError("TemplateModel::SetTemplatePrior : Did not find parameter."); 
		return 0; 
	}

	// get parameter
	BCParameter * par = this->GetParameter(parindex); 	

	// create new histogram
	TH1D hist("", "", fPriorNBins, par->GetLowerLimit(), par->GetUpperLimit()); 

	// loop over bins and fill histogram
	for (int i = 1; i < fPriorNBins; ++i) {
		double x = hist.GetBinCenter(i); 
		double fx = TMath::Gaus(x, mean, sigma);
		hist.SetBinContent(i, fx);
	}

	// set template prior
	int err = SetTemplatePrior(name, hist);

	// return error code
	return err; 
}

// ---------------------------------------------------------
int TemplateModel::AddSystError(const char* errorname, const char* errtype)
{
	// define parameter range
	double dx = 1.0;

	// check error type
	if (std::string(errtype) == std::string("gauss"))
		dx = 5.0; 
	else if (std::string(errtype) == std::string("flat"))
		dx = 1.0; 
	else {
		BCLog::OutError("TemplateModel::AddSystError : Unknown error type.");
		return 0;
	}

	// add a parameter for the expectation value 
	AddParameter(Form("systerr_%i", GetNSystErrors()), -dx, dx);

	// add name and index to containers
	fSystErrorNameContainer.push_back(errorname);
	fSystErrorParIndexContainer.push_back(GetNParameters()-1); 

	// create histogram
	TH1D hist = TH1D("", "", fNBins, fXmin, fXmax); 
	
	// fill histograms
	for (int i = 1; i <= fNBins; ++i) {
		hist.SetBinContent(i, 0.); 
	}

	// get number of templates
	int n = GetNTemplates(); 

	// define vector of histograms
	std::vector<TH1D> histvector; 

	// add histograms
	for (int i = 0; i < n; ++i) {
		histvector.push_back(hist); 
	}

	// add histogram vector
	fSystErrorHistogramContainer.push_back(histvector);

	// add error type to container
	fSystErrorTypeContainer.push_back(std::string(errtype));

	// no error 
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::SetTemplateSystError(const char* errorname, const char* templatename, TH1D parerror)
{
	// get error index
	int errindex = GetIndexSystError(errorname); 

	// check parameter index
	if (errindex < 0) {
		BCLog::OutError("TemplateModel::SetTemplateSystError : Did not find parameter."); 
		return 0; 
	}

	// get template index
	int tempindex = GetIndexTemplate(templatename); 

	// check index
	if (tempindex < 0) {
		BCLog::OutError("TemplateModel::SetTemplateSystError : Could not find template.");
		return 0;
	}

	// set style
	parerror.SetStats(kFALSE);

	// set histogram
	(fSystErrorHistogramContainer.at(errindex))[tempindex] = parerror;
	
	// no error 
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::Initialize()
{
	// check data integral
	if (fHistData.Integral() <= 0) {
		BCLog::OutError("TemplateModel::Initialize : Normalization of data histogram is zero or less than that."); 
		return 0;
	}

	// create histograms for uncertainty determination
	double maximum = 1.5 * fHistData.GetMaximum();

	fUncertaintyHistogramExp = new TH2D(Form("UncertaintyExp_%i", BCLog::GetHIndex()), "",
																			fHistData.GetNbinsX(),
																			fHistData.GetXaxis()->GetXmin(),
																			fHistData.GetXaxis()->GetXmax(),
																			100,
																			0.0,
																			maximum);
	
	fUncertaintyHistogramObsPosterior = new TH2D(Form("UncertaintyObsPosterior_%i", BCLog::GetHIndex()), "",
																							 fHistData.GetNbinsX(),
																							 fHistData.GetXaxis()->GetXmin(),
																							 fHistData.GetXaxis()->GetXmax(),
																							 int(maximum) + 1,
																							 -0.5,
																							 double(int(maximum))+0.5);

	// create histogram containing the normalization
	double xmin = 0; 
	double xmax = 0; 
	int ntemplates = int(fTemplateParIndexContainer.size());

	// calculate the limits on the norm from the sum of all parameter
	// limits
	for (int i = 0; i < ntemplates; ++i) {
		// get parameter index
		int parindex = fTemplateParIndexContainer.at(i);
		int effindex = fEffParIndexContainer.at(i);

		// get parameter 
		BCParameter * par = this->GetParameter(parindex); 
		BCParameter * eff = this->GetParameter(effindex); 

		// increate limits
		xmin += par->GetLowerLimit() * eff->GetLowerLimit(); 
		xmax += par->GetUpperLimit() * eff->GetUpperLimit(); 
	}

	// create new histogram for norm
	fHistNorm = TH1D("", ";N_{norm};dN/dN_{norm}", 100, xmin, xmax); 

	// no error 
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::CalculateRatio(int index, std::vector<int> indices, double rmin, double rmax)
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
	if (fFlagPhysicalLimits) {
		fmin = TMath::Max(rmin, 0.0); 
		fmax = TMath::Min(1.0, rmax); 
	}

	TH1D hist_ratio1d(Form("ratio %i", nratios), ";;", 100, fmin, fmax); 
	hist_ratio1d.SetXTitle(Form("r", nratios)); 
	hist_ratio1d.SetYTitle(Form("p(r|data)", nratios)); 
	fHistRatios1D.push_back(hist_ratio1d); 

	// no error
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::CompareHistogramProperties(TH1D hist1, TH1D hist2)
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
void TemplateModel::PrintStack(const char * filename, const char * options)
{
	int nbins = fHistData.GetNbinsX();
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
	double ymax = 1.1 * (fHistData.GetMaximum() + sqrt(fHistData.GetMaximum()));
	fHistData.GetYaxis()->SetRangeUser(ymin, ymax);
	fHistData.GetXaxis()->SetNdivisions(505); 
	if (flag_diff) {
		fHistData.GetXaxis()->SetLabelSize(fHistData.GetXaxis()->GetLabelSize()/(1.0-fraction_pads)); 
		fHistData.GetXaxis()->SetLabelOffset(fHistData.GetXaxis()->GetLabelOffset()*(1.0-fraction_pads)); 
		fHistData.GetXaxis()->SetTitleSize(fHistData.GetXaxis()->GetTitleSize()/(1.0-fraction_pads)); 
		fHistData.GetXaxis()->SetTitleOffset(fHistData.GetXaxis()->GetTitleOffset()*(1.0-fraction_pads)); 
		fHistData.GetYaxis()->SetLabelSize(fHistData.GetYaxis()->GetLabelSize()/(1.0-fraction_pads)); 
		fHistData.GetYaxis()->SetLabelOffset(fHistData.GetYaxis()->GetLabelOffset()/(fraction_pads)); 
		fHistData.GetYaxis()->SetTitleSize(fHistData.GetYaxis()->GetTitleSize()/(1.0-fraction_pads)); 
		fHistData.GetYaxis()->SetTitleOffset(fHistData.GetYaxis()->GetTitleOffset()*(1.0-fraction_pads)); 
	}
	fHistData.Draw("P");	

	// create a histogram with the sum of all contributions
 	TH1D * histsum = (TH1D*) fHistData.Clone("temp");
	
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
	legend1->AddEntry(&fHistData, "Data", "LEP");
	legend1->AddEntry(&fHistData, "Total expected uncertainty", "LE");

	double y = 0.99; 
	if (ntemplates > 2 && ntemplates <7)
		y -= 0.11 / 4. * double(ntemplates - 2); 
	legend2 = new TLegend(0.50,(y-fraction_pads)/(1-fraction_pads) , 0.85, 0.99);
	legend2->SetBorderSize(0);
	legend2->SetFillColor(kWhite);

	// scale histograms and add to stack and legend
	for (int itemp = 0; itemp < ntemplates; ++itemp)
	{
		int tempindex = fTemplateParIndexContainer.at(itemp); 
		int effindex = fEffParIndexContainer.at(itemp); 

		// scale histogram
		fTemplateHistogramContainer.at(itemp).Scale(GetBestFitParameter(tempindex) 
																										/ fTemplateHistogramContainer.at(itemp).Integral());

		// loop over bins and scale these
		for (int ibin = 1; ibin <= fNBins; ++ibin) {
			// get efficiency for the bin
			double efficiency = fEffHistogramContainer.at(itemp).GetBinContent(ibin); 

			// modify efficiency by uncertainty
			double efferr = fEffErrHistogramContainer.at(itemp).GetBinContent(ibin); 

			// check efficiency error 
			if (efferr > 0) 
				efficiency = TMath::Max(0., efficiency + GetBestFitParameter(effindex) * efferr); 

			fTemplateHistogramContainer.at(itemp).SetBinContent(ibin, 
																															fTemplateHistogramContainer.at(itemp).GetBinContent(ibin) * efficiency);
		}

		// add histogram to stack
		stack.Add(&(fTemplateHistogramContainer.at(itemp)));
		if (itemp < 2)
			legend1->AddEntry(&(fTemplateHistogramContainer.at(itemp)), fTemplateNameContainer.at(itemp).data(), "F");
		else if (itemp < 6)
			legend2->AddEntry(&(fTemplateHistogramContainer.at(itemp)), fTemplateNameContainer.at(itemp).data(), "F");
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
		hist_diff->GetXaxis()->SetTitle(fHistData.GetXaxis()->GetTitle()); 
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
			hist_diff->SetBinContent(i+1, fHistData.GetBinContent(i+1)-histsum->GetBinContent(i+1)); 
			hist_diff->SetBinError(i+1, fHistData.GetBinError(i+1)); 
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
	fHistData.Draw("SAMEP");

	if (flag_error0)
		fHistData.Draw("SAMEPE");

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
double TemplateModel::CalculateChi2()
{
	int nbins = fHistData.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	std::vector <double> parameters = GetBestFitParameters();

	double chi2 = 0;

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double nexp = 0;
		double ndata = fHistData.GetBinContent(ibin);

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp) {
			int tempindex = fTemplateParIndexContainer.at(itemp); 
			int effindex = fEffParIndexContainer.at(itemp);

			// get efficiency for the bin
			double efficiency = fEffHistogramContainer.at(itemp).GetBinContent(ibin); 

			// modify efficiency by uncertainty
			double efferr = fEffErrHistogramContainer.at(itemp).GetBinContent(ibin); 

			// check efficiency error 
			if (efferr > 0) 
				efficiency = TMath::Max(0., efficiency + parameters.at(effindex) * efferr); 

			// add expectation from bin
			nexp += parameters.at(tempindex) * efficiency * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);
		}

		// add to chi2
		chi2 += (nexp - ndata) * (nexp - ndata) / nexp;
	}

	// return chi2
	return chi2;
}

// ---------------------------------------------------------
double TemplateModel::CalculateChi2Prob()
{
	double chi2 = CalculateChi2();
	int ndf = GetNDF();

	// return chi2 probability
	return TMath::Prob(chi2, ndf);
}

// ---------------------------------------------------------
double TemplateModel::CalculateMaxLike()
{
	// return maximum likelihood
	return Eval( GetBestFitParameters() );
}

// ---------------------------------------------------------
double TemplateModel::CalculateKSProb()
{
	// create a histogram with the sum of all contributions
	TH1 * histsum = (TH1D*)(fTemplateHistogramContainer.at(0)).Clone("temp");

	int nbins = fHistData.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	std::vector <double> parameters = GetBestFitParameters();

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double bincontent = 0;

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp) {
			int tempindex = fTemplateParIndexContainer.at(itemp); 
			int effindex = fEffParIndexContainer.at(itemp);
			
			// get efficiency for the bin
			double efficiency = fEffHistogramContainer.at(itemp).GetBinContent(ibin); 
			
			// modify efficiency by uncertainty
			double efferr = fEffErrHistogramContainer.at(itemp).GetBinContent(ibin); 
			
			// check efficiency error 
			if (efferr > 0) 
				efficiency = TMath::Max(0., efficiency + parameters.at(effindex) * efferr); 
			
			// add expectation from bin
			bincontent += parameters.at(tempindex) * efficiency * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);
		}

		// set bin content
		histsum->SetBinContent(ibin, bincontent);
	}

	// perform KS test
	double ksprob = histsum->KolmogorovTest(&fHistData);

	// delete histogram
	delete histsum;

	return ksprob;
}

// ---------------------------------------------------------
double TemplateModel::CalculatePValue()
{
	// get best fit parameters
	std::vector<double> par = GetBestFitParameters();

	// check size of parameter vector
	if (par.size() != GetNParameters())
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCTemplateModel::CalculatePValueFast() : Number of parameters is inconsistent.");
		return -1;
	}

	// define temporary variables
	int nbins = fHistData.GetNbinsX();
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

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp) {
			int tempindex = fTemplateParIndexContainer.at(itemp); 
			int effindex = fEffParIndexContainer.at(itemp);
			
			// get efficiency for the bin
			double efficiency = fEffHistogramContainer.at(itemp).GetBinContent(ibin); 
			
			// modify efficiency by uncertainty
			double efferr = fEffErrHistogramContainer.at(itemp).GetBinContent(ibin); 
			
			// check efficiency error 
			if (efferr > 0) 
				efficiency = TMath::Max(0., efficiency + par.at(effindex) * efferr); 
			
			// add expectation from bin
			nexp += par.at(tempindex) * efficiency * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);
		}

		histogram[ibin-1]   = int(nexp);
		expectation[ibin-1] = nexp;

		// calculate p;
		logp += BCMath::LogPoisson(double(int(nexp)), nexp);
		logp_start += BCMath::LogPoisson(fHistData.GetBinContent(ibin), nexp);
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
void TemplateModel::MCMCUserIterationInterface()
{
	int nbins      = fHistData.GetNbinsX();
	int ntemplates = int(fTemplateHistogramContainer.size());

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		double bincontent = 0;

		// loop over all templates
		for (int itemp = 0; itemp < ntemplates; ++itemp) {
			int tempindex = fTemplateParIndexContainer.at(itemp); 
			int effindex = fEffParIndexContainer.at(itemp);

			// get efficiency for the bin
			double efficiency = fEffHistogramContainer.at(itemp).GetBinContent(ibin); 

			// modify efficiency by uncertainty
			double efferr = fEffErrHistogramContainer.at(itemp).GetBinContent(ibin); 

			// check efficiency error 
			if (efferr > 0) 
				efficiency = TMath::Max(0., efficiency + fMCMCx.at(effindex) * efferr); 

			bincontent += fMCMCx.at(tempindex) * efficiency * fTemplateHistogramContainer.at(itemp).GetBinContent(ibin);
		}

		// set bin content
		fUncertaintyHistogramExp->Fill(fHistData.GetBinCenter(ibin), bincontent);

		// loop over bins in the other direction
		int nbinsy = fUncertaintyHistogramObsPosterior->GetNbinsY();
		for (int jbin = 1; jbin <= nbinsy; ++jbin)
		{
			int n = jbin - 1;
			if (fabs(n - bincontent) < 2*sqrt(bincontent))
				fUncertaintyHistogramObsPosterior->Fill(fHistData.GetBinCenter(ibin), n, TMath::Poisson(bincontent, n));
		}
	}

	// fill normalization
	fHistNorm.Fill(fNorm);

	// fill ratios
	int nratios = int( fIndicesRatios1D.size() ); 

	// loop over fractions to fill 
	for (int i = 0; i < nratios; ++i) {
		int nsum = int( (fIndicesRatios1D.at(i)).size() ) - 1; 
		double sum = 0; 
		for (int j = 1; j <= nsum; ++j) {
			int indexsum = fIndicesRatios1D.at(i).at(j); 
			sum += fMCMCx.at(fTemplateParIndexContainer.at(indexsum)); 
		}
		
		fHistRatios1D.at(i).Fill(fMCMCx.at(fTemplateParIndexContainer.at(fIndicesRatios1D.at(i).at(0)))/sum); 
	}

}

// ---------------------------------------------------------
void TemplateModel::PrintRatios(const char * filename, int options, double ovalue)
{
	int nratios = int(fHistRatios1D.size());

	TCanvas* c1 = new TCanvas("c1");

	TPostScript * ps = new TPostScript(filename, 112);
	ps->NewPage(); 

	c1->cd(); 
	BCH1D* h1temp = new BCH1D(&fHistNorm); 
	h1temp->Draw(); 
	c1->Update(); 
	ps->NewPage(); 
	for (int i = 0; i < nratios; ++i) {
		c1->Update(); 
		ps->NewPage();
		c1->cd(); 
		BCH1D* h1temp = new BCH1D(&fHistRatios1D.at(i));
		h1temp->Draw(options, ovalue); 
	}
	c1->Update();
	ps->Close(); 

	delete c1; 
	delete ps; 
}

// ---------------------------------------------------------
int TemplateModel::SetTemplateEfficiency(const char* name, TH1D eff, TH1D efferr)
{
	// get index
	int index = GetIndexTemplate(name); 

	// check index
	if (index < 0) {
		BCLog::OutError("TemplateModel::SetTemplateEfficiency : Could not find template.");
		return 0;
	}

	// check efficiency histogram
	if (CompareHistogramProperties(fTemplateHistogramContainer.at(index), eff) != 1) {
		BCLog::OutError("TemplateModel::SetTemplate efficiency : Template and efficiency histogram properties are incompatible."); 
		return 0;
	}

	// set histogram style
	eff.SetXTitle((fHistData.GetXaxis())->GetTitle());
	eff.SetYTitle("Efficiency");

	efferr.SetXTitle((fHistData.GetXaxis())->GetTitle());
	efferr.SetYTitle("Efficiency uncertainty");

	// set efficiency histogram
	fEffHistogramContainer[index] = eff;
	fEffErrHistogramContainer[index] = efferr;

	// no error
	return 1;
}

// ---------------------------------------------------------
int TemplateModel::SetTemplateEfficiency(const char* name, double effmean, double effsigma)
{
	// get index
	int index = GetIndexTemplate(name); 

	// check index
	if (index < 0) {
		BCLog::OutError("TemplateModel::SetTemplateEfficiency : Could not find template.");
		return 0;
	}

	// create histograms
	TH1D histeff = TH1D("", "", fNBins, fXmin, fXmax); 
 	TH1D histefferr = TH1D("", "", fNBins, fXmin, fXmax); 
	
	// fill histograms
	for (int i = 1; i <= fNBins; ++i) {
		histeff.SetBinContent(i, effmean); 
		histefferr.SetBinContent(i, effsigma);
	}

	// set histograms
	int err = SetTemplateEfficiency(name, histeff, histefferr);

	// return error code
	return err;
}

// ---------------------------------------------------------
int TemplateModel::ConstrainSum(std::vector <int> indices, double mean, double rms)
{
	// add contraint to container(s)
	fConstraintSumIndices.push_back(indices); 
	fConstraintSumMean.push_back(mean); 
	fConstraintSumRMS.push_back(rms); 

	// no error
	return 1;
}

// ---------------------------------------------------------
void TemplateModel::PrintTemp()
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
int TemplateModel::PrintTemplate(const char* name, const char* filename)
{
	// get number of sources of systematic uncertainty
	int nsyst = GetNSystErrors();

	// get index
	int index = GetIndexTemplate(name); 

	// check index
	if (index < 0) {
		BCLog::OutError("TemplateModel::PrintTemplate : Could not find template.");
		return 0;
	}

	// create postscript 
	TPostScript* ps = new TPostScript(filename); 

	// create new canvas
	TCanvas* c1 = new TCanvas("c1", "", 700, 700); 

	c1->Update();
	ps->NewPage();
	c1->cd();

	// create legend
	TLegend l1(0.18, 0.75, 0.85, 0.85); 
	l1.SetBorderSize(0);
	l1.SetFillColor(kWhite);

	// draw histogram and uncertainties
	TH1D hist_template = fTemplateHistogramContainer.at(index); 
	hist_template.SetFillColor(kWhite);
	hist_template.SetFillStyle(0); 
	hist_template.SetMarkerSize(0); 
	hist_template.SetLineWidth(0); 
	l1.AddEntry(&hist_template, name, "L");
	TH1D hist_totalerr = CombineUncertainties(name); 
	TH1D hist_template_totalerr = CreateErrorHist(hist_template, hist_totalerr); 	
	hist_template_totalerr.SetFillColor(kYellow); 
	hist_template_totalerr.SetFillStyle(1001);
	hist_template_totalerr.SetMarkerSize(0);
	l1.AddEntry(&hist_template_totalerr, "Systematic uncertainties", "F");
	TH1D hist_efferr = fEffErrHistogramContainer.at(index); 
	TH1D hist_template_efferr = CreateErrorHist(hist_template, hist_efferr); 
	hist_template_totalerr.SetFillColor(kRed); 
	hist_template_totalerr.SetFillStyle(1001);
	hist_template_efferr.SetMarkerSize(0);
	l1.AddEntry(&hist_template_efferr, "Efficiency uncertainties", "F");
	int binmax = hist_template.GetMaximumBin();
	double ymax = hist_template.GetBinContent(binmax) + 2.0 * hist_template_totalerr.GetBinError(binmax);
	hist_template_totalerr.GetYaxis()->SetRangeUser(0.0, 1.25 * ymax);
	hist_template_totalerr.Draw("E2"); 
	hist_template_efferr.Draw("SAMEE2");
	hist_template.Draw("SAME");

	// draw legend
	l1.Draw();

	// update ps
	c1->Update();
	ps->NewPage();
	c1->cd();

	// create legend
	TLegend l2(0.18, 0.75, 0.85, 0.85); 
	l2.SetBorderSize(0);
	l2.SetFillColor(kWhite);

	// print uncertainties
	c1->cd(2); 
	hist_efferr = fEffErrHistogramContainer.at(index);
	double ymin = hist_efferr.GetMinimum(); 
	ymax = hist_efferr.GetMaximum(); 
	l2.AddEntry(&hist_efferr, "Efficiency", "L"); 
	hist_efferr.SetStats(kFALSE);
	hist_efferr.Draw(); 

	// loop over all uncertainties
	for (int i = 0; i < nsyst; ++i) {
		TH1D* hist = new TH1D(fSystErrorHistogramContainer.at(i).at(index));
		hist->SetLineColor(2 + i); 
		if (hist->GetMaximum()>ymax)
			ymax = hist->GetMaximum();
		if (hist->GetMinimum()<ymin)
			ymin = hist->GetMinimum(); 
		l2.AddEntry(hist, fSystErrorNameContainer.at(i).c_str(), "L");
		hist->Draw("SAME"); 
	}
	if (ymin < 0)
		ymin = 1.25*ymin;
	else 
		ymin = 0.8*ymin;

	if (ymax > 0)
		ymax = 1.25*ymax; 
	else
		ymax = 0.8*ymax; 

	hist_efferr.GetYaxis()->SetRangeUser(ymin, ymax);

	// draw legend
	l2.Draw();
	
	// close ps
	c1->Update(); 
	ps->Close();

	// print canvas
	//	c1->Print(filename); 

	// free memory
	delete c1; 
	delete ps;
	

	// no error 
	return 1;
}

// ---------------------------------------------------------
TH1D TemplateModel::CreateErrorHist(TH1D hist, TH1D histerr)
{
	// check histogram properties
	if (CompareHistogramProperties(fHistData, hist) != 1) {
		BCLog::OutError("TemplateModel::CreateErrorHist : Histograms are incompatible."); 
		return hist;
	}

	// copy histogram
	TH1D h = hist;

	// set style
	h.SetStats(kFALSE);
	h.SetFillColor(kYellow); 
	h.SetFillStyle(1001);

	// get number of bins
	int nbins = hist.GetNbinsX(); 

	// loop over bins
	for (int i = 1; i <= nbins; ++i) {
		h.SetBinError(i, histerr.GetBinContent(i) * hist.GetBinContent(i));
	}

	// return histogram
	return h;
}

// ---------------------------------------------------------
TH1D TemplateModel::CombineUncertainties(const char* name)
{
	// get number of sources of systematic uncertainty
	int nsyst = GetNSystErrors();

	// get template index
	int tempindex = GetIndexTemplate(name);

	// create new histogram
 	TH1D hist = TH1D("", "", fNBins, fXmin, fXmax); 
	
	// fill histogram
	for (int ibin = 1; ibin <= fNBins; ++ibin) {

		// define total uncertainty
		double err = 0; 

		// add efficiency uncertainty squared
		double erreff = fEffErrHistogramContainer.at(tempindex).GetBinContent(ibin); 
		err += erreff * erreff; 

		// add systematic uncertainty squared
		for (int isyst = 0; isyst < nsyst; ++isyst) {
			double errsyst = fSystErrorHistogramContainer.at(isyst).at(tempindex).GetBinContent(ibin);
			err += errsyst*errsyst; 
		}

		// take square root
		err = sqrt(err);

		// set bin content
		hist.SetBinContent(ibin, err); 
	}

	// return histogram
	return hist;
}

// ---------------------------------------------------------
	int  TemplateModel::GetIndexTemplate(const char* name)
	{
   int index = -1;
	 int n = GetNTemplates();

   for (int i = 0; i < n; i++)
		 if (name == fTemplateNameContainer.at(i))
			 index = i;

   if (index < 0) {
		 BCLog::OutWarning("TemplateModel::GetIndexTemplate : Template does not exist."); 
		 return 0;
   }

	 // return index
	 return index;		
	}

// ---------------------------------------------------------
int TemplateModel::GetIndexSystError(const char* name)
	{
   int index = -1;
	 int n = GetNSystErrors();

   for (int i = 0; i < n; i++)
		 if (name == fSystErrorNameContainer.at(i))
			 index = i;

   if (index < 0) {
		 BCLog::OutWarning("TemplateModel::GetIndexSystError : Template does not exist."); 
		 return 0;
   }

	 // return index
	 return index;		
	}

// ---------------------------------------------------------
	int  TemplateModel::GetParIndexTemplate(const char* name)
	{ 
   int index = -1;
	 int n = GetNTemplates();

   for (int i = 0; i < n; i++)
		 if (name == fTemplateNameContainer.at(i))
			 index = fTemplateParIndexContainer.at(i);

   if (index < 0) {
		 BCLog::OutWarning("TemplateModel::GetParIndexTemplate : Template does not exist."); 
		 return 0;
   }

	 // return index
	 return index;
	}
	
// ---------------------------------------------------------
int  TemplateModel::GetParIndexTemplate(int index)
{
	// get number of templates
	int n = GetNTemplates(); 
	
	if (index < 0 || index > n) {
		BCLog::OutError("TemplateModel::GetParIndexTemplate : Index out of range.");
		return -1;
	}

	// return index
	return fTemplateParIndexContainer.at(index);
}

// ---------------------------------------------------------
	int  TemplateModel::GetParIndexEff(const char* name)
	{
   int index = -1;
	 int n = GetNTemplates();

   for (int i = 0; i < n; i++)
		 if (name == fTemplateNameContainer.at(i))
			 index = fTemplateParIndexContainer.at(i);

   if (index < 0) {
		 BCLog::OutWarning("TemplateModel::GetParIndexEff : Template does not exist."); 
		 return 0;
   }

	 // return index
	 return index;

	}

// ---------------------------------------------------------
	int  TemplateModel::GetParIndexSystError(const char* name)
	{
   int index = -1;
	 int n = GetNTemplates();

   for (int i = 0; i < n; i++)
		 if (name == fSystErrorNameContainer.at(i))
			 index = fSystErrorParIndexContainer.at(i);

   if (index < 0) {
		 BCLog::OutWarning("TemplateModel::GetParIndexStatError : Systematic error does not exist."); 
		 return 0;
   }

	 // return index
	 return index;

	}

// ---------------------------------------------------------
int  TemplateModel::PerformFit()
{
	// initialize
	if (!Initialize()) {
		BCLog::OutError("TemplateModel::PerformFit : Could not initialize template fitter."); 
		return 0;
	}

	// run Markov Chains
	MarginalizeAll();
	
	// find global mode
	FindMode();

	// no error 
	return 1;
}
// ---------------------------------------------------------
