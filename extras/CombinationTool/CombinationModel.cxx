// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project CombinationTool
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "CombinationModel.h"

#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>

#include <TF1.h>
#include <TCanvas.h> 
#include <TGraphAsymmErrors.h>
#include <TH2D.h> 
#include <TPostScript.h> 
#include <TLatex.h> 
#include <TLine.h> 
#include <TBox.h> 
#include <TMath.h> 

double SplitGaussian(double* x, double* par);

// ---------------------------------------------------------
CombinationModel::CombinationModel(const char * name, double xmin, double xmax) : BCModel()
																																								, fFlagSystErrors(true)
{
	AddParameter(name, xmin, xmax);
};

// ---------------------------------------------------------
CombinationModel::~CombinationModel()
{
	for (int i = 0; i < int(fFunctionContainer.size()); ++i) {
		TF1* f = fFunctionContainer.at(i);
		delete f;
		f = 0;
	}
	fFunctionContainer.clear(); 

};  // default destructor

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob = 0.;

	return logprob;
}

// ---------------------------------------------------------
double CombinationModel::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	// For flat prior it's very easy.
//	for(unsigned int i=0; i < this -> GetNParameters(); i++)
//		logprob -= log(this -> GetParameter(i) -> GetRangeWidth());

	return logprob;
}

// ---------------------------------------------------------
int CombinationModel::GetContIndexSystError(const char* systerrorname)
{
	// prepare index counter 
	int index = -1;
	int n = GetNSystErrors();

	for (int i = 0; i < n; i++)
		if (systerrorname == fSystErrorNameContainer.at(i))
			index = i;

   if (index < 0) {
		 BCLog::OutWarning(Form("CombinationModel::GetContIndexSystError : Systematic error \"%s\" does not exist.", systerrorname));
		return -1;
   }

	// return index
	return index;
}

// ---------------------------------------------------------
int CombinationModel::GetContIndexChannel(const char* channelname)
{
	// prepare index counter 
	int index = -1;
	int n = GetNChannels();

	for (int i = 0; i < n; i++)
		if (channelname == fChannelNameContainer.at(i))
			index = i;

   if (index < 0) {
		 BCLog::OutWarning(Form("CombinationModel::GetContIndexChannel : Channel \"%s\" does not exist.", channelname));
		return -1;
   }

	// return index
	return index;
}

// ---------------------------------------------------------
int CombinationModel::GetContIndexChannelBackground(const char* channelname, const char* backgroundname)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname);

	// check channel index
	if (channelindex < 0) {
		return -1;
	}

	// prepare index counter 
	int index = -1;
	int n = GetNChannelBackgrounds(channelindex);

	for (int i = 0; i < n; i++)
		if (backgroundname == fChannelBackgroundNameContainer.at(channelindex).at(i))
			index = i;

   if (index < 0) {
		 BCLog::OutWarning(Form("CombinationModel::GetContIndexChannelBackground : Background \"%s\" does not exist for channel \"%s\".", backgroundname, channelname));
		return -1;
   }

	// return index
	return index;
}

// ---------------------------------------------------------
int CombinationModel::GetParIndexChannelBackground(const char* channelname, const char* backgroundname)
{
	// get channel container index
	int channelindex = GetContIndexChannel(channelname); 

	// check index
	if (channelindex < 0){ 
		return 1;
	}

	// get channel container index
	int backgroundindex = GetContIndexChannelBackground(channelname, backgroundname); 

	// check index
	if (backgroundindex < 0){ 
		return 1;
	}

	// return index
	return fParIndexChannelBackground.at(channelindex).at(backgroundindex); 
}

// ---------------------------------------------------------
int CombinationModel::GetParIndexChannelBackground(int channelindex, int backgroundindex)
{
	// no checks for CPU-time reasons
	return fParIndexChannelBackground.at(channelindex).at(backgroundindex); 
}

// ---------------------------------------------------------
int CombinationModel::AddChannel(const char* channelname)
{
	// add channel name to container
	fChannelNameContainer.push_back(channelname); 

	// add background container
	fChannelBackgroundNameContainer.push_back(std::vector<std::string>(0));
	fChannelBackgroundPriorContainer.push_back(std::vector<TF1*>(0)); 
	fParIndexChannelBackground.push_back(std::vector<int>(0));

	// add signal prior to container
	fChannelSignalPriorContainer.push_back(0); 

	// add observation container
	fChannelObservation.push_back(0);

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::AddChannelBackground(const char* channelname, const char* backgroundname, double xmin, double xmax)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::AddChannelBackground : Can not find channel index."); 
		return 0; 
	}

	// add background name
	fChannelBackgroundNameContainer.at(channelindex).push_back(backgroundname);

	// add prior
	fChannelBackgroundPriorContainer.at(channelindex).push_back(0); 

	// add parameter
	AddParameter(Form("%s_%s", channelname, backgroundname), xmin, xmax);
	
	// add parameter index
	fParIndexChannelBackground.at(channelindex).push_back(GetNParameters()-1);

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::AddSystError(const char* systerrorname)
{
	// add syst error name to container
	fSystErrorNameContainer.push_back(systerrorname); 

	// add parameter
	AddParameter(systerrorname, -5.0, 5.0);

	// add parameter index
	fParIndexSystErrorContainer.push_back(GetNParameters()-1);

	// get number of channels
	int nchannels = GetNChannels(); 

	std::vector< std::vector<double> > sdown;
	std::vector< std::vector<double> > sup;

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		// get number of background sources in this channel
		int nbackground = GetNChannelBackgrounds(i);
		
		std::vector<double> schanneldown(nbackground);
		std::vector<double> schannelup(nbackground);

		// loop over all background sources
		for (int j = 0; j < nbackground; ++j) {
			schanneldown[j] = 0.0;
			schannelup[j] = 0.0;
		}

		sdown.push_back(schanneldown);
		sup.push_back(schannelup);
	}	

	// add vectors to container
	fSystErrorSigmaDownContainer.push_back(sdown);
	fSystErrorSigmaUpContainer.push_back(sup);

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::SetSystErrorChannelSignal(const char* channelname, double sigmadown, double sigmaup)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::SetSystErrorChannelBackground : Can not find channel index."); 
		return 0; 
	}

	// set systematic error sigmas
	// ...

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::SetSystErrorChannelBackground(const char* systerrorname, const char* channelname, const char* backgroundname, double sigmadown, double sigmaup)
{
	// get syst error index
	int systerrorindex = GetContIndexSystError(systerrorname); 

	// check channel index
	if (systerrorindex < 0) {
		BCLog::OutError("CombinationModel::SetSystErrorChannelBackground : Can not find systematic error index."); 
		return 0; 
	}

	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::SetSystErrorChannelBackground : Can not find channel index."); 
		return 0; 
	}

	// get background index
	int backgroundindex = GetContIndexChannelBackground(channelname, backgroundname); 

	// check background index
	if (backgroundindex < 0) {
		BCLog::OutError("CombinationModel::SetSystErrorChannelBackground : Can not find background index."); 
		return 0; 
	}

	// set systematic error sigmas
	fSystErrorSigmaUpContainer.at(systerrorindex).at(channelindex)[backgroundindex] = sigmaup;
	fSystErrorSigmaDownContainer.at(systerrorindex).at(channelindex)[backgroundindex] = sigmadown;
	
	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::SetChannelObservation(const char* channelname, double observation)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::SetChannelObservation : Can not find channel index."); 
		return 0; 
	}

	// set observation
	fChannelObservation[channelindex] = observation;

	// no error 
	return 1; 
}

// ---------------------------------------------------------
int CombinationModel::SetChannelSignalPrior(const char* channelname, TF1* prior)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::SetChannelSignalPrior : Can not find channel index."); 
		return 0; 
	}

	// check if function exists
	if (!prior) {
		BCLog::OutError("CombinationModel::SetChannelSignalPrior : Function does not exist."); 
		return 0; 
	}

	// set prior
	fChannelSignalPriorContainer[channelindex] = prior;

	// no error
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::SetChannelBackgroundPrior(const char* channelname, const char* backgroundname, TF1* prior)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::SetChannelBackgroundPrior : Can not find channel index."); 
		return 0; 
	}

	// get background index
	int backgroundindex = GetContIndexChannelBackground(channelname, backgroundname); 

	// check background index
	if (backgroundindex < 0) {
		BCLog::OutError("CombinationModel::SetChannelBackgroundPrior : Can not find background index."); 
		return 0; 
	}

	// check if function exists
	if (!prior) {
		BCLog::OutError("CombinationModel::SetChannelBackgroundPrior : Function does not exist."); 
		return 0; 
	}

	// set prior
	(fChannelBackgroundPriorContainer.at(channelindex))[backgroundindex] = prior;

	// no error 
	return 1; 
}

// ---------------------------------------------------------
int CombinationModel::SetChannelSignalPriorGauss(const char* channelname, double mean, double sigma)
{
	// create new function
	TF1* f = new TF1(Form("f_%s", channelname), "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(- (x-[0])*(x-[0])/2/[1]/[1] )");
	f->SetParameter(0, mean); 
	f->SetParameter(1, sigma); 

	// add function to container
	fFunctionContainer.push_back(f); 

	// set channel background prior
	return SetChannelSignalPrior(channelname, f); 
}

// ---------------------------------------------------------
int CombinationModel::SetChannelSignalPriorGauss(const char* channelname, double mean, double sigmadown, double sigmaup)
{
	// create new function
	TF1* f = new TF1(Form("f_%s", channelname), SplitGaussian, mean - 5.0*sigmadown, mean + 5.0*sigmaup, 3);
	f->SetParameter(0, mean); 
	f->SetParameter(1, sigmadown); 
	f->SetParameter(2, sigmaup); 

	// add function to container
	fFunctionContainer.push_back(f); 

	// set channel background prior
	return SetChannelSignalPrior(channelname, f); 
}

// ---------------------------------------------------------
int CombinationModel::SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigma)
{
	// create new function
	TF1* f = new TF1(Form("f_%s_%s", channelname, backgroundname), "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(- (x-[0])*(x-[0])/2/[1]/[1] )");
	f->SetParameter(0, mean); 
	f->SetParameter(1, sigma); 

	// add function to container
	fFunctionContainer.push_back(f); 

	// set channel background prior
	return SetChannelBackgroundPrior(channelname, backgroundname, f); 
}

// ---------------------------------------------------------
int CombinationModel::SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigmadown, double sigmaup)
{
	// get parameter index
	int parindex = GetParIndexChannelBackground(channelname, backgroundname);

	// create new function
	TF1* f = new TF1(Form("f_%s_%s", channelname, backgroundname), SplitGaussian, GetParameter(parindex)->GetLowerLimit(), GetParameter(parindex)->GetUpperLimit(), 3);
	f->SetParameter(0, mean); 
	f->SetParameter(1, sigmadown); 
	f->SetParameter(2, sigmaup); 

	// add function to container
	fFunctionContainer.push_back(f); 

	// set channel background prior
	return SetChannelBackgroundPrior(channelname, backgroundname, f); 
}

// ---------------------------------------------------------
int CombinationModel::PerformAnalysis()
{
	// perform mcmc
	MarginalizeAll(); 

	// find mode using Minuit
	FindMode( GetBestFitParameters() );

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::GetParIndexSystError(int systerrorindex)
{
	return fParIndexSystErrorContainer.at(systerrorindex);
}

// ---------------------------------------------------------
void CombinationModel::PrintChannelOverview(const char* filename)
{
	// get number of channels
	int nchannels = GetNChannels(); 

	// define vector of histograms 
	std::vector<BCH1D*> histos_with_syst_error; 
	std::vector<BCH1D*> histos_without_syst_error; 

	// perform studies without systematics
	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		TH1D* h = new TH1D( PerformSingleChannelAnalysis( fChannelNameContainer.at(i).c_str(), false ));
		BCH1D* hist = new BCH1D( h ); 		
		histos_without_syst_error.push_back( hist );
	}

	// perform studies with systematics
	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		TH1D* h = new TH1D( PerformSingleChannelAnalysis( fChannelNameContainer.at(i).c_str(), true ));
		BCH1D* hist = new BCH1D( h ); 		
		histos_with_syst_error.push_back( hist );
	}

	// create graph
	TGraphAsymmErrors* graph_without_syst_error = new TGraphAsymmErrors(GetNChannels()+1); 
	graph_without_syst_error->SetMarkerStyle(20); 
	graph_without_syst_error->SetMarkerSize(0); 

	TGraphAsymmErrors* graph_with_syst_error = new TGraphAsymmErrors(GetNChannels()+1); 
	graph_with_syst_error->SetMarkerStyle(20); 
	graph_with_syst_error->SetMarkerSize(1); 

	// coordinate system
	double xmin = 0.0; 
	double xmax = 0.0; 
	double xwidth = 0.0;
	double ymin = -0.5;
	double ymax = double(nchannels)+0.5;

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {

		// get histogram
		BCH1D* h = histos_without_syst_error.at(i);

		// get summary information
		double median = h->GetMedian(); 
		double q16 = h->GetQuantile(0.16);
		double q84 = h->GetQuantile(0.84);
		double errlow = median - q16;
		double errhigh = q84 - median; 

		// update coordinate system
		if (q16 < xmin || i == 0)
			xmin = q16;
		if (q84 > xmax || i == 0)
			xmax = q84; 
		xwidth = xmax - xmin; 

		// set point and error 
		graph_without_syst_error->SetPoint(i, median, double(nchannels-i));
		graph_without_syst_error->SetPointError(i, errlow, errhigh, 0, 0);
	}

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {

		// get histogram
		BCH1D* h = histos_with_syst_error.at(i);

		// get summary information
		double median = h->GetMedian(); 
		double q16 = h->GetQuantile(0.16);
		double q84 = h->GetQuantile(0.84);
		double errlow = median - q16;
		double errhigh = q84 - median; 

		// update coordinate system
		if (q16 < xmin || i == 0)
			xmin = q16;
		if (q84 > xmax || i == 0)
			xmax = q84; 
		xwidth = xmax - xmin; 

		// set point and error 
		graph_with_syst_error->SetPoint(i, median, double(nchannels-i));
		graph_with_syst_error->SetPointError(i, errlow, errhigh, 0, 0);
	}

	// set point from combination
	BCH1D* h = GetMarginalized( GetParameter(0)->GetName().c_str() );
	double median_with_syst_error = h->GetMedian(); 
	double q16_with_syst_error = h->GetQuantile(0.16);
	double q84_with_syst_error = h->GetQuantile(0.84);
	double errlow_with_syst_error = median_with_syst_error - q16_with_syst_error;
	double errhigh_with_syst_error = q84_with_syst_error - median_with_syst_error; 
	
	// set point and error 
	graph_with_syst_error->SetPoint(nchannels, median_with_syst_error, 0.);
	graph_with_syst_error->SetPointError(nchannels, errlow_with_syst_error, errhigh_with_syst_error, 0, 0);

	// repeat studies without systematic errors
	SetFlagSystErrors(false);
	PerformAnalysis(); 
	h = GetMarginalized( GetParameter(0)->GetName().c_str() );
	SetFlagSystErrors(true);
	
	// set point from combination
	double median_without_syst_error = h->GetMedian(); 
	double q16_without_syst_error = h->GetQuantile(0.16);
	double q84_without_syst_error = h->GetQuantile(0.84);
	double errlow_without_syst_error = median_without_syst_error - q16_without_syst_error;
	double errhigh_without_syst_error = q84_without_syst_error - median_without_syst_error;
	
	// set point and error 
	graph_without_syst_error->SetPoint(nchannels, median_without_syst_error, 0.);
	graph_without_syst_error->SetPointError(nchannels, errlow_without_syst_error, errhigh_without_syst_error, 0, 0);

	// create histogram for axes
	TH2D* hist_axes = new TH2D("", Form(";%s;Channel", GetParameter(0)->GetName().c_str()), 1, xmin - 0.25*xwidth, xmax + 1.75 * xwidth, nchannels+1, ymin, ymax); 
	hist_axes->SetStats(kFALSE);
	hist_axes->GetYaxis()->SetNdivisions(0);
	hist_axes->GetYaxis()->SetTitleOffset(1.0);

	// create canvas
	TCanvas* c1 = new TCanvas();
	c1->cd(); 

	// create latex
	TLatex* latex = new TLatex(); 
	latex->SetTextSize(0.04);
	if (nchannels>=10)
		latex->SetTextSize(0.02);
	latex->SetTextAlign(12);

	// create lines
	TLine* line_mode = new TLine(); 
	line_mode->SetLineWidth(1);
	line_mode->SetLineStyle(1);
	line_mode->SetLineColor(kRed);
	TBox* box_total = new TBox();
	box_total->SetLineWidth(0); 
	box_total->SetFillColor(kYellow);
	box_total->SetFillStyle(1001);
	TBox* box_stat = new TBox();
	box_stat->SetLineWidth(0); 
	box_stat->SetFillColor(kYellow-7);
	box_stat->SetFillStyle(1001);
	TBox* box_help = new TBox();
	box_help->SetLineWidth(1); 
	box_help->SetFillColor(kBlack);
	box_help->SetFillStyle(0);

	// draw
	hist_axes->Draw();
	box_total->DrawBox(q16_with_syst_error, ymin, q84_with_syst_error, ymax);
	box_stat->DrawBox(q16_without_syst_error, ymin, q84_without_syst_error, ymax);
	line_mode->DrawLine(median_with_syst_error, ymin, median_with_syst_error, ymax);
	graph_with_syst_error->Draw("SAMEPZ"); 
	graph_without_syst_error->Draw("SAMEP"); 

	// loop over all channels and draw labels
	for (int i = 0; i < nchannels; ++i) {
		latex->DrawLatex(xmax + 0.25*xwidth, double(nchannels-i), fChannelNameContainer.at(i).c_str());
	}
	latex->DrawLatex(xmax + 0.25*xwidth, 0., "combination");
	hist_axes->Draw("SAMEAXIS");
	box_help->DrawBox(xmin - 0.25*xwidth, ymin, xmax + 1.75 * xwidth, ymax);

	// print to file 
	c1->Print(filename);

	// create postscript
	TPostScript* ps = new TPostScript("summary_signal.ps");

	// create canvas and prepare postscript
	TCanvas * c2 = new TCanvas("", "", 1000, 500);
	c2->Divide(2, 1);

	// loop over all channels
	for (int i = 0; i < nchannels; ++i) {
		c2->Update();
		ps->NewPage();
		c2->cd(1);
		histos_without_syst_error.at(i)->Draw(); 
		c2->cd(2);
		histos_with_syst_error.at(i)->Draw(); 
	}

	// close ps
	c2->Update();
	ps->Close();

	// free memory 
	delete c1; 
	delete c2;
	delete line_mode; 
	delete box_total; 
	delete box_stat;
	delete box_help; 
	delete latex;
	delete hist_axes;
	delete graph_without_syst_error;
	delete graph_with_syst_error;
	delete ps;
	for (int i = 0; i < nchannels; ++i) {
		delete histos_with_syst_error[i];
		delete histos_without_syst_error[i];
	} 
	histos_with_syst_error.clear();
	histos_without_syst_error.clear();
}

// ---------------------------------------------------------
double SplitGaussian(double* x, double* par)
{
	double mean = par[0]; 
	double sigmadown = par[1]; 
	double sigmaup = par[2];

	double sigma = sigmadown;

	if (x[0] > mean)
		sigma = sigmaup; 

	return 1.0/sqrt(2.0*TMath::Pi())/sigma * exp(- (x[0]-mean)*(x[0]-mean)/2./sigma/sigma);
}

// ---------------------------------------------------------

