#include "CombinationModel.h"
#include "ParameterSummary.h"

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

#include <fstream>
#include <iostream>
#include <iomanip>

double SplitGaussian(double* x, double* par);

// ---------------------------------------------------------
CombinationModel::CombinationModel(const char * name, double xmin, double xmax) : BCModel()
																																								, fFlagSystErrors(true)
																																								, fSummaryCombinationNoSyst(0)
																																								, fSummaryCombinationSyst(0)
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

	// delete summary container
	for (int i = 0; i < int(fSummaryChannelNoSyst.size()); ++i) 
		delete fSummaryChannelNoSyst[i]; 
	fSummaryChannelNoSyst.clear(); 

	for (int i = 0; i < int(fSummaryChannelSyst.size()); ++i) 
		delete fSummaryChannelSyst[i]; 
	fSummaryChannelSyst.clear(); 

	delete fSummaryCombinationNoSyst;
	delete fSummaryCombinationSyst;
}; 

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

	// add systematic uncertainties for signal
	std::vector< double > sigma_signal_down;
	sigma_signal_down.assign(nchannels, 0.0);
	std::vector< double > sigma_signal_up(nchannels);
	sigma_signal_up.assign(nchannels, 0.0);

	// add to container
	fSystErrorChannelSigmaDownContainer.push_back(sigma_signal_down);
	fSystErrorChannelSigmaUpContainer.push_back(sigma_signal_up);

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::SetSystErrorChannelSignal(const char* systerrorname, const char* channelname, double sigmadown, double sigmaup)
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

	// set systematic error sigmas
	fSystErrorChannelSigmaDownContainer.at(systerrorindex)[channelindex] = sigmadown;
	fSystErrorChannelSigmaUpContainer.at(systerrorindex)[channelindex] = sigmaup;

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
		//		BCLog::OutError("CombinationModel::SetChannelSignalPrior : Function does not exist."); 
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
		//		BCLog::OutError("CombinationModel::SetChannelBackgroundPrior : Function does not exist."); 
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
int CombinationModel::PerformFullAnalysis()
{
	// ---- perform combination analysis without systematics ---- //

	// print to screen
	BCLog::OutSummary("Perform combination without systematics"); 

	// set flags
	SetFlagSystErrors(false);

	// perform analysis
	PerformAnalysis();

	// delete old summary information
	delete fSummaryCombinationNoSyst; 

	// copy summary information
	fSummaryCombinationNoSyst = new ParameterSummary("combination");
	fSummaryCombinationNoSyst->Summarize( GetMarginalized( GetParameter(0)->GetName().c_str() ));

	// ---- perform combination analysis with systematics ---- //

	// print to screen
	BCLog::OutSummary("Perform combination with systematics"); 

	// set flags
	SetFlagSystErrors(true);

	// perform analysis
	PerformAnalysis();

	// delete old summary information
	delete fSummaryCombinationSyst; 

	// copy summary information
	fSummaryCombinationSyst = new ParameterSummary("combination");
	fSummaryCombinationSyst->Summarize( GetMarginalized( GetParameter(0)->GetName().c_str() ));

	// ---- perform single channel analysis without systematics ---- //

	// get number of channels
	int nchannels = GetNChannels();

	// print to screen
	BCLog::OutSummary("Perform single channel analysis without systematics"); 

	// delete old summaries
	for (int i = 0; i < int(fSummaryChannelNoSyst.size()); ++i) 
		delete fSummaryChannelNoSyst[i]; 
	fSummaryChannelNoSyst.clear(); 

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		ParameterSummary* ps = new ParameterSummary(fChannelNameContainer.at(i).c_str());
		*ps = PerformSingleChannelAnalysis( fChannelNameContainer.at(i).c_str(), false );
		fSummaryChannelNoSyst.push_back(ps);
	}

	// ---- perform single channel analysis with systematics ---- //

	// print to screen
	BCLog::OutSummary("Perform single channel analysis with systematics"); 

	// delete old summaries
	for (int i = 0; i < int(fSummaryChannelSyst.size()); ++i) 
		delete fSummaryChannelSyst[i]; 
	fSummaryChannelSyst.clear(); 

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		ParameterSummary* ps = new ParameterSummary(fChannelNameContainer.at(i).c_str());
		*ps = PerformSingleChannelAnalysis( fChannelNameContainer.at(i).c_str(), true );
		fSummaryChannelSyst.push_back(ps);
	}

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::GetParIndexSystError(int systerrorindex)
{
	return fParIndexSystErrorContainer.at(systerrorindex);
}

// ---------------------------------------------------------
int CombinationModel::PrintChannelOverview(const char* filename)
{
	// check if summary information is available
	if (!fSummaryCombinationNoSyst || !fSummaryCombinationSyst) {
		BCLog::OutError("CombinationModel::PrintChannelOverview : Summary information not available."); 
		return -1; 
	}

	// get number of channels
	int nchannels = GetNChannels(); 

	// create graph
	TGraphAsymmErrors* graph_nosyst = new TGraphAsymmErrors(GetNChannels()+1); 
	graph_nosyst->SetMarkerStyle(20); 
	graph_nosyst->SetMarkerSize(0); 

	TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(GetNChannels()+1); 
	graph_syst->SetMarkerStyle(20); 
	graph_syst->SetMarkerSize(1); 

	// coordinate system
	double xmin = 0.0; 
	double xmax = 0.0; 
	double xwidth = 0.0;
	double ymin = -0.5;
	double ymax = double(nchannels)+0.5;

	// ---- single channel analysis without systematics ---- // 

	// loop over all channels
	for (int i = 0; i < nchannels; ++i) {

		// get summary
		ParameterSummary* ps = fSummaryChannelNoSyst.at(i); 

		if (!ps) {
			BCLog::OutError("CombinationModel::PrintChannelOverview : Single channel summary information not available."); 
			return -1; 
		}

		// get summary information
		double median = ps->GetMedian(); 
		double q16 = ps->GetQuantile16();
		double q84 = ps->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 

		// update coordinate system
		if (q16 < xmin || i == 0)
			xmin = q16;
		if (q84 > xmax || i == 0)
			xmax = q84; 
		xwidth = xmax - xmin; 

		// set point and error 
		graph_nosyst->SetPoint(i, median, double(nchannels-i));
		graph_nosyst->SetPointError(i, errlow, errhigh, 0, 0);
	}

	// ---- single channel analysis with systematics ---- // 

	// loop over all channels
	for (int i = 0; i < nchannels; ++i) {

		// get summary
		ParameterSummary* ps = fSummaryChannelSyst.at(i); 

		if (!ps) {
			BCLog::OutError("CombinationModel::PrintChannelOverview : Single channel summary information not available."); 
			return -1; 
		}

		// get summary information
		double median = ps->GetMedian(); 
		double q16 = ps->GetQuantile16();
		double q84 = ps->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 

		// update coordinate system
		if (q16 < xmin || i == 0)
			xmin = q16;
		if (q84 > xmax || i == 0)
			xmax = q84; 
		xwidth = xmax - xmin; 

		// set point and error 
		graph_syst->SetPoint(i, median, double(nchannels-i));
		graph_syst->SetPointError(i, errlow, errhigh, 0, 0);
	}

	// ---- combined analysis without systematics

	// set point from combination
	double median_nosyst = fSummaryCombinationNoSyst->GetMedian(); 
	double q16_nosyst = fSummaryCombinationNoSyst->GetQuantile16();
	double q84_nosyst = fSummaryCombinationNoSyst->GetQuantile84();
	double errlow_nosyst = median_nosyst - q16_nosyst;
	double errhigh_nosyst = q84_nosyst - median_nosyst; 
	
	// set point and error 
	graph_nosyst->SetPoint(nchannels, median_nosyst, 0.);
	graph_nosyst->SetPointError(nchannels, errlow_nosyst, errhigh_nosyst, 0, 0);

	// ---- combined analysis with systematics

	// set point from combination
	double median_syst = fSummaryCombinationSyst->GetMedian(); 
	double q16_syst = fSummaryCombinationSyst->GetQuantile16();
	double q84_syst = fSummaryCombinationSyst->GetQuantile84();
	double errlow_syst = median_syst - q16_syst;
	double errhigh_syst = q84_syst - median_syst; 
	
	// set point and error 
	graph_syst->SetPoint(nchannels, median_syst, 0.);
	graph_syst->SetPointError(nchannels, errlow_syst, errhigh_syst, 0, 0);

	// ---- do the plotting ---- // 

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
	TLine* line_median = new TLine(); 
	line_median->SetLineWidth(1);
	line_median->SetLineStyle(1);
	line_median->SetLineColor(kRed);
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
	box_total->DrawBox(q16_syst, ymin, q84_syst, ymax);
	box_stat->DrawBox(q16_nosyst, ymin, q84_nosyst, ymax);
	line_median->DrawLine(median_syst, ymin, median_syst, ymax);
	graph_syst->Draw("SAMEPZ"); 
	graph_nosyst->Draw("SAMEP"); 

	// loop over all channels and draw labels
	for (int i = 0; i < nchannels; ++i) {
		latex->DrawLatex(xmax + 0.25*xwidth, double(nchannels-i), fChannelNameContainer.at(i).c_str());
	}
	latex->DrawLatex(xmax + 0.25*xwidth, 0., "combination");
	hist_axes->Draw("SAMEAXIS");
	box_help->DrawBox(xmin - 0.25*xwidth, ymin, xmax + 1.75 * xwidth, ymax);

	// print to file 
	c1->Print(filename);

	/*
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
		histos_nosyst.at(i)->Draw(); 
		c2->cd(2);
		histos_syst.at(i)->Draw(); 
	}

	// close ps
	c2->Update();
	ps->Close();
	*/

	// free memory 
	delete c1; 
	delete line_median; 
	delete box_total; 
	delete box_stat;
	delete box_help; 
	delete latex;
	delete hist_axes;
	delete graph_nosyst;
	delete graph_syst;
	/*
	delete c2;
	delete ps;
	for (int i = 0; i < nchannels; ++i) {
		delete histos_syst[i];
		delete histos_nosyst[i];
	} 
	histos_syst.clear();
	histos_nosyst.clear();
	*/

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::PrintChannelSummary(const char* filename)
{
	// check if summary information is available
	if (!fSummaryCombinationNoSyst || !fSummaryCombinationSyst) {
		BCLog::OutError("CombinationModel::PrintChannelOverview : Summary information not available."); 
		return -1; 
	}

	// create stream
	fstream output;
	output.open(filename, std::fstream::out);

	output << std::endl;
	output << " --------------------------------------------- " << std::endl;
	output << " Combination summary                           " << std::endl;
	output << " --------------------------------------------- " << std::endl;
	output << std::endl;

	// get number of channels
	int nchannels = GetNChannels(); 

	// ---- single channel analysis without systematics ---- // 

	output << " Single channels without systematics : " << std::endl;
	output << std::endl;
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Channel";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Median";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (low)"; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (high)"; 
	output << std::endl;

	// loop over all channels
	for (int i = 0; i < nchannels; ++i) {

		// get summary
		ParameterSummary* ps = fSummaryChannelNoSyst.at(i); 

		if (!ps) {
			BCLog::OutError("CombinationModel::PrintChannelSummary : Single channel summary information not available."); 
			return -1; 
		}

		// get summary information
		double median = ps->GetMedian(); 
		double q16 = ps->GetQuantile16();
		double q84 = ps->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 

		output << " ";
		output << std::setw(15) << std::setiosflags(std::ios::left) << ps->GetName();
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << median; 
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh << std::endl;
	}
	output << std::endl;

	// ---- single channel analysis with systematics ---- // 

	output << " Single channels with systematics : " << std::endl;
	output << std::endl;
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Channel";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Median";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (low)"; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (high)"; 
	output << std::endl;

	// loop over all channels
	for (int i = 0; i < nchannels; ++i) {

		// get summary
		ParameterSummary* ps = fSummaryChannelSyst.at(i); 

		if (!ps) {
			BCLog::OutError("CombinationModel::PrintChannelSummary : Single channel summary information not available."); 
			return -1; 
		}

		// get summary information
		double median = ps->GetMedian(); 
		double q16 = ps->GetQuantile16();
		double q84 = ps->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 

		output << " ";
		output << std::setw(15) << std::setiosflags(std::ios::left) << ps->GetName();
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << median; 
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh << std::endl;
	}
	output << std::endl;

	// ---- combined analysis with and without systematics ---- // 

	output << " Combined channels : " << std::endl;
	output << std::endl;
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Systematics";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Median";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (low)"; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (high)"; 
	output << std::endl;

	// get summaries
	ParameterSummary* ps_nosyst = fSummaryCombinationNoSyst;
	ParameterSummary* ps_syst = fSummaryCombinationSyst;

	// get summary information
	double median_nosyst = ps_nosyst->GetMedian(); 
	double q16_nosyst = ps_nosyst->GetQuantile16();
	double q84_nosyst = ps_nosyst->GetQuantile84();
	double errlow_nosyst = median_nosyst - q16_nosyst;
	double errhigh_nosyst = q84_nosyst - median_nosyst; 

	output << " ";
	output << std::setw(15) << std::setiosflags(std::ios::left) << "off";
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << median_nosyst; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow_nosyst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh_nosyst << std::endl;

	// get summary information
	double median_syst = ps_syst->GetMedian(); 
	double q16_syst = ps_syst->GetQuantile16();
	double q84_syst = ps_syst->GetQuantile84();
	double errlow_syst = median_syst - q16_syst;
	double errhigh_syst = q84_syst - median_syst; 

	output << " ";
	output << std::setw(15) << std::setiosflags(std::ios::left) << "on";
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << median_syst; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow_syst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh_syst << std::endl;

	output << std::endl;
	output << " --------------------------------------------- " << std::endl;
	output << std::endl;

	// close stream
	output.close();

	// no error
	return 1;
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

