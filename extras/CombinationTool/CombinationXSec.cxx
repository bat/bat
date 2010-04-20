#include "CombinationXSec.h"

#include <BAT/BCLog.h>
#include <BAT/BCMath.h>
#include <BAT/BCH1D.h> 

#include <TF1.h>
#include <TCanvas.h> 
#include <TGraphAsymmErrors.h>
#include <TH2D.h> 
#include <TPostScript.h> 

// ---------------------------------------------------------
CombinationXSec::CombinationXSec(const char * name, double xmin, double xmax) : CombinationModel(name, xmin, xmax)
{};

// ---------------------------------------------------------
CombinationXSec::~CombinationXSec()
{};

// ---------------------------------------------------------
double CombinationXSec::LogLikelihood(std::vector <double> parameters)
{
	double logprob = 0.;

	// get number of channels
	int nchannels = GetNChannels(); 

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		// calculate expectation for this channel
		double expectation = parameters.at(0) * fChannelLuminosity.at(i) * fChannelEfficiency.at(i); 
		
		// get number of background sources in this channel
		int nbackground = GetNChannelBackgrounds(i);

		// loop over all background sources
		for (int j = 0; j < nbackground; ++j) {
			// get parameter index
			int parindex = GetParIndexChannelBackground(i, j); 
			expectation += parameters.at(parindex);
		}
		
		// calculate Poisson term for this channel
		logprob += BCMath::LogPoisson( fChannelObservation.at(i), expectation); 
	}

	return logprob;
}

// ---------------------------------------------------------
double CombinationXSec::LogAPrioriProbability(std::vector <double> parameters)
{
	double logprob = 0.;

	// get number of channels
	int nchannels = GetNChannels(); 

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {

		// add channel signal prior
		if (fChannelSignalPriorContainer.at(i))
			logprob += log( fChannelSignalPriorContainer.at(i)->Eval(parameters.at(0) * fChannelLuminosity.at(i) * fChannelEfficiency.at(i)) );

		// get number of background sources in this channel
		int nbackground = GetNChannelBackgrounds(i);

		// loop over all background sources
		for (int j = 0; j < nbackground; ++j) {
			// get parameter index
			int parindex = GetParIndexChannelBackground(i, j); 

			// add channel background prior
			if (fChannelBackgroundPriorContainer.at(i).at(j))
				logprob += log( fChannelBackgroundPriorContainer.at(i).at(j)->Eval(parameters.at(parindex)) ); 
		}

	}

	return logprob;
}

// ---------------------------------------------------------
int CombinationXSec::AddChannel(const char* channelname)
{
	int errcode = 1; 

	// add channel 
	errcode = CombinationModel::AddChannel(channelname);

	if (errcode != 1)
		return errcode; 

	// add efficiency
	fChannelEfficiency.push_back(1.0); 
	fChannelEfficiencyPriorContainer.push_back(0); 

	// add luminosity
	fChannelLuminosity.push_back(1.0); 
	fChannelLuminosityPriorContainer.push_back(0); 

	// add histogram
	TH1D* hist = new TH1D(Form("sigma_%s", channelname), Form(";#sigma_{%s};p(#sigma_{%s})", channelname, channelname),
												100, 
												GetParameter(0)->GetLowerLimit(), 
												GetParameter(0)->GetUpperLimit());
	BCH1D* bchist = new BCH1D(hist); 
	fChannelSignal.push_back(bchist);

	// no error 
	return 1; 
}

// ---------------------------------------------------------
int CombinationXSec::SetChannelEfficiency(const char* channelname, double efficiency)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationXSec::SetChannelEfficiency : Can not find channel index."); 
		return 0; 
	}

	// set efficiency
	fChannelEfficiency[channelindex] = efficiency;

	// no error 
	return 1; 
}

// ---------------------------------------------------------
int CombinationXSec::SetChannelLuminosity(const char* channelname, double luminosity)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationXSec::SetChannelEfficiency : Can not find channel index."); 
		return 0; 
	}

	// set efficiency
	fChannelLuminosity[channelindex] = luminosity;

	// no error 
	return 1; 
}

// ---------------------------------------------------------
void CombinationXSec::MCMCUserIterationInterface()
{
	// debugKK
	return;

	// get number of channels
	int nchannels = GetNChannels(); 

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		// calculate expectation for this channel
		double expectation = fMCMCx.at(0) * fChannelLuminosity.at(i) * fChannelEfficiency.at(i); 

		// fill histogram
		fChannelSignal.at(i)->GetHistogram()->Fill(expectation); 
	}
}

// ---------------------------------------------------------
void CombinationXSec::PrintChannelOverview(const char* filename)
{
	// create graph
	TGraphAsymmErrors* g = new TGraphAsymmErrors(GetNChannels()+1); 
	g->SetMarkerStyle(20); 
	g->SetMarkerSize(1); 

	// get number of channels
	int nchannels = GetNChannels(); 

	TH2D* hist_axes = new TH2D("", ";#sigma;", 1, GetParameter(0)->GetLowerLimit(), GetParameter(0)->GetUpperLimit(), nchannels+1, -0.5, double(nchannels)+1.5); 
	hist_axes->SetStats(kFALSE);

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		// get histogram
		BCH1D* h = fChannelSignal.at(i); 
		g->SetPoint(i, h->GetMode(), double(nchannels-i+1));
	}

	TCanvas* c1 = new TCanvas();
	c1->cd(); 
	
	hist_axes->Draw();
	g->Draw("SAME"); 
	c1->Print(filename);
}

// ---------------------------------------------------------
void CombinationXSec::PrintChannels(const char* filename)
{
	TCanvas * c1 = new TCanvas("");

	TPostScript * ps = new TPostScript(filename, 112);
	ps->NewPage();

	c1->cd();
	for (int i = 0; i < GetNChannels(); ++i) {
		c1->Update();
		ps->NewPage();
		c1->cd();
		BCH1D* h = fChannelSignal.at(i); 
		h->GetHistogram()->Draw();
	}
	c1->Update();
	ps->Close();

	delete c1;
	delete ps;
}

// ---------------------------------------------------------

// to do: 
// add efficiency and luminosity as free parameter
