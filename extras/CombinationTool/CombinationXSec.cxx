#include "CombinationXSec.h"

#include <BAT/BCLog.h>
#include <BAT/BCMath.h>
#include <BAT/BCH1D.h> 

#include <TF1.h>
#include <TCanvas.h> 
#include <TGraphAsymmErrors.h>
#include <TH2D.h> 
#include <TPostScript.h> 
#include <TLatex.h> 
#include <TLine.h> 
#include <TBox.h> 

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
		double luminosity = parameters.at( fParIndexChannelLuminosity.at(i) ); 
		double efficiency = parameters.at( fParIndexChannelEfficiency.at(i) ); 
		double br = fChannelBR.at(i);
		double expectation = parameters.at(0) * luminosity * efficiency * br; 
		
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
		if (fChannelSignalPriorContainer.at(i)) {
			double luminosity = parameters.at( fParIndexChannelLuminosity.at(i) ); 
			double efficiency = parameters.at( fParIndexChannelEfficiency.at(i) ); 
			double br = fChannelBR.at(i);
			logprob += log( fChannelSignalPriorContainer.at(i)->Eval(parameters.at(0) * luminosity * efficiency * br) );
		}

		// add channel luminosty prior
		if (fChannelLuminosityPriorContainer.at(i))
			logprob += log( fChannelLuminosityPriorContainer.at(i)->Eval(parameters.at(fParIndexChannelLuminosity.at(i))) ); 

		// add channel efficiency prior
		if (fChannelEfficiencyPriorContainer.at(i))
			logprob += log( fChannelEfficiencyPriorContainer.at(i)->Eval(parameters.at(fParIndexChannelEfficiency.at(i))) ); 

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

	// add efficiency parameter
	AddParameter(Form("efficiency_%s", channelname), 0.0, 1.0);

	// add efficiency prior 
	fChannelEfficiencyPriorContainer.push_back(0); 

	// add parameter index
	fParIndexChannelEfficiency.push_back(GetNParameters()-1);

	// add Luminosity parameter
	AddParameter(Form("luminosity_%s", channelname), 0.0, 1.0);

	// add efficiency prior 
	fChannelLuminosityPriorContainer.push_back(0); 

	// add parameter index
	fParIndexChannelLuminosity.push_back(GetNParameters()-1);

	// add BR
	fChannelBR.push_back(1.0); 

	// no error 
	return 1; 
}

// ---------------------------------------------------------
int CombinationXSec::GetParIndexChannelEfficiency(const char* channelname)
{
 	// get channel container index
	int channelindex = GetContIndexChannel(channelname); 

	// check index
	if (channelindex < 0){ 
		return 1;
	}

	// return index
	return fParIndexChannelEfficiency.at(channelindex); 
}

// ---------------------------------------------------------
int CombinationXSec::GetParIndexChannelLuminosity(const char* channelname)
{
 	// get channel container index
	int channelindex = GetContIndexChannel(channelname); 

	// check index
	if (channelindex < 0){ 
		return 1;
	}

	// return index
	return fParIndexChannelLuminosity.at(channelindex); 
}

// ---------------------------------------------------------
int CombinationXSec::SetChannelEfficiencyPrior(const char* channelname, TF1* prior)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::SetChannelEfficiencyPrior : Can not find channel index."); 
		return 0; 
	}

	// set parameter range
	double effmin = 0.;
	double effmax = 1.;
	double x, y;
	prior->GetRange(x, y);
	if (x > effmin)
		effmin = x;
	if (y < effmax)
		effmax = y; 

	// rescale parameter ranges
	int parindex = GetParIndexChannelEfficiency(channelname); 
	GetParameter(parindex)->SetLowerLimit(effmin); 
	GetParameter(parindex)->SetUpperLimit(effmax); 
	fMCMCBoundaryMin[parindex] = effmin; 
	fMCMCBoundaryMax[parindex] = effmax; 

	// set prior
	fChannelEfficiencyPriorContainer[channelindex] = prior;

	// no error
	return 1;
}

// ---------------------------------------------------------
int CombinationXSec::SetChannelEfficiencyPriorGauss(const char* channelname, double mean, double sigma)
{
	// create new function
	TF1* f = new TF1(Form("f_eff_%s", channelname), "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(- (x-[0])*(x-[0])/2/[1]/[1] )");
	f->SetParameter(0, mean); 
	f->SetParameter(1, sigma); 

	// rescale parameter ranges
	int parindex = GetParIndexChannelEfficiency(channelname); 
	double effmin = TMath::Max(mean - 5.0*sigma, 0.);
	double effmax = TMath::Min(mean + 5.0*sigma, 1.);

	GetParameter(parindex)->SetLowerLimit(effmin); 
	GetParameter(parindex)->SetUpperLimit(effmax); 
	fMCMCBoundaryMin[parindex] = effmin; 
	fMCMCBoundaryMax[parindex] = effmax; 
	f->SetRange(effmin, effmax);

	// add function to container
	fFunctionContainer.push_back(f); 

	// set channel background prior
	return SetChannelEfficiencyPrior(channelname, f); 
}

// ---------------------------------------------------------
int CombinationXSec::SetChannelLuminosityPrior(const char* channelname, TF1* prior)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationModel::SetChannelLuminosityPrior : Can not find channel index."); 
		return 0; 
	}

	// set parameter range
	double lumimin = 0.;
	double lumimax = 1.;
	prior->GetRange(lumimin, lumimax);
	if (lumimin < 0.)
			lumimin = 0.;

	// rescale parameter ranges
	int parindex = GetParIndexChannelLuminosity(channelname); 
	GetParameter(parindex)->SetLowerLimit(lumimin); 
	GetParameter(parindex)->SetUpperLimit(lumimax); 
	fMCMCBoundaryMin[parindex] = lumimin; 
	fMCMCBoundaryMax[parindex] = lumimax; 

	// set prior
	fChannelLuminosityPriorContainer[channelindex] = prior;

	// no error
	return 1;
}

// ---------------------------------------------------------
int CombinationXSec::SetChannelLuminosityPriorGauss(const char* channelname, double mean, double sigma)
{
	// create new function
	TF1* f = new TF1(Form("f_eff_%s", channelname), "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(- (x-[0])*(x-[0])/2/[1]/[1] )");
	f->SetParameter(0, mean); 
	f->SetParameter(1, sigma); 

	// rescale parameter ranges
	int parindex = GetParIndexChannelLuminosity(channelname); 
	double effmin = TMath::Max(mean - 5.0*sigma, 0.);
	double effmax = mean + 5.0*sigma;

	GetParameter(parindex)->SetLowerLimit(effmin); 
	GetParameter(parindex)->SetUpperLimit(effmax); 
	fMCMCBoundaryMin[parindex] = effmin; 
	fMCMCBoundaryMax[parindex] = effmax; 
	f->SetRange(effmin, effmax);

	// add function to container
	fFunctionContainer.push_back(f); 

	// set channel background prior
	return SetChannelLuminosityPrior(channelname, f); 
}

// ---------------------------------------------------------
int CombinationXSec::SetChannelBR(const char* channelname, double BR)
{
	// get channel index
	int channelindex = GetContIndexChannel(channelname); 

	// check channel index
	if (channelindex < 0) {
		BCLog::OutError("CombinationXSec::SetChannelBR : Can not find channel index."); 
		return 0; 
	}

	// set efficiency
	fChannelBR[channelindex] = BR;

	// no error 
	return 1; 
}

// ---------------------------------------------------------
void CombinationXSec::MCMCUserIterationInterface()
{
	// debugKK
	return;

// 	// get number of channels
// 	int nchannels = GetNChannels(); 

// 	// loop over all channels 
// 	for (int i = 0; i < nchannels; ++i) {
// 		// calculate expectation for this channel
// 		double expectation = fMCMCx.at(0) * fChannelLuminosity.at(i) * fChannelEfficiency.at(i); 

// 		// fill histogram
// 		fChannelSignal.at(i)->GetHistogram()->Fill(expectation); 
// 	}
}

// ---------------------------------------------------------
void CombinationXSec::PrintChannelOverview(const char* filename)
{
	// get number of channels
	int nchannels = GetNChannels(); 

	// define vector of histograms 
	std::vector<BCH1D*> histos; 

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		
		// create analyses for each channel 
		CombinationXSec* model = new CombinationXSec( GetParameter(0)->GetName().c_str(), 
																									GetParameter(0)->GetLowerLimit(), 
																									GetParameter(0)->GetUpperLimit());
	
		// copy settings
		model->MCMCSetNLag( MCMCGetNLag() );
		model->MCMCSetNIterationsRun( MCMCGetNIterationsRun() );
		model->AddChannel( fChannelNameContainer.at(i).c_str() );
		model->SetChannelObservation( fChannelNameContainer.at(i).c_str(), 
																	fChannelObservation.at(i) );
		model->SetChannelEfficiencyPrior( fChannelNameContainer.at(i).c_str(),
																			fChannelEfficiencyPriorContainer.at(i) );
		model->SetChannelLuminosityPrior( fChannelNameContainer.at(i).c_str(),
																			fChannelLuminosityPriorContainer.at(i) );
		model->SetChannelBR( fChannelNameContainer.at(i).c_str(), 
												 fChannelBR.at(i) ); 

		int nbkg = int(fChannelBackgroundNameContainer.at(i).size());
		for (int j = 0; j < nbkg; ++j) {
			int parindex = GetParIndexChannelBackground(i, j);
			double xmin = GetParameter(parindex)->GetLowerLimit();
			double xmax = GetParameter(parindex)->GetUpperLimit();
			model->AddChannelBackground( fChannelNameContainer.at(i).c_str(),
																	 fChannelBackgroundNameContainer.at(i).at(j).c_str(),
																	 xmin, 
																	 xmax );
			model->SetChannelBackgroundPrior( fChannelNameContainer.at(i).c_str(), 
																				fChannelBackgroundNameContainer.at(i).at(j).c_str(),
																				fChannelBackgroundPriorContainer.at(i).at(j)); 
		}

		// debugKK
		// add systematics and all other settings here

		// perform analysis
		model->MarginalizeAll();
		model->FindMode( model->GetBestFitParameters() );

		// debugKK
		//		model->PrintAllMarginalized(Form("channel_%i_plots.ps", i)); 

		// save signal histogram for summary
		histos.push_back( model->GetMarginalized( model->GetParameter(0)->GetName().c_str() ) );
	}

	// create graph
	TGraphAsymmErrors* g = new TGraphAsymmErrors(GetNChannels()+1); 
	g->SetMarkerStyle(20); 
	g->SetMarkerSize(1); 

	// coordinate system
	double xmin = 0.0; 
	double xmax = 0.0; 
	double xwidth = 0.0;
	double ymin = -0.5;
	double ymax = double(nchannels)+0.5;

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		// get histogram
		BCH1D* h = histos.at(i);
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
		g->SetPoint(i, median, double(nchannels-i));
		g->SetPointError(i, errlow, errhigh, 0, 0);
	}

	// set point from combination
	BCH1D* h = GetMarginalized( GetParameter(0)->GetName().c_str() );
	double median = h->GetMedian(); 
	double q16 = h->GetQuantile(0.16);
	double q84 = h->GetQuantile(0.84);
	double errlow = median - q16;
	double errhigh = q84 - median; 
	
	// set point and error 
	g->SetPoint(nchannels, median, 0.);
	g->SetPointError(nchannels, errlow, errhigh, 0, 0);
	
	// create histogram for axes
	TH2D* hist_axes = new TH2D("", Form(";%s;Channel", GetParameter(0)->GetName().c_str()), 1, xmin - 0.25*xwidth, xmax + 1.75 * xwidth, nchannels+1, ymin, ymax); 
	hist_axes->SetStats(kFALSE);
	hist_axes->GetYaxis()->SetNdivisions(0);
	hist_axes->GetYaxis()->SetTitleOffset(1.0);

	// create canvas
	TCanvas* c1 = new TCanvas();
	c1->cd(); 

	// create latex
	TLatex* l = new TLatex(); 
	l->SetTextSize(0.04);
	if (nchannels>=10)
		l->SetTextSize(0.02);
	l->SetTextAlign(12);

	// create lines
	TLine* line_mode = new TLine(); 
	line_mode->SetLineWidth(1);
	line_mode->SetLineStyle(1);
	line_mode->SetLineColor(kRed);
	TBox* box = new TBox();
	box->SetLineWidth(0); 
	box->SetFillColor(kYellow);
	box->SetFillStyle(1001);
	TBox* box_help = new TBox();
	box_help->SetLineWidth(1); 
	box_help->SetFillColor(kBlack);
	box_help->SetFillStyle(0);

	// draw
	hist_axes->Draw();
	box->DrawBox(q16, ymin, q84, ymax);
	line_mode->DrawLine(median, ymin, median, ymax);
	g->Draw("SAMEP"); 
	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		l->DrawLatex(xmax + 0.25*xwidth, double(nchannels-i), fChannelNameContainer.at(i).c_str());
	}
	l->DrawLatex(xmax + 0.25*xwidth, 0., "combination");
	hist_axes->Draw("SAMEAXIS");
	box_help->DrawBox(xmin - 0.25*xwidth, ymin, xmax + 1.75 * xwidth, ymax);

	// print to file 
	c1->Print(filename);

}

// ---------------------------------------------------------

