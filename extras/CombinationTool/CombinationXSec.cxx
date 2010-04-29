#include "CombinationXSec.h"

#include <BAT/BCLog.h>
#include <BAT/BCMath.h>
#include <BAT/BCH1D.h> 

#include <TF1.h>
#include <TCanvas.h> 
#include <TGraphAsymmErrors.h>
#include <TH2D.h> 
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

	// get number of systematic uncertainties
	int nsysterrors = GetNSystErrors();

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

			if (fFlagSystErrors) {
				// loop over all systematic uncertainties
				for (int k = 0; k < nsysterrors; ++k) {
					int systparindex = GetParIndexSystError(k);
					double par = parameters.at(systparindex);
					if (par < 0)
						expectation += fSystErrorSigmaDownContainer.at(k).at(i).at(j) * par;
					else
						expectation += fSystErrorSigmaUpContainer.at(k).at(i).at(j) * par;
				}
			}
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
			if (fChannelLuminosityPriorContainer.at(i)->GetParameter(1) > 0)
				logprob += log( fChannelLuminosityPriorContainer.at(i)->Eval(parameters.at(fParIndexChannelLuminosity.at(i))) ); 

		// add channel efficiency prior
		if (fChannelEfficiencyPriorContainer.at(i))
			if (fChannelEfficiencyPriorContainer.at(i)->GetParameter(1) > 0)
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

	// loop over all systematic uncertainties
	int nsysterrors = GetNSystErrors();

	if (fFlagSystErrors) {
		for (int i = 0; i < nsysterrors; ++i) {
			logprob += BCMath::LogGaus(parameters.at( GetParIndexSystError(i) ) ); 
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
ParameterSummary CombinationXSec::PerformSingleChannelAnalysis(const char* channelname, bool flag_syst)
{
 	// get channel container index
	int channelindex = GetContIndexChannel(channelname); 

	// define summary
	ParameterSummary ps(channelname); 

	// create analyses for each channel 
	CombinationXSec* model = new CombinationXSec( GetParameter(0)->GetName().c_str(), 
																								GetParameter(0)->GetLowerLimit(), 
																								GetParameter(0)->GetUpperLimit());

	// copy settings
	model->MCMCSetNLag( MCMCGetNLag() );
	model->MCMCSetNIterationsRun( MCMCGetNIterationsRun() );
	model->SetFlagSystErrors(flag_syst);
	model->AddChannel( fChannelNameContainer.at(channelindex).c_str() );
	model->SetChannelObservation( fChannelNameContainer.at(channelindex).c_str(), 
																fChannelObservation.at(channelindex) );
	model->SetChannelEfficiencyPrior( fChannelNameContainer.at(channelindex).c_str(),
																		fChannelEfficiencyPriorContainer.at(channelindex) );
	model->SetChannelLuminosityPrior( fChannelNameContainer.at(channelindex).c_str(),
																		fChannelLuminosityPriorContainer.at(channelindex) );
	model->SetChannelBR( fChannelNameContainer.at(channelindex).c_str(), 
											 fChannelBR.at(channelindex) ); 
	
	int nbkg = int(fChannelBackgroundNameContainer.at(channelindex).size());
	for (int j = 0; j < nbkg; ++j) {
		int parindex = GetParIndexChannelBackground(channelindex, j);
		double xmin = GetParameter(parindex)->GetLowerLimit();
		double xmax = GetParameter(parindex)->GetUpperLimit();
		model->AddChannelBackground( fChannelNameContainer.at(channelindex).c_str(),
																 fChannelBackgroundNameContainer.at(channelindex).at(j).c_str(),
																 xmin, 
																 xmax );
		model->SetChannelBackgroundPrior( fChannelNameContainer.at(channelindex).c_str(), 
																			fChannelBackgroundNameContainer.at(channelindex).at(j).c_str(),
																			fChannelBackgroundPriorContainer.at(channelindex).at(j)); 
	}
	
	if (flag_syst) {
		int nsyst = GetNSystErrors();
		for (int i = 0; i < nsyst; ++i) {
			model->AddSystError( fSystErrorNameContainer.at(i).c_str() ); 
			
			for (int j = 0; j < nbkg; ++j) {
				model->SetSystErrorChannelBackground( fSystErrorNameContainer.at(i).c_str(), 
																							fChannelNameContainer.at(channelindex).c_str(),
																							fChannelBackgroundNameContainer.at(channelindex).at(j).c_str(), 
																							fSystErrorSigmaDownContainer.at(i).at(channelindex).at(j), 
																							fSystErrorSigmaUpContainer.at(i).at(channelindex).at(j) );
			}
		}
	}
	
	// perform analysis
	model->PerformAnalysis();

	// get histogram
	BCH1D* hist = model->GetMarginalized( model->GetParameter(0)->GetName().c_str() ); 

	// copy summary information
	ps.Summarize(hist);

	// free memory
	delete model;

	// return summary 
	return ps; 
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

