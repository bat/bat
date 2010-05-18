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
		double luminosity = fChannelLuminosity.at(i);
		double efficiency = fChannelEfficiency.at(i);  
		double br = fChannelBR.at(i);
		double exp_cross = parameters.at(0) * luminosity * efficiency * br;
		double expectation = exp_cross; 

		if (fFlagSystErrors) {
			// loop over all systematic uncertainties
			for (int k = 0; k < nsysterrors; ++k) {
				int systparindex = GetParIndexSystError(k);
				double par = parameters.at(systparindex);
				
				if (par < 0)
					expectation += fSystErrorChannelSigmaDownContainer.at(k).at(i) * par * exp_cross;
				else
					expectation += fSystErrorChannelSigmaUpContainer.at(k).at(i) * par * exp_cross;
			}
		}
		
		// get number of background sources in this channel
		int nbackground = GetNChannelBackgrounds(i);

		// loop over all background sources
		for (int j = 0; j < nbackground; ++j) {
			// get parameter index
			double exp_bkg = fChannelBackground.at(i).at(j);
			expectation += exp_bkg;

			if (fFlagSystErrors) {
				// loop over all systematic uncertainties
				for (int k = 0; k < nsysterrors; ++k) {
					int systparindex = GetParIndexSystError(k);
					double par = parameters.at(systparindex);
					if (par < 0)
						expectation += fSystErrorSigmaDownContainer.at(k).at(i).at(j) * par * exp_bkg;
					else
						expectation += fSystErrorSigmaUpContainer.at(k).at(i).at(j) * par * exp_bkg;
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
			double luminosity = fChannelLuminosity.at(i);
			double efficiency = fChannelEfficiency.at(i);
			double br = fChannelBR.at(i);
			logprob += log( fChannelSignalPriorContainer.at(i)->Eval(parameters.at(0) * luminosity * efficiency * br) );
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

	// set efficiency
	fChannelEfficiency.push_back(1.); 

	// set luminosity
	fChannelLuminosity.push_back(1.);

	// set branching ratio
	fChannelBR.push_back(1.); 

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
		BCLog::OutError("CombinationModel::SetChannelEfficiency : Can not find channel index."); 
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
		BCLog::OutError("CombinationModel::SetChannelLuminosity : Can not find channel index."); 
		return 0; 
	}

	// set luminosity
	fChannelLuminosity[channelindex] = luminosity;

	// no error
	return 1;
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
	model->SetChannelEfficiency( fChannelNameContainer.at(channelindex).c_str(),
															 fChannelEfficiency.at(channelindex) );
	model->SetChannelLuminosity( fChannelNameContainer.at(channelindex).c_str(),
															 fChannelLuminosity.at(channelindex) );
	model->SetChannelBR( fChannelNameContainer.at(channelindex).c_str(), 
											 fChannelBR.at(channelindex) ); 
	
	int nbkg = int(fChannelBackgroundNameContainer.at(channelindex).size());
	for (int j = 0; j < nbkg; ++j) {
		model->AddChannelBackground( fChannelNameContainer.at(channelindex).c_str(),
																 fChannelBackgroundNameContainer.at(channelindex).at(j).c_str(),
																 fChannelBackground.at(channelindex).at(j) );
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
	ps.SetGlobalMode( model->GetBestFitParameter(0) );

	// free memory
	delete model;

	// return summary 
	return ps; 
}

// ---------------------------------------------------------

