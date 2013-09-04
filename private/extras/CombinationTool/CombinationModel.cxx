#include "CombinationModel.h"
#include "ParameterSummary.h"

#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCMath.h>

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
double CombinationModel::LogAPrioriProbability(std::vector<double> parameters)
{
	double logprob = 0.;

	// get number of channels
	int nchannels = GetNChannels(); 

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {

		// get parameter
		double x = parameters.at(0);
		
		// systematics
		if (fFlagSystErrors) {
			// loop over all systematic uncertainties
			int nsysterrors = GetNSystErrors();
			
			// loop over all systematic uncertainties
			for (int k = 0; k < nsysterrors; ++k) {

				// check status
				if (fSystErrorStatusContainer.at(k)) {
					int systparindex = GetParIndexSystError(k);
					double par = parameters.at(systparindex);
				
					x -= fSystErrorChannelShiftContainer.at(k).at(i);
					
					if (par < 0)
						x -= fSystErrorChannelSigmaDownContainer.at(k).at(i) * par;
					else
						x -= fSystErrorChannelSigmaUpContainer.at(k).at(i) * par;
				} // end of check status
			} // end of loop over systematics
		} // end if systematics 

		// add measurements
		logprob += log( fChannelSignalPriorContainer.at(i)->Eval(x) );
	}

	return logprob;
}

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(std::vector<double> parameters)
{
	double logprob = 0.;

	// loop over all systematic uncertainties
	int nsysterrors = GetNSystErrors();

	if (fFlagSystErrors) {
	  for (int j = 0; j < nsysterrors; ++j) {
	    logprob += BCMath::LogGaus(parameters.at( GetParIndexSystError(j) ), 0., 1., true ); 
	  }
	}
	
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
int CombinationModel::AddChannel(const char* channelname)
{
	// add channel name to container
	fChannelNameContainer.push_back(channelname); 

	// add background container
	fChannelBackgroundNameContainer.push_back(std::vector<std::string>(0));
	fChannelBackground.push_back(std::vector<double>(0)); 

	// add signal prior to container
	fChannelSignalPriorContainer.push_back(0); 

	// add observation container
	fChannelObservation.push_back(0);

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::AddChannelBackground(const char* channelname, const char* backgroundname, double bkg)
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

	// add contribution
	fChannelBackground.at(channelindex).push_back(bkg); 

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::AddSystError(const char* systerrorname)
{
	// add syst error name to container
	fSystErrorNameContainer.push_back(systerrorname); 
	
	// add status "on"
	fSystErrorStatusContainer.push_back(true);

	// add parameter
	AddParameter(systerrorname, -4.0, 4.0);

	// add parameter index
	fParIndexSystErrorContainer.push_back(GetNParameters()-1);

	// get number of channels
	int nchannels = GetNChannels(); 

	std::vector< std::vector<double> > sdown;
	std::vector< std::vector<double> > sup;
	std::vector< std::vector<double> > sshift;

	// loop over all channels 
	for (int i = 0; i < nchannels; ++i) {
		// get number of background sources in this channel
		int nbackground = GetNChannelBackgrounds(i);
		
		std::vector<double> schanneldown(nbackground);
		std::vector<double> schannelup(nbackground);
		std::vector<double> schannelshift(nbackground);

		// loop over all background sources
		for (int j = 0; j < nbackground; ++j) {
			schanneldown[j] = 0.0;
			schannelup[j] = 0.0;
			schannelshift[j] = 0.0;
		}

		sdown.push_back(schanneldown);
		sup.push_back(schannelup);
		sshift.push_back(schannelup);
	}	

	// add vectors to container
	fSystErrorSigmaDownContainer.push_back(sdown);
	fSystErrorSigmaUpContainer.push_back(sup);
	fSystErrorShiftContainer.push_back(sshift);

	// add systematic uncertainties for signal
	std::vector< double > sigma_signal_down;
	sigma_signal_down.assign(nchannels, 0.0);
	std::vector< double > sigma_signal_up(nchannels);
	sigma_signal_up.assign(nchannels, 0.0);
	std::vector< double > sigma_signal_shift(nchannels);
	sigma_signal_shift.assign(nchannels, 0.0);

	// add to container
	fSystErrorChannelSigmaDownContainer.push_back(sigma_signal_down);
	fSystErrorChannelSigmaUpContainer.push_back(sigma_signal_up);
	fSystErrorChannelShiftContainer.push_back(sigma_signal_shift);

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::SetSystErrorChannelSignal(const char* systerrorname, const char* channelname, double sigmadown, double sigmaup, double shift)
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
	fSystErrorChannelShiftContainer.at(systerrorindex)[channelindex] = shift;

	// no error 
	return 1;
}

// ---------------------------------------------------------
int CombinationModel::SetSystErrorChannelBackground(const char* systerrorname, const char* channelname, const char* backgroundname, double sigmadown, double sigmaup, double shift)
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
	fSystErrorShiftContainer.at(systerrorindex).at(channelindex)[backgroundindex] = shift;
	
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
int CombinationModel::PerformFullAnalysis(int index_syst)
{
	// ---- perform combination analysis without systematics ---- //

	// print to screen
	BCLog::OutSummary("Perform combination without systematics"); 

	// copy flags
	bool flagsysterrors = GetFlagSystErrors();
	int nsyst = int(fSystErrorStatusContainer.size());
	std::vector<bool> systerrorstatuscontainer(nsyst);
	for (int i = 0; i < nsyst; ++i) {
		systerrorstatuscontainer[i]=fSystErrorStatusContainer.at(i);
	}

	// set flags
	for (int i = 0; i < nsyst; ++i) {
		fSystErrorStatusContainer[i]=false;
	}


	// perform analysis
	PerformAnalysis();

	// delete old summary information
	delete fSummaryCombinationNoSyst; 

	// copy summary information
	fSummaryCombinationNoSyst = new ParameterSummary("combination");
	fSummaryCombinationNoSyst->Summarize( GetMarginalized( GetParameter(0)->GetName().c_str() ));
	fSummaryCombinationNoSyst->SetGlobalMode( GetBestFitParameter(0) );

	// ---- perform combination analysis with systematics ---- //

	// print to screen
	BCLog::OutSummary("Perform combination with systematics"); 

	// set flags
	for (int i = 0; i < nsyst; ++i) {
		fSystErrorStatusContainer[i]=true;
	}

	// perform analysis
	PerformAnalysis();

	// delete old summary information
	delete fSummaryCombinationSyst; 

	// copy summary information
	fSummaryCombinationSyst = new ParameterSummary("combination");
	fSummaryCombinationSyst->Summarize( GetMarginalized( GetParameter(0)->GetName().c_str() ));
	fSummaryCombinationSyst->SetGlobalMode( GetBestFitParameter(0) );

	// ---- perform combination analysis with each systematic separately ---- //

	// print to screen
	BCLog::OutSummary("Perform combination with each systematic separately"); 

	// set flags
	SetFlagSystErrors(true);

	// delete old summaries
	for (int i = 0; i < int(fSummaryCombinationSingleSyst.size()); ++i) 
		delete fSummaryCombinationSingleSyst[i]; 
	fSummaryCombinationSingleSyst.clear(); 

	// loop over all systematic uncertainties
	int nsysterrors = GetNSystErrors();
	for (int k = 0; k < nsysterrors; ++k) {

		// loop over all systematic uncertainties and set status of syst. errors to false
		for (int j = 0; j < nsysterrors; ++j) {
			fSystErrorStatusContainer[j] = false;
		}

		// set status to true
		fSystErrorStatusContainer[k] = true;

		// perform analysis
		PerformAnalysis();

		// copy summary information
		ParameterSummary* ps = new ParameterSummary(fSystErrorNameContainer.at(k).c_str());
		ps->Summarize( GetMarginalized( GetParameter(0)->GetName().c_str() ));
		ps->SetGlobalMode( GetBestFitParameter(0) );
		ps->SetBuffer( GetMarginalized( GetParameter(0)->GetName().c_str(), GetParameter(fParIndexSystErrorContainer.at(k))->GetName().c_str())->GetHistogram()->GetCorrelationFactor() );
		fSummaryCombinationSingleSyst.push_back(ps);
	}

	// loop over all systematic uncertainties and set status of syst. errors to false
	for (int j = 0; j < nsysterrors; ++j) {
		fSystErrorStatusContainer[j] = false;
	}

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

	// reset flags
	SetFlagSystErrors(flagsysterrors);
	for (int i = 0; i < nsyst; ++i) {
		fSystErrorStatusContainer[i] = systerrorstatuscontainer.at(i);
	}

	// no error 
	return 1;
}

// ---------------------------------------------------------
ParameterSummary CombinationModel::PerformSingleChannelAnalysis(const char* channelname, bool flag_syst, int index_syst)
{
 	// get channel container index
	int channelindex = GetContIndexChannel(channelname); 

	// define summary
	ParameterSummary ps(channelname); 

	// create analyses for each channel 
	CombinationModel* model = new CombinationModel( GetParameter(0)->GetName().c_str(), 
																									GetParameter(0)->GetLowerLimit(), 
																									GetParameter(0)->GetUpperLimit());

	// copy settings
	model->MCMCSetNChains( MCMCGetNChains() );
	model->MCMCSetNLag( MCMCGetNLag() );
	model->MCMCSetNIterationsMax( MCMCGetNIterationsMax() );
	model->MCMCSetNIterationsRun( MCMCGetNIterationsRun() );
	model->MCMCSetNIterationsPreRunMin( MCMCGetNIterationsPreRunMin() ); 
	model->MCMCSetNIterationsUpdate( MCMCGetNIterationsUpdate() );
	model->MCMCSetNIterationsUpdateMax( MCMCGetNIterationsUpdateMax() );
	model->MCMCSetRValueCriterion( MCMCGetRValueCriterion() );
	model->MCMCSetRValueParametersCriterion( MCMCGetRValueParametersCriterion() );

	model->SetFlagSystErrors(flag_syst);
	model->AddChannel( fChannelNameContainer.at(channelindex).c_str() );
	model->SetChannelSignalPrior( channelname, fChannelSignalPriorContainer.at(channelindex) ); 
	
	int nbkg = int(fChannelBackgroundNameContainer.at(channelindex).size());
	for (int j = 0; j < nbkg; ++j) {
		model->AddChannelBackground( fChannelNameContainer.at(channelindex).c_str(),
																 fChannelBackgroundNameContainer.at(channelindex).at(j).c_str(),
																 fChannelBackground.at(channelindex).at(j) );
	}
	
	if (flag_syst) {
		int nsyst = GetNSystErrors();
		for (int i = 0; i < nsyst; ++i) {
			if (index_syst < 0 || index_syst == i) {
				model->AddSystError( fSystErrorNameContainer.at(i).c_str() ); 
				
				model->SetSystErrorChannelSignal( fSystErrorNameContainer.at(i).c_str(),
																					fChannelNameContainer.at(channelindex).c_str(),
																					fSystErrorChannelSigmaDownContainer.at(i).at(channelindex), 
																					fSystErrorChannelSigmaUpContainer.at(i).at(channelindex),
																					fSystErrorChannelShiftContainer.at(i).at(channelindex) );
				
				
				for (int j = 0; j < nbkg; ++j) {
					model->SetSystErrorChannelBackground( fSystErrorNameContainer.at(i).c_str(), 
																								fChannelNameContainer.at(channelindex).c_str(),
																								fChannelBackgroundNameContainer.at(channelindex).at(j).c_str(), 
																								fSystErrorSigmaDownContainer.at(i).at(channelindex).at(j), 
																								fSystErrorSigmaUpContainer.at(i).at(channelindex).at(j),
																								fSystErrorShiftContainer.at(i).at(channelindex).at(j) );
				}
			}
		}
	}
	
	// set default histogram binning to the one of the original model
	for (int i = 0; i < int(model->GetNParameters()); ++i) {
		// this construct has to go here, because otherwise there is a
		// warning from BCEngineMCMC:: MCMCGetH1Marginalized
		if (model->GetBestFitParameters().size() > 0){
			BCH1D* hist = model->GetMarginalized( model->GetParameter(i) );
			if (hist) {
				int nbins = hist->GetHistogram()->GetNbinsX();
				SetNbins( (model->GetParameter(i)->GetName()).c_str(), nbins);
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
int CombinationModel::GetParIndexSystError(int systerrorindex)
{
	return fParIndexSystErrorContainer.at(systerrorindex);
}

// ---------------------------------------------------------
int CombinationModel::PrintChannelOverview(const char* filename1, const char* filename2)
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

	TGraphAsymmErrors* graph_singlesyst = new TGraphAsymmErrors(GetNSystErrors()+2); 
	graph_singlesyst->SetMarkerStyle(20); 
	graph_singlesyst->SetMarkerSize(1); 

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

	// ---- combined analysis with each systematic separately

	// loop over all systematic uncertainties
	int nsysterrors = GetNSystErrors();
			
	// loop over all systematic uncertainties
	for (int k = 0; k < nsysterrors; ++k) {

		// get summary
		ParameterSummary* ps = fSummaryCombinationSingleSyst.at(k); 

		if (!ps) {
			BCLog::OutError("CombinationModel::PrintChannelOverview : Combined channel summary information not available."); 
			return -1; 
		}

		// set point from combination
		double median = ps ->GetMedian(); 
		double q16 = ps ->GetQuantile16();
		double q84 = ps ->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 
	
		// set point and error 
		graph_singlesyst->SetPoint(k, median, k);
		graph_singlesyst->SetPointError(k, errlow, errhigh, 0, 0);
	}
	graph_singlesyst->SetPoint(nsysterrors, median_nosyst, nsysterrors);
	graph_singlesyst->SetPointError(nsysterrors, errlow_nosyst, errhigh_nosyst, 0, 0);
	graph_singlesyst->SetPoint(nsysterrors+1, median_syst, nsysterrors+1.);
	graph_singlesyst->SetPointError(nsysterrors+1, errlow_syst, errhigh_syst, 0, 0);

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
	c1->Print(filename1);

	// create histogram for axes
	TH2D* hist_axes2 = new TH2D("", Form(";%s;Systematic uncertainty", GetParameter(0)->GetName().c_str()), 1, xmin - 0.25*xwidth, xmax + 1.75 * xwidth, GetNSystErrors() + 2, -0.5, GetNSystErrors() - 0.5 + 2.0); 
	hist_axes2->SetStats(kFALSE);
	hist_axes2->GetYaxis()->SetNdivisions(0);
	hist_axes2->GetYaxis()->SetTitleOffset(1.0);

	// create canvas
	if (GetNSystErrors() > 0) {
		TCanvas* c2 = new TCanvas();
		c2->cd(); 

		hist_axes2->Draw();
		graph_singlesyst->Draw("SAMEP");

		// loop over all systematics and draw labels
		for (int i = 0; i < GetNSystErrors(); ++i) {
			latex->DrawLatex(xmax + 0.25*xwidth, i, fSystErrorNameContainer.at(i).c_str());
		}
		latex->DrawLatex(xmax + 0.25*xwidth, GetNSystErrors(), "no systematics");
		latex->DrawLatex(xmax + 0.25*xwidth, GetNSystErrors()+1., "all systematics");
		line_median->DrawLine(median_syst, -0.5, median_syst, GetNSystErrors() - 0.5 + 2.0);
		line_median->DrawLine(median_nosyst, -0.5, median_nosyst, GetNSystErrors() - 0.5 + 2.0);

		c2->Print(filename2);
	}
	
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
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Gl. mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Median";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " RMS";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 16% quantile";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 84% quantile";
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
		double glmode = ps->GetGlobalMode(); 
		double mode = ps->GetMode(); 
		double rms = ps->GetRMS(); 
		double q16 = ps->GetQuantile16();
		double q84 = ps->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 

		output << " ";
		output << std::setw(15) << std::setiosflags(std::ios::left) << ps->GetName();
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << glmode; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << median; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << mode; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << rms; 
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q16;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q84;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh << std::endl;
	}
	output << std::endl;

	// ---- single channel analysis with systematics ---- // 

	output << " Single channels with systematics : " << std::endl;
	output << std::endl;
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Channel";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Gl. mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Median";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " RMS";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 16% quantile";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 84% quantile";
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
		double glmode = ps->GetGlobalMode(); 
		double mode = ps->GetMode(); 
		double rms = ps->GetRMS(); 
		double q16 = ps->GetQuantile16();
		double q84 = ps->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 

		output << " ";
		output << std::setw(15) << std::setiosflags(std::ios::left) << ps->GetName();
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << glmode; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << median; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << mode; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << rms; 
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q16;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q84;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh << std::endl;
	}
	output << std::endl;

	// ---- combined analysis with and without systematics ---- // 

	output << " Combined channels : " << std::endl;
	output << std::endl;
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Systematics";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Gl. mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Median";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " RMS";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 16% quantile";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 84% quantile";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (low)"; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (high)"; 
	output << std::endl;

	// get summaries
	ParameterSummary* ps_nosyst = fSummaryCombinationNoSyst;
	ParameterSummary* ps_syst = fSummaryCombinationSyst;

	// get summary information
	double median_nosyst = ps_nosyst->GetMedian(); 
	double glmode_nosyst = ps_nosyst->GetGlobalMode(); 
	double mode_nosyst = ps_nosyst->GetMode(); 
	double rms_nosyst = ps_nosyst->GetRMS(); 
	double q16_nosyst = ps_nosyst->GetQuantile16();
	double q84_nosyst = ps_nosyst->GetQuantile84();
	double errlow_nosyst = median_nosyst - q16_nosyst;
	double errhigh_nosyst = q84_nosyst - median_nosyst; 

	output << " ";
	output << std::setw(15) << std::setiosflags(std::ios::left) << "off";
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << glmode_nosyst; 
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << median_nosyst; 
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << mode_nosyst; 
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << rms_nosyst; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q16_nosyst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q84_nosyst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow_nosyst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh_nosyst << std::endl;

	// get summary information
	double glmode_syst = ps_syst->GetGlobalMode(); 
	double median_syst = ps_syst->GetMedian(); 
	double mode_syst = ps_syst->GetMode(); 
	double rms_syst = ps_syst->GetRMS(); 
	double q16_syst = ps_syst->GetQuantile16();
	double q84_syst = ps_syst->GetQuantile84();
	double errlow_syst = median_syst - q16_syst;
	double errhigh_syst = q84_syst - median_syst; 

	output << " ";
	output << std::setw(15) << std::setiosflags(std::ios::left) << "on";
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << glmode_syst; 
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << median_syst; 
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << mode_syst; 
	output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << rms_syst; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q16_syst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q84_syst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow_syst;
	output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh_syst << std::endl;
	output << std::endl;

	// ---- combined analysis with each systematic separately ---- // 

	output << " Combined analysis with each systematic separately : " << std::endl;
	output << std::endl;
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Channel";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Gl. mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Median";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " Mode";
	output << std::setw(10) << std::setiosflags(std::ios::left) << " RMS";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 16% quantile";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " 84% quantile";
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (low)"; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Uncert. (high)"; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Contribution"; 
	output << std::setw(15) << std::setiosflags(std::ios::left) << " Corr." << std::endl;; 
	output << std::endl;

	// loop over all channels
	for (int i = 0; i < GetNSystErrors(); ++i) {

		// get summary
		ParameterSummary* ps = fSummaryCombinationSingleSyst.at(i); 

		if (!ps) {
			BCLog::OutError("CombinationModel::PrintChannelSummary : Single channel summary information not available."); 
			return -1; 
		}

		// get summary information
		double median = ps->GetMedian(); 
		double glmode = ps->GetGlobalMode(); 
		double mode = ps->GetMode(); 
		double rms = ps->GetRMS(); 
		double q16 = ps->GetQuantile16();
		double q84 = ps->GetQuantile84();
		double errlow = median - q16;
		double errhigh = q84 - median; 
		double corr = ps->GetBuffer();

		output << " ";
		output << std::setw(15) << std::setiosflags(std::ios::left) << ps->GetName();
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << glmode; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << median; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << mode; 
		output << std::setw(10) << std::setiosflags(std::ios::left) << std::setprecision(4) << rms; 
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q16;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << q84;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errlow;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << errhigh;
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << sqrt(rms*rms-rms_nosyst*rms_nosyst);
		output << std::setw(15) << std::setiosflags(std::ios::left) << std::setprecision(4) << corr << std::endl;
	}
	output << std::endl << std::endl;

	output << " Contribution   = sqrt(rms*rms - rms(no syst)*rms(nosyst))" << std::endl;
	output << " Uncert. (low)  = median - 16% quantile" << std::endl;
	output << " Uncert. (high) = 84% quantile - median" << std::endl;

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

