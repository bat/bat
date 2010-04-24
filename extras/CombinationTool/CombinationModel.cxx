// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project CombinationTool
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "CombinationModel.h"

#include <BAT/BCLog.h>

#include <TF1.h>
#include <TMath.h> 

double SplitGaussian(double* x, double* par);

// ---------------------------------------------------------
CombinationModel::CombinationModel(const char * name, double xmin, double xmax) : BCModel()
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

