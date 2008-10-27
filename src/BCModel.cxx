/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCModelTest.h"
#include "BCModel.h"
#include "BCLog.h"
#include "BCErrorCodes.h"
#include "BCMath.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include <fstream>
#include <iomanip>

// ---------------------------------------------------------

BCModel::BCModel(const char * name) : BCIntegrate()
{
	fNormalization = -1.0;
	fDataSet = 0;
	fParameterSet = new BCParameterSet;

	fIndex = -1;
	fPValue = -1;

	fName = (char *) name;
	flag_ConditionalProbabilityEntry = true;

	fDataPointUpperBoundaries = 0;
	fDataPointLowerBoundaries = 0;

	fErrorBandXY = 0;
}

// ---------------------------------------------------------

BCModel::BCModel() : BCIntegrate()
{
	fNormalization = -1.0;
	fDataSet = 0;
	fParameterSet = new BCParameterSet();

	fIndex = -1;
	fPValue = -1;

	fName = "model";
	fDataPointUpperBoundaries = 0;
	fDataPointLowerBoundaries = 0;

	flag_ConditionalProbabilityEntry = true;
}

// ---------------------------------------------------------

BCModel::~BCModel()
{
	delete fParameterSet;

	if (fDataPointLowerBoundaries)
		delete fDataPointLowerBoundaries;

	if (fDataPointUpperBoundaries)
		delete fDataPointUpperBoundaries;
}

// ---------------------------------------------------------

int BCModel::GetNDataPoints()
{
	int npoints = 0;
	if (fDataSet)
		npoints = fDataSet -> GetNDataPoints();
	else
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCModel::GetNDataPoints(). No data set defined.");
		return ERROR_NOEVENTS;
	}

	return npoints;
}

// ---------------------------------------------------------

BCDataPoint * BCModel::GetDataPoint(int index)
{
	if (fDataSet)
		return fDataSet -> GetDataPoint(index);

	BCLog::Out(BCLog::warning, BCLog::warning,"BCModel::GetDataPoint. No data set defined.");
	return 0;
}

// ---------------------------------------------------------

BCParameter * BCModel::GetParameter(int index)
{
	if (!fParameterSet)
		return 0;

	if (index < 0 || index >= this -> GetNParameters())
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				Form("BCModel::GetParameter. Parameter index %d not within range.", index));
		return 0;
	}

	return fParameterSet -> at(index);
}

// ---------------------------------------------------------

BCParameter * BCModel::GetParameter(const char * name)
{
	if (!fParameterSet)
		return 0;

	int index = -1;
	for (int i = 0; i < this->GetNParameters(); i++)
		if (name == this -> GetParameter(i) -> GetName())
			index = i;

	if (index<0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				Form("BCModel::GetParameter. Model %s has no parameter named %s.", (this -> GetName()).data(), name));
		return 0;
	}

	return this->GetParameter(index);
}

// ---------------------------------------------------------

std::vector <double> BCModel::GetErrorBand(double level)
{
	std::vector <double> errorband;

	if (!fErrorBandXY)
		return errorband;

	int nx = fErrorBandXY -> GetNbinsX();
	errorband.assign(nx - 1, 0.0);

	// loop over x and y bins
	for (int ix = 1; ix < nx; ix++)
	{
		TH1D * temphist = fErrorBandXY -> ProjectionY("temphist", ix, ix);

		int nprobSum = 1;
		double q[1];
		double probSum[1];
		probSum[0] = level;

		temphist -> GetQuantiles(nprobSum, q, probSum);

		errorband[ix-1] = q[0];
	}

	return errorband;
}

// ---------------------------------------------------------

TGraph * BCModel::GetErrorBandGraph(double level1, double level2)
{
	if (!fErrorBandXY)
		return 0;

	// define new graph
	int nx = fErrorBandXY -> GetNbinsX() - 1;

	TGraph * graph = new TGraph(2 * nx);
	graph -> SetFillStyle(1001);
	graph -> SetFillColor(kYellow);

	// get error bands
	std::vector <double> ymin = this -> GetErrorBand(level1);
	std::vector <double> ymax = this -> GetErrorBand(level2);

	for (int i = 0; i < nx; i++)
	{
		graph -> SetPoint(i, fErrorBandXY -> GetXaxis() -> GetBinCenter(i + 1), ymin.at(i));
		graph -> SetPoint(nx + i, fErrorBandXY -> GetXaxis() -> GetBinCenter(nx - i), ymax.at(nx - i - 1));
	}

	return graph;
}

// ---------------------------------------------------------

TGraph * BCModel::GetFitFunctionGraph(std::vector <double> parameters)
{
	if (!fErrorBandXY)
		return 0;

	// define new graph
	int nx = fErrorBandXY -> GetNbinsX();
	TGraph * graph = new TGraph(nx);

	// loop over x values
	for (int i = 0; i < nx; i++)
	{
		double x = fErrorBandXY -> GetXaxis() -> GetBinCenter(i + 1);

		std::vector <double> xvec;
		xvec.push_back(x);
		double y = this -> FitFunction(xvec, parameters);
		xvec.clear();

		graph -> SetPoint(i, x, y);
	}

	return graph;
}

// ---------------------------------------------------------

TGraph * BCModel::GetFitFunctionGraph(std::vector <double> parameters, double xmin, double xmax, int n)
{
	// define new graph
	TGraph * graph = new TGraph(n+1);

	double dx = (xmax-xmin)/(double)n;

	// loop over x values
	for (int i = 0; i <= n; i++)
	{
		double x = (double)i*dx+xmin;
		std::vector <double> xvec;
		xvec.push_back(x);
		double y = this -> FitFunction(xvec, parameters);

		xvec.clear();

		graph -> SetPoint(i, x, y);
	}

	return graph;
}

// ---------------------------------------------------------

bool BCModel::GetFlagBoundaries()
{
	if (!fDataPointLowerBoundaries)
		return false;

	if (!fDataPointUpperBoundaries)
		return false;

	if (int(fDataPointLowerBoundaries -> GetNValues()) != fDataSet -> GetDataPoint(0) -> GetNValues())
		return false;

	if (int(fDataPointUpperBoundaries -> GetNValues()) != fDataSet -> GetDataPoint(0) -> GetNValues())
		return false;

	return true;
}

// ---------------------------------------------------------

void BCModel::SetSingleDataPoint(BCDataPoint * datapoint)
{
	// create new data set consisting of a single data point
	BCDataSet * dataset = new BCDataSet();

	// add the data point
	dataset -> AddDataPoint(datapoint);

	// set this new data set
	this -> SetDataSet(dataset);

}

// ---------------------------------------------------------

void BCModel::SetSingleDataPoint(BCDataSet * dataset, unsigned int index)
{
	if (index < 0 || index > dataset -> GetNDataPoints())
		return;

	this -> SetSingleDataPoint(dataset -> GetDataPoint(index));
}

// ---------------------------------------------------------

void BCModel::SetDataBoundaries(int index, double lowerboundary, double upperboundary, bool fixed)
{
	// check if data set exists
	if (!fDataSet)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				 "BCModel::SetDataBoundaries. Need to define data set first.");
		return;
	}

	// check if index is within range
	if (index < 0 || index > fDataSet -> GetDataPoint(0) -> GetNValues())
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				"BCModel::SetDataBoundaries. Index out of range.");
		return;
	}

	// check if boundary data points exist
	if (!fDataPointLowerBoundaries)
		fDataPointLowerBoundaries = new BCDataPoint(fDataSet -> GetDataPoint(0) -> GetNValues());

	if (!fDataPointUpperBoundaries)
		fDataPointUpperBoundaries = new BCDataPoint(fDataSet -> GetDataPoint(0) -> GetNValues());

	if (fDataFixedValues.size() == 0)
		fDataFixedValues.assign(fDataSet -> GetDataPoint(0) -> GetNValues(), false);

	// set boundaries
	fDataPointLowerBoundaries -> SetValue(index, lowerboundary);
	fDataPointUpperBoundaries -> SetValue(index, upperboundary);
	fDataFixedValues[index] = fixed;
}

// ---------------------------------------------------------

void BCModel::SetErrorBandContinuous(bool flag)
{
	fErrorBandContinuous = flag;

	if (flag)
		return;

	// clear x-values
	fErrorBandX.clear();

	// copy data x-values
	for (unsigned int i = 0; i < fDataSet -> GetNDataPoints(); ++i)
		fErrorBandX.push_back(fDataSet -> GetDataPoint(i) -> GetValue(fFitFunctionIndexX));
}

// ---------------------------------------------------------

int BCModel::AddParameter(const char * name, double lowerlimit, double upperlimit)
{
	// create new parameter
	BCParameter * parameter = new BCParameter(name, lowerlimit, upperlimit);

	int flag_ok = this -> AddParameter(parameter);
	if (flag_ok != 0)
		delete parameter;

	return flag_ok;
}

// ---------------------------------------------------------

int BCModel::AddParameter(BCParameter * parameter)
{
	// check if parameter set exists
	if (!fParameterSet)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				"BCModel::AddParameter. Parameter set does not exist");
			return ERROR_PARAMETERSETDOESNOTEXIST;
	}

	// check if parameter with same name exists
	int flag_exists = 0;
	for (int i = 0; i < this -> GetNParameters(); i++)
		if (this -> CompareStrings(parameter -> GetName().data(), this -> GetParameter(i) -> GetName().data()) == 0)
			flag_exists = -1;

	if (flag_exists < 0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				Form("BCModel::AddParameter. Parameter with name %s exists already. ", parameter -> GetName().data()));
		return ERROR_PARAMETEREXISTSALREADY;
	}

	// define index of new parameter
	int index = int(fParameterSet -> size());
	parameter -> SetIndex(index);

	// add parameter to parameter container
	fParameterSet -> push_back(parameter);

	// add parameters to integation methods
	this -> SetParameters(fParameterSet);

	// add parameter to MCMC
	//	this -> MCMCAddParameter(parameter -> GetLowerLimit(), parameter -> GetUpperLimit());

	// re-initialize
	this -> MCMCInitialize();

	return 0;
}

// ---------------------------------------------------------

double BCModel::LogProbabilityNN(std::vector <double> parameters)
{
	// add log of conditional probability
	double logprob = this -> LogLikelihood(parameters);

	// add log of prior probability
	logprob += this -> LogAPrioriProbability(parameters);

	return logprob;
}

// ---------------------------------------------------------

double BCModel::LogProbability(std::vector <double> parameters)
{
	// check if normalized
	if (fNormalization<0. || fNormalization==0.)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				 "BCModel::LogProbability. Normalization not available or zero.");
		return 0.;
	}

	return this -> LogProbabilityNN(parameters) - log(fNormalization);
}

// ---------------------------------------------------------

double BCModel::LogLikelihood(std::vector <double> parameters)
{
	int ndatapoints = fDataSet -> GetNDataPoints();
	double logprob = 0.;

	// add log of conditional probabilities event-by-event
	for (int i=0;i<ndatapoints;i++)
	{
		BCDataPoint * datapoint = this -> GetDataPoint(i);
		logprob += this -> LogConditionalProbabilityEntry(datapoint, parameters);
	}

	return logprob;
}

// ---------------------------------------------------------

double BCModel::LogEval(std::vector <double> parameters)
{
	return this -> LogProbabilityNN(parameters);
}

// ---------------------------------------------------------

double BCModel::EvalSampling(std::vector <double> parameters)
{
	return this -> SamplingFunction(parameters);
}

// ---------------------------------------------------------

double BCModel::SamplingFunction(std::vector <double> parameters)
{
	double probability = 1.0;
	for (std::vector<BCParameter*>::const_iterator it = fParameterSet -> begin(); it != fParameterSet -> end(); ++it)
		probability *= 1.0 / ((*it) -> GetUpperLimit() - (*it) -> GetLowerLimit());
	return probability;
}

// ---------------------------------------------------------

double BCModel::Normalize()
{
	BCLog::Out(BCLog::summary, BCLog::summary, Form("Model \'%s\': Normalizing probability",this->GetName().data()));

	int n = this -> GetNvar();

	// initialize BCIntegrate if not done already
	if (n == 0)
	{
		this->SetParameters(fParameterSet);
		n = this->GetNvar();
	}

	// integrate and get best fit parameters
	// maybe we have to remove the mode finding from here in the future
	fNormalization = this -> Integrate();

	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Normalization : %.2lf", fNormalization));

	return fNormalization;
}

// ---------------------------------------------------------

int BCModel::CheckParameters(std::vector <double> parameters)
{
	// check if vectors are of equal size
	if (!parameters.size() == fParameterSet -> size())
		return  ERROR_INVALIDNUMBEROFPARAMETERS;

	// check if parameters are within limits
	for (int i = 0; i < int(fParameterSet -> size()); i++)
	{
		BCParameter * modelparameter = fParameterSet -> at(i);

		if (modelparameter -> GetLowerLimit() > parameters.at(i) ||
				modelparameter -> GetUpperLimit() < parameters.at(i))
		{
			BCLog::Out(BCLog::warning, BCLog::warning,
					 Form("BCModel::CheckParameters. Parameter %s not within limits.", fParameterSet -> at(i) -> GetName().data()));
			return ERROR_PARAMETERNOTWITHINRANGE;
				}
		}

	return 0;
}

// ---------------------------------------------------------

void BCModel::FindMode()
{
// this implementation is CLEARLY not good we have top work on this.

	BCLog::Out(BCLog::summary, BCLog::summary,
		Form("Model \'%s\': Finding mode", this -> GetName().data()));

	// synchronize parameters in BCIntegrate
	this -> SetParameters(fParameterSet);

	switch(this -> GetOptimizationMethod())
	{
		case BCIntegrate::kOptSimulatedAnnealing:
			this -> FindModeSA();
			return;
		case BCIntegrate::kOptMinuit:
			this -> FindModeMinuit();
			return;
		case BCIntegrate::kOptMetropolis:
			this -> MarginalizeAll();
			return;
		}

	BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::FindMode. Invalid mode finding method: %d. Return.",
			this->GetOptimizationMethod()));

	return;
}

// ---------------------------------------------------------

void BCModel::WriteMode(const char * file)
{
	ofstream ofi(file);
	if(!ofi.is_open())
	{
		std::cerr<<"Couldn't open file "<<file<<std::endl;
		return;
	}

	int npar = fParameterSet -> size();
	for (int i=0; i<npar; i++)
		ofi<<fBestFitParameters.at(i)<<std::endl;

	ofi<<std::endl;
	ofi<<"#######################################################################"<<std::endl;
	ofi<<"#"<<std::endl;
	ofi<<"#  This file was created automatically by BCModel::WriteMode() call."<<std::endl;
	ofi<<"#  It can be read in by call to BCModel::ReadMode()."<<std::endl;
	ofi<<"#  Do not modify it unless you know what you're doing."<<std::endl;
	ofi<<"#"<<std::endl;
	ofi<<"#######################################################################"<<std::endl;
	ofi<<"#"<<std::endl;
	ofi<<"#  Best fit parameters (mode) for model:"<<std::endl;
	ofi<<"#  \'"<<fName.data()<<"\'"<<std::endl;
	ofi<<"#"<<std::endl;
	ofi<<"#  Number of parameters: "<<npar<<std::endl;
	ofi<<"#  Parameters ordered as above:"<<std::endl;

	for (int i=0; i<npar; i++)
	{
		ofi<<"#     "<<i<<": ";
		ofi<<fParameterSet->at(i)->GetName().data()<<" = ";
		ofi<<fBestFitParameters.at(i)<<std::endl;
	}

	ofi<<"#"<<std::endl;
	ofi<<"########################################################################"<<std::endl;
}

// ---------------------------------------------------------

int BCModel::ReadMode(const char * file)
{
	ifstream ifi(file);
	if(!ifi.is_open())
	{
		std::cerr<<"Couldn't open file "<<file<<std::endl;
		return 0;
	}

	int npar = fParameterSet -> size();
	std::vector <double> mode;

	int i=0;
	while (i<npar && !ifi.eof())
	{
		double a;
		ifi>>a;
		mode.push_back(a);
		i++;
	}

	if(i<npar)
	{
		std::cerr<<"Couldn't read mode from file "<<file<<std::endl;
		std::cerr<<"Expected "<<npar<<" parameters, found "<<i<<std::endl;
		return 0;
	}

	BCLog::Out(BCLog::summary,BCLog::summary,
			Form("#  Read in best fit parameters (mode) for model \'%s\' from file %s:",fName.data(),file));
	this->SetMode(mode);
	for(int j=0 ; j<npar; j++)
		BCLog::Out(BCLog::summary,BCLog::summary,
				Form("#    -> Parameter %d : %s = %e", j, fParameterSet->at(j)->GetName().data(), fBestFitParameters[j]));

	BCLog::Out(BCLog::warning,BCLog::warning,
			"#  ! Best fit values obtained before this call will be overwritten !");

	return npar;
}

// ---------------------------------------------------------

int BCModel::MarginalizeAll()
{
	// prepare function fitting
	double dx = 0.0;
	double dy = 0.0;

	if (fFitFunctionIndexX >= 0)
	{
		dx = (fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexX) - fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexX))
				/ (double)fErrorBandNbinsX;

		dy = (fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexY) - fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexY))
				/ (double)fErrorBandNbinsY;

		fErrorBandXY = new TH2D("errorbandxy", "",
				fErrorBandNbinsX + 1,
				fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexX) - .5 * dx,
				fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexX) + .5 * dx,
				fErrorBandNbinsY + 1,
				fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexY) - .5 * dy,
				fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexY) + .5 * dy);
		fErrorBandXY -> SetStats(kFALSE);

		for (int ix = 1; ix <= fErrorBandNbinsX; ++ix)
			for (int iy = 1; iy <= fErrorBandNbinsX; ++iy)
				fErrorBandXY -> SetBinContent(ix, iy, 0.0);
	}

	this -> MCMCMetropolis();
	this -> FindModeMCMC();

	//	this -> PrintResults(Form("%s.txt", this -> GetName().data()));

	return 1;
}


// ---------------------------------------------------------

BCH1D * BCModel::GetMarginalized(BCParameter * parameter)
{
//	if(fMCMCH1Marginalized.size()==0)
//	{
//		BCLog::Out(BCLog::warning, BCLog::warning,
//				"BCModel::GetMarginalized. MarginalizeAll() has to be run prior to this.");
//		return 0;
//	}

	int index = parameter -> GetIndex();

	// get histogram
	TH1D * hist = this -> MCMCGetH1Marginalized(index);
	if(!hist)
		return 0;

	BCH1D * hprob = new BCH1D();

	// set axis labels
	hist -> SetName(Form("hist_%s_%s", this -> GetName().data(), parameter -> GetName().data()));
	hist -> SetXTitle(parameter -> GetName().data());
	hist -> SetYTitle(Form("p(%s|data)", parameter -> GetName().data()));
	hist -> SetStats(kFALSE);

	// set histogram
	hprob -> SetHistogram(hist);

	// set best fit parameter
	double bestfit = hprob -> GetMode();

	if (fBestFitParametersMarginalized.size() == 0)
		for (int i = 0; i < this -> GetNParameters(); i++)
			fBestFitParametersMarginalized.push_back(0.0);

	fBestFitParametersMarginalized[index] = bestfit;

	hprob->SetGlobalMode(fBestFitParameters.at(index));

	return hprob;
}

// ---------------------------------------------------------

int BCModel::ReadMarginalizedFromFile(const char * file)
{
	TFile * froot = new TFile(file);
	if(!froot->IsOpen())
	{
		BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::ReadMarginalizedFromFile. Couldn't open file %s.",file));
		return 0;
	}

	int k=0;
	int n=this->GetNParameters();
	for(int i=0;i<n;i++)
	{
		BCParameter * a = this->GetParameter(i);
		TH1D * h1 = (TH1D*)froot -> Get(Form("hist_%s_%s", this -> GetName().data(), a -> GetName().data()));
		if(h1)
		{
			h1->SetDirectory(0);
			if(this->SetMarginalized(i,h1))
				k++;
		}
		else
			BCLog::Out(BCLog::warning, BCLog::warning,
				Form("BCModel::ReadMarginalizedFromFile. Couldn't read histogram \"hist_%s_%s\" from file %s.",
					this -> GetName().data(), a -> GetName().data(), file));
	}

	for(int i=0;i<n-1;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			BCParameter * a = this->GetParameter(i);
			BCParameter * b = this->GetParameter(j);
			TH2D * h2 = (TH2D*)froot -> Get(Form("hist_%s_%s_%s", this -> GetName().data(), a -> GetName().data(), b -> GetName().data()));
			if(h2)
			{
				h2->SetDirectory(0);
				if(this->SetMarginalized(i,j,h2))
					k++;
			}
			else
				BCLog::Out(BCLog::warning, BCLog::warning,
					Form("BCModel::ReadMarginalizedFromFile. Couldn't read histogram \"hist_%s_%s_%s\" from file %s.",
						this -> GetName().data(), a -> GetName().data(), b -> GetName().data(), file));
		}
	}

	froot->Close();

	return k;
}

// ---------------------------------------------------------

int BCModel::ReadErrorBandFromFile(const char * file)
{
	TFile * froot = new TFile(file);
	if(!froot->IsOpen())
	{
		BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::ReadErrorBandFromFile. Couldn't open file %s.",file));
		return 0;
	}

	int r=0;

	TH2D * h2 = (TH2D*)froot->Get("errorbandxy");
	if(h2)
	{
		h2->SetDirectory(0);
		this->SetErrorBandHisto(h2);
		r=1;
	}
	else
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCModel::ReadErrorBandFromFile. Couldn't read histogram \"errorbandxy\" from file %s.", file));

	froot->Close();

	return r;
}

// ---------------------------------------------------------

int BCModel::PrintAllMarginalized1D(const char * filebase)
{
	if(fMCMCH1Marginalized.size()==0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				"BCModel::PrintAllMarginalized2D. MarginalizeAll() has to be run prior to this to fill the distributions.");
		return 0;
	}

	int n=this->GetNParameters();
	for(int i=0;i<n;i++)
	{
		BCParameter * a = this->GetParameter(i);
		this -> GetMarginalized(a) -> Print(Form("%s_1D_%s.ps",filebase,a->GetName().data()));
	}

	return n;
}

// ---------------------------------------------------------

int BCModel::PrintAllMarginalized2D(const char * filebase)
{
	if(fMCMCH2Marginalized.size()==0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				"BCModel::PrintAllMarginalized2D. MarginalizeAll() has to be run prior to this to fill the distributions.");
		return 0;
	}

	int k=0;
	int n=this->GetNParameters();
	for(int i=0;i<n-1;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			BCParameter * a = this->GetParameter(i);
			BCParameter * b = this->GetParameter(j);

			double meana = (a -> GetLowerLimit() + a -> GetUpperLimit()) / 2.0;
			double deltaa = (a -> GetUpperLimit() - a -> GetLowerLimit());
			if (deltaa <= 1e-7 * meana)
				continue;

			double meanb = (b -> GetLowerLimit() + b -> GetUpperLimit()) / 2.0;
			double deltab = (b -> GetUpperLimit() - b -> GetLowerLimit());
			if (deltab <= 1e-7 * meanb)
				continue;

			this -> GetMarginalized(a,b) -> Print(Form("%s_2D_%s_%s.ps",filebase,a->GetName().data(),b->GetName().data()));
			k++;
		}
	}

	return k;
}

// ---------------------------------------------------------

int BCModel::PrintAllMarginalized(const char * file, int hdiv, int vdiv)
{
	if (fMCMCH1Marginalized.size()==1 && fMCMCH2Marginalized.size()==0)
	{
		BCParameter * a = this->GetParameter(0);
		this -> GetMarginalized(a) -> Print(file);
		return 1;
	}

	if(fMCMCH1Marginalized.size()==0 || fMCMCH2Marginalized.size()==0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				"BCModel::PrintAllMarginalized. MarginalizeAll() has to be run prior to this to fill the distributions.");
		return 0;
	}

	int n=0;

	// setup the canvas and postscript file
	TCanvas * c = new TCanvas("c","canvas");

	int type = 112;
	TPostScript * ps = new TPostScript(file,type);
	ps->Range(24,18);

	// draw all 1D distributions
	ps->NewPage();
	c->cd();
	c->Clear();
	c->Divide(hdiv,vdiv);

	int n1d=this->GetNParameters();
	for(int i=0;i<n1d;i++)
	{
		// if current page is full, swith to new page
		if(i!=0 && i%(hdiv*vdiv)==0)
		{
			c->Update();
			ps->NewPage();
			c->cd();
			c->Clear();
			c->Divide(hdiv,vdiv);
		}

		// go to next pad
		c->cd(i%(hdiv*vdiv)+1);

		BCParameter * a = this->GetParameter(i);
		this -> GetMarginalized(a) -> Draw();
		n++;
	}

	c->Update();

	// draw all the 2D distributions
	ps->NewPage();
	c->cd();
	c->Clear();
	c->Divide(hdiv,vdiv);

	int k=0;
	int n2d=this->GetNParameters();
	for(int i=0;i<n2d-1;i++)
	{
		for(int j=i+1;j<n2d;j++)
		{
			// if current page is full, switch to new page
			if(k!=0 && k%(hdiv*vdiv)==0)
			{
				c->Update();
				ps->NewPage();
				c->cd();
				c->Clear();
				c->Divide(hdiv,vdiv);
			}

			// go to next pad
			c->cd(k%(hdiv*vdiv)+1);

			BCParameter * a = this->GetParameter(i);
			BCParameter * b = this->GetParameter(j);

			double meana = (a -> GetLowerLimit() + a -> GetUpperLimit()) / 2.0;
			double deltaa = (a -> GetUpperLimit() - a -> GetLowerLimit());
			if (deltaa <= 1e-7 * meana)
				continue;

			double meanb = (b -> GetLowerLimit() + b -> GetUpperLimit()) / 2.0;
			double deltab = (b -> GetUpperLimit() - b -> GetLowerLimit());
			if (deltab <= 1e-7 * meanb)
				continue;

			this -> GetMarginalized(a,b) -> Draw(52);
			k++;
		}
	}

	c->Update();
	ps->Close();

	delete c;
	delete ps;

	// return total number of drawn histograms
	return n+k;
}

// ---------------------------------------------------------

BCH2D * BCModel::GetMarginalized(BCParameter * parameter1, BCParameter * parameter2)
{
//	if(fMCMCH2Marginalized.size()==0)
//	{
//		BCLog::Out(BCLog::warning, BCLog::warning,
//				 "BCModel::GetMarginalized. MarginalizeAll() has to be run prior to this.");
//		return 0;
//	}

	int index1 = parameter1 -> GetIndex();
	int index2 = parameter2 -> GetIndex();

	if (index1 > index2)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				 "BCModel::GetMarginalized. Index 1 has to be smaller than index 2.");
		return 0;
	}

	// get histogram
	TH2D * hist = this -> MCMCGetH2Marginalized(index1, index2);

	if(hist==0)
		return 0;

	BCH2D * hprob = new BCH2D();

	// set axis labels
	hist -> SetName(Form("hist_%s_%s_%s", this -> GetName().data(), parameter1 -> GetName().data(), parameter2 -> GetName().data()));
	hist -> SetXTitle(Form("%s", parameter1 -> GetName().data()));
	hist -> SetYTitle(Form("%s", parameter2 -> GetName().data()));
	hist -> SetStats(kFALSE);

	double gmode[] = { fBestFitParameters.at(index1), fBestFitParameters.at(index2) };
	hprob->SetGlobalMode(gmode);

	// set histogram
	hprob -> SetHistogram(hist);

	return hprob;
}

// ---------------------------------------------------------

double BCModel::GetPvalueFromChi2(std::vector<double> par, int sigma_index)
{
	double ll = this -> LogLikelihood(par);
	int n = this -> GetNDataPoints();

	double sum_sigma=0;
	for (int i=0;i<n;i++)
		sum_sigma += log(this -> GetDataPoint(i) -> GetValue(sigma_index));

	double chi2 = -2.*(ll + (double)n/2. * log(2.*M_PI) + sum_sigma);

	return TMath::Prob(chi2,n);
}

// ---------------------------------------------------------

BCH1D * BCModel::CalculatePValue(std::vector<double> par, bool flag_histogram)
{
	BCH1D * hist = 0;

	// print log
	BCLog::Out(BCLog::summary, BCLog::summary, "Do goodness-of-fit-test");

	// create model test
	BCModelTest * modeltest = new BCModelTest("modeltest");

	// set this model as the model to be tested
	modeltest -> SetTestModel(this);

	// set the point in parameter space which is tested an initialize
	// the model testing
	if (!modeltest -> SetTestPoint(par))
		return 0;

	// get p-value
	fPValue = modeltest -> GetCalculatedPValue(flag_histogram);

	// get histogram
	if (flag_histogram)
	{
		hist = new BCH1D();
		hist -> SetHistogram(modeltest -> GetHistogramLogProb());
	}

	// delete model test
	delete modeltest;

	// return histogram
	return hist;
}

// ---------------------------------------------------------

void BCModel::CorrelateDataPointValues(std::vector<double> &x)
{
	// ...
}

// ---------------------------------------------------------

double BCModel::HessianMatrixElement(BCParameter * parameter1, BCParameter * parameter2, std::vector<double> point)
{
	// check number of entries in vector
	if (int(point.size()) != this -> GetNParameters())
	{
		BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::HessianMatrixElement. Invalid number of entries in the vector."));
		return -1;
	}

	// define steps
	double nsteps = 1e5;
	double dx1 = (parameter1 -> GetUpperLimit() - parameter1 -> GetLowerLimit()) / nsteps;
	double dx2 = (parameter2 -> GetUpperLimit() - parameter2 -> GetLowerLimit()) / nsteps ;

	// define points at which to evaluate
	std::vector<double> xpp = point;
	std::vector<double> xpm = point;
	std::vector<double> xmp = point;
	std::vector<double> xmm = point;

	xpp.at(parameter1 -> GetIndex()) = xpp.at(parameter1 -> GetIndex()) + dx1;
	xpp.at(parameter2 -> GetIndex()) = xpp.at(parameter2 -> GetIndex()) + dx2;

	xpm.at(parameter1 -> GetIndex()) = xpm.at(parameter1 -> GetIndex()) + dx1;
	xpm.at(parameter2 -> GetIndex()) = xpm.at(parameter2 -> GetIndex()) - dx2;

	xmp.at(parameter1 -> GetIndex()) = xmp.at(parameter1 -> GetIndex()) - dx1;
	xmp.at(parameter2 -> GetIndex()) = xmp.at(parameter2 -> GetIndex()) + dx2;

	xmm.at(parameter1 -> GetIndex()) = xmm.at(parameter1 -> GetIndex()) - dx1;
	xmm.at(parameter2 -> GetIndex()) = xmm.at(parameter2 -> GetIndex()) - dx2;

	// calculate probability at these points
	double ppp = this -> Likelihood(xpp);
	double ppm = this -> Likelihood(xpm);
	double pmp = this -> Likelihood(xmp);
	double pmm = this -> Likelihood(xmm);

	// calculate derivative
	double derivative = (ppp + pmm - ppm - pmp) / (4.0 * dx1 * dx2);

	return derivative;
}

// ---------------------------------------------------------

void BCModel::FixDataAxis(int index, bool fixed)
{
	// check if index is within range
	if (index < 0 || index > fDataSet -> GetDataPoint(0) -> GetNValues())
	{
		BCLog::Out(BCLog::warning, BCLog::warning, "BCModel::FixDataAxis. Index out of range.");
		return;
	}

	if (fDataFixedValues.size() == 0)
		fDataFixedValues.assign(fDataSet -> GetDataPoint(0) -> GetNValues(), false);

	fDataFixedValues[index] = fixed;
}

// ---------------------------------------------------------

bool BCModel::GetFixedDataAxis(int index)
{
	// check if index is within range
	if (index < 0 || index > fDataSet -> GetDataPoint(0) -> GetNValues())
	{
		BCLog::Out(BCLog::warning, BCLog::warning, "BCModel::GetFixedDataAxis. Index out of range.");
		return false;
	}

	return fDataFixedValues.at(index);
}

// ---------------------------------------------------------

void BCModel::PrintSummary()
{
	int nparameters = this -> GetNParameters();

	// model summary
	std::cout
		<< std::endl
		<< "   ---------------------------------" << std::endl
		<< "    Model : " << fName.data() << std::endl
		<< "   ---------------------------------"<< std::endl
		<< "     Index                : " << fIndex << std::endl
		<< "     Number of parameters : " << nparameters << std::endl
		<< std::endl
		<< "     - Parameters : " << std::endl
		<< std::endl;

	// parameter summary
	for (int i=0; i<nparameters; i++)
		fParameterSet -> at(i) -> PrintSummary();

	// best fit parameters
	if (this -> GetBestFitParameters().size() > 0)
	{
		std::cout
			<< std::endl
			<< "     - Best fit parameters :" << std::endl
			<< std::endl;

		for (int i=0; i<nparameters; i++)
		{
			std::cout
				<< "       " << fParameterSet -> at(i) -> GetName().data()
				<< " = " << this -> GetBestFitParameter(i)
				<< " (overall)" << std::endl;
			if ((int)fBestFitParametersMarginalized.size() == nparameters)
				std::cout
					<< "       " << fParameterSet -> at(i) -> GetName().data()
					<< " = " << this -> GetBestFitParameterMarginalized(i)
					<< " (marginalized)" << std::endl;
		}
	}

	std::cout << std::endl;

	// model testing
	if (fPValue >= 0)
	{
		double likelihood = this -> Likelihood(this -> GetBestFitParameters());
		std::cout
			<< "   - Model testing:" << std::endl
			<< std::endl
			<< "       p(data|lambda*) = " << likelihood << std::endl
			<< "       p-value         = " << fPValue << std::endl
			<< std::endl;
	}

	// normalization
	if (fNormalization > 0)
		std::cout << "     Normalization        : " << fNormalization << std::endl;
}

// ---------------------------------------------------------

void BCModel::PrintResults(const char * file)
{
	// print summary of Markov Chain Monte Carlo

	// open file
	ofstream ofi(file);

	// check if file is open
	if(!ofi.is_open())
	{
		std::cerr << "Couldn't open file " << file <<std::endl;
		return;
	}

	// number of parameters
	int npar = fParameterSet -> size();

	// check convergence
	bool flag_conv = ((this -> MCMCGetNIterationsConvergenceGlobal() > 0)?1:0);

	ofi
		<< std::endl
		<< " -----------------------------------------------------" << std::endl
		<< " Summary of the Markov Chain Monte Carlo run" << std::endl
		<< " -----------------------------------------------------" << std::endl
		<< std::endl;

	if (!flag_conv)
	{
		ofi
			<< " WARNING: the Markov Chain did not converge! Be" << std::endl
			<< " cautious using the following results!" << std::endl
			<< std::endl;
	}

	ofi
		<< " Model summary" << std::endl
		<< " =============" << std::endl
		<< " Model: " << fName.data() << std::endl
		<< " Number of parameters: " << npar << std::endl
		<< " List of Parameters and ranges:" << std::endl;
	for (int i = 0; i < npar; ++i)
	{
		ofi
			<< "  (" << i << ") Parameter \""
			<< fParameterSet -> at(i) -> GetName().data() << "\""
			<< ": " << fParameterSet -> at(i) -> GetLowerLimit()
			<< " - "
			<< fParameterSet -> at(i) -> GetUpperLimit() << std::endl;
	}
	ofi << std::endl;

	if (flag_conv)
	{
		ofi
			<< " Results of the marginalization" << std::endl
			<< " ==============================" << std::endl
			<< " List of parameters and properties of the marginalized" << std::endl
			<< " distributions:" << std::endl;
		for (int i = 0; i < npar; ++i)
		{
			BCH1D * bch1d = this -> GetMarginalized(fParameterSet -> at(i));

			ofi
				<< "  (" << i << ") Parameter \""
					<< fParameterSet -> at(i) -> GetName().data() << "\"" << std::endl
				<< "      Mean +- RMS:         "
					<< std::setprecision(4) << bch1d -> GetMean()
					<< " +- "
					<< std::setprecision(4) << bch1d -> GetRMS() << std::endl
				<< "      Median +- sigma:     "
					<< std::setprecision(4) << bch1d -> GetMedian()
					<< " +  " << std::setprecision(4) << bch1d -> GetQuantile(0.84) - bch1d -> GetMedian()
					<< " - " << std::setprecision(4) << bch1d -> GetMedian() - bch1d -> GetQuantile(0.16) << std::endl
				<< "      (Marginalized) mode: " << bch1d -> GetMode() << std::endl
				<< "      Smallest interval(s) containing 68% and local modes:" << std::endl;

			std::vector <double> v;
			v = bch1d -> GetSmallestIntervals(0.68);
			int ninter = int(v.size());

			for (int j = 0; j < ninter; j+=5)
				ofi << "       " << v.at(j) << " - " << v.at(j+1) << " (local mode at " << v.at(j+3) << " with rel. height " << v.at(j+2) << "; rel. area " << v.at(j+4) << ")" << std::endl;
		}
		ofi << std::endl;
	}

	ofi
		<< " Results of the optimization" << std::endl
		<< " ===========================" << std::endl
		<< " Optimization algorithm used: ";
	switch(this -> GetOptimizationMethod())
	{
		case BCIntegrate::kOptSimulatedAnnealing:
			ofi << " Simulated Annealing" << std::endl;
			break;
		case BCIntegrate::kOptMinuit:
			ofi << " Minuit" << std::endl;
			break;
		case BCIntegrate::kOptMetropolis:
			ofi << " MCMC " << std::endl;
			break;
	}

	ofi << " List of parameters and global mode:" << std::endl;
	for (int i = 0; i < npar; ++i)
		ofi
			<< "  (" << i << ") Parameter \""
			<< fParameterSet -> at(i) -> GetName().data() << "\": "
			<< fBestFitParameters.at(i) << std::endl;
	ofi << std::endl;

	if (fPValue >= 0.)
	{
		ofi
			<< " Results of the model test" << std::endl
			<< " =========================" << std::endl
			<< " p-value at global mode: " << fPValue << std::endl
			<< std::endl;
	}

	ofi
		<< " Status of the MCMC" << std::endl
		<< " ==================" << std::endl
		<< " Convergence reached: " << ((flag_conv)?"yes":"no") << std::endl;

	if (flag_conv)
		ofi << " Number of iterations until convergence: " << this -> MCMCGetNIterationsConvergenceGlobal() << std::endl;
	else
		ofi
			<< " WARNING: the Markov Chain did not converge! Be\n"
			<< " cautious using the following results!" << std::endl
			<< std::endl;
	ofi
		<< " Number of chains:                       " << this -> MCMCGetNChains() << std::endl
		<< " Number of iterations of each chain:     " << this -> MCMCGetNIterationsMax() << std::endl
		<< std::endl;

	ofi
		<< " -----------------------------------------------------" << std::endl
		<< std::endl
		<< " Notes" << std::endl
		<< " =====" << std::endl
		<< " (i) Median +- sigma denotes the median, m, of the" << std::endl
		<< "     marginalized distribution and the intervals from" << std::endl
		<< "     m to the 16% and 84% quantiles." << std::endl
		<< " -----------------------------------------------------" << std::endl;

	// close file
//	ofi.close;
}

// ---------------------------------------------------------

void BCModel::PrintHessianMatrix(std::vector<double> parameters)
{
	// check number of entries in vector
	if (int(parameters.size()) != this -> GetNParameters())
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				Form("BCModel::PrintHessianMatrix. Invalid number of entries in the vector."));
		return;
	}

	// print to screen
	std::cout
		<< std::endl
		<< " Hessian matrix elements: " << std::endl
		<< " Point: ";

	for (int i = 0; i < int(parameters.size()); i++)
		std::cout << parameters.at(i) << " ";
	std::cout << std::endl;

	// loop over all parameter pairs
	for (int i = 0; i < this -> GetNParameters(); i++)
		for (int j = 0; j < i; j++)
		{
			// calculate Hessian matrix element
			double hessianmatrixelement = this -> HessianMatrixElement(fParameterSet -> at(i),
					fParameterSet -> at(j), parameters);

			// print to screen
			std::cout << " " << i << " " << j << " : " << hessianmatrixelement << std::endl;
		}
}

// ---------------------------------------------------------

BCDataPoint * BCModel::VectorToDataPoint(std::vector<double> data)
{
	int sizeofvector = int(data.size());
	BCDataPoint * datapoint = new BCDataPoint(sizeofvector);
	datapoint -> SetValues(data);
	return datapoint;
}

// ---------------------------------------------------------

int BCModel::CompareStrings(const char * string1, const char * string2)
{
	int flag_same = 0;

	if (strlen(string1) != strlen(string2))
		return -1;

	for (int i = 0; i < int(strlen(string1)); i++)
		if (string1[i] != string2[i])
			flag_same = -1;

	return flag_same;
}

// ---------------------------------------------------------

