#include "BCIntegrate.h"
#include "BCLog.h" 

#include "cuba.h" 

// --------------------------------------------------------- 

class BCIntegrate * global_this; 

// *********************************************
BCIntegrate::BCIntegrate()
{
	fNvar=0;
	fNiterPerDimension = 100;
	fNSamplesPer2DBin = 100;
	fRandom = new TRandom3(0);

	fNIterationsMax    = 1000000; 
	fRelativePrecision = 1e-3; 

	fIntegrateMethod   = BCIntegrate::kIMonteCarlo;
	fMarginalizeMethod = BCIntegrate::kMMetropolis;

	fNbins=100;
}

// *********************************************
BCIntegrate::BCIntegrate(BCParameterSet * par)
{
	fNvar=0;
	fNiterPerDimension = 100;
	fNSamplesPer2DBin = 100;
	fRandom = new TRandom3(0);

	fNbins=100;

	this->SetParameters(par);
}

// *********************************************
BCIntegrate::~BCIntegrate()
{
	DeleteVarList();

	fx=0;

	delete fRandom;
	fRandom=0;

	int n1 = fHProb1D.size();
	if(n1>0)
	{
		for (int i=0;i<n1;i++)
			delete fHProb1D.at(i);
	}

	int n2 = fHProb2D.size();
	if(n2>0)
	{
		for (int i=0;i<n2;i++)
			delete fHProb2D.at(i);
	}
}

// *********************************************
void BCIntegrate::SetParameters(BCParameterSet * par)
{
	DeleteVarList();

	fx = par;
	fNvar = fx->size();
	fMin = new double[fNvar];
	fMax = new double[fNvar];
	fVarlist = new int[fNvar];

	this->ResetVarlist(1);

	for(int i=0;i<fNvar;i++)
	{
		fMin[i]=(fx->at(i))->GetLowerLimit();
	 	fMax[i]=(fx->at(i))->GetUpperLimit();
	}

	fXmetro0.assign(fNvar,0.);
	fXmetro1.assign(fNvar,0.);
}

// *********************************************
void BCIntegrate::DeleteVarList()
{
	if(fNvar)
	{
		delete[] fVarlist;
		fVarlist=0;

		delete[] fMin;
		fMin=0;

		delete[] fMax;
		fMax=0;

		fx=0;
		fNvar=0;
	}
}

// *********************************************
void BCIntegrate::SetNbins(int n)
{
	if(n<2)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCIntegrate::SetNbins. Number of bins less than 2 makes no sense. Setting to 2.");
		n=2;
	}
	fNbins=n;
}

// *********************************************
void BCIntegrate::SetVarList(int * varlist)
{
	for(int i=0;i<fNvar;i++)
		fVarlist[i]=varlist[i];
}

// *********************************************
void BCIntegrate::ResetVarlist(int v)
{
	for(int i=0;i<fNvar;i++)
		fVarlist[i]=v;
}

// *********************************************
double BCIntegrate::Eval(std::vector <double> x)
{
	BCLog::Out(BCLog::warning, BCLog::warning, "BCIntegrate::Eval. No function. Function needs to be overloaded."); 
	return 0;
}

// *********************************************
double BCIntegrate::LogEval(std::vector <double> x)
{
	// this method should better also be overloaded
	return TMath::Log(this->Eval(x));
}

// *********************************************
double BCIntegrate::EvalSampling(std::vector <double> x)
{
	BCLog::Out(BCLog::warning, BCLog::warning, "BCIntegrate::EvalSampling. No function. Function needs to be overloaded."); 
	return 0;
}

// *********************************************
double BCIntegrate::LogEvalSampling(std::vector <double> x)
{
	return TMath::Log(this->EvalSampling(x));
}

// *********************************************
double BCIntegrate::EvalPrint(std::vector <double> x)
{
	double val=this->Eval(x);

	BCLog::Out(BCLog::detail, BCLog::detail, Form("BCIntegrate::EvalPrint. Value: %d.", val)); 

	return val;
}

// *********************************************
double BCIntegrate::Integrate()
{
	std::vector <double> parameter;
	parameter.assign(fNvar, 0.0);

	switch(fIntegrateMethod)
	{
		case BCIntegrate::kIMonteCarlo:
			return IntegralMC(parameter);

		case BCIntegrate::kIMetropolis: 
			return this -> IntegralMetro(parameter); 

		case BCIntegrate::kIImportance: 
			return this -> IntegralImportance(parameter); 

		case BCIntegrate::kICuba:
			return this -> CubaIntegrate();
	}

	BCLog::Out(BCLog::warning, BCLog::warning, Form("BCIntegrate::Integrate. Invalid integration method: %d. Return 0.", 
							fIntegrateMethod)); 

	return 0;
}

// *********************************************
double BCIntegrate::IntegralMC(std::vector <double> x, int * varlist)
{
	this->SetVarList(varlist);
	return IntegralMC(x);
}
  
// *********************************************
double BCIntegrate::IntegralMC(std::vector <double> x)
{
	// count the variables to integrate over
	int NvarNow=0;

	for(int i = 0; i < fNvar; i++)
		if(fVarlist[i])
			NvarNow++;

	// print to log
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Running MC integation over %i dimensions.", NvarNow));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Maximum number of iterations: %i", this->GetNIterationsMax()));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Aimed relative precision:     %e", this->GetRelativePrecision()));

	// calculate D (the integrated volume)
	double D = 1.0;
	for(int i = 0; i < fNvar; i++)
		if (fVarlist[i])
			D = D * (fMax[i] - fMin[i]);

	// reset variables
	double pmax = 0.0;
	double sumW  = 0.0;
	double sumW2 = 0.0;
	double integral = 0.0;
	double variance = 0.0;
	double relprecision = 1.0;

	std::vector <double> randx;
	randx.assign(fNvar, 0.0);

	// reset number of iterations
	fNIterations = 0;

	// iterate while precision is not reached and the number of iterations is lower than maximum number of iterations
	while ((fRelativePrecision < relprecision && fNIterationsMax > fNIterations) || fNIterations < 10)
	{
		// increase number of iterations
		fNIterations++;

		// get random numbers
		this -> GetRndmVector(randx);

		// scale random numbers
		for(int i = 0; i < fNvar; i++)
		{
			if(fVarlist[i])
				randx[i]=fMin[i]+randx[i]*(fMax[i]-fMin[i]);
			else
				randx[i]=x[i];
		}

		// evaluate function at sampled point
		double value = this->Eval(randx);

		// add value to sum and sum of squares
		sumW  += value;
		sumW2 += value * value;

		// search for maximum probability
		if (value > pmax)
		{
			// set new maximum value
			pmax = value;

			// delete old best fit parameter values
			fBestFitParameters.clear();

			// write best fit parameters
			for(int i = 0; i < fNvar; i++)
				fBestFitParameters.push_back(randx.at(i));
		}

		// calculate integral and variance
		integral = D * sumW / fNIterations;

		if (fNIterations%10000 == 0)
		{
			variance = (1.0 / double(fNIterations)) * (D * D * sumW2 / double(fNIterations) - integral * integral);
			double error = TMath::Sqrt(variance);
			relprecision = error / integral;

			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::IntegralMC. Iteration %i, integral: %e +- %e.", fNIterations, integral, error));
		}
	}

	fError = variance / fNIterations;

	// print to log
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Result of integration:        %e +- %e   in %i iterations.", integral, TMath::Sqrt(variance), fNIterations));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Obtained relative precision:  %e. ", TMath::Sqrt(variance) / integral));

	return integral;
}


// *********************************************
double BCIntegrate::IntegralMetro(std::vector <double> x)
{
	// print debug information
	BCLog::Out(BCLog::debug, BCLog::debug, Form("BCIntegrate::IntegralMetro. Integate over %i dimensions.", fNvar));

	// get total number of iterations
	double Niter = pow(fNiterPerDimension, fNvar);

	// print if more than 100,000 iterations
	if(Niter>1e5)
		BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Total number of iterations in Metropolis: %d.", Niter));

	// reset sum
	double sumI = 0;

	// prepare Metropolis algorithm
	std::vector <double> randx;
	randx.assign(fNvar,0.);
	InitMetro();

	// prepare maximum value
	double pmax = 0.0;

	// loop over iterations
	for(int i = 0; i <= Niter; i++)
	{
		// get random point from sampling distribution
		this -> GetRandomPointSamplingMetro(randx);

		// calculate probability at random point
		double val_f = this -> Eval(randx);

		// calculate sampling distributions at that point
		double val_g = this -> EvalSampling(randx);

		// add ratio to sum
		if (val_g > 0)
			sumI += val_f / val_g;

		// search for maximum probability
		if (val_f > pmax)
		{
			// set new maximum value
			pmax = val_f;

			// delete old best fit parameter values
			fBestFitParameters.clear();

			// write best fit parameters
			for(int i = 0; i < fNvar; i++)
				fBestFitParameters.push_back(randx.at(i));
		}

		// write intermediate results
		if((i+1)%100000 == 0)
			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::IntegralMetro. Iteration %d, integral: %d.", i+1, sumI/(i+1)));
	}

	// calculate integral
	double result=sumI/Niter;

	// print debug information
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Integral is %f in %i iterations. ", result, Niter));

	return result;
}

// *********************************************
double BCIntegrate::IntegralImportance(std::vector <double> x)
{
	// print debug information
	BCLog::Out(BCLog::debug, BCLog::debug, Form("BCIntegrate::IntegralImportance. Integate over %i dimensions.", fNvar));

	// get total number of iterations
	double Niter = pow(fNiterPerDimension, fNvar);

	// print if more than 100,000 iterations
	if(Niter>1e5)
		BCLog::Out(BCLog::detail, BCLog::detail, Form("BCIntegrate::IntegralImportance. Total number of iterations: %d.", Niter));

	// reset sum
	double sumI = 0;

	std::vector <double> randx;
	randx.assign(fNvar,0.);

	// prepare maximum value
	double pmax = 0.0;

	// loop over iterations
	for(int i = 0; i <= Niter; i++)
	{
		// get random point from sampling distribution
		this -> GetRandomPointImportance(randx);

		// calculate probability at random point
		double val_f = this -> Eval(randx);

		// calculate sampling distributions at that point
		double val_g = this -> EvalSampling(randx);

		// add ratio to sum
		if (val_g > 0)
			sumI += val_f / val_g;

		// search for maximum probability
		if (val_f > pmax)
		{
			// set new maximum value
			pmax = val_f;

			// delete old best fit parameter values
			fBestFitParameters.clear();

			// write best fit parameters
			for(int i = 0; i < fNvar; i++)
				fBestFitParameters.push_back(randx.at(i));
		}

		// write intermediate results
		if((i+1)%100000 == 0)
			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::IntegralImportance. Iteration %d, integral: %d.", i+1, sumI/(i+1)));
	}

	// calculate integral
	double result=sumI/Niter;

	// print debug information
	BCLog::Out(BCLog::debug, BCLog::debug, Form("BCIntegrate::IntegralImportance. Integral %f in %i iterations. ", result, Niter));

	return result;  
}

// *********************************************
TH1D* BCIntegrate::Marginalize(BCParameter * parameter)
{
	BCLog::Out(BCLog::detail, BCLog::detail,
		Form(" --> Marginalizing model wrt. parameter %s using method %d.", parameter->GetName(), fMarginalizeMethod));

	switch(fMarginalizeMethod)
	{
		case BCIntegrate::kIMonteCarlo:
			return MarginalizeByIntegrate(parameter);

		case BCIntegrate::kMMetropolis:
			return MarginalizeByMetro(parameter);
	}

	BCLog::Out(BCLog::warning, BCLog::warning,
		Form("BCIntegrate::Marginalize. Invalid marginalization method: %d. Return 0.", fMarginalizeMethod));

	return 0;
}

// *********************************************
TH2D * BCIntegrate::Marginalize(BCParameter * parameter1, BCParameter * parameter2)
{
	switch(fMarginalizeMethod)
	{
		case BCIntegrate::kIMonteCarlo:
			return MarginalizeByIntegrate(parameter1,parameter2);

		case BCIntegrate::kMMetropolis:
			return MarginalizeByMetro(parameter1,parameter2);
	}

	BCLog::Out(BCLog::warning, BCLog::warning,
		Form("BCIntegrate::Marginalize. Invalid marginalization method: %d. Return 0.", fMarginalizeMethod));

	return 0;
}

// *********************************************
TH1D* BCIntegrate::MarginalizeByIntegrate(BCParameter * parameter)
{
	// set parameter to marginalize
	this->ResetVarlist(1);
	int index = parameter->GetIndex();
	this->UnsetVar(index);

	// define histogram
	double hmin = parameter -> GetLowerLimit();
	double hmax = parameter -> GetUpperLimit();
	double hdx  = (hmax - hmin) / double(fNbins);
	TH1D * hist = new TH1D("hist","", fNbins, hmin, hmax);

	// integrate
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);

	for(int i=0;i<=fNbins;i++)
	{
		randx[index] = hmin + (double)i * hdx;

		double val = IntegralMC(randx);
		hist->Fill(randx[index], val);
	}

	// normalize
	hist -> Scale( 1./hist->Integral("width") );

	return hist;
}

// *********************************************
TH2D * BCIntegrate::MarginalizeByIntegrate(BCParameter * parameter1, BCParameter * parameter2)
{
	// set parameter to marginalize
	this->ResetVarlist(1);
	int index1 = parameter1->GetIndex();
	this->UnsetVar(index1);
	int index2 = parameter2->GetIndex();
	this->UnsetVar(index2);

	// define histogram
	double hmin1 = parameter1 -> GetLowerLimit(); 
	double hmax1 = parameter1 -> GetUpperLimit(); 
	double hdx1  = (hmax1 - hmin1) / double(fNbins); 

	double hmin2 = parameter2 -> GetLowerLimit(); 
	double hmax2 = parameter2 -> GetUpperLimit(); 
	double hdx2  = (hmax2 - hmin2) / double(fNbins); 

	TH2D * hist = new TH2D(Form("hist_%s_%s", parameter1 -> GetName(), parameter2 -> GetName()),"",
			fNbins, hmin1, hmax1,
			fNbins, hmin2, hmax2); 

	// integrate
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);

	for(int i=0;i<=fNbins;i++)
	{
		randx[index1] = hmin1 + (double)i * hdx1;
		for(int j=0;j<=fNbins;j++)
		{
			randx[index2] = hmin2 + (double)j * hdx2;

			double val = IntegralMC(randx);
			hist->Fill(randx[index1],randx[index2], val);
		}
	}

	// normalize
	hist -> Scale(1.0/hist->Integral("width"));

	return hist; 
}

// *********************************************
TH1D * BCIntegrate::MarginalizeByMetro(BCParameter * parameter)
{
	int niter=fNbins*fNbins*fNSamplesPer2DBin*fNvar;

	BCLog::Out(BCLog::detail, BCLog::detail,
		Form(" --> Number of samples in Metropolis marginalization: %d.", niter));

	// set parameter to marginalize
	int index = parameter->GetIndex();

	// define histogram
	double hmin = parameter -> GetLowerLimit();
	double hmax = parameter -> GetUpperLimit();
	TH1D * hist = new TH1D("hist","", fNbins, hmin, hmax);

	// prepare Metro
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);
	InitMetro();

	for(int i=0;i<=niter;i++)
	{
		GetRandomPointMetro(randx);
		hist->Fill(randx[index]);
	}

	// normalize
	hist -> Scale(1.0/hist->Integral("width"));

	return hist;
}

// *********************************************
TH2D * BCIntegrate::MarginalizeByMetro(BCParameter * parameter1, BCParameter * parameter2)
{
	int niter=fNbins*fNbins*fNSamplesPer2DBin*fNvar;

	// set parameter to marginalize
	int index1 = parameter1->GetIndex();
	int index2 = parameter2->GetIndex();

	// define histogram
	double hmin1 = parameter1 -> GetLowerLimit();
	double hmax1 = parameter1 -> GetUpperLimit();

	double hmin2 = parameter2 -> GetLowerLimit();
	double hmax2 = parameter2 -> GetUpperLimit();

	TH2D * hist = new TH2D(Form("hist_%s_%s", parameter1 -> GetName(), parameter2 -> GetName()),"",
			fNbins, hmin1, hmax1,
			fNbins, hmin2, hmax2);

	// prepare Metro
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);
	InitMetro();

	for(int i=0;i<=niter;i++)
	{
		GetRandomPointMetro(randx);
		hist->Fill(randx[index1],randx[index2]);
	}

	// normalize
	hist -> Scale(1.0/hist->Integral("width"));

	return hist; 
}

// *********************************************
int BCIntegrate::MarginalizeAllByMetro(const char * name="")
{
	int niter=fNbins*fNbins*fNSamplesPer2DBin*fNvar;

	BCLog::Out(BCLog::detail, BCLog::detail,
		Form(" --> Number of samples in Metropolis marginalization: %d.", niter));

	// define 1D histograms
	for(int i=0;i<fNvar;i++)
	{
		double hmin1 = fx->at(i) -> GetLowerLimit();
		double hmax1 = fx->at(i) -> GetUpperLimit();

		TH1D * h1 = new TH1D(Form("h%s_%s", name, fx->at(i) -> GetName()),"",
			fNbins, hmin1, hmax1);

		fHProb1D.push_back(h1);
	}

	// define 2D histograms
	for(int i=0;i<fNvar-1;i++)
		for(int j=i+1;j<fNvar;j++)
		{
			double hmin1 = fx->at(i) -> GetLowerLimit();
			double hmax1 = fx->at(i) -> GetUpperLimit();

			double hmin2 = fx->at(j) -> GetLowerLimit();
			double hmax2 = fx->at(j) -> GetUpperLimit();

			TH2D * h2 = new TH2D(Form("h%s_%s_%s", name, fx->at(i) -> GetName(), fx->at(j) -> GetName()),"",
				fNbins, hmin1, hmax1,
				fNbins, hmin2, hmax2);

			fHProb2D.push_back(h2);
		}

	// get number of 2d distributions
	int nh2d = fHProb2D.size();

	BCLog::Out(BCLog::detail, BCLog::detail,
		Form(" --> Marginalizing %d 1D distributions and %d 2D distributions.", fNvar, nh2d));

	// prepare Metro
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);
	InitMetro();

	// run Metro and fill histograms
	for(int i=0;i<=niter;i++)
	{
		GetRandomPointMetro(randx);

		for(int j=0;j<fNvar;j++)
			fHProb1D[j] -> Fill( randx[j] );

		int ih=0;
		for(int j=0;j<fNvar-1;j++)
			for(int k=j+1;k<fNvar;k++)
			{
				fHProb2D[ih] -> Fill(randx[j],randx[k]);
				ih++;
			}

		if((i+1)%100000==0)
			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::MarginalizeAllByMetro. %d samples done.",i+1));
	}

	// normalize histograms
	for(int i=0;i<fNvar;i++)
		fHProb1D[i] -> Scale( 1./fHProb1D[i]->Integral("width") );
	for (int i=0;i<nh2d;i++)
		fHProb2D[i] -> Scale( 1./fHProb2D[i]->Integral("width") );

	return fNvar+nh2d;
}

// *********************************************
TH1D * BCIntegrate::GetH1D(int parIndex)
{
	if(fHProb1D.size()==0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCModel::GetH1D. MarginalizeAll() has to be run prior to this to fill the distributions.");
		return 0;
	}

	if(parIndex<0 || parIndex>fNvar)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH1D. Parameter index %d is invalid.",parIndex));
		return 0;
	}

	return fHProb1D[parIndex];
}

// *********************************************
int BCIntegrate::GetH2DIndex(int parIndex1, int parIndex2)
{
	if(parIndex1>fNvar-1 || parIndex1<0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH2DIndex. Parameter index %d is invalid", parIndex1));
		return -1;
	}

	if(parIndex2>fNvar-1 || parIndex2<0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH2DIndex. Parameter index %d is invalid", parIndex2));
		return -1;
	}

	if(parIndex1==parIndex2)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH2DIndex. Parameters have equal indeces: %d , %d", parIndex1, parIndex2));
		return -1;
	}

	if(parIndex1>parIndex2)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCIntegrate::GetH2DIndex. First parameters must be smaller than second (sorry).");
		return -1;
	}

	int index=0;
	for(int i=0;i<fNvar-1;i++)
		for(int j=i+1;j<fNvar;j++)
		{
			if(i==parIndex1 && j==parIndex2)
				return index;
			index++;
		}

	BCLog::Out(BCLog::warning, BCLog::warning,
		"BCIntegrate::GetH2DIndex. Invalid index combination.");

	return -1;
}

// *********************************************
TH2D * BCIntegrate::GetH2D(int parIndex1, int parIndex2)
{
	if(fHProb2D.size()==0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCModel::GetH2D. MarginalizeAll() has to be run prior to this to fill the distributions.");
		return 0;
	}

	int hindex = this -> GetH2DIndex(parIndex1, parIndex2);
	if(hindex==-1)
		return 0;

	if(hindex>(int)fHProb2D.size()-1)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCIntegrate::GetH2D. Got invalid index from GetH2DIndex(). Something went wrong.");
		return 0;
	}

	return fHProb2D[hindex];
}

// *********************************************
double BCIntegrate::GetRandomPoint(std::vector <double> &x)
{
	this->GetRndmVector(x);

	for(int i=0;i<fNvar;i++)
		x[i]=fMin[i]+x[i]*(fMax[i]-fMin[i]);

	return this->Eval(x);
}

// *********************************************
double BCIntegrate::GetRandomPointImportance(std::vector <double> &x)
{
	double p = 1.1; 
	double g = 1.0; 

	while (p>g)
	{
		this->GetRndmVector(x);

		for(int i=0;i<fNvar;i++)
			x[i] = fMin[i] + x[i] * (fMax[i]-fMin[i]);

		p = fRandom->Rndm();

		g = this -> EvalSampling(x);
	}

	return this->Eval(x);
}

// *********************************************
void BCIntegrate::InitMetro()
{
	fNmetro=0;

	// start in the center of the phase space
	for(int i=0;i<fNvar;i++)
		fXmetro0[i]=(fMin[i]+fMax[i])/2.;

	// run metropolis for a few times and dump the result... (to forget the initial position)
	std::vector <double> x;
	x.assign(fNvar,0.); 

	for(int i=0;i<1000;i++)
		GetRandomPointMetro(x);
}

// *********************************************
void BCIntegrate::GetRandomPointMetro(std::vector <double> &x)
{
	// define the stepsize relative to max-min
	// note: we need to make this setable in the future
	double stepsize=0.1;

	// get new point
	this->GetRndmVector(fXmetro1);

	// scale the point to the allowed region and stepsize
	int in=1;
	for(int i=0;i<fNvar;i++)
	{
		fXmetro1[i] = fXmetro0[i] + 2 * (fXmetro1[i]-0.5) * stepsize * (fMax[i]-fMin[i]);

		// check whether the generated point is inside the allowed region
		if( fXmetro1[i]<fMin[i] || fXmetro1[i]>fMax[i] )
			in=0;
	}

	// calculate the log probabilities and compare old and new point
	double p0 = this -> LogEval(fXmetro0); // old point
	double p1=0; // new point
	if(in)
		p1= this -> LogEval(fXmetro1);

	// compare log probabilities
	int accept=0;
	if(in)
	{
		if(p1>=p0)
			accept=1;
		else
		{
			double r=TMath::Log(fRandom->Rndm());
			if(r<p1-p0)
				accept=1;
		}
	}

	// fill the return point after the decision
	if(accept)
		for(int i=0;i<fNvar;i++)
		{
			fXmetro0[i]=fXmetro1[i];
			x[i]=fXmetro1[i];
		}
	else
		for(int i=0;i<fNvar;i++)
			x[i]=fXmetro0[i];
	
	fNmetro++;
}

// *********************************************
void BCIntegrate::GetRandomPointSamplingMetro(std::vector <double> &x)
{
	// define the stepsize relative to max-min
	// note: we need to make this setable in the future
	double stepsize=0.1;

	// get new point
	this->GetRndmVector(fXmetro1);

	// scale the point to the allowed region and stepsize
	int in=1;
	for(int i=0;i<fNvar;i++)
	{
		fXmetro1[i] = fXmetro0[i] + 2 * (fXmetro1[i]-0.5) * stepsize * (fMax[i]-fMin[i]);

		// check whether the generated point is inside the allowed region
		if( fXmetro1[i]<fMin[i] || fXmetro1[i]>fMax[i] )
			in=0;
	}

	// calculate the log probabilities and compare old and new point
	double p0 = this -> LogEvalSampling(fXmetro0); // old point
	double p1=0; // new point
	if(in)
		p1= this -> LogEvalSampling(fXmetro1);

	// compare log probabilities
	int accept=0;
	if(in)
	{
		if(p1>=p0)
			accept=1;
		else
		{
			double r=TMath::Log(fRandom->Rndm());
			if(r<p1-p0)
				accept=1;
		}
	}

	// fill the return point after the decision
	if(accept)
		for(int i=0;i<fNvar;i++)
		{
			fXmetro0[i]=fXmetro1[i];
			x[i]=fXmetro1[i];
		}
	else
		for(int i=0;i<fNvar;i++)
			x[i]=fXmetro0[i];
	
	fNmetro++;
}

// *********************************************
void BCIntegrate::GetRndmVector(std::vector <double> &x)
{
	double * randx = new double[fNvar];

	fRandom -> RndmArray(fNvar, randx);

	for(int i=0;i<fNvar;i++)
		x[i] = randx[i];

	delete[] randx;
	randx = 0;
}

// *********************************************
void BCIntegrate::FindModeSA()
{
	// number of initial samples to determine starting temperature
	int npresamples = fNvar*10000;

	// point with maximum -log(Eval())
	vector <double> xmax;
	xmax.assign(fNvar, 0.0);

	// initial maximum value and location
	double lastmax = -TMath::Log(this->GetRandomPoint(xmax));

	vector <double> x;
	x.assign(fNvar, 0.0);

	// determine the mean of -log(Eval())
	double mean=0;
	for(int i=0;i<npresamples;i++)
	{
		this->GetRandomPoint(x);
		double val = -this->LogEval(x);
		mean += val;

		// store new location of maximum if found
		if(val < lastmax)
		{
			lastmax=val;
			xmax.clear();
			for(int j=0;j<fNvar;j++)
				xmax.push_back(x[j]);
		}
	}
	mean /= (double)npresamples;

	// start in the starting point
	for(int i=0;i<fNvar;i++)
		fXmetro0[i]=xmax[i];

	double T = mean/2.; // set starting temperature to mean/2
	double xstep = .1; // define steps in individual parameters relative to range
	int nSmin = fNvar*10000;  // minimum number of samples per stage
	int nSmax = fNvar*100000; // maximum number of samples per stage

	// factor defining the cooling schedule
	double factorT = 1.2; // lower T in next stage by factor

	BCLog::Out(BCLog::debug, BCLog::debug,
		Form("BCIntegrate::FindModeSA. Starting SA mode finding after %d initial iterations with T=%f",npresamples,T));

	int nsave=1000; // number of points to store
	// array to store last nsave points for RMS calculation
	double * ysave = new double[nsave];

	// clear mode
	xmax.clear();
//	int nomode=1;
	double maxval=0.;

	// minimum temperature
	// when reached, stop mode finding
	double Tmin=.1;

	// initiate rms
	double lastrms=2.*T;

	int j=0;

	// count total number of temperature stages
	int stage=0;

	// count total number of iterations
	int niter=0;

	// stop if Tmin is reached
	while(T>Tmin)
	{
		stage++;

		int i=0;
		j=0;

		// do minimum nSmin samples and stop if RMS of last nsave points is small enough
		// or if nSmax samples is reached
		while(  i<nSmax && (i<nSmin || (i>=nSmin && lastrms>T) ) )
		{
			niter++;

			// get random point from the SA algorithm
			this->GetRandomPointSA(x,T,xstep);

			// get -log of the value at generated point
			double val = - this->LogEval(x);

			// keep last nsave points if close to nSmin
			if(i>nSmin-nsave-1)
			{
				// save current point
				ysave[j%nsave] = val;
				j++;
			}

			if(i>=nSmin-1)
				// calculate RMS of last nsave points
				lastrms = TMath::RMS(nsave,ysave);

			// save the location of maximum
			// as we're looking at -log Eval(), maximum is actually minimum
			if(val<maxval)
			{
				maxval=val;
				for(int k=0;k<fNvar;k++)
					xmax.push_back(x[k]);
			}

			i++;
		}

		BCLog::Out(BCLog::debug, BCLog::debug,
			Form("BCIntegrate::FindModeSA. Stage %d finished at T=%f after %d samples with RMS = %f",stage,T,i,lastrms));

		if(i>=nSmax && lastrms>T)
			BCLog::Out(BCLog::detail, BCLog::detail,
				Form(" --> SA not converged in stage %d for T=%f in %d iterations (RMS reached is %f)",stage,T,nSmax,lastrms));

		// adjust temperature for next stage
		T /= factorT;
	}

	// fill the mode
	fBestFitParameters.clear();
	for(int i=0;i<fNvar;i++)
		fBestFitParameters.push_back(xmax[i]);

	BCLog::Out(BCLog::detail, BCLog::detail,
		Form(" --> SA mode finding finished in %d stages after %d iterations",stage,niter));

	delete[] ysave;
}

// *********************************************
void BCIntegrate::GetRandomPointSA(std::vector <double> &x, double T, double step)
{
	// get new point
	this->GetRndmVector(fXmetro1);

	// scale the point to the allowed region and step
	int in=1;
	for(int i=0;i<fNvar;i++)
	{
		fXmetro1[i] = fXmetro0[i] + 2 * (fXmetro1[i]-0.5) * step * (fMax[i]-fMin[i]);

		// check whether the generated point is inside the allowed region
		if( fXmetro1[i]<fMin[i] || fXmetro1[i]>fMax[i] )
			in=0;
	}

	// calculate the log probabilities and compare old and new point
	double p0 = - this->LogEval(fXmetro0); // old point
	double p1 = 0; // new point

	// compare
	int accept=0;
	if(in)
	{
		p1 = - this->LogEval(fXmetro1);

		if(p1<=p0)
			accept=1;
		else
		{
			double pval = TMath::Exp(-(p1-p0)/T);
			double r=fRandom->Rndm();
			if(r<pval)
				accept=1;
		}
	}

	// fill the return point after the decision
	if(accept)
		for(int i=0;i<fNvar;i++)
		{
			fXmetro0[i]=fXmetro1[i];
			x[i]=fXmetro1[i];
		}
	else
		for(int i=0;i<fNvar;i++)
			x[i]=fXmetro0[i];
	
	fNmetro++;
}

// *********************************************
void BCIntegrate::CubaIntegrand(const int *ndim, const double xx[], 
				const int *ncomp, double ff[]) 
{
	// scale variables 
	double jacobian = 1.0; 

	std::vector<double> scaled_parameters; 

	for (int i = 0; i < *ndim; i++)
	{
		double range = global_this -> fx -> at(i) -> GetUpperLimit() -  global_this -> fx -> at(i) -> GetLowerLimit(); 

		// multiply range to jacobian 
		jacobian *= range; 

		// get the scaled parameter value 
		scaled_parameters.push_back(global_this -> fx -> at(i) -> GetLowerLimit() + xx[i] * range); 
	}

	// call function to integrate 
	ff[0] = global_this -> Eval(scaled_parameters); 

	// multiply jacobian 
	ff[0] *= jacobian; 

	// multiply fudge factor 
	ff[0] *= 1e99; 

	// remove parameter vector 
	scaled_parameters.clear(); 
}

// *********************************************
double BCIntegrate::CubaIntegrate() 
{
	double EPSREL = 1e-3; 
	double EPSABS = 1e-12; 
	int VERBOSE   = 0; 
	int MINEVAL   = 0; 
	int MAXEVAL   = 2000000; 
	int NSTART    = 25000; 
	int NINCREASE = 25000; 

	std::vector<double> parameters_double; 
	std::vector<int>    parameters_int; 

	parameters_double.push_back(EPSREL); 
	parameters_double.push_back(EPSABS); 

	parameters_int.push_back(VERBOSE); 
	parameters_int.push_back(MINEVAL); 
	parameters_int.push_back(MAXEVAL); 
	parameters_int.push_back(NSTART); 
	parameters_int.push_back(NINCREASE); 

	return this -> CubaIntegrate(0, parameters_double, parameters_int); 
}

// *********************************************
double BCIntegrate::CubaIntegrate(int method, std::vector<double> parameters_double, std::vector<int> parameters_int) 
{
	const int NDIM      = int(fx ->size()); 
	const int NCOMP     = 1; 

	const double EPSREL = parameters_double[0]; 
	const double EPSABS = parameters_double[1]; 
	const int VERBOSE   = parameters_int[0]; 
	const int MINEVAL   = parameters_int[1]; 
	const int MAXEVAL   = parameters_int[2]; 

	int neval;
	int fail;
	int nregions; 
	double integral[NCOMP];
	double error[NCOMP];
	double prob[NCOMP];

	global_this = this; 

	integrand_t an_integrand = &BCIntegrate::CubaIntegrand; 

	if (method == 0)
	{
		// set VEGAS specific parameters 
		const int NSTART    = parameters_int[3]; 
		const int NINCREASE = parameters_int[4]; 

		// call VEGAS integration method  
		Vegas(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
			NSTART, NINCREASE,
			&neval, &fail, integral, error, prob);
	}
	else if (method == 1)
	{
		// set SUAVE specific parameters 
		const int LAST     = parameters_int[3]; 
		const int NNEW     = parameters_int[4]; 
		const int FLATNESS = parameters_int[5]; 

		// call SUAVE integration method 
		Suave(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL,
			NNEW, FLATNESS,
			&nregions, &neval, &fail, integral, error, prob);
	}
	else if (method == 2)
	{
		// set DIVONNE specific parameters 
		const int KEY1         = parameters_int[3]; 
		const int KEY2         = parameters_int[4]; 
		const int KEY3         = parameters_int[5]; 
		const int MAXPASS      = parameters_int[6]; 
		const int BORDER       = parameters_int[7]; 
		const int MAXCHISQ     = parameters_int[8]; 
		const int MINDEVIATION = parameters_int[9]; 
		const int NGIVEN       = parameters_int[10]; 
		const int LDXGIVEN     = parameters_int[11]; 
		const int NEXTRA       = parameters_int[12]; 

		// call DIVONNE integration method 
		Divonne(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
			KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION,
			NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
			&nregions, &neval, &fail, integral, error, prob);
	}
	else if (method == 3) 
	{
		// set CUHRE specific parameters 
		const int LAST = parameters_int[3]; 
		const int KEY  = parameters_int[4]; 

		// call CUHRE integration method 
		Cuhre(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL, KEY,
			&nregions, &neval, &fail, integral, error, prob);
	}
	else
	{
		std::cout << " Integration method not available. " << std::endl; 
		integral[0] = -1e99; 
	}

	if (fail != 0) 
	{
		std::cout << " Warning, integral did not converge with the given set of parameters. "<< std::endl; 
		std::cout << " neval    = " << neval       << std::endl; 
		std::cout << " fail     = " << fail        << std::endl; 
		std::cout << " integral = " << integral[0] << std::endl; 
		std::cout << " error    = " << error[0]    << std::endl; 
		std::cout << " prob     = " << prob[0]     << std::endl; 
	}

	return integral[0] / 1e99; 
}
