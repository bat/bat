// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project TestProject
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "TestModel.h"

#include <TMath.h>

// ---------------------------------------------------------
TestModel::TestModel() : BCModel()
{  // default constructor
	DefineParameters();
};

// ---------------------------------------------------------
TestModel::TestModel(const char * name) : BCModel(name)
{  // constructor
	DefineParameters();
};

// ---------------------------------------------------------
TestModel::~TestModel()
{};  // default destructor

// ---------------------------------------------------------
void TestModel::DefineParameters()
{
	// Add parameters to your model here.
	// You can then use them in the methods below by calling the
	// parameters.at(i) or parameters[i], where i is the index
	// of the parameter. The indices increase from 0 according to the
	// order of adding the parameters.

	this -> AddParameter("par0", -5.0, 5.0);
 	this -> AddParameter("par1", -5.0, 5.0);
 	this -> AddParameter("par2", -5.0, 0.0);
}

// ---------------------------------------------------------
double TestModel::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob = 0.;

 	logprob += log ( 2.0 * TMath::Gaus(parameters.at(0), -3.0, 0.5) + 
 									 1.0 * TMath::Gaus(parameters.at(0), -1.0, 0.5) + 
									 1.0 * TMath::Gaus(parameters.at(0),  2.0, 0.7) );

	double rho = -0.8; 
 	double a = (parameters.at(1) - 1.0)/1.5; 
	double b = (parameters.at(2) + 2.0)/0.5; 
 	logprob += - 1.0/(1-rho*rho)*(a*a + b*b - 2*rho*a*b); 

	return logprob;
}

// ---------------------------------------------------------
double TestModel::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	return logprob;
}
// ---------------------------------------------------------

