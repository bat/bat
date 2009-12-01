// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project TestProject
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "TestModel.h"

#include <TMath.h>
#include <BAT/BCMath.h> 

// ---------------------------------------------------------
TestModel::TestModel() : BCModel()
{  // default constructor
	DefineParameters();

	// set data
	double x[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}; 
	double y[10] = {1.1, 2.1, 2.9, 4.2, 5.1, 6.2, 7.3, 7.8, 9.0, 10.1}; 
	
	for (int i = 0; i < 10; ++i) {
		fDataX[i] = x[i]; 
		fDataY[i] = y[i]; 
	}

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

	this -> AddParameter("offset",  -1.0, 2.0);
 	this -> AddParameter("slope",  0.8, 1.5);
}

// ---------------------------------------------------------
double TestModel::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob = 0.;

	for (int i = 0; i < 10; ++i) {
		double x = fDataX[i]; 
		double yexp = x * parameters.at(1) + parameters.at(0); 
		logprob += BCMath::LogGaus(yexp, fDataY[i], 0.1); 
	}

	return logprob;
}

// ---------------------------------------------------------
double TestModel::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	logprob += BCMath::LogGaus(parameters.at(0), -0.2, 0.2); 
	logprob += BCMath::LogGaus(parameters.at(1),  1.2, 0.2); 

	return logprob;
}
// ---------------------------------------------------------

