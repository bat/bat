// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project TestProject
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __TESTMODEL__H
#define __TESTMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------

class TestModel : public BCModel
{
	public:

		// Constructors and destructor
		TestModel();
		TestModel(const char * name);
		~TestModel();

		// Methods to overload, see file TestModel.cxx
		void DefineParameters();
		double LogAPrioriProbability(std::vector <double> parameters);
		double LogLikelihood(std::vector <double> parameters);

 private: 

		double fDataX[10]; 
		double fDataY[10]; 
};
// ---------------------------------------------------------

#endif

