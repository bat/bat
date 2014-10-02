// ***************************************************************
// This file was created using the |:PROGRAM:| script.
// |:PROGRAM:| is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__|:UP_MODEL:|__H
#define __BAT__|:UP_MODEL:|__H

#include <BAT/BCModel.h>

// This is a |:Model:| header file.
// Model source code is located in file |:Project:|/|:Model:|.cxx

// ---------------------------------------------------------
class |:Model:| : public BCModel {
 public:

	// Constructor and destructor
	|:Model:|(const char * name = "|:Model:|");
	~|:Model:|();

	// Methods to overload, see file |:Model:|.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	// double LogAPrioriProbability(const std::vector<double> & parameters);

};
// ---------------------------------------------------------

#endif
