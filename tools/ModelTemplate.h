// ***************************************************************
// This file was created using the ((PROGRAM)) script.
// ((PROGRAM)) is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__((UPMODEL))__H
#define __BAT__((UPMODEL))__H

#include <BAT/BCModel.h>

#include <string>
#include <vector>

// This is a ((MODEL)) header file.
// Model source code is located in file ((PROJECT))/((MODEL)).cxx

// ---------------------------------------------------------
class ((MODEL)) : public BCModel
{

public:

    // Constructor
    ((MODEL))(const std::string& name);

    // Destructor
    ~((MODEL))();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& parameters);

    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & parameters);

};
// ---------------------------------------------------------

#endif
