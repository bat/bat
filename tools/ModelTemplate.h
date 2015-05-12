// ***************************************************************
// This file was created using the |:PROGRAM:| script.
// |:PROGRAM:| is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__((UPMODEL))__H
#define __BAT__((UPMODEL))__H

#include <BAT/BCModel.h>

// This is a ((MODEL)) header file.
// Model source code is located in file ((PROJECT))/((MODEL)).cxx

// ---------------------------------------------------------
class ((MODEL)) : public BCModel
{
public:

    // Constructor and destructor
    ((MODEL))(std::string name);
    ~((MODEL))();

    // Methods tox overload, see file |:Model:|.cxx
    double LogLikelihood(const std::vector<double> & parameters);
    // double LogAPrioriProbability(const std::vector<double> & parameters);

};
// ---------------------------------------------------------

#endif
