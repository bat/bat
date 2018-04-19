// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__MYMOD__H
#define __BAT__MYMOD__H

#include <BAT/BCModel.h>

#include <TH1D.h>

#include <string>
#include <vector>

// This is a MyMod header file.
// Model source code is located in file MyTut/MyMod.cxx

// ---------------------------------------------------------
class MyMod : public BCModel
{

public:

    MyMod(const std::string& name);

    // Destructor
    ~MyMod();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);

    // Overload CalculateObservables if using observables
    void CalculateObservables(const std::vector<double>& pars);

private:

    TH1D fDataHistogram;

};
// ---------------------------------------------------------

#endif
