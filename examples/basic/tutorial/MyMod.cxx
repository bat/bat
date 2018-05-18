// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "MyMod.h"

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCMath.h>

#include <TF1.h>
#include <TMath.h>

// ---------------------------------------------------------
MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", "mass [GeV];count", 100, 5.0, 5.6)
{
    // create function to fill data according to
    TF1 data_func("data_func", "exp(-0.5*((x - 5.27926) / 0.04)^2)", 5.0, 5.6);

    // fill data histogram randomly from data_func 1,000 times
    fDataHistogram.FillRandom("data_func", 1e3);

    // add parameters for Gaussian distribution
    AddParameter("mu",    5.27, 5.29, "#mu", "[GeV]");
    GetParameters().Back().SetPrior(new BCGaussianPrior(5.28, 2e-3));

    AddParameter("sigma", 25e-3, 45e-3, "#sigma", "[GeV]");
    GetParameters().Back().SetPrior(new BCGaussianPrior(35e-3, 3e-3));

    AddParameter("height", 0, 10, "", "[events]");
    GetParameters().Back().SetPriorConstant();

    AddObservable("SignalYield", 900, 1100, "Y_{S}", "[events]");
    AddObservable("Resolution",
                  100. * GetParameter("sigma").GetLowerLimit() / GetParameter("mu").GetUpperLimit(),
                  100. * GetParameter("sigma").GetUpperLimit() / GetParameter("mu").GetLowerLimit(),
                  "#sigma / #mu", "[%]");
}


// ---------------------------------------------------------
MyMod::~MyMod()
{
    // destructor
}

// ---------------------------------------------------------
double MyMod::LogLikelihood(const std::vector<double>& pars)
{
    // store our log-likelihood as we loop through bins
    double LL = 0.;

    // loop over bins of our data
    for (int i = 1; i <= fDataHistogram.GetNbinsX(); ++i) {

        // retrieve observed number of events
        double x = fDataHistogram.GetBinContent(i);

        // retrieve bin center
        double m = fDataHistogram.GetBinCenter(i);

        // calculate expected number of events, using ROOT Gaus function * height
        double nu = pars[2] * TMath::Gaus(m, pars[0], pars[1], true);

        // add to log-likelihood sum
        LL += BCMath::LogPoisson(x, nu);

    }

    // return log-likelihood
    return LL;
}

// ---------------------------------------------------------
// double MyMod::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// ---------------------------------------------------------
void MyMod::CalculateObservables(const std::vector<double>& pars)
{
    // store total of number events expected
    double nu = 0;

    // loop over bins of our data
    for (int i = 1; i <= fDataHistogram.GetNbinsX(); ++i)
        // calculate expected number of events in that bin
        // and add to total expectation
        nu += pars[2] * TMath::Gaus(fDataHistogram.GetBinCenter(i), pars[0], pars[1], true);

    // store in the observable
    GetObservable(0) = nu;
    GetObservable(1) = 100. * pars[1] / pars[0];
}
