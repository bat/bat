#include "GaussModel.h"

#include <test.h>
#include <BAT/BCMath.h>

// ---------------------------------------------------------
GaussModel::GaussModel(const char * name, const unsigned & nParameters, long loopIterations) :
    BCModel(name),
    fLoopIterations(loopIterations)
{
   // add identical, independent parameters
   for (unsigned i = 0; i < nParameters ; i++)
   {
      std::string parName("par");
       AddParameter((parName + test::stringify(i)).c_str(), -15.0, 15.0);
   }
}

// ---------------------------------------------------------
GaussModel::~GaussModel()
{
}

// ---------------------------------------------------------
double GaussModel::LogLikelihood(const std::vector<double> & parameters)
{
    // run extra loop to make likelihood evaluation slower(!)
    {
        std::vector<double> chain_means(3);
        chain_means[0] = 4.2;
        chain_means[1] = 4.25;
        chain_means[2] = 4.22;
        std::vector<double> chain_variances(3);
        chain_variances[0] = 0.1;
        chain_variances[1] = 0.15;
        chain_variances[2] = 0.19;

        const bool relaxed = false;
        const unsigned points = 500;

        for(unsigned long i = 0 ; i < fLoopIterations ; ++i)
            BCMath::Rvalue(chain_means, chain_variances, points, relaxed);
    }

    // assume a simple Gaussian Likelihood with N independent variables
    double logprob = 0.;
    for (unsigned i=0; i < parameters.size(); i++)
    {
      // Gaussian Likelihood
      logprob += BCMath::LogGaus(parameters.at(i), 0.0, 2.0);
    }
    return logprob;
}

// ---------------------------------------------------------
double GaussModel::LogAPrioriProbability(const std::vector<double> &)
{
    // assume flat prior in all variables
    return 0.0;
}
