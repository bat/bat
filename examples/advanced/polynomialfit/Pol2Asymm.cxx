#include "Pol2Asymm.h"

#include <BAT/BCLog.h>
#include <BAT/BCParameter.h>

#include <TString.h>

// ---------------------------------------------------------
Pol2Asymm::Pol2Asymm() : Pol1Asymm()
{ // default constructor
}

// ---------------------------------------------------------
Pol2Asymm::Pol2Asymm(const char * name) : Pol1Asymm(name)
{ // constructor
}

// ---------------------------------------------------------
Pol2Asymm::~Pol2Asymm()
{ // destructor
}

// ---------------------------------------------------------
void Pol2Asymm::DefineParameters()
{
   // Adding all the parameters of the model and their allowed ranges
   // Keep track of index of the parameter as you need it to be able
   // to use the parameter lateron. Usage by using the parameter
   // name is also possible but due to the array searching it is very
   // CPU time expensive and therefore not recommended in routines
   // called for every iteration, like LogLikelihood() or
   // LogAPrioriProbability()
   this -> AddParameter("p0",  0.  ,  2.2);    // index 0
   this -> AddParameter("p1", -0.1 ,  0.05);  // index 1
   this -> AddParameter("p2",  0.0 ,  0.001); // index 2

   // Print parameter summary
   BCLog::OutSummary(
      Form("Model \'%s\' has %d parameters:",this->GetName().data(),this -> GetNParameters()));
   for(unsigned int i=0; i< this -> GetNParameters(); i++)
      BCLog::OutSummary(Form("   %d. %s    range: %g - %g",
                             i,
                             this -> GetParameter(i) -> GetName().data(),
                             this -> GetParameter(i) -> GetLowerLimit(),
                             this -> GetParameter(i) -> GetUpperLimit() ) );
}

// ---------------------------------------------------------
double Pol2Asymm::FitFunction(const std::vector<double> & x, const std::vector<double> & par)
{
   // function to fit with
   // get the parameters of the function
   double p0 = par[0];
   double p1 = par[1];
   double p2 = par[2];

   // calculate the value of the function at x[0] for a given set of parameters
   return p0 + p1*x[0] + p2*x[0]*x[0];
}

// ---------------------------------------------------------
