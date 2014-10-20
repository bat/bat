
#include <IntegrationModel.h>

#include <TString.h>
#include <TMath.h>
#include <Math/SpecFuncMathMore.h>
#include <TRandom3.h>

#include <cmath>
#include <iostream>

using namespace std;

IntegrationModel::IntegrationModel() : BCModel() {
   SetDimensionality(1);
   SetModality(0);
   SetComplexity(0);
}

IntegrationModel::~IntegrationModel() {
}

bool IntegrationModel::DefineParameters() {
   ClearParameters(true);
   for (unsigned int k = 0; k < GetDimensionality(); k++)
      AddParameter(Form("x%0d",k),0,1);
   SetPriorConstantAll();
   return true;
}

double IntegrationModel::Phi(unsigned int m, std::vector<double> x) const {
   double output = 1.;
   for (unsigned int k = 0; k < fDimensionality; k++)
      output *= S(fPolynomialDegrees[m][k],x[k]);
   return output;
}

double IntegrationModel::S(unsigned int i, double x) const
{
   return sqrt(2. * i + 1) * ROOT::Math::legendre(i,2*x-1);
}

double IntegrationModel::Likelihood(const std::vector<double> &parameters) {
   double output = 0;
   for (unsigned int m = 0; m < fComplexity; m++)
      output += Phi(m,parameters);
   return output*output;
}

double IntegrationModel::LogLikelihood(const std::vector<double> &parameters) {
   return log(Likelihood(parameters));
}

double IntegrationModel::Integral() {
   double sum = 0;
   for (unsigned int m = 0; m < GetComplexity(); m++)
      for (unsigned int n = 0; n < GetComplexity(); n++) {
         double product = 1;
         for (unsigned int k = 0; k < GetDimensionality(); k++)
            product *= (fPolynomialDegrees[m][k] == fPolynomialDegrees[n][k]);
         sum += product;
      }
   return sum;
}

std::vector<unsigned int> IntegrationModel::Divisors(unsigned int d) {
   std::vector<unsigned int> ds;
   for (unsigned int i=1; i <= d; i++)
      if (d % i == 0)
         ds.push_back(i);
   return ds;
}

void IntegrationModel::PopulatePolynomialDegrees() {
   fPolynomialDegrees.clear();
   fPolynomialModality.clear();
   unsigned int d = GetModality();
   for (unsigned int m = 0; m < GetComplexity(); m++) {
      unsigned int dm = (m == GetComplexity()-1) ? d :
            fRandom->Integer(d-(GetComplexity()-1-m)-1)+1;
      fPolynomialModality.push_back(dm);
      d -= dm;
      std::vector<unsigned int> temp1;
      for (unsigned int k = 0; k < GetDimensionality()-1; k++) {
         std::vector<unsigned int> ds = Divisors(dm);
         temp1.push_back(ds[fRandom->Integer(ds.size())]);
         dm = dm/temp1.back();
      }
      temp1.push_back(dm);
      std::vector<unsigned int> temp2;
      for (unsigned int k = 0; k < GetDimensionality(); k++) {
         int i = fRandom->Integer(temp1.size());
         temp2.push_back(Degree(temp1[i]));
         temp1.erase(temp1.begin()+i);
      }
      fPolynomialDegrees.push_back(temp2);
   }
}

void IntegrationModel::SetDimensionality(unsigned int k){
   fDimensionality = k;
   DefineParameters();
}
