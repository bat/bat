#include "ReferenceCounting.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>
#include <BAT/BCH1D.h>

#include <TMath.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>

#include <iostream>
#include <iomanip>

// ---------------------------------------------------------
ReferenceCounting::ReferenceCounting() : BCModel()
                                       , fEvalOption(kHistogram)
                                       , fNObs(0)
                                       , fEps(1e-6)
                                       , fAlpha(10)
                                       , fBeta(5)
                                       , logs(0)
                                       , maxn(0)
                                       , maxk(10000)
                                       , fHistPriorS(0)
                                       , fHistPriorB(0)
{
   // default constructor
   DefineParameters();
}

// ---------------------------------------------------------
ReferenceCounting::ReferenceCounting(const char * name) : BCModel(name)
                                                        , fEvalOption(kHistogram)
                                                        , fNObs(0)
                                                        , fEps(1e-6)
                                                        , fAlpha(10)
                                                        , fBeta(5)
                                                        , logs(0)
                                                        , maxn(0)
                                                        , maxk(10000)
                                                        , fHistPriorS(0)
                                                        , fHistPriorB(0)
{
   // constructor
   DefineParameters();
}

// ---------------------------------------------------------
ReferenceCounting::~ReferenceCounting()
// default destructor
{
   if (fHistPriorS)
      delete fHistPriorS;

   if (fHistPriorB)
      delete fHistPriorB; 
}

// ---------------------------------------------------------
void ReferenceCounting::SetAlphaBeta(double alpha, double beta)
{
   fAlpha = alpha; 
   fBeta  = beta;

   // clear vectors
   helper_a.clear();
   helper_b.clear();
   helper_c.clear();
}

// ---------------------------------------------------------
void ReferenceCounting::DefineParameters()
{
   AddParameter("s", 0, 50.); // index 0
   AddParameter("b", 0, 50.); // index 1
}

// ---------------------------------------------------------
double ReferenceCounting::LogLikelihood(const std::vector<double> &parameters)
{
   if (!fHistPriorS)
      FillPriorS();

   if (!fHistPriorB)
      FillPriorB();

   double logprob = 0.;

   double s = parameters.at(0);
   double b = parameters.at(1);
   double nu = s+b;

   logprob += BCMath::LogPoisson( double(fNObs), nu);

   return logprob;
}

// ---------------------------------------------------------
double ReferenceCounting::LogAPrioriProbability(const std::vector<double> &parameters)
{
   double logprob = 0.;

   double s = parameters.at(0);
   double b = parameters.at(1);

   if (fEvalOption == kAnalytic) {
      logprob += log(RefPriorS(s));
      logprob += log(ConjPriorPoisson(b, fAlpha, fBeta));
   }
   else if (fEvalOption == kHistogram) {
      int bins = fHistPriorS->FindBin(s);
      int binb = fHistPriorB->FindBin(b);
      logprob += log( fHistPriorS->GetBinContent(bins) );
      logprob += log( fHistPriorB->GetBinContent(binb) );
   }
   else if (fEvalOption == kApprox) {
      logprob += log(fFuncPriorS->Eval(s));
      logprob += log(ConjPriorPoisson(b, fAlpha, fBeta));
   }

   return logprob;
}

// ---------------------------------------------------------
double ReferenceCounting::ConjPriorPoisson(double nu, double alpha, double beta)
{
   return TMath::GammaDist(nu, alpha, 0., 1/beta);
}

// ---------------------------------------------------------
double ReferenceCounting::RefPriorS(double s)
{
   return sqrtI(s) / sqrtI(0.);
}

// ---------------------------------------------------------
double ReferenceCounting::sqrtI(double s)
{
   double sum = 0;

   double inc = 1.;

   for (int n=0; n<=maxk; ++n) {
      double t1 = f(s, n);
      double t2 = f(s, n+1);
      double dsum = t1*t1/t2;
      sum += dsum;
      inc = fabs( dsum/sum );

      // check if sum does not change by more than predefined precision
      if (inc < fEps)
         break;
   }

   return sqrt( fabs( TMath::Power(fBeta/(1+fBeta), fAlpha) * exp(-s) * sum - 1 ) );
}

// ---------------------------------------------------------
double ReferenceCounting::f(double s, int k)
{
   if (s > 0)
      logs = log(s);
   else {
      if (k==0) 
         return 1;
      else 
         return (fAlpha+k-1)/(k*(1+fBeta)) * f(0., k-1);
   }

   double ftemp = 0;

   extend_abc(k+2);

   for (int n = 0; n <= k; ++n) {
      ftemp += df(s, k, n);
   }

   return ftemp;
}

// ---------------------------------------------------------
double ReferenceCounting::df(double s, int k, int n)
{
   return exp(logdf(s, k, n));
}

// ---------------------------------------------------------
double ReferenceCounting::logdf(double s, int k, int n)
{
   double logdf = 0;

   logdf = helper_a[n] + (k-n) * logs - helper_b[n] - helper_c[k-n];

   return logdf;
}

// ---------------------------------------------------------
void ReferenceCounting::extend_abc(int n)
{
   if (n <= maxn)
      return;
	
   int nfrom = maxn;
	
   if (nfrom == 0) {
      helper_a.push_back( log(1) );
      helper_b.push_back( 0. );
      helper_c.push_back( 0. );
      ++nfrom;
   }
   if (nfrom == 1) {
      helper_a.push_back( log(fAlpha) );
      helper_b.push_back( log( 1. + fBeta) );
      helper_c.push_back( helper_c[0] );
      ++nfrom;
   }
	
   if (maxn > 1)
      nfrom++;

   for (int i = nfrom; i <= n; ++i) {
      double temp_a = helper_a[i-1] + log( 1. + (fAlpha-1.)/double(i)); 
      helper_a.push_back(temp_a);
      helper_b.push_back( int(i) * helper_b[1] );
      helper_c.push_back( log( double(i) ) + helper_c[i-1]); 
   }

   maxn = int( helper_a.size() - 1);
}

// ---------------------------------------------------------
void ReferenceCounting::FillPriorS()
{
   // remove old prior histogram
   if (fHistPriorS)
      delete fHistPriorS;

   // create new histogram
   BCParameter* p = GetParameter(0);
   fHistPriorS = new TH1D("hist_prior_s", ";s;p(s)", p->GetNbins(), p->GetLowerLimit(), p->GetUpperLimit());
	
   // fill histogram
   for (unsigned i = 1; i <= p->GetNbins(); ++i){
      double s = fHistPriorS->GetBinCenter(i);
      double p = RefPriorS(s);
      fHistPriorS->SetBinContent(i, p);
   }
   
   // fit histogram

   // remove old prior function
   if (fFuncPriorS)
      delete fFuncPriorS;

   //   fFuncPriorS = new TF1("func_prior_s", "sqrt([0] / ([0] + x)) * exp([1]*TMath::Power(x, 0.25))", p->GetLowerLimit(), p->GetUpperLimit());
   //   fFuncPriorS = new TF1("func_prior_s", "exp([0]*x^0.25+[1]*x^0.5+[2]*x)", p->GetLowerLimit(), p->GetUpperLimit());
   fFuncPriorS = new TF1("func_prior_s", "sqrt(([0]*exp([1]*x^0.125))/(x+[0]*exp([1]*x^0.125)))", p->GetLowerLimit(), p->GetUpperLimit());

   fHistPriorS->Fit(fFuncPriorS);
}

// ---------------------------------------------------------
void ReferenceCounting::FillPriorB()
{
   // remove old prior histogram
   if (fHistPriorB)
      delete fHistPriorB;

   // create new histogram
   BCParameter * p = GetParameter(1);
   fHistPriorB = new TH1D("hist_prior_b", ";b;p(b)", p->GetNbins(), p->GetLowerLimit(), p->GetUpperLimit());
	
   // fill histogram
   for (unsigned i = 1; i <= p->GetNbins(); ++i){
      double b = fHistPriorB->GetBinCenter(i);
      double p = ConjPriorPoisson(b, fAlpha, fBeta);
      fHistPriorB->SetBinContent(i, p);
   }
}

// ---------------------------------------------------------
void ReferenceCounting::PrintPriors(std::string filename)
{
   if (!fHistPriorS || !fHistPriorB) 
      return;

   TCanvas* c1 = new TCanvas("c1");
   c1->cd();
   fHistPriorS->Draw();
   c1->Print(std::string( filename + "(").c_str());

   fHistPriorB->Draw();
   c1->Print(std::string( filename + ")").c_str());

   return;   
}

// ---------------------------------------------------------
