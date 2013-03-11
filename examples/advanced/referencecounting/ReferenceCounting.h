// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __BAT__REFERENCECOUNTING__H
#define __BAT__REFERENCECOUNTING__H

#include <BAT/BCModel.h>

// This is a ReferenceCounting header file.
// Model source code is located in file ReferenceCounting/ReferenceCounting.cxx

// ---------------------------------------------------------
class ReferenceCounting : public BCModel
{
   public:

	    // An enumerator for the prior evaluation options.
      // kAnalytic: use analytic approximation
      // kHistogram: fill a histogram with the prior before running
      // kExpo: use an expotential approximation for the signal prior
	    enum EvalOption{ kAnalytic, kHistogram, kExpo };

      // Constructors and destructor
      ReferenceCounting();
      ReferenceCounting(const char * name);
      ~ReferenceCounting();

      // Methods to overload, see file ReferenceCounting.cxx
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);
      // void MCMCIterationInterface();

			// set option of how to evaluate the reference prior
			void SetPriorEvalOption(ReferenceCounting::EvalOption option)
			{ fEvalOption = option; }; 

			void SetAlphaBeta(double alpha, double beta);
			void SetNObs(int n)
			{ fNObs = n; };

			// the reference prior for the signal strength
			double RefPriorS(double s);

			// the conjugate prior for the background
			double ConjPriorPoisson(double nu, double alpha, double beta);

			// option of how to evaluate the signal prior
			EvalOption fEvalOption;

			// helper functions for calculating priors
			double sqrtI(double s);

			void extend_abc(int n);
			double f(double s, int k);
			double df(double s, int k, int n) {
				return exp(logdf(s, k, n)); }
			double logdf(double s, int k, int n);

			// helper functions for histogram
			void FillPriorS();
			void FillPriorB();

			// number of observed events
			int fNObs;

			// the precision of the reference prior sum
			double fEps;

			// helper variables
			double fAlpha; 			// the shape variable for the conjugate priors
			double fBeta;	  		// the shape variable for the conjugate priors
			double logs; // lograrithm of s
			int maxn; // upper limit of the helpers
			int maxk; // maximum number of sums

			std::vector<double> helper_a;
			std::vector<double> helper_b;
			std::vector<double> helper_c;

			// the histogrammed priors
			TH1D* fHistPriorS;
			TH1D* fHistPriorB;

			
};
// ---------------------------------------------------------

#endif

