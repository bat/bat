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
        // kApprox: use approximation
        enum EvalOption{ kAnalytic, kHistogram, kApprox };

        // Constructors and destructor
        ReferenceCounting();
        ReferenceCounting(const char * name);
        ~ReferenceCounting();

        // Methods to overload, see file ReferenceCounting.cxx
        void DefineParameters();
        double LogAPrioriProbability(const std::vector<double> &parameters);
        double LogLikelihood(const std::vector<double> &parameters);

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
        double df(double s, int k, int n);
        double logdf(double s, int k, int n);

        // helper functions for histogram
        void FillPriorS();
        void FillPriorB();

        // print prior histograms
        void PrintPriors(std::string filename);

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

        // the approximate priors as functions
        TF1* fFuncPriorS;
        TF1* fFuncPriorB;


};
// ---------------------------------------------------------

#endif

