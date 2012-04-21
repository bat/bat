#ifndef __BCMODELHISTOGRAMFITTER__H
#define __BCMODELHISTOGRAMFITTER__H

/*!
 * \class BCHistogramFitter
 * \brief A class for fitting histograms with functions
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 11.2008
 * \detail This class allows fitting of a TH1D histogram using
 * a TF1 function.
 */

/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>

#include "../../BAT/BCModel.h"

class BCHistogramFitterToyDataInterface;

// ROOT classes
class TH1D;
class TF1;

// ---------------------------------------------------------

class BCHistogramFitter : public BCModel
{
   public:

      /** \name Constructors and destructors */
      /* @{ */

      /**
       * The default constructor. */
      BCHistogramFitter();

      /**
       * Constructor
       * @param name name of the model */
      BCHistogramFitter(const char * name);

      /**
       * Constructor.
       * @param hist histogram to fit
       * @param func fit function */
      BCHistogramFitter(TH1D * hist, TF1 * func);

      /**
       * Constructor
       * @param name name of the model
       * @param hist histogram to fit
       * @param func fit function */
      BCHistogramFitter(const char * name, TH1D * hist, TF1 * func);

      /**
       * The default destructor. */
      ~BCHistogramFitter();

      /* @} */

      /** \name Member functions (get) */
      /* @{ */

      /**
       * @return The data histogram */
      TH1D * GetHistogram()
         { return fHistogram; };

      /**
       * @return The  histogram of expected counts*/
      TH1D * GetHistogramExpected()
         { return fHistogramExpected; };


      /**
       * @return The fit function */
      TF1 * GetFitFunction()
         { return fFitFunction; };

      /**
       * @return pointer to the error band */
      TGraph * GetErrorBand()
         { return fErrorBand; };

      /**
       * @return pointer to a graph for the fit function */
      TGraph * GetGraphFitFunction()
         { return fGraphFitFunction; };

      /* @} */

      /** \name Member functions (set) */
      /* @{ */

      /**
       * @param hist The histogram containing the data
       * @return An error code (1:pass, 0:fail).
       */
      int SetHistogram(TH1D * hist);

      /**
       * @param hist The histogram with the expected counts (typically non-integer values!)
       * @return An error code (1:pass, 0:fail).
       */
      int SetHistogramExpected(const std::vector<double>& parameters);

      /**
       * @param func The fit function
       * @return An error code (1:pass, 0:fail).
       */
      int SetFitFunction(TF1 * func);

      /**
       * Sets the flag for integration. \n
       * true: use ROOT's TH1D::Integrate() \n
       * false: use linear interpolation */
      void SetFlagIntegration(bool flag)
         { fFlagIntegration = flag; };

      /* @} */
      /** \name Member functions (miscellaneous methods) */
      /* @{ */

      /**
       * The log of the prior probability. Overloaded from BCModel.
       * @param parameters A vector of doubles containing the parameter values. */
//      virtual double LogAPrioriProbability(const std::vector<double> & parameters);

      /**
       * The log of the conditional probability. Overloaded from BCModel.
       * @param parameters A vector of doubles containing the parameter values. */
      virtual double LogLikelihood(const std::vector<double> & parameters);

      /**
       * Plots the histogram
       * @param options Options for plotting.
       * @param filename Name of the file which the histogram is printed into.
       * The following options are available:\n
       * F : plots the fit function on top of the data
       * E0 : plots the fit function and the 68% prob. uncertainty band of the fit function on top of the data
       * E1 : plots the expectation from the fit function and the uncertainty bin-by-bin as error bars. */
//      void PrintHistogram(const char * options = "", const char * filename = "");

      /**
       * Returns the y-value of the 1-dimensional fit function at an x and
       * for a set of parameters.
       * @param x A vector with the x-value.
       * @param parameters A set of parameters. */
      double FitFunction(const std::vector<double> & x, const std::vector<double> & parameters);

      /**
       * Performs the fit.
       * @return An error code. */
      int Fit();

      /**
       * Performs the fit.
       * @param hist The histogram (TH1D).
       * @param func The fit function.
       * @return An error code. */
      int Fit(TH1D * hist, TF1 * func);

      /**
       * Draw the fit in the current pad. */
      void DrawFit(const char * options = "HIST", bool flaglegend = false);

      /**
       * Calculate the p-value using fast-MCMC and the likelihood as test statistic. The method is explained in the appendix of
       * http://arxiv.org/abs/1011.1674
       * @param par A set of parameter values
       * @param  pvalue The p-value
       * @param nIterations number of pseudo experiments generated by the Markov chain
       * @return An error code */
      int CalculatePValueFast(const std::vector<double> & par, double &pvalue, int nIterations = 100000);

       /**
       * Calculate the p-value using fast-MCMC. In every iteration, a new toy data set is created.
       * By providing a suitable implementation of BCHistogramFitterToyDataInterface, the user can
       * calculate the distribution of an arbitrary statistic. Each toy data set as well as the
       * expected values for the parameter values are passed on to the interface.
       * @param par A set of parameter values
       * @param  pvalue The p-value for the default likelihood statistic
       * @param callback requires class with operator(...) defined.
       * @param nIterations number of toy data sets generated by the Markov chain
       * @return An error code */
      int CalculatePValueFast(const std::vector<double> & par, BCHistogramFitterToyDataInterface* callback,
            double &pvalue,  int nIterations = 100000);

      /**
       * Calculate the p-value using approximate chi^2 distribution of scaled likelihood.
       * Approximation is valid for bin contents >5, see eq. (32.12),
       * PDG: Statistics, Monte Carlo, Group Theory. Physics Letters B 667, 316-339(2008).
       * @param par The set of parameter values used in the model, usually the best fit parameters
       * @param  pvalue The pvalue
       * @return An error code */
      int CalculatePValueLikelihood(const std::vector<double> &par, double &pvalue);

      /**
       * Calculate the p-value using approximate chi^2 distribution of squared difference
       * for conventional weights.
       * Approximation is valid for bin contents >5 and not as as good for little data as
       * CalculatePValueLikelihood, see eq. (32.13),
       * PDG: Statistics, Monte Carlo, Group Theory. Physics Letters B 667, 316-339(2008).
       * @param par The set of parameter values used in the model, usually the best fit parameters
       * @param  pvalue The pvalue
       * @param  weight use the variance from the expected #counts (true) or the measured counts (false)
       * @return An error code */
      int CalculatePValueLeastSquares(const std::vector<double> &par, double &pvalue, bool weightExpect=true);

      /**
       * Calculate the p-value using Kolmogorov-Smirnov test. Note that the
       * reference distribution is known only asymptotically. Some explanation is given in
       * http://root.cern.ch/root/htmldoc/TMath.html
       * @param par The set of parameter values used in the model, usually the best fit parameters
       * @param  pvalue The pvalue
       * @return An error code */
      int CalculatePValueKolmogorov(const std::vector<double> &par, double &pvalue);


      double CDF(const std::vector<double>& parameters, int index, bool lower=false);

      /* @} */

   private:

      /**
       * The histogram containing the data.
       */
      TH1D * fHistogram;

      /**
       * The fit function */
      TF1 * fFitFunction;

      /**
       * Flag for using the ROOT TH1D::Integral method (true), or linear
       * interpolation (false) */
      bool fFlagIntegration;

      /**
       * Pointer to the error band (for legend) */
      TGraph * fErrorBand;

      /**
       * Pointer to a graph for displaying the fit function */
      TGraph * fGraphFitFunction;

      /**
       * The histogram containing the expected data.
       */
      TH1D * fHistogramExpected;

};

// ---------------------------------------------------------

/**Abstract class which doesn't do anything
 * but offers the right interface
 * to allow calculation the distribution of any statistic.
 * User has to create a subclass and implement the operator().
 */
class BCHistogramFitterToyDataInterface
{
public:
   /**
    * operator() is called for each generated toy data set of the fast p-value calculation.
    * @param expectation the expected number of events for the parameter values
    * chosen in the call to CalculatePValueFast
    * @param toyData one toy data set
    */
   virtual void operator()(const std::vector<double>& expectation, const std::vector<int>& toyData)=0;

   /**pure abstract */
   virtual ~BCHistogramFitterToyDataInterface(){}
};

#endif
