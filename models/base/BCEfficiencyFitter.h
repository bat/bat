#ifndef __BCEFFICIENCYFITTER__H
#define __BCEFFICIENCYFITTER__H

/*!
 * \class BCEfficiencyFitter
 * \brief A class for fitting histograms with functions
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 11.2008
 * \detail This class allows fitting of efficiencies defined as
 * a ratio of two TH1D histograms using a TF1 function. It uses
 * binomial probabilities calculated based on the number of entries
 * in histograms. This is only applicable if the numerator is
 * a subset of the denominator.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>

#include "BCFitter.h"

// ROOT classes
class TH1D;
class TF1;
class TGraphAsymmErrors;

// ---------------------------------------------------------

class BCEfficiencyFitter : public BCFitter
{
   public:

      /**Abstract class which doesn't do anything
       * but offers the right interface
       * to allow calculation the distribution of any statistic.
       * User has to create a subclass and implement the operator().
       */
      class ToyDataInterface
      {
      public:
         /**
          * operator() is called for each generated toy data set of the fast p-value calculation.
          * @param expectation the expected number of events for the parameter values
          * chosen in the call to CalculatePValueFast
          * @param toyData one toy data set
          */
         virtual void operator()(const std::vector<double>& expectation, const std::vector<unsigned>& toyData) = 0;

         /**pure abstract */
         virtual ~ToyDataInterface(){}
      };

      /** \name Constructors and destructors */
      /* @{ */

      /**
       * The default constructor. */
      BCEfficiencyFitter();

      /**
       * Constructor
       * @param name name fo the model */
      BCEfficiencyFitter(const char * name);

      /**
       * A constructor.
       * @param hist1 The histogram with the larger numbers
       * @param hist2 The histogram with the smaller numbers
       * @param func The fit function. */
      BCEfficiencyFitter(TH1D * hist1, TH1D * hist2, TF1 * func);

      /**
       * Constructor.
       * @param name name fo the model
       * @param hist1 The histogram with the larger numbers
       * @param hist2 The histogram with the smaller numbers
       * @param func The fit function. */
      BCEfficiencyFitter(const char * name, TH1D * hist1, TH1D * hist2, TF1 * func);

      /**
       * The default destructor. */
      ~BCEfficiencyFitter();

      /* @} */
      /** \name Member functions (get) */
      /* @{ */

      /**
       * @return The histogram 1 */
      TH1D * GetHistogram1()
         { return fHistogram1; };

      /**
       * @return The histogram 2 */
      TH1D * GetHistogram2()
         { return fHistogram2; };

      /**
       * @return The histogram ratio */
      TGraphAsymmErrors * GetHistogramRatio()
         { return fHistogramRatio; };

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

      /**
       * Calculates the central value and the lower and upper limits for a given probability.
       * @param n n for the binomial.
       * @param k k for the binomial.
       * @param p The central probability defining the limits.
       * @param xmin The central value.
       * @param xmin The lower limit.
       * @param xmax The upper limit.
       * @return A flag (=1 plot point, !=1 do not plot point). */
      int GetUncertainties(int n, int k, double p, double &xexp, double &xmin, double &xmax);

      /* @} */
      /** \name Member functions (set) */
      /* @{ */

      /**
       * @param hist The histogram 1
       * @param hist The histogram 2
       * @return An error code (1:pass, 0:fail). */
      int SetHistograms(TH1D * hist1, TH1D * hist2);

      /**
       * @param func The fit function
       * @return An error code (1:pass, 0:fail). */
      int SetFitFunction(TF1 * func);

      /**
       * Sets the flag for integration. \n
       * true: use ROOT's TH1D::Integrate() \n
       * false: use linear interpolation */
      void SetFlagIntegration(bool flag)
         { fFlagIntegration = flag; };

      /** Set type of point to be used to plot the efficiency data
        * 0 - mean + rms
        * 1 - mode + smallest 68%
        * 2 - median + central 68% */
      void SetDataPointType(int type);

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
       * @param hist1 The histogram with the larger number.
       * @param hist2 The histogram with the smaller number.
       * @param func The fit function.
       * @return An error code. */
      int Fit(TH1D * hist1, TH1D * hist2, TF1 * func);

      /**
       * Draw the fit in the current pad. */
      void DrawFit(const char * options = "", bool flaglegend = false);

      /**
      * Calculate the p-value using fast-MCMC. In every iteration, a new toy data set is created.
      * By providing a suitable implementation of ToyDataInterface, the user can
      * calculate the distribution of an arbitrary statistic. Each toy data set as well as the
      * expected values for the parameter values are passed on to the interface.
      * @param par A set of parameter values
      * @param pvalue The p-value for the default likelihood statistic
      * @param pvalueCorrected The p-value for the default likelihood statistic
      *                        with corrections for DoF to reduce fitting bias, in all nonpathological
      *                        cases the correction lowers the p-value.
      * @param callback requires class with operator(...) defined.
      * @param nIterations number of toy data sets generated by the Markov chain
      * @return An error code */
      int CalculatePValueFast(const std::vector<double> & par,
                              BCEfficiencyFitter::ToyDataInterface * callback, double & pvalue,
                              double & pvalueCorrected, unsigned nIterations = 100000);

      /**
       * Calculate the p-value using fast-MCMC.
       * @param par A set of parameter values
       * @param  pvalue The pvalue
       * @param pvalueCorrected The p-value for the default likelihood statistic
       *                        with corrections for DoF to reduce fitting bias
       * @param nIterations number of pseudo experiments generated by the Markov chain
       * @return An error code */
      int CalculatePValueFast(const std::vector<double> & par, double & pvalue,
                              double & pvalueCorrected, unsigned nIterations = 100000);

      /* @} */

   private:

      /**
       * The histogram containing the larger numbers. */
      TH1D * fHistogram1;

      /**
       * The histogram containing the smaller numbers. */
      TH1D * fHistogram2;

      /**
       * The efficiency histogram. */
      TGraphAsymmErrors * fHistogramRatio;

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
       * Temporary histogram for calculating the binomial qunatiles */
      TH1D * fHistogramBinomial;

      /** Type of point to plot for efficiency data
        * 0 - mean + rms
        * 1 - mode + smallest 68%
        * 2 - median + central 68% */
      unsigned int fDataPointType;

};

// ---------------------------------------------------------

#endif
