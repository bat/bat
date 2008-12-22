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
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>

#include "BCModel.h"
#include "BCH1D.h"

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
		 * A constructor.
		 * @param hist The histogram (TH1D).
		 * @param func The fit function. */
		BCHistogramFitter(TH1D * hist, TF1 * func);

		/**
		 * The default destructor. */
		~BCHistogramFitter();

		/* @} */

		/** \name Member functions (get) */
		/* @{ */

		/**
		 * @return The histogram */
		TH1D * GetHistogram()
			{ return fHistogram; };

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
		 * @param hist The histogram
		 * @ return An error code (1:pass, 0:fail).
		 */
		int SetHistogram(TH1D * hist);

		/**
		 * @param func The fit function
		 * @ return An error code (1:pass, 0:fail).
		 */
		int SetFitFunction(TF1 * func);

		/**
		 * Sets the flag for integration. \n
		 * true: use ROOT's TH1D::Integrate() \n
		 * false: use linear interpolation 
		 */ 
		void SetFlagIntegration(bool flag)
		{ fFlagIntegration = flag; }; 

		/* @} */
		/** \name Member functions (miscellaneous methods) */
		/* @{ */

		/**
		 * The log of the prior probability. Overloaded from BCModel.
		 * @param parameters A vector of doubles containing the parameter values. */
		double LogAPrioriProbability(std::vector <double> parameters);

		/**
		 * The log of the conditional probability. Overloaded from BCModel.
		 * @param parameters A vector of doubles containing the parameter values. */
		double LogLikelihood(std::vector <double> parameters);

		/**
		 * Plots the histogram
		 * @param options Options for plotting.
		 * @param filename Name of the file which the histogram is printed into.
		 * The following options are available:\n
		 * F : plots the fit function on top of the data
		 * E0 : plots the fit function and the 68% prob. uncertainty band of the fit function on top of the data
		 * E1 : plots the expectation from the fit function and the uncertainty bin-by-bin as error bars. */
//		void PrintHistogram(const char * options = "", const char * filename = "");

		/**
		 * Returns the y-value of the 1-dimensional fit function at an x and
		 * for a set of parameters.
		 * @param x A vector with the x-value.
		 * @param parameters A set of parameters. */
		double FitFunction(std::vector <double> x, std::vector <double> parameters);

		/**
		 * Performs the fit.
		 * @return An error code. */
		int Fit()
			{ return this -> Fit(fHistogram, fFitFunction); };

		/**
		 * Performs the fit.
		 * @param hist The histogram (TH1D).
		 * @param func The fit function.
		 * @return An error code. */
		int Fit(TH1D * hist, TF1 * func);

		/**
		 * Draw the fit in the current pad. */
		void DrawFit(const char * options = "", bool flaglegend = false);

		/**
		 * Print a summary of the fit to the screen
		 */
		void PrintFitSummary();

		/**
		 * Calculate the p-value using fast-MCMC.
		 * @param par A set of parameter values 
		 * @param  pvalue The pvalue
		 * @return An error code 
		 */ 
		int CalculatePValueFast(std::vector<double> par, double &pvalue); 

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
};

// ---------------------------------------------------------

#endif
