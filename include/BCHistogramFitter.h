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

#ifndef __BCMODELHISTOGRAMFITTER__H
#define __BCMODELHISTOGRAMFITTER__H

#include "BCModel.h"

#include <TH1D.h>
#include <TF1.h>

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

		/* @} */

		/** \name Member functions (set) */
		/* @{ */

		/**
		 * @param hist The histogram */
		void SetHistogram(TH1D * hist);

		/**
		 * @param func The fit function */
		void SetFitFunction(TF1 * func);

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
		void PrintHistogram(const char * options = "", const char * filename = "");

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
		void DrawFit(const char * options = "");

		/* @} */

	private:

		/**
		 * The histogram containing the data.
		 */
		TH1D * fHistogram;

		/**
		 * The fit function */
		TF1 * fFitFunction;

};

// ---------------------------------------------------------

#endif
