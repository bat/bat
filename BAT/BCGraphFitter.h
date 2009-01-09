#ifndef __BCMODELGRAPHFITTER__H
#define __BCMODELGRAPHFITTER__H

/*!
 * \class BCGraphFitter
 * \brief A class for fitting graphs with functions
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 2008
 * \detail This class allows fitting of a TGraphErrors using
 * a TF1 function. It doeasn't take the x uncertainties into account.
 * For that look at BCGraphXFitter (not yet implemented).
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>

#include "BAT/BCModel.h"

class TGraphErrors;
class TF1;

// ---------------------------------------------------------

class BCGraphFitter : public BCModel
{
	public:

		/** \name Constructors and destructors */
		/* @{ */

		/**
		 * Default constructor */
		BCGraphFitter();

		/**
		 * Constructor
		 * @param graph pointer to TGraphErrors
		 * @param func pointer to TF1 */
		BCGraphFitter(TGraphErrors * graph, TF1 * func);

		/**
		 * The default destructor. */
		~BCGraphFitter();

		/* @} */

		/** \name Member functions (get) */
		/* @{ */

		/**
		 * @return pointer to TGraphErrors */
		TGraphErrors * GetGraph()
			{ return fGraph; };

		/**
		 * @return pointer to TF1 */
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
		 * @param graph pointer to TGraphErrors object */
		int SetGraph(TGraphErrors * graph);

		/**
		 * @param func pointer to TF1 object */
		int SetFitFunction(TF1 * func);

		/* @} */
		/** \name Member functions (miscellaneous methods) */
		/* @{ */

		/**
		 * The log of the prior probability. It is set to be flat in all parameters.
		 * @param parameters vector containing the parameter values */
		double LogAPrioriProbability(std::vector <double> parameters);

		/**
		 * The log of the conditional probability.
		 * @param parameters vector containing the parameter values */
		double LogLikelihood(std::vector <double> parameters);

		/**
		 * Returns the value of the 1D fit function for a given set of parameters
		 * at a given x.
		 * @param x point to calculate the function value at
		 * @param parameters parameters of the function */
		double FitFunction(std::vector <double> x, std::vector <double> parameters);

		/**
		 * Performs the fit. The graph and the function has to beset beforehand.
		 * @return An error code. */
		int Fit()
			{ return this -> Fit(fGraph, fFitFunction); };

		/**
		 * Performs the fit of the graph with the function.
		 * @param graph pointer to TGraphErrors object
		 * @param func pointer to TF1 object
		 * @return An error code. */
		int Fit(TGraphErrors * graph, TF1 * func);

		/**
		 * Draw the fit in the current pad. */
		void DrawFit(const char * options = "", bool flaglegend = false);

		/* @} */

	private:

		/**
		 * The graph containing the data. */
		TGraphErrors * fGraph;

		/**
		 * The fit function */
		TF1 * fFitFunction;

		/**
		 * Pointer to the error band (for legend) */
		TGraph * fErrorBand;

		/**
		 * Pointer to a graph for displaying the fit function */
		TGraph * fGraphFitFunction;

};

// ---------------------------------------------------------

#endif
