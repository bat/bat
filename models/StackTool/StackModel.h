#ifndef __STACKMODEL__H
#define __STACKMODEL__H

/*!
 * \class StackModel
 * This class can be used for fitting several template
 * histograms to a data histogram. The templates are assumed to have
 * no statistical uncertainty whereas the data are assumed to have
 * Poissonian fluctuations in each bin. Several methods to judge the
 * validity of the model are available.
 * \brief A class for fitting several templates to a data set.
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 10.2009
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <BAT/BCModel.h>

#include <TH1D.h>
#include <TH2D.h>

// ---------------------------------------------------------

class StackModel : public BCModel
{
 public:

	/**
	 * The default constructor.
	 */
	StackModel();

	/**
	 * A constructor.
	 */
	StackModel(const char * name);

	/**
	 * The default destructor.
	 */
	~StackModel();

	/**
	 * Calculates and returns the log of the prior probability at a
	 * given point in parameter space.
	 */
	double LogAPrioriProbability(std::vector <double> parameters);

	/**
	 * Calculates and returns the log of the Likelihood at a given point
	 * in parameter space.
	 */
	double LogLikelihood(std::vector <double> parameters);

	/**
	 * Returns the number of degrees of freedom of the model. These
	 * are the number of bins - the number of templates.
	 */
	int GetNDF()
		{ return fDataHistogram.GetNbinsX() - GetNParameters(); };

	/**
	 * Set the histogram containing the data.
	 */
	int SetDataHistogram(TH1D hist);

	/**
	 * Add a template histogram. The templates do not have to be
	 * normalized. The histogram has to have the same number of bins
	 * and cover the same region as the data histogram.
	 * @param hist The template histogram
	 * @param name The name of the template
	 * @param Nmin The lower limit of the normalization.
	 * @param Nmax The upper limit of the normalization.
	 */
	int AddTemplateHistogram(TH1D hist, const char * name="", double Nmin=0, double Nmax=0);

	/**
	 * Checks if two histograms have the same properties.
	 * @return 0 if not, 1 if they have the same properties.
	 */
	int CompareHistogramProperties(TH1D hist1, TH1D hist2);

	/**
	 * Prints the stack of templates scaled with the global mode. The
	 * data is plotted on top. The following options are available:\n\n
	 * "L"  : adds a legend\n\n
	 * "E0" : symmetric errorbars on the data points with a length of sqrt(obs).\n\n
	 * "E1" : symmetric errorbars on the sum of templates with a length of sqrt(exp).\n\n
	 * "E2" : asymmetric errorbars on the sum of templates from error
	 * propagation of the parameters.\n\n
	 * "E3" : asymmetric errorbars on the data points. The errorbars
	 * mark the 16% and 85% quantiles of the Poisson distribution
	 * rounded to the lower or upper integer, respectively.
	 * @param filename The name of the file the output is printed to.
	 * @param options Plotting options
	 */
	void PrintStack(const char * filename, const char * options="");

	/**
	 * Calculates and returns the chi2 value. The chi2 is calculated using
	 * the expectation value for the uncertainty.
	 */
	double CalculateChi2();

	/**
	 * Calculates and returns the chi2-probability.
	 */
	double CalculateChi2Prob();

	/**
	 * Calculates and returns the Likelihood at the global mode.
	 */
	double CalculateMaxLike();

	/**
	 * Calculates and returns the p-value for the global mode.
	 */
	double CalculatePValue();

	/**
	 * Calculates and returns the Kolmogorov-Smirnov-probability.
	 */
	double CalculateKSProb();

	/**
	 * Overloaded from BCIntegrate.
	 */
	void MCMCUserIterationInterface();

	/**
	 * Temporary entry. Used for debugging.
	 */
	void PrintTemp();

 protected:

	/**
	 * The data histogram.
	 */
	TH1D fDataHistogram;

	/**
	 * A container for the template histograms.
	 */
	std::vector <TH1D> fTemplateHistogramContainer;

	/**
	 * A container for the template names.
	 */
	std::vector <std::string> fTemplateNameContainer;

	/**
	 * A 2-d histogram for calculating the error bars
	 */
	TH2D * fUncertaintyHistogramExp;

	/**
	 * A 2-d histogram for calculating the error bars
	 */
	TH2D * fUncertaintyHistogramObsPosterior;

};

// ---------------------------------------------------------

#endif

