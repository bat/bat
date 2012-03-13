#ifndef __BINLIKEMODEL__H
#define __BINLIKEMODEL__H

/*!
 * \class BinLikeModel
 * A class for fitting several histograms with one function. The
 * function depends on one particular "leading parameter", e.g. the a
 * mass parameter when parameterizing templates.
 * \brief A class for fitting several histograms with one function 
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 10.2009
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <BAT/BCModel.h>

#include <TH1D.h>
#include <TH2D.h>

// ---------------------------------------------------------

class BinLikeModel : public BCModel
{
 public:

	/**
	 * The default constructor.
	 */
	BinLikeModel();

	/**
	 * The default destructor.
	 */
	~BinLikeModel();

	/**
	 * Defines the parameters of the fit. */ 
	virtual void DefineParameters() = 0; 

	/**
	 * Calculates the expectation value. 
	 * @param parameters The current parameter values.
	 * @param parvalue the leading parameter value. 
	 * @param x the x-value 
	 * @return The expectation value. */ 
	virtual double Expectation(std::vector<double> parameters, double parvalue, double x) = 0; 

	/**
	 * Calculates and returns the log of the prior probability at a
	 * given point in parameter space.
	 */
	virtual double LogAPrioriProbability(std::vector<double> parameters);

	/**
	 * Calculates and returns the log of the Likelihood at a given point
	 * in parameter space.
	 */
	virtual double LogLikelihood(std::vector<double> parameters);

	/**
	 * Add a histogram for a certain value of the leading parameter. 
	 * @param hist a pointer to the histogram
	 * @param parvalue The value of the leading parameter corresponding
	 * to the histogram.
	 * @param weight The weight of the histogram. 
	 */
	int AddHistogram(TH1D* hist, double parvalue, double weight = 1);

	/**
	 * Return the number of histograms.
	 * @return the number of histograms. */ 
	int GetNHistograms()
	{ return int(fHistogramContainer.size()); }; 

	/**
	 * Return the number of bins in a histogram. 
	 * @param index the index of the histogram. 
	 * @return the number of bins. */ 
	int GetNBins(int index); 

	/**
	 * Return a histogram. 
	 * @param index the index of the histogram. 
	 * @return the histogram. */
	TH1D* GetHistogram(int index) 
	{ return fHistogramContainer.at(index); }; 

	/** 
	 * Return the value of the leading parameter. 
	 * @param index the index of the histogram
	 * @return the value of the leading parameter. */ 
	double GetLeadParValue(int index); 

	/**
	 * Print all histograms. */ 
	void PrintHistograms(const char* path = "./"); 

	/**
	 * Print summary. */
	void PrintChi2Summary(); 

	/**
	 * Overloaded from BCIntegrate.
	 */
	void MCMCUserIterationInterface();

	/**
	 * Returns the number of degrees of freedom of the model. These
	 * are the number of bins - the number of parameters.
	 */
	int GetNDF(); 

	/**
	 * Calculates and returns the chi2 value for all histograms. The
	 * chi2 is calculated using the expectation value for the
	 * uncertainty.  
	 * @return the chi2-value.
	 */
	double CalculateChi2(); 

	/**
	 * Calculates and returns the chi2 value for one histogram. The chi2
	 * is calculated using the expectation value for the uncertainty.
	 * @param index the index of the histogram.
	 * @return the chi2-value.
	 */
	double CalculateChi2Hist(int index); 

	/**
	 * Calculates and returns the chi2-probability.
	 */
	double CalculateChi2Prob();

 protected:

	/**
	 * A container for the template histograms.
	 */
	std::vector<TH1D*> fHistogramContainer;

	/**
	 * A container for the leading parameter values. 
	 */
	std::vector<double> fLeadParContainer;

	/**
	 * A container for the histogram weights. 
	 */
	std::vector<double> fHistWeightContainer;

	/**
	 * A container of 2-d histograms for calculating the error bars
	 */
	std::vector<TH2D*> fErrHistContainer;
};

// ---------------------------------------------------------

#endif

