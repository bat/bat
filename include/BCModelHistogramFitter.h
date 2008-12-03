/*!
 * \class BCModelHistogramFitter
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

class BCModelHistogramFitter : public BCModel
{

 public:

	/** \name Constructors and destructors */ 
	/* @{ */ 

	/** 
	 * The default constructor. 
	 */ 
	BCModelHistogramFitter(); 

	/** 
	 * A constructor. 
	 * @param hist The histogram (TH1D). 
	 * @param func The fit function. 
	 */ 
	BCModelHistogramFitter(TH1D * hist, TF1 * func); 

	/**
	 * The default destructor. 
	 */ 
	~BCModelHistogramFitter(); 

	/* @} */ 

	/** \name Member functions (get) */ 
	/* @{ */ 

	/**
	 * @return The histogram 
	 */ 
	TH1D * GetHistogram()
	{ return fHistogram; }; 

	/**
	 * @return The fit function
	 */ 
	TF1 * GetFitFunction()
		{ return fFitFunction; }; 

	/* @} */ 

	/** \name Member functions (set) */ 
	/* @{ */ 

	/** 
	 * @param hist The histogram 
	 */ 
	void SetHistogram(TH1D * hist); 
	/**
	 * @param func The fit function
	 */ 
	void SetFitFunction(TF1 * func);

	/* @} */ 

	/** \name Member functions (miscellaneous methods) */ 
	/* @{ */ 

	double LogAPrioriProbability(std::vector <double> parameters); 
  
  double LogLikelihood(std::vector <double> parameters); 

	/**
	 * Plots the histogram 
	 * @param options Options for plotting.
	 * @param filename Name of the file which the histogram is printed into. 
	 * The following options are available:\n
	 * F : plots the fit function on top of the data 
	 * E0 : plots the fit function and the 68% prob. uncertainty band of the fit function on top of the data
	 * E1 : plots the expectation from the fit function and the uncertainty bin-by-bin as error bars.
	 */ 
	void PrintHistogram(const char * options = "", const char * filename = ""); 

	double FitFunction(std::vector <double> x, std::vector <double> parameters); 

	/* @} */ 

 private: 

	/**
	 * The histogram containing the data. 
	 */
	TH1D * fHistogram; 

	/**
	 * The fit function 
	 */ 
	TF1 * fFitFunction; 

}; 

// --------------------------------------------------------- 

#endif 
