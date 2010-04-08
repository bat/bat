// ***************************************************************
// This file was created using the ./CreateFitModel.sh script
// ./CreateFitModel.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __TEMPLATEFIT__H
#define __TEMPLATEFIT__H

#include <../BinLikeModel.h>

// This is a TemplateFit header file.
// Model source code is located in file TemplateFit/TemplateFit.cxx

// ---------------------------------------------------------
class TemplateFit : public BinLikeModel
{
 public:
	
  /**
   * The default constructor.
   */
  TemplateFit();

  /**
   * The default destructor.
   */
  ~TemplateFit();

  /**
   * Defines the parameters of the fit. */ 
  void DefineParameters(); 

	/**
	 * Calculates and returns the log of the prior probability at a
	 * given point in parameter space.
	 */
	double LogAPrioriProbability(std::vector <double> parameters); 

  /**
   * Calculates the expectation value. 
   * @param parameters The current parameter values.
   * @param parvalue the leading parameter value. 
   * @param x the x-value
   * @return The expectation value. */ 
  double Expectation(std::vector <double> parameters, double parvalue, double x); 
};
// ---------------------------------------------------------

#endif

