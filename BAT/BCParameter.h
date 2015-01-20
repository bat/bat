#ifndef __BCPARAMETER__H
#define __BCPARAMETER__H

/*!
 * \class BCParameter
 * \brief A class representing a parameter of a model.
 * \author Daniel Greenwald
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a parameter of a model. It contains
 * information about the name and the range of the parameter.
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCVariable.h"

#include <TNamed.h>
#include <TH1.h>
#include <TF1.h>

#include <limits>

class TRandom;

// ---------------------------------------------------------

class BCParameter : public BCVariable {

public:

	/** \name Enumerators */
	/** @{ */

	enum BCPriorType {
		kPriorUnset         = -1,
		kPriorConstant      =  0,
		kPriorTF1           =  1,
		kPriorTH1           =  2,
		kPriorGaussian      =  3,
		kPriorSplitGaussian =  4
	};

	/** @} */

	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * The default constructor. */
	BCParameter();

	/**
	 * Copy constructor. */
	BCParameter(const BCParameter & other);

	/**
	 * A constructor.
	 * @param name The name of the parameter.
	 * @param lowerlimit The lower limit of the parameter values.
	 * @param upperlimit The upper limit of the parameter values.
	 * @param latexname The latex name of the parameter used in axis labeling.
	 */
	BCParameter(const char* name, double lowerlimit, double upperlimit, const char* latexname = "");

	/**
	 * Destructor. */
	virtual ~BCParameter();

	/** \name Member functions (get) */
	/** @{ */

	/**
	 * @return Whether parameter is fixed to a value. */
	bool Fixed() const
	{ return fFixed; }

	/**
	 * @return vector of parameters used for priors. */
	std::vector<double> & GetPriorParameters()
	{ return fPriorParameters; }

	/**
	 * @return Value parameter may be fixed to. */
	double GetFixedValue() const
	{ return fFixedValue; }

	/**
	 * @return type of prior for parameter. */
	BCParameter::BCPriorType GetPriorType()
	{ return fPriorType; }

	/**
	 * @return prior container object. */
	TNamed * GetPriorContainer() const
	{ return fPriorContainer; }

	/**
	 * @return prior container object as TF1. */
	TF1 * GetPriorTF1() const
	{ return (fPriorContainer) ? dynamic_cast<TF1*>(fPriorContainer) : 0 ; }

	/**
	 * @return prior container object as TH1. */
	TH1 * GetPriorTH1() const
	{ return (fPriorContainer) ? dynamic_cast<TH1*>(fPriorContainer) : 0 ; }
 
	/**
	 * @return whether to interpolate prior (if TH1). */
	bool GetInterpolatePrior()
	{ return fInterpolatePrior; }

	/**
	 * Get log of value of prior at parameter value.
	 * @param value of parameter to return prior of.
	 * @return log of prior value at parameter value. */
	virtual double GetLogPrior(double x) const;

	/**
	 * @return a random value distributed according to the prior.
	 * @param rnd Pointer to the random generator to be used, if needed.
	 * @param N Maximum number of tries to make to generate value within parameter range. */
	 virtual double GetRandomValueAccordingToPrior(TRandom * const R, unsigned N=1000000) const;

	/**
	 * Get exact or estimated mean of prior.
	 * @return mean of prior. */
	virtual double GetPriorMean() const;

	/**
	 * Get exact or estimated variance of prior.
	 * @return variance of prior. */
	virtual double GetPriorVariance() const;

	/**
	 * Get exact or estimated standard deviation of prior.
	 * @return standard deviation of prior. */
	virtual double GetPriorStandardDeviation() const;

	/**
	 * Get random value distributed according to normal distribution
	 * with mean of prior distribution and
	 * standard deviation of prior distribution multiplied by the expansion_factor.
	 * @param R Random number generator to use.
	 * @param expansion_factor Constant to multiple standard deviation by.
	 * @param N Maximum number of tries to make to generate value within parameter range.
	 * @return random value of normal distribution approximation to prior. */
	virtual double GetRandomValueAccordingToGaussianOfPrior(TRandom * const R, double expansion_factor=1., unsigned N=1000000) const;

	/** @} */
	
	/** \name Member functions (set) */
	/** @{ */
	
	/**
	 * Fix parameter to value (set prior to delta).
	 * @param value value to fix parameter to. */
	bool Fix(double value)
	{	fFixed = true; fFixedValue = value; return true;}

	/**
	 * Unfix parameter. */
	bool Unfix()
	{ fFixed = false; return true;}

	/**
	 * Set prior to constant;
	 * @return success of action. */
	bool SetPriorConstant();
	
	/**
	 * Set prior to a ROOT TF1.
	 * @param f pointer to ROOT TF1.
	 * @return success of action. */
	bool SetPrior(const TF1 * const f);

	/**
	 * Set prior to a ROOT TH1.
	 * @param h pointer to ROOT TH1.
	 * @param interpolate flag for controlling whether to interpolate histogram for prior.
	 * @return success of action. */
	bool SetPrior(const TH1 * const h, bool interpolate=false);

	/**
	 * Copy prior from other parameter.
	 * @param other Other parameter.
	 * @return success of action. */
	bool CopyPrior(const BCParameter & other);

	/**
	 * Set prior to ROOT TF1 holding a normal distribution.
	 * @param mean mean of normal distribution.
	 * @param sigma standard deviation of normal distribution.
	 * @return success of action. */
	bool SetPriorGauss(double mean, double sigma);

	/**
	 * Set prior to ROOT TF1 holding a normal distribution with different width
	 * for above and below the mean.
	 * @param mean mean of normal distribution.
	 * @param sigma_below standard deviation of normal distribution above mean.
	 * @param sigma_above standard deviation of normal distribution below mean.
	 * @return success of action. */
	bool SetPriorGauss(double mean, double sigma_below, double sigma_above);

	/** @} */
	
	/** \name Member functions (miscellaneous methods) */
	/** @{ */

	std::string OneLineSummary() const;
	
	/** @} */

private:
	/// Flag to fix parameter; useful for example, for integration.
	bool fFixed;

	/// The fixed value of the parameter.
	double fFixedValue;

	/// type of prior
	BCParameter::BCPriorType fPriorType;
	
	/// parameters needed for built in prior types (Flat, Gaus, etc)
	std::vector<double> fPriorParameters;

	/// ROOT prior object (function, histogram, graph, etc.)
	TNamed * fPriorContainer;

	/// Flag for controlling interpolation of prior
	bool fInterpolatePrior;

};
#endif
