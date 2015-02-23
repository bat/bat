#ifndef __BCPARAMETERSET__H
#define __BCPARAMETERSET__H

/**
 * @class BCParameterSet Wrapper to allow access by name into list of BCParameter.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Parameters are not owned, and will not be deleted by BCParameterSet.
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>
#include <string>

#include "BCVariableSet.h"

// ---------------------------------------------------------

class BCParameterSet : public BCVariableSet {
public:

	/**
	 * Constructor */
	BCParameterSet() : BCVariableSet() 
	{}

	/*
	 * Destructor */
	virtual ~BCParameterSet()
	{}

	/*
	 * Assignment operator. */
	BCParameterSet & operator=(const BCParameterSet & rhs);

	/**
	 * @return The number of fixed parameters. */
	virtual unsigned int GetNFixedParameters() const;

	/**
	 * @return The number of free parameters. */
	virtual unsigned int GetNFreeParameters() const
	{ return Size() - GetNFixedParameters(); }

	/**
	 * @return volume of the set. */
	virtual double Volume() const ;

	/**
	 * Check whether all parameters have factorized priors set.
	 * @param ignore_fixed Whether to ignore fixed parameters.
	 * @return Whether all parameters have factorized priors set. */
	virtual bool ArePriorsSet(bool ignore_fixed=true) const;

	using BCVariableSet::IsWithinLimits;

	/**
	 * Check if vector of values is within limits.
	 * @param x Values to check
	 * @param ignore_fixed Flag for ignoring fixed values.
	 * @return Whether values are within limits of variables. */
	virtual bool IsWithinLimits(const std::vector<double> & x, bool ignore_fixed) const;

	/**
	 * Check if vector of values is at fixed values
	 * @param x Values to check. 
	 * @param ignore_unfixed Ignore unfixed parameters.
	 * @return Whether values are at fixed values. */
	virtual bool IsAtFixedValues(const std::vector<double> & x, bool ignore_unfixed=true) const;

	using BCVariableSet::ValueFromPositionInRange;

	/**
	 * Translate from unit interval to values in variable ranges.
	 * @param p vector of positions in the unit interval (0 = lower limit, 1 = upper limit).
	 * @param fix If true, fix fixed parameters to fixed value regardless of value in argument p.*/
	virtual void ValueFromPositionInRange(std::vector<double> &p, bool fix) const;

	using BCVariableSet::GetRangeCenters;
	
	/**
	 * Get range centers.
	 * @param fix If true, return fixed value for fixed parameters.
	 * @return vector fo range centers. */
	virtual std::vector<double> GetRangeCenters(bool fix) const;

	using BCVariableSet::GetUniformRandomValues;
	
	/**
	 * Get vector of values uniformly distributed in parameter ranges.
	 * @param fix If true, return fixed value for fixed parameters.
	 * @return vector of uniformly distributed random values. */
	virtual std::vector<double> GetUniformRandomValues(TRandom * const R, bool fix) const;

	/**
	 * Get vector values distributed randomly by the parameter priors.
	 * Parameters with unset priors will have infinite values.
	 * One should first call BCParameterSet::ArePriorsSet to be safe.
	 * @param R Random number generator to use.
	 * @param fix Whether to fix fixed parameters to their fixed values.
	 * @return vector of random values distributed according to priors. */
	virtual std::vector<double> GetRandomValuesAccordingToPriors(TRandom * const R, bool fix) const;

	/**
	 * Get random values distributed according to normal distributions
	 * with means of prior distributions and
	 * standard deviations of prior distributions multiplied by the expansion_factor.
	 * @param R Random number generator to use.
	 * @param fix Whether to fix fixed parameters to their fixed values.
	 * @param expansion_factor Constant to multiple standard deviations by.
	 * @param N Maximum number of tries to make to generate value within each parameter range.
	 * @param over_range Flag for whether to calculate means and std dev's only in parameter ranges
	 * @return vector of random values of normal distribution approximations to priors. */
	virtual std::vector<double> GetRandomValuesAccordingToGaussiansOfPriors(TRandom * const R, bool fix, double expansion_factor=1., unsigned N=1000000, bool over_range=true) const;

	/**
	 * Set all priors to constant. */
	virtual void SetPriorConstantAll();

	/**
	 * Get log of prior;
	 * assumes independent priors given for all parameters in set.
	 * @param parameters vector of parameters to return prior at.
	 * @return log of prior at parameter set value. */
	virtual double GetLogPrior(const std::vector<double> & parameters) const;

	/**
	 * Get vector of fixed values.
	 * @param include_unfixed Flag for whether to return fixed values (true) or infinity (false) for unfixed parameters.
	 * @return vector of fixed values for all parameters. */
	virtual std::vector<double> GetFixedValues(bool include_unfixed=true) const;
	
	/**
	 * Change values to fixed values for fixed parameters.
	 * @param x Vector of parameter values to adjust. */
	virtual bool ApplyFixedValues(std::vector<double> &x) const;

};
#endif
