#ifndef __BCPARAMETERSET__H
#define __BCPARAMETERSET__H

/**
 * @class BCParameterSet
 * @brief Wrapper to allow access by name into list of BCParameter.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Parameters are not owned, and will not be deleted by BCParameterSet.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>
#include <string>

#include "BCParameter.h"
#include "BCVariableSet.h"

// ---------------------------------------------------------

class BCParameterSet : public BCVariableSet<BCParameter>
{
public:

    /**
     * Constructor */
    BCParameterSet();

    /**
     * Destructor */
    virtual ~BCParameterSet() {};

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
    virtual bool ArePriorsSet(bool ignore_fixed = true) const;

    /**
     * Check if vector of values is within limits. Check if fixed parameters are at fixed values.
     * @param x Values to check
     * @return Whether values are within limits of variables. */
    virtual bool IsWithinLimits(const std::vector<double>& x) const;

    /**
     * Check if fixed parameters in vector of values are at fixed values
     * @param x Values to check.
     * @return Whether values are at fixed values. */
    virtual bool IsAtFixedValues(const std::vector<double>& x) const;

    /**
     * Translate from unit interval to values in variable ranges, fixing fixed parameters along the way.
     * @param p vector of positions in the unit interval (0 = lower limit, 1 = upper limit). */
    virtual void ValueFromPositionInRange(std::vector<double>& p) const;

    /**
     * Get range centers, leaving fixed parameters at fixed values
     * @return vector of range centers & fixed values. */
    virtual std::vector<double> GetRangeCenters() const;

    /**
     * Get vector of values uniformly distributed in parameter ranges (or at fixed values, if fixed)
     * @return vector of uniformly distributed random values. */
    virtual std::vector<double> GetUniformRandomValues(TRandom* const R) const;

    /**
     * Get vector values distributed randomly by the parameter priors.
     * Parameters with unset priors will have infinite values.
     * Fixed parameters will be at fixed values.
     * One should first call BCParameterSet::ArePriorsSet to be safe.
     * @param R Random number generator to use.
     * @return vector of random values distributed according to priors. */
    virtual std::vector<double> GetRandomValuesAccordingToPriors(TRandom* const R) const;

    /**
     * Set all priors to constant. */
    virtual void SetPriorConstantAll();

    /**
     * Get log of prior;
     * assumes independent priors given for all parameters in set.
     * @param parameters vector of parameters to return prior at.
     * @return log of prior at parameter set value. */
    virtual double GetLogPrior(const std::vector<double>& parameters) const;

    /**
     * Get vector of fixed values.
     * @param include_unfixed Flag for whether to return fixed values (true) or infinity (false) for unfixed parameters.
     * @return vector of fixed values for all parameters. */
    virtual std::vector<double> GetFixedValues(bool include_unfixed = true) const;

    /**
     * Change values to fixed values for fixed parameters.
     * @param x Vector of parameter values to adjust. */
    virtual bool ApplyFixedValues(std::vector<double>& x) const;

};
#endif
