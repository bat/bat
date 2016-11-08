#ifndef __BCPARAMETER__H
#define __BCPARAMETER__H

/**
 * @class BCParameter
 * @brief A class representing a parameter of a model.
 * @author Daniel Greenwald
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class represents a parameter of a model. It contains
 * information about the name and the range of the parameter.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCVariable.h"

#include <limits>

class TRandom;
class BCPrior;

// ---------------------------------------------------------

class BCParameter : public BCVariable
{

public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * The default constructor. */
    BCParameter();

    /**
     * Copy constructor. */
    BCParameter(const BCParameter& other);

    /**
     * A constructor.
     * @param name The name of the parameter.
     * @param lowerlimit The lower limit of the parameter values.
     * @param upperlimit The upper limit of the parameter values.
     * @param latexname The latex name of the parameter used in axis labeling.
     * @param unitstring Unit string to be printed for parameter. */
    BCParameter(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname = "", const std::string& unitstring = "");

    /**
     * Destructor. */
    virtual ~BCParameter();

    /** @} */

    /** \name operators and swap */
    /** @{ */

    /** Copy operator */
    BCParameter& operator=(BCParameter other);

    /** swap */
    friend void swap(BCParameter& A, BCParameter& B);

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return Whether parameter is fixed to a value. */
    virtual bool Fixed() const
    { return fFixed; }

    /**
     * @return Value parameter may be fixed to. */
    virtual double GetFixedValue() const
    { return fFixedValue; }

    /**
     * @return prior object*/
    virtual BCPrior* GetPrior()
    { return fPrior;}

    /**
     * @return prior object*/
    virtual const BCPrior* GetPrior() const
    { return fPrior;}

    /**
     * @return prior evaluated from prior object */
    virtual double GetPrior(double x) const
    { double lp = GetLogPrior(x); return (std::isfinite(lp)) ? exp(lp) : ((lp < 0) ? 0 : std::numeric_limits<double>::infinity()); }

    /**
     * Get log of value of prior at parameter value.
     * @param value of parameter to return prior of.
     * @return log of prior value at parameter value. */
    virtual double GetLogPrior(double x) const;

    /**
     * @return prior's mode in parameter range. (For absolute mode, get prior object.) */
    virtual double GetPriorMode() const;

    /**
     * @return prior's mean in parameter range. (For absolute mean, get prior object.) */
    virtual double GetPriorMean() const;

    /**
     * @return prior's variance in parameter range. (For absolute variance, get prior object.) */
    virtual double GetPriorVariance() const;

    /**
     * @return a random value distributed according to the prior.
     * @param rng Pointer to the random generator to be used, if needed. */
    virtual double GetRandomValueAccordingToPrior(TRandom* const rng) const;

    /** @} */

    /** \name Member functions (set) */
    /** @{ */

    /**
     * Set the limits of the parameter values.
     * @param lowerlimit The lower limit of the variable values.
     * @param upperlimit The upper limit of the variable values. */
    virtual void SetLimits(double lowerlimit = 0, double upperlimit = 1);

    /**
     * Fix parameter to value (set prior to delta).
     * @param value value to fix parameter to. */
    virtual bool Fix(double value)
    {	fFixed = true; fFixedValue = value; return true;}

    /**
     * Unfix parameter. */
    virtual bool Unfix()
    { fFixed = false; return true;}

    /**
     * Set prior. Parameter will own prior! */
    virtual void SetPrior(BCPrior* const prior);

    /**
     * Set constant prior. */
    virtual void SetPriorConstant();

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    std::string OneLineSummary(bool print_prefix = true, int name_length = -1) const;

    /** @} */

private:
    /// Flag to fix parameter; useful for example, for integration.
    bool fFixed;

    /// The fixed value of the parameter.
    double fFixedValue;

    /// prior
    BCPrior* fPrior;

};
#endif
