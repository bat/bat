#ifndef __BCDATAPOINT__H
#define __BCDATAPOINT__H

/**
 * @class BCDataPoint
 * @brief A class representing a data point.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class represents a data point which is the basic unit
 * of information. A data point can be an event, a bin content, etc.
 * Each data point can store several variables of type  double.
 * The variables are organized in a vector.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <string>
#include <vector>

#include "BCLog.h"

// ---------------------------------------------------------

class BCDataPoint
{
public:

    /** \name Constructors and destructor */
    /** @{ */

    /**
     * A constructor.
     * @param nvariables The number of variables stored in a data object. */
    BCDataPoint(int nvariables = 0);

    /**
     * A constructor.
     * @param x The vector containing the data. */
    BCDataPoint(const std::vector<double>& x);
    /** @} */

    /** \name operators */
    /** @{ */

    /**
     * Raw and fast access. */
    double& operator[](unsigned index)
    {	return fData[index]; }

    /**
     * Raw and fast access. */
    const double& operator[](unsigned index) const
    {	return fData[index]; }

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * Safer access, but slower
     * @param index The index of the variable.
     * @return The value of the variable. */
    double GetValue(unsigned index) const
    { return fData.at(index); }

    /**
     * @return A vector of values. */
    std::vector<double>& GetValues()
    { return fData; };

    /**
     * @return A vector of values. */
    const std::vector<double>& GetValues() const
    { return fData; };

    /**
     * Returns the number of values. */
    unsigned int GetNValues() const
    { return fData.size(); };

    /** @} */

    /** \name Member functions (set) */
    /** @{ */

    /**
     * Safer, but slower, value setting of a variable.
     * @param index The index of the variable
     * @param value The value of the variable */
    void SetValue(unsigned index, double value)
    { fData.at(index) = value; }

    /**
     * Set the values of all variables.
     * @param values A vector of values */
    void SetValues(const std::vector<double>& values);

    /**
     * Set the number of variables.
     * Use with caution!
     * @param n New dimensionality to set for data point.
     * @param val Value to fill into new data values if enlarging data point size. */
    void SetNValues(unsigned n, double val = 0.)
    { fData.resize(n, val); }

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Print summary of data point to the string handler.
     * @param output String handler (default = BCLog::OutSummary. */
    void PrintSummary(void (*output)(const std::string&) = BCLog::OutSummary) const;

    /** @} */

private:

    /**
     * The vector containing the values of the variables. */
    std::vector<double> fData;

};

// ---------------------------------------------------------

#endif
