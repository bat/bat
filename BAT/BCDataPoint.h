#ifndef __BCDATAPOINT__H
#define __BCDATAPOINT__H

/*!
 * \class BCDataPoint
 * \brief A class representing a data point.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a data point which is the basic unit
 * of information. A data point can be an event, a bin content, etc.
 * Each data point can store several variables of type  double.
 * The variables are organized in a vector.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>
#include "BCLog.h"

// ---------------------------------------------------------

class BCDataPoint
{
public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * A constructor.
     * @param nvariables The number of variables stored in a data object. */
    BCDataPoint(int nvariables = 0);

    /**
     * A constructor.
     * @param x The vector containing the data. */
    BCDataPoint(const std::vector<double>& x);

    /**
     * Copy constructor. */
    BCDataPoint(const BCDataPoint& other);

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * @param index The index of the variable.
     * @return The value of the variable. */
    double GetValue(unsigned index) const;

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
     * Set the value of a variable.
     * @param index The index of the variable
     * @param value The value of the variable */
    void SetValue(unsigned index, double value);

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
     * Copy from other data point.
     * @param other Data point to copy*/
    void Copy(const BCDataPoint& other)
    { fData = other.fData; }

    /**
     * Assignment operator. */
    BCDataPoint& operator=(const BCDataPoint& rhs)
    { Copy(rhs); return *this; }

    /**
     * Dump the data to the standard output */
    void Dump(void (*output)(const char*) = BCLog::OutSummary) const;

    /** @} */

private:

    /**
     * The vector containing the values of the variables. */
    std::vector<double> fData;

};

// ---------------------------------------------------------

#endif

