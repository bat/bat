#ifndef __BCOBSERVABLE__H
#define __BCOBSERVABLE__H

/**
 * @class BCObservable
 * @brief A class representing a variable of a model.
 * @author Daniel Greenwald
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class represents a variable of a model. It contains
 * information about the name and the range of the variable.
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include <string>

#include "BCVariable.h"
#include "BCVariableSet.h"

// ---------------------------------------------------------

class BCObservable : public BCVariable
{

public:

    /** \name Constructors and destructor */
    /** @{ */

    /**
     * The default constructor. */
    BCObservable();

    /**
     * Constructor with names and limits.
     * @param name The name of the variable.
     * @param lowerlimit The lower limit of the variable values.
     * @param upperlimit The upper limit of the variable values.
     * @param latexname The latex name of the variable used in axis labeling.
     * @param unitstring Unit string to be printed for variable. */
    BCObservable(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname = "", const std::string& unitstring = "");

    /** @} */

    /** \name operators and swap */
    /** @{ */

    /** assignment to a double operator */
    BCObservable& operator=(const double& value)
    { Value(value); return *this; }

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return Value of the observable. */
    double Value() const
    { return fObservableValue; }

    /** @} */


    /** \name Member functions (set) */
    /** @{ */

    /**
     * Set value of observable. */
    virtual void Value(double val)
    { fObservableValue = val; }

    /** @} */

private:

    double fObservableValue;

};

typedef BCVariableSet<BCObservable> BCObservableSet;

#endif
