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
#include <limits>
#include <vector>

#include "BCVariable.h"

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
     * Copy constructor.
     * PLEASE NOTE: if flag is true, the pointer for the internal value of the observable
     * will be shared by the copy. A change of value to any copy propagates to all others!
     * @param share_pointer Flag for sharing pointer of value to new observable */
    BCObservable(const BCObservable& other, bool share_pointer = false);

    /**
     * Function-pointer constructor.
     * @param name The name of the variable.
     * @param lowerlimit The lower limit of the variable values.
     * @param upperlimit The upper limit of the variable values.
     * @param obs Pointer to double which stores value to be plotted (the value must be set by model during calculation of likelihood).
     * @param latexname The latex name of the variable used in axis labeling.
     * @param unitstring Unit string to be printed for variable. */
    BCObservable(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname = "", const std::string& unitstring = "");

    /**
     * Destructor */
    virtual ~BCObservable();

    /** @} */

    /** \name operators and swap */
    /** @{ */

    /** Copy operator (Creates new pointer for value)*/
    BCObservable& operator=(BCObservable other);

    /** swap */
    friend void swap(BCObservable& A, BCObservable& B);

    /** assignment to a double operator */
    BCObservable& operator=(const double& value)
    { Value(value); return *this; }

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return Value of the observable. */
    virtual double Value()
    { return (fObservableValue) ? *fObservableValue : std::numeric_limits<double>::quiet_NaN(); }

    /** @} */


    /** \name Member functions (set) */
    /** @{ */

    /**
     * Set value of observable. */
    virtual void Value(double val)
    { if (fObservableValue) *fObservableValue = val; }

    /**
     * Set value location. */
    virtual void SetValueLocation(double* location)
    { fObservableValue = location; }

    /** @} */

private:

    double* fObservableValue;

};
#endif
