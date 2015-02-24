#ifndef __BCOBSERVABLE__H
#define __BCOBSERVABLE__H

/*!
 * \class BCObservable
 * \brief A class representing a variable of a model.
 * \author Daniel Greenwald
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a variable of a model. It contains
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
#include <vector>
#include "BCVariable.h"

// ---------------------------------------------------------

class BCObservable : public BCVariable {

public:

	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * The default constructor. */
	BCObservable();

	/**
	 * Copy constructor.
	 * PLEASE NOTE: the pointer for the internal value of the observable
	 * will be shared by the copy. A change of value to any copy propegates to all others!*/
	BCObservable(const BCObservable & other);

	/**
	 * Function-pointer constructor.
	 * @param name The name of the variable.
	 * @param lowerlimit The lower limit of the variable values.
	 * @param upperlimit The upper limit of the variable values.
	 * @param obs Pointer to double which stores value to be plotted (the value must be set by model during calculation of likelihood).
	 * @param latexname The latex name of the variable used in axis labeling.
	 * @param unitstring Unit string to be printed for variable. */
	BCObservable(const char* name, double lowerlimit, double upperlimit, const char* latexname = "", const char * unitstring="");

	/**
	 * Destructor */
	virtual ~BCObservable();


	/** \name Member functions (get) */
	/** @{ */

	/**
	 * @return Value of the observable. */
	virtual double Value()
	{ return *fObservableValue; }

	/** @} */
	

	/** \name Member functions (set) */
	/** @{ */
	
	/**
	 * Set value of observable. */
	virtual void Value(double val)
	{ *fObservableValue = val; }

	/** @} */
	

	/** \name Member functions (miscellaneous methods) */
	/** @{ */

	/** @} */

private:

	double * fObservableValue;
	
};
#endif
