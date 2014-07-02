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

typedef double (*ObservableFunction)(std::vector<double> const &);

class BCObservable : public BCVariable {

public:

	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * The default constructor. */
	BCObservable();

	/**
	 * Function-pointer constructor.
	 * @param name The name of the variable.
	 * @param lowerlimit The lower limit of the variable values.
	 * @param upperlimit The upper limit of the variable values.
	 * @param obs Pointer to double which stores value to be plotted (the value must be set by model during calculation of likelihood).
	 * @param latexname The latex name of the variable used in axis labeling.
	 */
	BCObservable(const char* name, double lowerlimit, double upperlimit, double * obs, const char* latexname = "");

	/**
	 * Function-pointer constructor.
	 * @param name The name of the variable.
	 * @param lowerlimit The lower limit of the variable values.
	 * @param upperlimit The upper limit of the variable values.
	 * @param fn Pointer to function returning double evaluating on parameter set std::vector<double>.
	 * @param latexname The latex name of the variable used in axis labeling.
	 */
	BCObservable(const char* name, double lowerlimit, double upperlimit, ObservableFunction fn, const char* latexname = "");

	/**
	 * Destructor */
	virtual ~BCObservable();


	/** \name Member functions (get) */
	/** @{ */

	/**
	 * Evaluate function over parameter set. Can be overloaded by user
	 * to make more complicated calculations.
	 */
	virtual double Evaluate(std::vector<double> const & parameters);

	/** @} */
	

	/** \name Member functions (set) */
	/** @{ */

	/** @} */
	

	/** \name Member functions (miscellaneous methods) */
	/** @{ */

	/** @} */

private:

	ObservableFunction fFunctionPointer;

	double * fObservablePointer;

};
#endif
