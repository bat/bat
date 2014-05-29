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
#include "BCObservable.h"

// ---------------------------------------------------------

class BCParameter : public BCObservable {

public:

	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * The default constructor. */
	BCParameter();

	/**
	 * A constructor.
	 * @param name The name of the parameter.
	 * @param lowerlimit The lower limit of the parameter values.
	 * @param upperlimit The upper limit of the parameter values.
	 * @param latexname The latex name of the parameter used in axis labeling.
	 */
	BCParameter(const char* name, double lowerlimit, double upperlimit, const char* latexname = "");

	/** \name Member functions (get) */
	/** @{ */

	bool Fixed() const
	{ return fFixed; }

	double GetFixedValue() const
	{ return fFixedValue; }

	/** @} */
	
	/** \name Member functions (set) */
	/** @{ */
	
	void Fix(double value)
	{
		fFixed = true;
		fFixedValue = value;
	}

	void Unfix()
	{ fFixed = false; }
	/** @} */
	
	/** \name Member functions (miscellaneous methods) */
	/** @{ */

	/**
	 * Returns true if the value is at a parameter limit.
	 * @return flag States if value is at parameter limit. */
	bool IsAtLimit(double value) const;
	
	bool IsValid(double value) const
	{ return (fLowerLimit <= value) && (value <= fUpperLimit) ? true : false; }
	
	/** @} */

private:
	/// Flag to fix parameter; useful for example, for integration.
	bool fFixed;

	/// The fixed value of the parameter.
	double fFixedValue;

};
#endif
