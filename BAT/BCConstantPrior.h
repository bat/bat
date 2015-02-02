#ifndef __BCCONSTANTPRIOR__H
#define __BCCONSTANTPRIOR__H

/*!
 * \class BCConstantPrior
 * \brief A class to represent a constant prior of a parameter
 * \author Daniel Greenwald
 * \version 1.0
 * \date 01.2015
 */


/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <limits>

class TF1;

// ---------------------------------------------------------

/** \class BCConstantPrior **/
class BCConstantPrior : public BCPrior {
public:

	/** Constructor */
	BCConstantPrior()
	{ }

	/** Copy constructor */
	BCConstantPrior(const BCConstantPrior & other)
	{ }

	/** Destructor */
	virtual ~BCConstantPrior()
	{ }

	/** Get as TF1 for drawing purposes.
	 * Parameter zero is 1/(range width).
	 * Parameter one is integral ( = 1)
	 * @param xmin lower limit of range for TF1
	 * @param xmax upper limit of range for TF1
	 * @param normalize whether to normalize TF1 over range*/
	virtual TF1 * GetAsTF1(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity(), bool normalize=true) const;

	/** @return constant log(prior) = 0 */
	virtual double GetLogPrior(double x) const
	{ return 0; }

	/** @return raw moments of uniform continuous distribtion. */
	virtual double GetRawMoment(unsigned n, double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const;

	/** @return integral = 1 */
	virtual double GetIntegral(double xmin=-std::numeric_limits<double>::infinity(), double xmax=std::numeric_limits<double>::infinity()) const
	{ return 1; }

};
