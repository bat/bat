#ifndef __BCUSEROBSERVABLE__H
#define __BCUSEROBSERVABLE__H

/*!
 * \class BCUserObservable
 * \brief A class representing a observable of a model.
 * \author Daniel Greenwald
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a observable of a model. It contains
 * information about the name and the range of the observable.
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
#include "BCObservable.h"

// ---------------------------------------------------------

typedef double (*ObservableFunction)(std::vector<double> const &);

class BCUserObservable : public BCObservable {

public:

	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * The default constructor. */
	BCUserObservable();

	/**
	 * A constructor.
	 * @param name The name of the observable.
	 * @param lowerlimit The lower limit of the observable values.
	 * @param upperlimit The upper limit of the observable values.
	 * @param latexname The latex name of the observable used in axis labeling.
	 */
	BCUserObservable(const char* name, double lowerlimit, double upperlimit, ObservableFunction * fn, const char* latexname = "");

	/** \name Member functions (get) */
	/** @{ */

	

	/** @} */
	
	/** \name Member functions (set) */
	/** @{ */
	

	/** @} */
	
	/** \name Member functions (miscellaneous methods) */
	/** @{ */

	/** @} */

private:

	ObservableFunction fFunction;

};
#endif
