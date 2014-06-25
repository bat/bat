/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCObservable.h"

#include <string>


// ---------------------------------------------------------

BCObservable::BCObservable() : BCVariable(), fFunction(0)
{
	fPrefix = "Observable";
}

// ---------------------------------------------------------

BCObservable::BCObservable(const char * name, double lowerlimit, double upperlimit, ObservableFunction fn, const char * latexname) :
	BCVariable(name,lowerlimit,upperlimit,latexname)
	, fFunction(fn)
{
	fPrefix = "Observable";
}

// ---------------------------------------------------------

BCObservable::~BCObservable() {
}

// ---------------------------------------------------------

double BCObservable::Evaluate(std::vector<double> const & parameters) {
	return (fFunction) ? (*fFunction)(parameters) : 0;
}
