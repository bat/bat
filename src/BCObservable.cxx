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

BCObservable::BCObservable()
	: BCVariable()
	, fFunctionPointer(0)
	, fObservablePointer(0)
{
	fPrefix = "Observable";
}

// ---------------------------------------------------------

BCObservable::BCObservable(const char * name, double lowerlimit, double upperlimit, double * obs, const char * latexname)
	:	BCVariable(name,lowerlimit,upperlimit,latexname)
	, fFunctionPointer(0)
	, fObservablePointer(obs)
{
	fPrefix = "Observable";
}

// ---------------------------------------------------------

BCObservable::BCObservable(const char * name, double lowerlimit, double upperlimit, ObservableFunction fn, const char * latexname)
	:	BCVariable(name,lowerlimit,upperlimit,latexname)
	, fFunctionPointer(fn)
	, fObservablePointer(0)
{
	fPrefix = "Observable";
}

// ---------------------------------------------------------

BCObservable::~BCObservable() {
}

// ---------------------------------------------------------

double BCObservable::Evaluate(std::vector<double> const & parameters) {
	if (fObservablePointer)
		return *fObservablePointer;
	if (fFunctionPointer)
		return (*fFunctionPointer)(parameters);
	return 0;
}
