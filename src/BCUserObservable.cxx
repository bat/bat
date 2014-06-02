/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCUserObservable.h"

#include <string>


// ---------------------------------------------------------

BCUserObservable::BCUserObservable() :
	BCObservable()
{
	fPrefix = "UserObservable";
}

// ---------------------------------------------------------

BCUserObservable::BCUserObservable(const char * name, double lowerlimit, double upperlimit, ObservableFunction * fn, const char * latexname) :
	BCObservable(name,lowerlimit,upperlimit,latexname)
{
	fPrefix = "UserObservable";
}


