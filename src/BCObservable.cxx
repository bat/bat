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
{
	fPrefix = "Observable";
}

// ---------------------------------------------------------

BCObservable::BCObservable(const char * name, double lowerlimit, double upperlimit, const char * latexname)
	:	BCVariable(name,lowerlimit,upperlimit,latexname)
{
	fPrefix = "Observable";
}

// ---------------------------------------------------------

BCObservable::~BCObservable() {
}

// ---------------------------------------------------------
