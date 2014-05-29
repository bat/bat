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
	BCObservable(),
	fFixed(false),
	fFixedValue(-1.e+111)
{
	fPrefix = "UserObservable";
}

// ---------------------------------------------------------

BCUserObservable::BCUserObservable(const char * name, double lowerlimit, double upperlimit, const char * latexname) :
	BCObservable(name,lowerlimit,upperlimit,latexname),
	fFixed(false),
	fFixedValue(-1.e+111)
{
	fPrefix = "UserObservable";
}

// ---------------------------------------------------------

bool BCUserObservable::IsAtLimit(double value) const
{
   if (fLowerLimit == fUpperLimit)
      return false;

   if ( ( (value-fLowerLimit)*(value-fLowerLimit)/fLowerLimit/fLowerLimit <= 1e-10) ||
         ( (value-fUpperLimit)*(value-fUpperLimit)/fUpperLimit/fUpperLimit <= 1e-10))
      return true;
   else
      return false;
}

