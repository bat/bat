/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCParameter.h"

#include <string>

#include <TString.h>


// ---------------------------------------------------------

BCParameter::BCParameter()
	:	BCVariable()
	, fFixed(false)
	, fFixedValue(-1.e+111)
{
	fPrefix = "Parameter";
}

// ---------------------------------------------------------

BCParameter::BCParameter(const BCParameter & other)
	: BCVariable(other)
	, fFixed(other.fFixed)
	, fFixedValue(other.fFixedValue)
{
}

// ---------------------------------------------------------

BCParameter::BCParameter(const char * name, double lowerlimit, double upperlimit, const char * latexname) :
	BCVariable(name,lowerlimit,upperlimit,latexname),
	fFixed(false),
	fFixedValue(-1.e+111)
{
	fPrefix = "Parameter";
}

// ---------------------------------------------------------

bool BCParameter::IsAtLimit(double value) const
{
   if (fLowerLimit == fUpperLimit)
      return false;

   if ( ( (value-fLowerLimit)*(value-fLowerLimit)/fLowerLimit/fLowerLimit <= 1e-10) ||
         ( (value-fUpperLimit)*(value-fUpperLimit)/fUpperLimit/fUpperLimit <= 1e-10))
      return true;
   else
      return false;
}

// ---------------------------------------------------------

std::string BCParameter::OneLineSummary() {
	if (!Fixed())
		return BCVariable::OneLineSummary();
	return std::string(Form("%s (fixed at %.*f)",BCVariable::OneLineSummary().data(),GetPrecision(),GetFixedValue()));
}
