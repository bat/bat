/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCParameter.h"
#include "BCLog.h"

#include <TString.h>

#include <algorithm>
#include <fstream>
#include <iostream>

// ---------------------------------------------------------

BCParameter::BCParameter():
    fName("parameter"),
    fLowerLimit(0),
    fUpperLimit(1),
    fFixed(false),
    fFixedValue(-1.e+111),
    fFillHistograms(true),
    fNbins(100)
{
}

// ---------------------------------------------------------

BCParameter::BCParameter(const char * name, double lowerlimit, double upperlimit, const char * latexname) :
    fName(name),
    fLowerLimit(lowerlimit),
    fUpperLimit(upperlimit),
    fLatexName(latexname),
    fFixed(false),
    fFixedValue(-1.e+111),
    fFillHistograms(true),
    fNbins(100)
{
}

// ---------------------------------------------------------

void BCParameter::PrintSummary() const
{
	BCLog::OutSummary("Parameter summary:");
	BCLog::OutSummary(Form("Parameter   : %s", fName.c_str()));
	BCLog::OutSummary(Form("Lower limit : %f", fLowerLimit));
	BCLog::OutSummary(Form("Upper limit : %f", fUpperLimit));
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

