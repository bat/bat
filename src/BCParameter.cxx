/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <iostream>
#include <fstream>

#include <TROOT.h>

#include "BCLog.h"
#include "BCParameter.h"

// ---------------------------------------------------------

BCParameter::BCParameter()
{
   fName       = "parameter";
   fLowerLimit = 0.;
   fUpperLimit = 1.;
   fNuisance   = 0;
}

// ---------------------------------------------------------

BCParameter::BCParameter(const char * name, double lowerlimit, double upperlimit)
{
   fName       = "parameter";
   fLowerLimit = 0.;
   fUpperLimit = 1.;
   fNuisance   = 0;

   fName       = name;
   fLowerLimit = lowerlimit;
   fUpperLimit = upperlimit;
}

// ---------------------------------------------------------

BCParameter::BCParameter(const BCParameter & parameter)
{
   fName       = parameter.fName;
   fIndex      = parameter.fIndex;
   fLowerLimit = parameter.fLowerLimit;
   fUpperLimit = parameter.fUpperLimit;
   fNuisance   = parameter.fNuisance;
}

// ---------------------------------------------------------

BCParameter & BCParameter::operator = (const BCParameter & parameter)
{
   fName       = parameter.fName;
   fIndex      = parameter.fIndex;
   fLowerLimit = parameter.fLowerLimit;
   fUpperLimit = parameter.fUpperLimit;
   fNuisance   = parameter.fNuisance;

   // return this
   return *this;
}

// ---------------------------------------------------------

BCParameter::~BCParameter()
{}

// ---------------------------------------------------------

void BCParameter::PrintSummary()
{
	BCLog::OutSummary("Parameter summary:");
	BCLog::OutSummary(Form("Parameter   : %s", fName.c_str()));
	BCLog::OutSummary(Form("Index       : %d", fIndex));
	BCLog::OutSummary(Form("Lower limit : %f", fLowerLimit));
	BCLog::OutSummary(Form("Upper limit : %f", fUpperLimit));
}

// ---------------------------------------------------------

bool BCParameter::IsAtLimit(double value)
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
