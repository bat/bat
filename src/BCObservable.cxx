/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCObservable.h"
#include "BCLog.h"

#include <TString.h>

// ---------------------------------------------------------

BCObservable::BCObservable():
	fPrefix("Observable"),
	fName("parameter"),
	fLowerLimit(0),
	fUpperLimit(1),
	fFillHistograms(true),
	fNbins(100)
{
}

// ---------------------------------------------------------

BCObservable::BCObservable(const char * name, double lowerlimit, double upperlimit, const char * latexname) :
	fPrefix("Observable"),
	fName(name),
	fLowerLimit(lowerlimit),
	fUpperLimit(upperlimit),
	fLatexName(latexname),
	fFillHistograms(true),
	fNbins(100)
{
}

// ---------------------------------------------------------

void BCObservable::PrintSummary() const {
	BCLog::OutSummary(Form("%s summary:",fPrefix.c_str()));
	BCLog::OutSummary(Form("%11s : %s", fPrefix.c_str(),fName.c_str()));
	BCLog::OutSummary(Form("Lower limit : %f", fLowerLimit));
	BCLog::OutSummary(Form("Upper limit : %f", fUpperLimit));
}
