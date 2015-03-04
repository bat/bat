/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCTF1LogPrior.h"

#include <cmath>

#include <TF1.h>

// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(TF1 const * const f)
	: BCPrior()
	, fLogPriorFunction(NULL)
{
	if (f)
		fLogPriorFunction = new TF1(*f);
}


// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(const char * formula, double xmin, double xmax)
	: BCPrior()
	, fLogPriorFunction(new TF1("f1_logprior",formula,xmin,xmax))
{
}

// ---------------------------------------------------------
BCTF1LogPrior::BCTF1LogPrior(BCTF1LogPrior const & other)
	: BCPrior(other)
{
	if (other.fLogPriorFunction)
		fLogPriorFunction = new TF1(*fLogPriorFunction);
}

// ---------------------------------------------------------
BCTF1LogPrior::~BCTF1LogPrior() {
	delete fLogPriorFunction;
}


