/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <TH1D.h>

#include "BCFBUBackground.h"

// ---------------------------------------------------------
BCFBUBackground::BCFBUBackground(std::string name)
	: fHistogram(0)
	, fName(name)
{
}

// ---------------------------------------------------------
BCFBUBackground::~BCFBUBackground()
{
	if (fHistogram)
		delete fHistogram;
}

// ---------------------------------------------------------
void BCFBUBackground::SetHistogram(TH1D * hist)
{
	fHistogram = hist;
}

// ---------------------------------------------------------
