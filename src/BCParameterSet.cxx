/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCParameterSet.h"
#include "BCParameter.h"

// ---------------------------------------------------------
unsigned int BCParameterSet::GetNFixedParameters() {
	unsigned int n = 0;
	for (unsigned int i = 0; i < Size(); ++i)
		if (((BCParameter*)fPars[i])->Fixed())
			++n;
	return n;
}
