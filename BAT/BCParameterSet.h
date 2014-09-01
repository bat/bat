#ifndef __BCPARAMETERSET__H
#define __BCPARAMETERSET__H

/**
 * @class BCParameterSet Wrapper to allow access by name into list of BCParameter.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Parameters are not owned, and will not be deleted by BCParameterSet.
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>
#include <string>

#include "BCVariableSet.h"

// ---------------------------------------------------------

class BCParameterSet : public BCVariableSet {
public:

	/**
	 * Constructor */
	BCParameterSet() : BCVariableSet() 
	{}

	/*
	 * Destructor */
	~BCParameterSet()
	{}

	/**
	 * @return The number of fixed parameters. */
	unsigned int GetNFixedParameters();

	/**
	 * @return The number of free parameters. */
	unsigned int GetNFreeParameters()
	{ return Size() - GetNFixedParameters(); }


};
#endif
