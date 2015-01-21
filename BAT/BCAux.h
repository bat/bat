#ifndef __BCAUX__H
#define __BCAUX__H

/*!
 * \namespace BCAux
 * \brief Some functions not fitting anywhere else
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 01.2009
 * \detail A namespace which encapsulates auxiliary functions
 * necessary for BAT.
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <string>

// ---------------------------------------------------------

namespace BCAux {

	/**
	 * Sets the default BAT style for drawing plots. */
	void SetStyle();
	
	/**
	 * Force file extension to be .pdf if not already .pdf or .ps
	 * @param filename Filename to be altered */
	void ForceToBePDF(std::string & filename);
}

// ---------------------------------------------------------

#endif

