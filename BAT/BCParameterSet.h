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

class BCParameter;

// ---------------------------------------------------------

class BCParameterSet {
public:

	/**
	 * Constructor */
	BCParameterSet();

	/*
	 * Destructor */
	~BCParameterSet()
	{}

   /**
    * Add a parameter if no parameter of same name exists yet.
    *
    * @param par Parameter
    * @return True if successful.
    */
	bool Add(BCParameter * par);

	void Clear(bool hard);

	/**
	 * Raw and fast access.
	 *
	 * @param index Index
	 * @return Parameter
	 */
	BCParameter * operator[](unsigned index) const
	{
		return fPars[index];
	}

	/**
	 * Safe access, but slightly less efficient access to parameter.
	 *
	 * @param index Index gets checked.
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCParameter * Get(unsigned index) const
	{
		return ValidIndex(index, "Get") ? fPars[index] : NULL;
	}
	
	/**
	 * Safe access, but slightly less efficient access to parameter.
	 *
	 * @param name Look up name in list
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCParameter * Get(const std::string & name) const
	{
		return Get(Index(name));
	}
	
	/**
	 * Find index of parameter identified by name
	 */
	unsigned Index(const std::string & name) const;
	
	/**
	 * Number of parameters contained
	 */
	unsigned Size() const
	{ return fPars.size(); }
	
	/**
	 * Check if indes is in range
	 * @param index Index
	 * @param caller Optional string to identify caller for debug output
	 * @return
	 */
	bool ValidIndex(unsigned index, const std::string caller="CheckIndex") const;
	
	/**
	 * Set number of bins for all parameters
	 * @param nbins Number of bins. */
	void SetNBins(unsigned nbins);

	/**
	 * Set precision for output of all parameters
	 * @param n Number of significant digits for printing.*/
	void SetPrecision(unsigned n);

	/**
	 * Set fill-histograms flag for all parameters.
	 * @parag flag Filling flag. */
	void FillHistograms(bool flag);

	/**
	 * @return Length of longest parameter name. */
	unsigned MaxNameLength()
	   { return fMaxNameLength; }

private:
	/// Don't own parameters
	std::vector<BCParameter*> fPars;

	unsigned fMaxNameLength;
};
#endif
