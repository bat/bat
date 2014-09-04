#ifndef __BCVARIABLESET__H
#define __BCVARIABLESET__H

/**
 * @class BCVariableSet Wrapper to allow access by name into list of BCVariable.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Variables are not owned, and will not be deleted by BCVariableSet.
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

class BCVariable;

// ---------------------------------------------------------

class BCVariableSet {
public:

	/**
	 * Constructor */
	BCVariableSet();

	/*
	 * Destructor */
	~BCVariableSet()
	{}

   /**
    * Add a parameter if no parameter of same name exists yet.
    *
    * @param par Variable
    * @return True if successful.
    */
	bool Add(BCVariable * par);

	void Clear(bool hard);

	/**
	 * Raw and fast access.
	 *
	 * @param index Index
	 * @return Variable
	 */
	BCVariable * operator[](unsigned index) const
	{
		return fPars[index];
	}

	/**
	 * Safe access, but slightly less efficient access to parameter.
	 *
	 * @param index Index gets checked.
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCVariable * Get(unsigned index) const
	{
		return ValidIndex(index, "Get") ? fPars[index] : NULL;
	}
	
	/**
	 * Safe access, but slightly less efficient access to parameter.
	 *
	 * @param name Look up name in list
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCVariable * Get(const std::string & name) const
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
	 * Whether container is empty
	 */
	bool Empty() const
	{ return fPars.empty(); }

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

protected:
	/// Don't own parameters
	std::vector<BCVariable*> fPars;

	unsigned fMaxNameLength;
};
#endif
