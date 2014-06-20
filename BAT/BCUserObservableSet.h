#ifndef __BCUSEROBSERVABLESET__H
#define __BCUSEROBSERVABLESET__H

/**
 * @class BCUserObservableSet Wrapper to allow access by name into list of BCUserObservable.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note UserObservables are not owned, and will not be deleted by BCUserObservableSet.
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

class BCUserObservable;

// ---------------------------------------------------------

class BCUserObservableSet {
public:

	/**
	 * Constructor */
	BCUserObservableSet();

	/*
	 * Destructor */
	~BCUserObservableSet()
	{}

   /**
    * Add a user-defined observable if no user-defined observable of same name exists yet.
    *
    * @param par UserObservable
    * @return True if successful.
    */
	bool Add(BCUserObservable * par);

	void Clear(bool);

	/**
	 * Raw and fast access.
	 *
	 * @param index Index
	 * @return UserObservable
	 */
	BCUserObservable * operator[](unsigned index) const
	{
		return fPars[index];
	}

	/**
	 * Safe access, but slightly less efficient access to user-defined observable.
	 *
	 * @param index Index gets checked.
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCUserObservable * Get(unsigned index) const
	{
		return ValidIndex(index, "Get") ? fPars[index] : NULL;
	}
	
	/**
	 * Safe access, but slightly less efficient access to user-defined observable.
	 *
	 * @param name Look up name in list
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCUserObservable * Get(const std::string & name) const
	{
		return Get(Index(name));
	}
	
	/**
	 * Find index of user-defined observable identified by name
	 */
	unsigned Index(const std::string & name) const;
	
	/**
	 * Number of user-defined observables contained
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
	 * Set number of bins for all observables.
	 * @param nbins Number of bins. */
	void SetNBins(unsigned nbins);

	/**
	 * Set precision for output of all observables.
	 * @param n Number of significant digits for printing.*/
	void SetPrecision(unsigned n);

	/**
	 * Set fill-histograms flag for all observables.
	 * @parag flag Filling flag. */
	void FillHistograms(bool flag);

	/**
	 * @return Length of longest parameter name. */
	unsigned MaxNameLength()
	   { return fMaxNameLength; }
	
private:
	/// Don't own user-defined observables
	std::vector<BCUserObservable*> fPars;

	unsigned fMaxNameLength;
};
#endif
