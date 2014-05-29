#ifndef __BCOBSERVABLESET__H
#define __BCOBSERVABLESET__H

/**
 * @class BCObservableSet Wrapper to allow access by name into list of BCObservable.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Observables are not owned, and will not be deleted by BCObservableSet.
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

class BCObservable;

// ---------------------------------------------------------

class BCObservableSet {
public:
   /**
    * Add a observable if no observable of same name exists yet.
    *
    * @param par Observable
    * @return True if successful.
    */
	bool Add(BCObservable * par);

	void Clear(bool);

	/**
	 * Raw and fast access.
	 *
	 * @param index Index
	 * @return Observable
	 */
	BCObservable * operator[](unsigned index) const
	{
		return fPars[index];
	}

	/**
	 * Safe access, but slightly less efficient access to observable.
	 *
	 * @param index Index gets checked.
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCObservable * Get(unsigned index) const
	{
		return ValidIndex(index, "Get") ? fPars[index] : NULL;
	}
	
	/**
	 * Safe access, but slightly less efficient access to observable.
	 *
	 * @param name Look up name in list
	 * @return The pointer at index position or NULL if invalid index.
	 */
	BCObservable * Get(const std::string & name) const
	{
		return Get(Index(name));
	}
	
	/**
	 * Find index of observable identified by name
	 */
	unsigned Index(const std::string & name) const;
	
	/**
	 * Number of observables contained
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
	
private:
	/// Don't own observables
	std::vector<BCObservable*> fPars;
};
#endif
