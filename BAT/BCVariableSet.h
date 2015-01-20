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
class TRandom;

// ---------------------------------------------------------

class BCVariableSet {
public:

	/**
	 * Constructor */
	BCVariableSet();

	/*
	 * Destructor */
	virtual ~BCVariableSet()
	{}

	/*
	 * Assignment operator. */
	BCVariableSet & operator=(const BCVariableSet & rhs);

   /**
    * Add a parameter if no parameter of same name exists yet.
    * @param par Variable
    * @return True if successful. */
	virtual bool Add(BCVariable * par);

	/**
	 * Clear vector of variables. */
	virtual void Clear(bool hard);

	/**
	 * Raw and fast access.
	 * @param index Index
	 * @return Variable */
	virtual BCVariable * operator[](unsigned index) const
	{	return fPars[index]; }

	/**
	 * Safe access, but slightly less efficient access to parameter.
	 * @param index Index gets checked.
	 * @return The pointer at index position or NULL if invalid index. */
	virtual BCVariable * Get(unsigned index) const
	{	return ValidIndex(index, "Get") ? fPars[index] : NULL; }
	
	/**
	 * Safe access, but slightly less efficient access to parameter.
	 * @param name Look up name in list
	 * @return The pointer at index position or NULL if invalid index. */
	virtual BCVariable * Get(const std::string & name) const
	{	return Get(Index(name)); }
	
	/**
	 * Find index of parameter identified by name;
	 * return Size() if name not found. */
	virtual unsigned Index(const std::string & name) const;
	
	/**
	 * Number of parameters contained */
	virtual unsigned Size() const
	{ return fPars.size(); }
	
	/**
	 * Whether container is empty */
	virtual bool Empty() const
	{ return fPars.empty(); }
	
	/**
	 * Check if indes is in range
	 * @param index Index
	 * @param caller Optional string to identify caller for debug output
	 * @return	 */
	virtual bool ValidIndex(unsigned index, const std::string caller="CheckIndex") const;
	
	/**
	 * Set number of bins for all parameters
	 * @param nbins Number of bins. */
	virtual void SetNBins(unsigned nbins);

	/**
	 * Set precision for output of all parameters
	 * @param n Number of significant digits for printing.*/
	virtual void SetPrecision(unsigned n);

	/**
	 * Set fill-histograms flag for 1D and 2D histograms for all parameters.
	 * @parag flag Filling flag for both 1D and 2D histograms. */
	virtual void FillHistograms(bool flag)
	{ FillHistograms(flag,flag); }

	/**
	 * Set fill-histograms flag for 1D and 2D histograms for all parameters.
	 * @parag flag_1d Filling flag for 1D histograms.
	 * @parag flag_2d Filling flag for 2D histograms. */
	virtual void FillHistograms(bool flag_1d, bool flag_2d)
	{ FillH1(flag_1d); FillH2(flag_2d); }

	/**
	 * Set fill-histograms flag for all 1D histograms for all parameters.
	 * @parag flag Filling flag. */
	virtual void FillH1(bool flag);

	/**
	 * Set fill-histograms flag for all 2D histograms for all parameters.
	 * @param flag Filling flag. */
	virtual void FillH2(bool flag);

	/**
	 * Set filling of individual H2 histogram
	 * @param x index of variable to be horizontal axis
	 * @param y index of variable to be vertical axis
	 * @param flag bool to toggle filling of H2*/
	virtual void FillH2(unsigned x, unsigned y, bool flag)
	{ if (x<fFillH2.size() and y<fFillH2[x].size()) fFillH2[x][y]=flag;}

	/**
	 * Set filling of individual H2 histogram
	 * @param x name of variable to be horizontal axis
	 * @param y name of variable to be vertical axis
	 * @param flag bool to toggle filling of H2*/
	virtual void FillH2(std::string x, std::string y, bool flag)
	{ FillH2(Index(x),Index(y)); }
	
	/**
	 * @return Whether to fill particular 2D marginalization.
	 * @param x index of variable on horizontal axis
	 * @param y index of variable on vertical axis
	 * @return whether to fill histogram */
	bool FillH2(unsigned x, unsigned y) const;

	/**
	 * @return Whether to fill particular 2D marginalization.
	 * @param x name of variable on horizontal axis
	 * @param y name of variable on vertical axis
	 * @return whether to fill histogram */
	bool FillH2(std::string x, std::string y) const
	{ return FillH2(Index(x),Index(y)); }

	/**
	 * @return Length of longest parameter name. */
	virtual unsigned MaxNameLength()
	{ return fMaxNameLength; }

	/**
	 * @return Volume of set. */
	virtual double Volume() const;

	/**
	 * @return Whether values are within limits of variables. */
	virtual bool IsWithinLimits(const std::vector<double> & x) const;

	/**
	 * return positions in ranges of given values
	 * from 0 (at lower limit) to 1 (at upper limit)
	 * for each variable in the set.
	 * @param x vector of values to report positions of.
	 * @return vector of positions of values in ranges. */
	virtual std::vector<double> PositionInRange(const std::vector<double> & x) const;

	/**
	 * Translate from unit interval to values in variable ranges.
	 * @param p vector of positions in the unit interval (0 = lower limit, 1 = upper limit). */
	virtual void ValueFromPositionInRange(std::vector<double> &p) const;

	/**
	 * @return vector of range centers. */
	virtual std::vector<double> GetRangeCenters() const;
	
	/**
	 * Get vector of uniformly distributed random values.
	 * @param R Random number generator to use.
	 * @return vector of random values uniformly distributed in variable ranges. */
	virtual std::vector<double> GetUniformRandomValues(TRandom * const R) const;

protected:
	/// Don't own parameters
	std::vector<BCVariable*> fPars;

	unsigned fMaxNameLength;

	std::vector<std::vector<bool> > fFillH2;
};
#endif
