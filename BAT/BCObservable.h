#ifndef __BCOBSERVABLE__H
#define __BCOBSERVABLE__H

/*!
 * \class BCObservable
 * \brief A class representing a observable of a model.
 * \author Daniel Greenwald
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a observable of a model. It contains
 * information about the name and the range of the observable.
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

class BCObservable {

public:

	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * The default constructor. */
	BCObservable();

	/**
	 * A constructor.
	 * @param name The name of the observable.
	 * @param lowerlimit The lower limit of the observable values.
	 * @param upperlimit The upper limit of the observable values.
	 * @param latexname The latex name of the observable used in axis labeling.
	 */
	BCObservable(const char* name, double lowerlimit, double upperlimit, const char* latexname = "");

	/** \name Member functions (get) */
	/** @{ */
	
	/**
	 * @return The name of the observable. */
	const std::string & GetName() const
	{ return fName; }

	/**
	 * Returns latex name if set, else identical to GetName().
	 */
	const std::string & GetLatexName() const
	{ return (fLatexName.empty()) ? fName : fLatexName; }
	
	/**
	 * @return The lower limit of the observable values. */
	double GetLowerLimit() const
	{ return fLowerLimit; }

	/**
	 * @return The upper limit of the observable values. */
	double GetUpperLimit() const
	{ return fUpperLimit; }

	/**
	 * Returns the range width of the observable values. It is
	 * always a positive value.
	 * @return The range width of the observable values. */
	double GetRangeWidth() const
	{ return (fUpperLimit > fLowerLimit) ? fUpperLimit - fLowerLimit : fLowerLimit - fUpperLimit; }
	
	bool FillHistograms() const
	{ return fFillHistograms; }
	
	unsigned GetNbins() const
	{ return fNbins; }

	/** @} */

	/** \name Member functions (set) */
	/** @{ */

	/**
	 * @param name The name of the observable. */
	void SetName(const char * name)
	{ fName = name; }

	void SetLatexName(const char * latex_name)
	{ fLatexName = latex_name; }

	/**
	 * Set the lower limit of the observable values.
	 * @param limit The lower limit of the observable values. */
	void SetLowerLimit(double limit = 0)
	{ fLowerLimit = limit; }

	/**
	 * Set the upper limit of the observable values.
	 * @param limit The upper limit of the observable values. */
	void SetUpperLimit(double limit = 1)
	{ fUpperLimit = limit; }

	/**
	 * Set the limits of the observable values.
	 * @param lowerlimit The lower limit of the observable values.
	 * @param upperlimit The upper limit of the observable values. */
	void SetLimits(double lowerlimit = 0, double upperlimit = 1)
	{ fLowerLimit = lowerlimit; fUpperLimit = upperlimit; }

	void FillHistograms(bool flag)
	{ fFillHistograms = flag; }

	void SetNbins(unsigned nbins)
	{ fNbins = nbins; }
	/** @} */

	/** \name Member functions (miscellaneous methods) */
	/** @{ */
	
	/**
	 * Prints a observable summary on the screen. */
	void PrintSummary() const;
	
	/** @} */
	
protected:
	/// prefix for output
	std::string fPrefix;

	/// The name of the observable.
	std::string fName;
	
	///  The lower limit of the observable value.
	double fLowerLimit;
	
	/// The upper limit of the observable value.
	double fUpperLimit;
	
	/// The latex name of the observable.
	std::string fLatexName;
	
	/// Flag to store MCMC samples in histograms.
	bool fFillHistograms;
	
	/// The number of equal-size bins used in histograms involving this observable.
	unsigned fNbins;

};
#endif