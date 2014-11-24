#ifndef __BCVARIABLE__H
#define __BCVARIABLE__H

/*!
 * \class BCVariable
 * \brief A class representing a variable of a model.
 * \author Daniel Greenwald
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a variable of a model. It contains
 * information about the name and the range of the variable.
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
#include <cmath>

class TH1D;
class TH2D;

// ---------------------------------------------------------

class BCVariable {

public:

	/** \name Constructors and destructors */
	/** @{ */
	
	/**
	 * The default constructor. */
	BCVariable();

	/**
	 * Copy constructor. */
	BCVariable(const BCVariable & other);

	/**
	 * A constructor.
	 * @param name The name of the variable.
	 * @param lowerlimit The lower limit of the variable values.
	 * @param upperlimit The upper limit of the variable values.
	 * @param latexname The latex name of the variable used in axis labeling.
	 */
	BCVariable(const char* name, double lowerlimit, double upperlimit, const char* latexname = "");

	/**
	 * Destructor */
	virtual ~BCVariable();

	/** \name Member functions (get) */
	/** @{ */
	
	/**
	 * @return The name of the variable. */
	const std::string & GetName() const
	{ return fName; }

	/**
	 * @return Safe name of the variable. */
	const std::string & GetSafeName() const
	{ return fSafeName; }

	/**
	 * Returns latex name if set, else identical to GetName().
	 */
	const std::string & GetLatexName() const
	{ return (fLatexName.empty()) ? fName : fLatexName; }
	
	/**
	 * @return The lower limit of the variable values. */
	double GetLowerLimit() const
	{ return fLowerLimit; }

	/**
	 * @return The upper limit of the variable values. */
	double GetUpperLimit() const
	{ return fUpperLimit; }

	/**
	 * Returns the range width of the variable values. It is
	 * always a positive value.
	 * @return The range width of the variable values. */
	double GetRangeWidth() const
	{ return (fUpperLimit > fLowerLimit) ? fUpperLimit - fLowerLimit : fLowerLimit - fUpperLimit; }

	/**
	 * Returns center of variable range.*/
	double GetRangeCenter() const
	{ return (fUpperLimit+fLowerLimit)/2.; }

	/**
	 * @return precision of output */
	unsigned GetPrecision()
	{ return fPrecision;}
	
	bool FillHistograms() const
	{ return fFillHistograms; }
	
	unsigned GetNbins() const
	{ return fNbins; }

	/** @} */

	/** \name Member functions (set) */
	/** @{ */

	/**
	 * @param name The name of the variable. */
	void SetName(const char * name);

	void SetLatexName(const char * latex_name)
	{ fLatexName = latex_name; }

	/**
	 * Set the lower limit of the variable values.
	 * @param limit The lower limit of the variable values. */
	void SetLowerLimit(double limit = 0)
	{ fLowerLimit = limit; CalculatePrecision();}

	/**
	 * Set the upper limit of the variable values.
	 * @param limit The upper limit of the variable values. */
	void SetUpperLimit(double limit = 1)
	{ fUpperLimit = limit; CalculatePrecision();}

	/**
	 * Set the limits of the variable values.
	 * @param lowerlimit The lower limit of the variable values.
	 * @param upperlimit The upper limit of the variable values. */
	void SetLimits(double lowerlimit = 0, double upperlimit = 1)
	{ fLowerLimit = lowerlimit; fUpperLimit = upperlimit; CalculatePrecision();}

	/**
	 * Set the precision of the output of variable
	 * @param precision The precision of the variable for output. */
	void SetPrecision(unsigned precision)
	{ fPrecision = precision; }

	void FillHistograms(bool flag)
	{ fFillHistograms = flag; }

	void SetNbins(unsigned nbins)
	{ fNbins = nbins; }
	/** @} */

	/** \name Member functions (miscellaneous methods) */
	/** @{ */

	/**
	 * Check if name is that provided by parameter.
	 * @param name Name to check against parameter name. */
	bool IsNamed(std::string name)
	{ return fName.compare(name) == 0; }

	/**
	 * return position in range of given value
	 * from 0 (at lower limit) to 1 (at upper limit)
	 * @param x Value to report position of.
	 * @return Position of value in range. */
	double PositionInRange(double x)
	{ return (x - fLowerLimit)/(fUpperLimit-fLowerLimit); }

	/**
	 * Calculate the necessary precision for outputting this
	 * parameter */
	void CalculatePrecision()
	{ SetPrecision(ceil(-log10(2.*fabs(fUpperLimit-fLowerLimit)/(fabs(fUpperLimit)+fabs(fLowerLimit))))); }

	/**
	 * Prints a variable summary on the screen. */
	void PrintSummary() const;

	/**
	 * @return A one line summary of the variable. */
	virtual std::string OneLineSummary();

	/**
	 * Creates a 1D Histogram for this variable.
	 * @param name Name of the histogram.
	 * @return pointer to histogram object. */
	virtual TH1D * CreateH1(const char * name);
	
	/**
	 * Creates a 2D Histogram for this variable as the abcissa
	 * and a second as the ordinate.
	 * @name name The name of the histogram.
	 * @param ordinate The variable to be used for the ordinate. */
	virtual TH2D * CreateH2(const char * name, BCVariable * ordinate);

	/** @} */
	
protected:
	/// prefix for output
	std::string fPrefix;

	/// The name of the variable.
	std::string fName;
	
	/// Safe name of the variable for use in ROOT object naming
	std::string fSafeName;

	///  The lower limit of the variable value.
	double fLowerLimit;
	
	/// The upper limit of the variable value.
	double fUpperLimit;
	
	/// Necessary precision for output
	unsigned fPrecision;

	/// The latex name of the variable.
	std::string fLatexName;
	
	/// Flag to store MCMC samples in histograms.
	bool fFillHistograms;
	
	/// The number of equal-size bins used in histograms involving this variable.
	unsigned fNbins;

};
#endif
