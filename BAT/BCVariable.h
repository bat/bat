#ifndef __BCVARIABLE__H
#define __BCVARIABLE__H

/**
 * @class BCVariable
 * @brief A class representing a variable of a model.
 * @author Daniel Greenwald
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class represents a variable of a model. It contains
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

class TH1;
class TH2;
class TH3;
class TRandom;

// ---------------------------------------------------------

class BCVariable
{

public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * The default constructor. */
    BCVariable();

    /**
     * A constructor.
     * @param name The name of the variable.
     * @param lowerlimit The lower limit of the variable values.
     * @param upperlimit The upper limit of the variable values.
     * @param latexname The latex name of the variable used in axis labeling.
     * @param unitstring Unit string to be printed for variable. */
    BCVariable(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname = "", const std::string& unitstring = "");

    /**
     * Destructor */
    virtual ~BCVariable() {};

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return Prefix for name of type of variable ("Parameter", "Observable") */
    virtual const std::string& GetPrefix() const
    { return fPrefix; }

    /**
     * @return The name of the variable. */
    virtual const std::string& GetName() const
    { return fName; }

    /**
     * @return Safe name of the variable. */
    virtual const std::string& GetSafeName() const
    { return fSafeName; }

    /**
     * @return latex name if set, else identical to GetName(). */
    virtual const std::string& GetLatexName() const
    { return (fLatexName.empty()) ? fName : fLatexName; }

    /**
     * @return unit string */
    virtual const std::string& GetUnitString() const
    { return fUnitString; }

    /**
     * @return latex name with unit string */
    virtual std::string GetLatexNameWithUnits() const
    { if (GetUnitString().empty()) return GetLatexName(); return GetLatexName() + " " + GetUnitString(); }

    /**
     * @return The lower limit of the variable values. */
    virtual double GetLowerLimit() const
    { return fLowerLimit; }

    /**
     * @return The upper limit of the variable values. */
    virtual double GetUpperLimit() const
    { return fUpperLimit; }

    /**
     * Returns the range width of the variable values. It is
     * always a positive value.
     * @return The range width of the variable values. */
    virtual double GetRangeWidth() const
    { return fUpperLimit - fLowerLimit; }

    /**
     * Returns center of variable range.*/
    virtual double GetRangeCenter() const
    { return (fUpperLimit + fLowerLimit) / 2.; }

    /**
     * @return precision of output */
    virtual unsigned GetPrecision() const
    { return fPrecision;}

    /**
     * @return whether to fill 1D histogram. */
    virtual bool FillH1() const
    { return fFillH1; }

    /**
     * @return whether to fill 2D histograms. */
    virtual bool FillH2() const
    { return fFillH2; }

    /**
     * @return Number of bins on axis in 1D and 2D histograms. */
    virtual unsigned GetNbins() const
    { return fNbins; }

    /** @} */

    /** \name Member functions (set) */
    /** @{ */

    /**
     * @param name The name of the variable. */
    virtual void SetName(const std::string& name);

    /**
     * @param latex_name Latex-formatted name of Variable. */
    virtual void SetLatexName(const std::string& latex_name)
    { fLatexName = latex_name; }

    /**
     * @param unit_string String to printed to mark units of variable when needed. */
    virtual void SetUnitString(const std::string& unit_string)
    { fUnitString = unit_string; }

    /**
     * Set the lower limit of the variable values.
     * @param limit The lower limit of the variable values. */
    virtual void SetLowerLimit(double limit)
    { SetLimits(limit, fUpperLimit); }

    /**
     * Set the upper limit of the variable values.
     * @param limit The upper limit of the variable values. */
    virtual void SetUpperLimit(double limit)
    { SetLimits(fLowerLimit, limit); }

    /**
     * Set the limits of the variable values.
     * @param lowerlimit The lower limit of the variable values.
     * @param upperlimit The upper limit of the variable values. */
    virtual void SetLimits(double lowerlimit = 0, double upperlimit = 1);

    /**
     * Set the precision of the output of variable
     * @param precision The precision of the variable for output. */
    virtual void SetPrecision(unsigned precision)
    { fPrecision = precision; }

    /**
     * Set the filling of 1D and 2D histograms.
     * @param flag Toggles filling of 1D and 2D histograms. */
    virtual void FillHistograms(bool flag)
    { FillHistograms(flag, flag); }

    /**
     * Set the filling of 1D and 2D histograms.
     * @param fill_1d Toggles filling of 1D histogram.
     * @param fill_2d Toggles filling of 2D histograms. */
    virtual void FillHistograms(bool fill_1d, bool fill_2d)
    { FillH1(fill_1d); FillH2(fill_2d); }

    /**
     * Set the filling of 1D histogram. */
    virtual void FillH1(bool flag)
    { fFillH1 = flag; }

    /**
     * Set the filling of 2D histograms. */
    virtual void FillH2(bool flag)
    { fFillH2 = flag; }


    virtual void SetNbins(unsigned nbins)
    { fNbins = nbins; }
    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Check if name is that of variable.
     * @param name Name to check against variable name. */
    virtual bool IsNamed(const std::string& name) const
    { return fName.compare(name) == 0; }

    /**
     * Check if safe name is that of variable.
     * @param safename Safe name to check against variable name. */
    virtual	bool IsSafeNamed(const std::string& safename) const
    { return fSafeName.compare(safename) == 0; }

    /**
     * return position in range of given value
     * from 0 (at lower limit) to 1 (at upper limit)
     * @param x Value to report position of.
     * @return Position of value in range. */
    virtual double PositionInRange(double x) const
    { return (x - fLowerLimit) / (fUpperLimit - fLowerLimit); }

    /**
     * Translate from unit interval to value in variable range.
     * @param p Position in unit interval (0 = lower limit, 1 = upper limit).
     * @return Value of variable at position in range. */
    virtual double ValueFromPositionInRange(double p) const
    { return fLowerLimit + p * (fUpperLimit - fLowerLimit); }

    /**
     * Calculate the necessary precision for outputting this parameter
     * and replace current precision is smaller or if force is set true
     * @param force replace current precision even if calculated precision is lower than current precision.*/
    virtual void CalculatePrecision(bool force = false);

    /**
     * @return Whether value is at upper or lower limit. */
    virtual bool IsAtLimit(double value) const;

    /**
     * @return Whether value is within limits. */
    virtual bool IsWithinLimits(double value) const
    { return (value >= fLowerLimit and value <= fUpperLimit); }

    /**
     * Prints a variable summary on the screen. */
    virtual void PrintSummary() const;

    /**
     * @return A one line summary of the variable.
     * @param print_prefix Whether to print prefix before variable name.
     * @param name_length Number of characters to pad name to (if positive) */
    virtual std::string OneLineSummary(bool print_prefix = true, int name_length = -1) const;

    /**
     * Creates a 1D Histogram for this variable.
     * @param name Name of the histogram.
     * @return pointer to histogram object. */
    virtual TH1* CreateH1(const std::string& name) const;

    /**
     * Creates a 2D Histogram for this variable as the abcissa
     * and a second as the ordinate.
     * @name name The name of the histogram.
     * @param ordinate The variable to be used for the ordinate. */
    virtual TH2* CreateH2(const std::string& name, const BCVariable& ordinate) const;

    /**
     * Creates a 3D Histogram for this variable as the abcissa
     * and a second as the ordinate.
     * @name name The name of the histogram.
     * @param ordinate_y The variable to be used for the y ordinate.
     * @param ordinate_z The variable to be used for the z ordinate. */
    virtual TH3* CreateH3(const std::string& name, const BCVariable& ordinate_y, const BCVariable& ordinate_z) const;

    /**
     * Get random value uniformly distributed in range.
     * @param R random number generator to use*/
    virtual double GetUniformRandomValue(TRandom* const R) const;

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

    /// Unit string for variable
    std::string fUnitString;

    /// Flag to store MCMC samples in 1D histogram.
    bool fFillH1;

    /// Flag to store MCMC samples in 2D histograms
    bool fFillH2;

    /// The number of equal-size bins used in histograms involving this variable.
    unsigned fNbins;

};
#endif
