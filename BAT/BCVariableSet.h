#ifndef __BCVARIABLESET__H
#define __BCVARIABLESET__H

/**
 * @class BCVariableSet
 * @brief Wrapper to allow access by name into list of BCVariable.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Variables are owned by and will be deleted by BCVariableSet.
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

#include "BCLog.h"
#include "BCAux.h"

#include <TString.h>

class TRandom;

// ---------------------------------------------------------

template<class T>
class BCVariableSet
{
public:

    /**
     * Constructor */
    BCVariableSet() :
        fMaxNameLength(0)
    {
    }

    /*
     * Destructor */
    virtual ~BCVariableSet()
    {
    }

    /**
     * Add a variable if no variable of same name exists yet.
     * @param var Variable
     * @return Success of action. */
    virtual bool Add(const T& var)
    {
        // check if variable with same name or same safe name exists
        for (unsigned int i = 0; i < fVars.size() ; ++i)
            if (var.IsNamed(fVars[i].GetName())) {
                BCLog::OutError("BCVariableSet::Add : Variable with name" + var.GetName() + "exists already.");
                return false;
            } else if ( var.IsSafeNamed(fVars[i].GetSafeName()) ) {
                BCLog::OutError("BCVariableSet::Add : Variable with safe name " + var.GetSafeName() + " exists already.");
                return false;
            }

        // add var to container
        fVars.push_back(var);
        fMaxNameLength = std::max(fMaxNameLength, (unsigned)var.GetName().length());

        return true;
    }

    /**
     * Add a variable
     * @param name Name of variable
     * @param min minimum value of the variable
     * @param max maximum value of the variable
     * @param latexname Optional latexname used for plotting
     * @param unitstring Unit string to be printed for variable
     * @return Success of action. */
    virtual bool Add(const std::string& name, double min, double max, const std::string& latexname = "", const std::string& unitstring = "")
    {
        // check if variable with same name or same safe name exists
        for (unsigned int i = 0; i < fVars.size() ; ++i)
            if (fVars[i].IsNamed(name)) {
                BCLog::OutError("BCVariableSet::Add : Variable with name " + name + " exists already.");
                return false;
            } else if ( fVars[i].IsSafeNamed(BCAux::SafeName(name)) ) {
                BCLog::OutError("BCVariableSet::Add : Variable with safe name " + fVars[i].GetSafeName() + "%s exists already.");
                return false;
            }

        // add var to container
        fVars.push_back(T(name, min, max, latexname, unitstring));
        fMaxNameLength = std::max(fMaxNameLength, (unsigned)name.length());

        return true;
    }

    /**
     * Raw and fast access.
     * @param index Index
     * @return Variable */
    virtual T& operator[](unsigned index)
    {	return fVars[index]; }

    /**
     * Raw and fast access.
     * @param index Index
     * @return Variable */
    virtual const T& operator[](unsigned index) const
    {	return fVars[index]; }

    /**
     * Safe access, but slightly less efficient access to parameter.
     * @param index Index gets checked.
     * @return The pointer at index position or NULL if invalid index. */
    virtual T& At(unsigned index)
    { return fVars.at(index); }

    /**
     * Safe access, but slightly less efficient access to parameter.
     * @param index Index gets checked.
     * @return The pointer at index position or NULL if invalid index. */
    virtual const T& At(unsigned index) const
    { return fVars.at(index); }

    /**
     * Safe access, but slightly less efficient access to parameter.
     * @param name Look up name in list
     * @return The pointer at index position or NULL if invalid index. */
    virtual T& Get(const std::string& name)
    {	return At(Index(name)); }

    /**
     * Safe access, but slightly less efficient access to parameter.
     * @param name Look up name in list
     * @return The pointer at index position or NULL if invalid index. */
    virtual const T& Get(const std::string& name) const
    {	return At(Index(name)); }

    /**
     * Access to last pushed variable. */
    virtual T& Back()
    { return fVars.back(); }

    /**
     * Access to last pushed variable. */
    virtual const T& Back() const
    { return fVars.back(); }

    /**
     * Find index of parameter identified by name;
     * return Size() if name not found. */
    virtual unsigned Index(const std::string& name) const
    {
        for (unsigned int i = 0; i < fVars.size() ; ++i)
            if ( fVars[i].IsNamed(name) )
                return i;
        BCLog::OutWarning("BCVariableSet::Index : no variable named '" + name + "'");
        return fVars.size();
    }

    /**
     * Number of variables contained */
    virtual unsigned Size() const
    { return fVars.size(); }

    /**
     * Whether container is empty */
    virtual bool Empty() const
    { return fVars.empty(); }

    /**
     * Set number of bins for all parameters
     * @param nbins Number of bins. */
    virtual void SetNBins(unsigned nbins)
    {
        for (unsigned i = 0 ; i < fVars.size() ; ++i )
            fVars[i].SetNbins(nbins);
    }

    /**
     * Set precision for output of all parameters
     * @param n Number of significant digits for printing.*/
    virtual void SetPrecision(unsigned n)
    {
        for (unsigned i = 0 ; i < fVars.size() ; ++i )
            fVars[i].SetPrecision(n);
    }

    /**
     * Set fill-histograms flag for 1D and 2D histograms for all parameters.
     * @param flag Filling flag for both 1D and 2D histograms. */
    virtual void FillHistograms(bool flag)
    { FillHistograms(flag, flag); }

    /**
     * Set fill-histograms flag for 1D and 2D histograms for all parameters.
     * @param flag_1d Filling flag for 1D histograms.
     * @param flag_2d Filling flag for 2D histograms. */
    virtual void FillHistograms(bool flag_1d, bool flag_2d)
    { FillH1(flag_1d); FillH2(flag_2d); }

    /**
     * Set fill-histograms flag for all 1D histograms for all parameters.
     * @param flag Filling flag. */
    virtual void FillH1(bool flag)
    {
        for (unsigned i = 0 ; i < fVars.size() ; ++i )
            fVars[i].FillH1(flag);
    }

    /**
     * Set fill-histograms flag for all 2D histograms for all parameters.
     * @param flag Filling flag. */
    virtual void FillH2(bool flag)
    {
        for (unsigned i = 0 ; i < fVars.size() ; ++i )
            fVars[i].FillH2(flag);
    }

    /**
     * @return Length of longest parameter name. */
    virtual unsigned MaxNameLength() const
    { return fMaxNameLength; }

    /**
     * @return Volume of set. */
    virtual double Volume() const
    {
        double volume = 1;
        for (unsigned i = 0; i < fVars.size(); ++i)
            volume *= fVars[i].GetRangeWidth();
        if (volume < 0)
            return 0;
        return volume;
    }

    /**
     * @return Whether values are within limits of variables. */
    virtual bool IsWithinLimits(const std::vector<double>& x) const
    {
        if (x.size() != fVars.size())
            return false;
        for (unsigned i = 0; i < fVars.size(); ++i)
            if (!fVars[i].IsWithinLimits(x[i]))
                return false;
        return true;
    }

    /**
     * return positions in ranges of given values
     * from 0 (at lower limit) to 1 (at upper limit)
     * for each variable in the set.
     * @param x vector of values to report positions of.
     * @return vector of positions of values in ranges. */
    virtual std::vector<double> PositionInRange(const std::vector<double>& x) const
    {
        std::vector<double> p;
        for (unsigned i = 0; i < fVars.size(); ++i)
            p.push_back(fVars[i].PositionInRange(x[i]));
        return p;
    }

    /**
     * Translate from unit interval to values in variable ranges.
     * @param p vector of positions in the unit interval (0 = lower limit, 1 = upper limit). */
    virtual void ValueFromPositionInRange(std::vector<double>& p) const
    {
        if ( p.size() != fVars.size() )
            return;
        for (unsigned i = 0; i < fVars.size(); ++i)
            p[i] = fVars[i].ValueFromPositionInRange(p[i]);
    }

    /**
     * @return vector of range centers. */
    virtual std::vector<double> GetRangeCenters() const
    {
        std::vector<double> p;
        for (unsigned i = 0; i < fVars.size(); ++i)
            p.push_back(fVars[i].GetRangeCenter());
        return p;
    }

    /**
     * Get vector of uniformly distributed random values.
     * @param R Random number generator to use.
     * @return vector of random values uniformly distributed in variable ranges. */
    virtual std::vector<double> GetUniformRandomValues(TRandom* const R) const
    {
        std::vector<double> p;
        for (unsigned i = 0; i < fVars.size(); ++i)
            p.push_back(fVars[i].GetUniformRandomValue(R));
        return p;
    }

    /**
     * Print summary of variable set to logs. */
    virtual void PrintSummary() const
    {
        unsigned n = (int)log10(fVars.size()) + 1;
        for (unsigned i = 0; i < fVars.size(); ++i)
            BCLog::OutSummary(Form(" %*u) ", n, i) + fVars[i].OneLineSummary(false, fMaxNameLength));
    }
protected:
    /**
     * Vector of BCVariables that forms the set.
     * BCVariables are not owned by set, and are not deleted upon deletion of set. */
    std::vector<T> fVars;

    /**
     * Maximum length (in characters) of variable names. */
    unsigned fMaxNameLength;

};
#endif
