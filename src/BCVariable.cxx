/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCVariable.h"
#include "BCLog.h"
#include "BCAux.h"

#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TRandom.h>

#include <limits>

// ---------------------------------------------------------
BCVariable::BCVariable()
    :	fPrefix("Variable")
    , fLowerLimit(-std::numeric_limits<double>::infinity())
    , fUpperLimit(+std::numeric_limits<double>::infinity())
    , fPrecision(3)
    , fFillH1(true)
    , fFillH2(true)
    , fNbins(100)
{
    SetName("parameter");
}

// ---------------------------------------------------------
BCVariable::BCVariable(const BCVariable& other)
    : fPrefix(other.fPrefix)
    , fName(other.fName)
    , fSafeName(other.fSafeName)
    , fLowerLimit(other.fLowerLimit)
    , fUpperLimit(other.fUpperLimit)
    , fPrecision(other.fPrecision)
    , fLatexName(other.fLatexName)
    , fUnitString(other.fUnitString)
    , fFillH1(other.fFillH1)
    , fFillH2(other.fFillH2)
    , fNbins(other.fNbins)
{
}

// ---------------------------------------------------------
BCVariable::BCVariable(std::string name, double lowerlimit, double upperlimit, std::string latexname, std::string unitstring)
    :	fPrefix("Variable")
    , fLowerLimit(lowerlimit)
    , fUpperLimit(upperlimit)
    , fPrecision(3)
    , fLatexName(latexname)
    , fUnitString(unitstring)
    , fFillH1(true)
    , fFillH2(true)
    , fNbins(100)
{
    SetName(name);
    CalculatePrecision();
}

// ---------------------------------------------------------
BCVariable::~BCVariable()
{
}

// ---------------------------------------------------------
void BCVariable::SetName(std::string name)
{
    fName = name;
    fSafeName = BCAux::SafeName(name);
}

// ---------------------------------------------------------
void BCVariable::SetLimits(double lowerlimit, double upperlimit)
{
    fLowerLimit = lowerlimit;
    fUpperLimit = upperlimit;
    if (BCAux::RangeType(fLowerLimit, fUpperLimit) == BCAux::kFiniteRange)
        CalculatePrecision();
}

// ---------------------------------------------------------
void BCVariable::CalculatePrecision(bool force)
{
    double new_precision = ceil(-log10(2.*fabs(fUpperLimit - fLowerLimit) / (fabs(fUpperLimit) + fabs(fLowerLimit))));
    if (force or new_precision > GetPrecision())
        SetPrecision(force);
}
// ---------------------------------------------------------
bool BCVariable::IsAtLimit(double value) const
{
    if ( (value - fLowerLimit) * (value - fLowerLimit) / fLowerLimit / fLowerLimit <= 1e-10
            or
            (value - fUpperLimit) * (value - fUpperLimit) / fUpperLimit / fUpperLimit <= 1e-10 )
        return true;

    return false;
}

// ---------------------------------------------------------
void BCVariable::PrintSummary() const
{
    BCLog::OutSummary(fPrefix + " summary:");
    BCLog::OutSummary(Form("%11s : %s", fPrefix.c_str(), fName.c_str()));
    BCLog::OutSummary(Form("Lower limit : % .*f", fPrecision, fLowerLimit));
    BCLog::OutSummary(Form("Upper limit : % .*f", fPrecision, fUpperLimit));
}

// ---------------------------------------------------------
std::string BCVariable::OneLineSummary(bool print_prefix, int name_length) const
{
    if (name_length < 0)
        name_length = fName.size();
    if (print_prefix)
        return std::string(Form("%s \"%*s\" : [%.*g, %.*g]", fPrefix.data(), name_length, fName.data(), fPrecision, fLowerLimit, fPrecision, fUpperLimit));
    else
        return std::string(Form("%-*s : [%.*g, %.*g]", name_length, fName.data(), fPrecision, fLowerLimit, fPrecision, fUpperLimit));
}

// ---------------------------------------------------------
TH1* BCVariable::CreateH1(std::string name) const
{
    return new TH1D(name.data(),
                    (";" + GetLatexNameWithUnits() +
                     ";P(" + GetLatexName() + " | Data)").data(),
                    fNbins, fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
TH2* BCVariable::CreateH2(std::string name, const BCVariable& ordinate) const
{
    return new TH2D(name.data(),
                    (";" + GetLatexNameWithUnits() +
                     ";" + ordinate.GetLatexNameWithUnits() +
                     ";P(" + GetLatexName() + ", " + ordinate.GetLatexName() + " | Data)").data(),
                    fNbins, fLowerLimit, fUpperLimit,
                    ordinate.GetNbins(), ordinate.GetLowerLimit(), ordinate.GetUpperLimit());
}

// ---------------------------------------------------------
TH3* BCVariable::CreateH3(std::string name, const BCVariable& ordinate_y, const BCVariable& ordinate_z) const
{
    return new TH3D(name.data(),
                    (";" + GetLatexNameWithUnits() +
                     ";" + ordinate_y.GetLatexNameWithUnits() +
                     ";" + ordinate_z.GetLatexNameWithUnits() +
                     ";P(" + GetLatexName() + ", " + ordinate_y.GetLatexName() + ", " + ordinate_z.GetLatexName() + " | Data)").data(),
                    fNbins, fLowerLimit, fUpperLimit,
                    ordinate_y.GetNbins(), ordinate_y.GetLowerLimit(), ordinate_y.GetUpperLimit(),
                    ordinate_z.GetNbins(), ordinate_z.GetLowerLimit(), ordinate_z.GetUpperLimit());
}

// ---------------------------------------------------------
double BCVariable::GetUniformRandomValue(TRandom* const R) const
{
    if (!R)
        return std::numeric_limits<double>::quiet_NaN();
    return ValueFromPositionInRange(R->Rndm());
}
