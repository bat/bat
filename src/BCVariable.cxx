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
BCVariable::BCVariable() :
    fPrefix("Variable"),
    fLowerLimit(-std::numeric_limits<double>::infinity()),
    fUpperLimit(+std::numeric_limits<double>::infinity()),
    fPrecision(4),
    fFillH1(true),
    fFillH2(true),
    fNbins(100)
{
    SetName("parameter");
}

// ---------------------------------------------------------
BCVariable::BCVariable(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname, const std::string& unitstring) :
    fPrefix("Variable"),
    fLowerLimit(-std::numeric_limits<double>::infinity()),
    fUpperLimit(+std::numeric_limits<double>::infinity()),
    fPrecision(3),
    fLatexName(latexname),
    fUnitString(unitstring),
    fFillH1(true),
    fFillH2(true),
    fNbins(100)
{
    SetName(name);
    SetLimits(lowerlimit, upperlimit);
}

// ---------------------------------------------------------
void BCVariable::SetName(const std::string& name)
{
    fName = name;
    fSafeName = BCAux::SafeName(name);
}

// ---------------------------------------------------------
void BCVariable::SetLimits(double lowerlimit, double upperlimit)
{
    fLowerLimit = lowerlimit;
    fUpperLimit = upperlimit;
    if (fLowerLimit > fUpperLimit)
        BCLog::OutError(Form("BCVariable:SetLimits : lower limit (%f) is greater than upper limit (%f) for variable %s", fLowerLimit, fUpperLimit, fName.data()));
    if (BCAux::RangeType(fLowerLimit, fUpperLimit) == BCAux::kFiniteRange)
        CalculatePrecision();
}

// ---------------------------------------------------------
void BCVariable::CalculatePrecision(bool force)
{
    // need some extra digits to see where things become insignificant
    unsigned new_precision = 4 + ceil(-log10(fabs(fUpperLimit - fLowerLimit) / std::max(fabs(fUpperLimit), fabs(fLowerLimit))));
    if (force or new_precision > GetPrecision())
        SetPrecision(new_precision);
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
TH1* BCVariable::CreateH1(const std::string& name) const
{
    return new TH1D(name.data(),
                    (";" + GetLatexNameWithUnits() +
                     ";P(" + GetLatexName() + " | Data)").data(),
                    fNbins, fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
TH2* BCVariable::CreateH2(const std::string& name, const BCVariable& ordinate) const
{
    return new TH2D(name.data(),
                    (";" + GetLatexNameWithUnits() +
                     ";" + ordinate.GetLatexNameWithUnits() +
                     ";P(" + GetLatexName() + ", " + ordinate.GetLatexName() + " | Data)").data(),
                    fNbins, fLowerLimit, fUpperLimit,
                    ordinate.GetNbins(), ordinate.GetLowerLimit(), ordinate.GetUpperLimit());
}

// ---------------------------------------------------------
TH3* BCVariable::CreateH3(const std::string& name, const BCVariable& ordinate_y, const BCVariable& ordinate_z) const
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
