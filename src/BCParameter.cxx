/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCParameter.h"
#include "BCPrior.h"
#include "BCConstantPrior.h"
#include "BCLog.h"
#include "BCAux.h"

#include <string>
#include <limits>

#include <TString.h>
#include <TRandom.h>
#include <TMath.h>
#include <TF1.h>


// ---------------------------------------------------------
BCParameter::BCParameter()
    :	BCVariable()
    , fFixed(false)
    , fFixedValue(std::numeric_limits<double>::infinity())
    , fPrior(NULL)
{
    fPrefix = "Parameter";
}

// ---------------------------------------------------------
BCParameter::BCParameter(const BCParameter& other)
    : BCVariable(other)
    , fFixed(other.fFixed)
    , fFixedValue(other.fFixedValue)
    , fPrior(NULL)
{
    if (other.fPrior)
        SetPrior(other.fPrior->Clone());
}

// ---------------------------------------------------------
BCParameter::BCParameter(std::string name, double lowerlimit, double upperlimit, std::string latexname, std::string unitstring)
    : BCVariable(name, lowerlimit, upperlimit, latexname, unitstring)
    , fFixed(false)
    , fFixedValue(std::numeric_limits<double>::infinity())
    , fPrior(NULL)
{
    fPrefix = "Parameter";
}

// ---------------------------------------------------------
BCParameter::~BCParameter()
{
    delete fPrior;
}

// ---------------------------------------------------------
void BCParameter::SetLimits(double lowerlimit, double upperlimit)
{
    BCVariable::SetLimits(lowerlimit, upperlimit);
    if (BCAux::RangeType(fLowerLimit, fUpperLimit) == BCAux::kFiniteRange and
            fPrior and fPrior->GetFunction())
        fPrior->GetFunction()->SetRange(fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
double BCParameter::GetLogPrior(double x) const
{
    if ( fFixed )
        return 0;
    if (fPrior)
        return fPrior->GetLogPrior(x);
    return -std::numeric_limits<double>::infinity();
}

// ---------------------------------------------------------
double BCParameter::GetPriorMode() const
{
    if (!fPrior)
        return std::numeric_limits<double>::quiet_NaN();
    return fPrior->GetMode(fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
double BCParameter::GetPriorMean() const
{
    if (!fPrior)
        return std::numeric_limits<double>::quiet_NaN();
    return fPrior->GetMean(fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
double BCParameter::GetPriorVariance() const
{
    if (!fPrior)
        return std::numeric_limits<double>::quiet_NaN();
    return fPrior->GetVariance(fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
double BCParameter::GetRandomValueAccordingToPrior(TRandom* const R) const
{
    if (!fPrior) {
        BCLog::OutError("BCParameter::GetRandomValueAccordingToPrior : no prior specified.");
        return std::numeric_limits<double>::quiet_NaN();
    }

    return fPrior->GetRandomValue(fLowerLimit, fUpperLimit, R);
}

// ---------------------------------------------------------
void BCParameter::SetPrior(BCPrior* const prior)
{
    delete fPrior;
    fPrior = prior;
    if (fPrior and fPrior->GetFunction())
        fPrior->GetFunction()->SetRange(fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
void BCParameter::SetPriorConstant()
{
    SetPrior(new BCConstantPrior(GetRangeWidth()));
}

// ---------------------------------------------------------
std::string BCParameter::OneLineSummary() const
{
    if (!Fixed())
        return BCVariable::OneLineSummary();
    return std::string(TString::Format("%s (fixed at %.*f)", BCVariable::OneLineSummary().data(), GetPrecision(), GetFixedValue()));
}

