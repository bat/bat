/*
 * Copyright (C) 2007-2018, the BAT core developer team
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
BCParameter::BCParameter() :
    BCVariable(),
    fFixed(false),
    fFixedValue(std::numeric_limits<double>::infinity()),
    fPrior(NULL)
{
    fPrefix = "Parameter";
}

// ---------------------------------------------------------
BCParameter::BCParameter(const BCParameter& other) :
    BCVariable(other),
    fFixed(other.fFixed),
    fFixedValue(other.fFixedValue),
    fPrior(NULL)
{
    if (other.fPrior)
        SetPrior(other.fPrior->Clone());
}

// ---------------------------------------------------------
BCParameter::BCParameter(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname, const std::string& unitstring) :
    BCVariable(name, lowerlimit, upperlimit, latexname, unitstring),
    fFixed(false),
    fFixedValue(std::numeric_limits<double>::infinity()),
    fPrior(NULL)
{
    fPrefix = "Parameter";
    if (lowerlimit == upperlimit)
        Fix(lowerlimit);
}

// ---------------------------------------------------------
BCParameter::~BCParameter()
{
    delete fPrior;
}

// ---------------------------------------------------------
BCParameter& BCParameter::operator=(BCParameter other)
{
    swap(*this, other);
    return *this;
}

// ---------------------------------------------------------
void swap(BCParameter& A, BCParameter& B)
{
    std::swap(static_cast<BCVariable&>(A), static_cast<BCVariable&>(B));
    std::swap(A.fFixed,      B.fFixed);
    std::swap(A.fFixedValue, B.fFixedValue);
    std::swap(A.fPrior,      B.fPrior);
}

// ---------------------------------------------------------
void BCParameter::SetLimits(double lowerlimit, double upperlimit)
{
    BCVariable::SetLimits(lowerlimit, upperlimit);
    if (BCAux::RangeType(fLowerLimit, fUpperLimit) == BCAux::kFiniteRange && fPrior)
        fPrior->GetFunction().SetRange(fLowerLimit, fUpperLimit);
    if (lowerlimit == upperlimit)
        fFixedValue = lowerlimit;
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
    if (fPrior)
        fPrior->SetFunctionRange(fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
void BCParameter::SetPriorConstant()
{
    SetPrior(new BCConstantPrior(GetRangeWidth()));
}

// ---------------------------------------------------------
std::string BCParameter::OneLineSummary(bool print_prefix, int name_length) const
{
    if (!Fixed())
        return BCVariable::OneLineSummary(print_prefix, name_length);
    return BCVariable::OneLineSummary(print_prefix, name_length) + Form(" (fixed at %.*f)", GetPrecision(), GetFixedValue());
}
