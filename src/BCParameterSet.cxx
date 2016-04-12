/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCParameterSet.h"
#include "BCLog.h"

#include <Math/Util.h>

// ---------------------------------------------------------
BCParameterSet::BCParameterSet()
    : BCVariableSet<BCParameter>()
{
}

// ---------------------------------------------------------
unsigned int BCParameterSet::GetNFixedParameters() const
{
    unsigned int n = 0;
    for (unsigned int i = 0; i < Size(); ++i)
        if (fVars[i].Fixed())
            ++n;
    return n;
}

// ---------------------------------------------------------
double BCParameterSet::Volume() const
{
    double volume = -1;
    for (unsigned i = 0; i < fVars.size(); ++i)
        if (!fVars[i].Fixed()) {
            if ( volume < 0 )
                volume = 1;
            volume *= fVars[i].GetRangeWidth();
        }
    if (volume < 0)
        return 0;
    return volume;
}

// ---------------------------------------------------------
void BCParameterSet::SetPriorConstantAll()
{
    for (unsigned i = 0; i < fVars.size(); ++i)
        fVars[i].SetPriorConstant();
}

// ---------------------------------------------------------
double BCParameterSet::GetLogPrior(const std::vector<double>& parameters) const
{
    if (parameters.size() < fVars.size()) {
        std::string msg("incorrect size of parameter set provided. Expected ");
        msg += ROOT::Math::Util::ToString(fVars.size()) + ", got ";
        msg += ROOT::Math::Util::ToString(parameters.size());
        BCLOG_ERROR(msg);
        return -std::numeric_limits<double>::infinity();
    }
    double logprob = 0;
    for (unsigned i = 0; i < fVars.size(); ++i)
        logprob += fVars[i].GetLogPrior(parameters[i]);
    return logprob;
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetFixedValues(bool include_unfixed) const
{
    std::vector<double> v (fVars.size(), std::numeric_limits<double>::infinity());
    for (unsigned i = 0; i < fVars.size(); ++i)
        if (include_unfixed or fVars[i].Fixed())
            v[i] = fVars[i].GetFixedValue();
    return v;
}

// ---------------------------------------------------------
bool BCParameterSet::ArePriorsSet(bool ignore_fixed) const
{
    for (unsigned i = 0; i < fVars.size(); ++i)
        if (ignore_fixed and fVars[i].Fixed())
            continue;
        else if (fVars[i].GetPrior() == NULL)
            return false;
    return true;
}

// ---------------------------------------------------------
bool BCParameterSet::IsWithinLimits(const std::vector<double>& x) const
{
    if (x.size() != fVars.size())
        return false;
    for (unsigned i = 0; i < fVars.size(); ++i)
        if (fVars[i].Fixed()) {
            if (fabs(x[i] - fVars[i].GetFixedValue()) > std::numeric_limits<double>::epsilon())
                return false;
        } else if (!fVars[i].IsWithinLimits(x[i]))
            return false;
    return true;
}

// ---------------------------------------------------------
bool BCParameterSet::IsAtFixedValues(const std::vector<double>& x) const
{
    if (x.size() < fVars.size())
        return false;
    for (unsigned i = 0; i < fVars.size(); ++i)
        if (fVars[i].Fixed() and (fVars[i].GetFixedValue() - x[i]) > std::numeric_limits<double>::epsilon())
            return false;
    return true;
}

// ---------------------------------------------------------
void BCParameterSet::ValueFromPositionInRange(std::vector<double>& p) const
{
    if ( p.size() != fVars.size() )
        return;
    for (unsigned i = 0; i < fVars.size(); ++i)
        p[i] = (fVars[i].Fixed()) ? fVars[i].GetFixedValue() : fVars[i].ValueFromPositionInRange(p[i]);
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetRangeCenters() const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fVars.size(); ++i)
        p.push_back((fVars[i].Fixed()) ? fVars[i].GetFixedValue() : fVars[i].GetRangeCenter());
    return p;
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetUniformRandomValues(TRandom* const R) const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fVars.size(); ++i)
        p.push_back((fVars[i].Fixed()) ? fVars[i].GetFixedValue() : fVars[i].GetUniformRandomValue(R));
    return p;
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetRandomValuesAccordingToPriors(TRandom* const R) const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fVars.size(); ++i)
        p.push_back((fVars[i].Fixed()) ? fVars[i].GetFixedValue() : fVars[i].GetRandomValueAccordingToPrior(R));
    return p;
}

// ---------------------------------------------------------
bool BCParameterSet::ApplyFixedValues(std::vector<double>& x) const
{
    if (x.size() != fVars.size())
        return false;
    for (unsigned i = 0; i < fVars.size(); ++i)
        if (fVars[i].Fixed())
            x[i] = fVars[i].GetFixedValue();
    return true;
}
