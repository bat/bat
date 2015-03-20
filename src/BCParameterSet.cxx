/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCParameterSet.h"
#include "BCParameter.h"
#include "BCLog.h"

// ---------------------------------------------------------
BCParameterSet& BCParameterSet::operator=(const BCParameterSet& rhs)
{
    for (unsigned i = 0; i < rhs.fPars.size(); ++i)
        Add(new BCParameter(*((BCParameter*)(fPars[i]))));
    return *this;
}

// ---------------------------------------------------------
bool BCParameterSet::Add(BCParameter* par)
{
    return BCVariableSet::Add(dynamic_cast<BCVariable*>(par));
}

// ---------------------------------------------------------
unsigned int BCParameterSet::GetNFixedParameters() const
{
    unsigned int n = 0;
    for (unsigned int i = 0; i < Size(); ++i)
        if (dynamic_cast<BCParameter*>(fPars[i])->Fixed())
            ++n;
    return n;
}

// ---------------------------------------------------------
double BCParameterSet::Volume() const
{
    double volume = -1;
    for (unsigned i = 0; i < fPars.size(); ++i)
        if ( ! ((BCParameter*)fPars[i])->Fixed() ) {
            if ( volume < 0 )
                volume = 1;
            volume *= fPars[i]->GetRangeWidth();
        }
    if (volume < 0)
        return 0;
    return volume;
}

// ---------------------------------------------------------
void BCParameterSet::SetPriorConstantAll()
{
    for (unsigned i = 0; i < fPars.size(); ++i)
        dynamic_cast<BCParameter*>(fPars[i])->SetPriorConstant();
}

// ---------------------------------------------------------
double BCParameterSet::GetLogPrior(const std::vector<double>& parameters) const
{
    if (parameters.size() != fPars.size()) {
        BCLog::OutError("BCParameterSet::GetLogPrior : incorrect size of parameter set provided.");
        return -std::numeric_limits<double>::infinity();
    }
    double logprob = 0;
    for (unsigned i = 0; i < fPars.size(); ++i)
        logprob += dynamic_cast<BCParameter*>(fPars[i])->GetLogPrior(parameters[i]);
    return logprob;
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetFixedValues(bool include_unfixed) const
{
    std::vector<double> v (fPars.size(), std::numeric_limits<double>::infinity());
    for (unsigned i = 0; i < fPars.size(); ++i)
        if (include_unfixed or dynamic_cast<BCParameter*>(fPars[i])->Fixed())
            v[i] = dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue();
    return v;
}

// ---------------------------------------------------------
bool BCParameterSet::ArePriorsSet(bool ignore_fixed) const
{
    for (unsigned i = 0; i < fPars.size(); ++i)
        if (ignore_fixed and dynamic_cast<BCParameter*>(fPars[i])->Fixed())
            continue;
        else if (dynamic_cast<BCParameter*>(fPars[i])->GetPrior() == NULL)
            return false;
    return true;
}

// ---------------------------------------------------------
bool BCParameterSet::IsWithinLimits(const std::vector<double>& x) const
{
    if (x.size() != fPars.size())
        return false;
    for (unsigned i = 0; i < fPars.size(); ++i) 
        if (dynamic_cast<BCParameter*>(fPars[i])->Fixed()) {
            if (fabs(x[i] - dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue()) > std::numeric_limits<double>::epsilon())
                return false;
        } else if (!fPars[i]->IsWithinLimits(x[i]))
            return false;
    return true;
}

// ---------------------------------------------------------
bool BCParameterSet::IsAtFixedValues(const std::vector<double>& x) const
{
    if (x.size() != fPars.size())
        return false;
    for (unsigned i = 0; i < fPars.size(); ++i)
        if (!dynamic_cast<BCParameter*>(fPars[i])->Fixed())
            continue;
        else if (dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue() - x[i] > std::numeric_limits<double>::epsilon())
            return false;
    return true;
}

// ---------------------------------------------------------
void BCParameterSet::ValueFromPositionInRange(std::vector<double>& p) const
{
    if ( p.size() != fPars.size() )
        return;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p[i] = (dynamic_cast<BCParameter*>(fPars[i])->Fixed()) ? dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue() : fPars[i]->ValueFromPositionInRange(p[i]);
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetRangeCenters() const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p.push_back((dynamic_cast<BCParameter*>(fPars[i])->Fixed()) ? dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue() : fPars[i]->GetRangeCenter());
    return p;
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetUniformRandomValues(TRandom* const R) const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p.push_back((dynamic_cast<BCParameter*>(fPars[i])->Fixed()) ? dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue() : fPars[i]->GetUniformRandomValue(R));
    return p;
}

// ---------------------------------------------------------
std::vector<double> BCParameterSet::GetRandomValuesAccordingToPriors(TRandom* const R) const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p.push_back((dynamic_cast<BCParameter*>(fPars[i])->Fixed()) ? dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue() : dynamic_cast<BCParameter*>(fPars[i])->GetRandomValueAccordingToPrior(R));
    return p;
}

// ---------------------------------------------------------
bool BCParameterSet::ApplyFixedValues(std::vector<double>& x) const
{
    if (x.size() != fPars.size())
        return false;
    for (unsigned i = 0; i < fPars.size(); ++i)
        if (dynamic_cast<BCParameter*>(fPars[i])->Fixed())
            x[i] = dynamic_cast<BCParameter*>(fPars[i])->GetFixedValue();
    return true;
}
