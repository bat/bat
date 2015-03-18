/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCVariableSet.h"

#include "BCLog.h"
#include "BCVariable.h"

#include <TString.h>

#include <algorithm>

// ---------------------------------------------------------
BCVariableSet::BCVariableSet()
    : fMaxNameLength(0)
    , fPartnerSet(0)
{
}

// ---------------------------------------------------------
BCVariableSet& BCVariableSet::operator=(const BCVariableSet& rhs)
{
    for (unsigned i = 0; i < rhs.fPars.size(); ++i)
        Add(new BCVariable(*fPars[i]));
    return *this;
}

// ---------------------------------------------------------
bool BCVariableSet::Add(BCVariable* parameter)
{
    if ( !parameter)
        return false;

    // check if variable with same name or same safe name exists
    for (unsigned int i = 0; i < fPars.size() ; ++i)
        if ( parameter->IsNamed(fPars[i]->GetName()) ) {
            BCLog::OutError(TString::Format("BCVariableSet::Add : Variable with name %s exists already.", parameter->GetName().data()));
            return false;
        } else if ( parameter->IsSafeNamed(fPars[i]->GetSafeName()) ) {
            BCLog::OutError(TString::Format("BCVariableSet::Add : Variable with safe name %s exists already.", parameter->GetSafeName().data()));
            return false;
        }

    // check if variable with same name or safe name exists in partnered set
    for (unsigned int i = 0; i < fPartnerSet->Size() ; ++i)
        if ( parameter->IsNamed(fPartnerSet->Get(i)->GetName()) ) {
            BCLog::OutError(TString::Format("BCVariableSet::Add : Variable with name %s exists already in partner set.", parameter->GetName().data()));
            return false;
        } else if ( parameter->IsSafeNamed(fPartnerSet->Get(i)->GetSafeName()) ) {
            BCLog::OutError(TString::Format("BCVariableSet::Add : Variable with safe name %s exists already in partner set.", parameter->GetSafeName().data()));
            return false;
        }

    // add parameter to parameter container
    fPars.push_back(parameter);
    fMaxNameLength = std::max(fMaxNameLength, (unsigned)parameter->GetName().length());

    // extend existing entries of fFillH2
    for (unsigned i = 0; i < fFillH2.size(); ++i)
        fFillH2[i].resize(Size(), false);

    // extend fFillH2
    fFillH2.resize(Size(), std::vector<bool>(Size(), false));

    if (fPartnerSet) {
        // extend fFillH2Partner
        fFillH2Partner.resize(Size(), std::vector<bool>(fPartnerSet->Size(), false));
        // extend partner's fFillH2Partner
        for (unsigned i = 0; i < fPartnerSet->fFillH2Partner.size(); ++i)
            fPartnerSet->fFillH2Partner[i].resize(Size(), false);
    }

    return true;
}

// ---------------------------------------------------------
void BCVariableSet::SetPartner(BCVariableSet* const set)
{
    fPartnerSet = set;
    if (!fPartnerSet)
        fFillH2Partner.clear();
    else
        fFillH2Partner.assign(Size(), std::vector<bool>(fPartnerSet->Size(), false));
}

// ---------------------------------------------------------
void BCVariableSet::Clear(bool hard)
{
    if (hard)
        for (unsigned int i = 0; i < fPars.size() ; ++i)
            delete fPars[i];
    fPars.clear();
}

// ---------------------------------------------------------
bool BCVariableSet::ValidIndex(unsigned index, const std::string caller) const
{
    if (index < fPars.size())
        return true;
    BCLog::OutError(TString::Format("BCVariableSet::%s : Variable index %u out of range", caller.c_str(), index));
    return false;
}

// ---------------------------------------------------------
unsigned BCVariableSet::Index(const std::string& name) const
{
    for (unsigned int i = 0; i < fPars.size() ; ++i)
        if ( fPars[i]->IsNamed(name) )
            return i;
    BCLog::OutWarning(TString::Format("BCVariableSet::Index: no parameter named '%s'", name.c_str()));
    return fPars.size();
}


// ---------------------------------------------------------
void BCVariableSet::SetNBins(unsigned nbins)
{
    for (unsigned i = 0 ; i < fPars.size() ; ++i )
        fPars[i]->SetNbins(nbins);
}

// ---------------------------------------------------------
void BCVariableSet::SetPrecision(unsigned n)
{
    for (unsigned i = 0 ; i < fPars.size() ; ++i )
        fPars[i]->SetPrecision(n);
}

// ---------------------------------------------------------
void BCVariableSet::FillH1(bool flag)
{
    for (unsigned i = 0 ; i < fPars.size() ; ++i )
        fPars[i]->FillH1(flag);
}

// ---------------------------------------------------------
void BCVariableSet::FillH2(bool flag)
{
    for (unsigned i = 0 ; i < fPars.size() ; ++i )
        fPars[i]->FillH2(flag);
}

// ---------------------------------------------------------
void BCVariableSet::FillH2(std::string x, std::string y, bool flag)
{
    unsigned index_x = Index(x);

    if (index_x >= Size()) {
        BCLog::OutError("BCVariableSet::FillH2 : Called with abscissa not in set. If abcissa belongs to partner set, call its function instead!");
        return;
    }

    unsigned index_y = Index(y);

    // referencing members of this set
    if (index_y < Size()) {
        FillH2(index_x, index_y, flag);
        return;
    }

    // otherwise possibly involving partner set's variables
    // if no partner set, return
    if (!fPartnerSet)
        return;
    // try setting filling of partner-involved H2:
    FillH2Partner(index_x, fPartnerSet->Index(y), flag);
}

// ---------------------------------------------------------
void BCVariableSet::FillAllH2(unsigned index, int axis, bool flag, bool include_partner_set)
{
    for (unsigned i = 0; i < fPars.size(); ++i)
        if (axis == 0)								// index as abscissa
            FillH2(index, i, flag);
        else												// index as ordinate
            FillH2(i, index, flag);

    // set for partner-involved H2's as well
    if (include_partner_set and fPartnerSet) {
        for (unsigned i = 0; i < fPartnerSet->Size(); ++i)
            if (axis == 0)							// index as abscissa
                FillH2Partner(index, i, flag);
            else
                fPartnerSet->FillH2Partner(i, index, flag);
    }
}

// ---------------------------------------------------------
bool BCVariableSet::FillH2(unsigned x, unsigned y) const
{
    return x != y																		 // don't fill H2's with identical axes
           and x < fFillH2.size() and y < fFillH2[x].size() // x and y within range
           and (fFillH2[x][y]					// explicitely designated to be filled
                or (y > x and fPars[x]->FillH2() and fPars[y]->FillH2())); // or implicitely
}

// ---------------------------------------------------------
bool BCVariableSet::FillH2Partner(unsigned x, unsigned y) const
{
    return x < fFillH2Partner.size() and y < fFillH2Partner[x].size() // x and y within range
           and (fFillH2Partner[x][y]		// explicitely designated to be filled
                or (fPars[x]->FillH2() and fPartnerSet and fPartnerSet->Get(y) and fPartnerSet->Get(y)->FillH2())); // or implicitely
}

// ---------------------------------------------------------
bool BCVariableSet::FillH2(std::string x, std::string y) const
{
    unsigned index_x = Index(x);

    if (index_x >= Size()) {
        BCLog::OutError("BCVariableSet::FillH2 : Called with abscissa not in set. If abcissa belongs to partner set, call its function instead!");
        return false;
    }

    unsigned index_y = Index(y);

    // y is in this set
    if (index_y < Size())
        return FillH2(index_x, index_y);

    // else try partner set
    return FillH2Partner(index_x, fPartnerSet->Index(y));
}

// ---------------------------------------------------------
double BCVariableSet::Volume() const
{
    double volume = 1;
    for (unsigned i = 0; i < fPars.size(); ++i)
        volume *= fPars[i]->GetRangeWidth();
    if (volume < 0)
        return 0;
    return volume;
}

// ---------------------------------------------------------
bool BCVariableSet::IsWithinLimits(const std::vector<double>& x) const
{
    if (x.size() != fPars.size())
        return false;
    for (unsigned i = 0; i < fPars.size(); ++i)
        if (!fPars[i]->IsWithinLimits(x[i]))
            return false;
    return true;
}

// ---------------------------------------------------------
std::vector<double> BCVariableSet::PositionInRange(const std::vector<double>& x) const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p.push_back(fPars[i]->PositionInRange(x[i]));
    return p;
}

// ---------------------------------------------------------
void BCVariableSet::ValueFromPositionInRange(std::vector<double>& p) const
{
    if ( p.size() != fPars.size() )
        return;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p[i] = fPars[i]->ValueFromPositionInRange(p[i]);
}

// ---------------------------------------------------------
std::vector<double> BCVariableSet::GetRangeCenters() const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p.push_back(fPars[i]->GetRangeCenter());
    return p;
}

// ---------------------------------------------------------
std::vector<double> BCVariableSet::GetUniformRandomValues(TRandom* const R) const
{
    std::vector<double> p;
    for (unsigned i = 0; i < fPars.size(); ++i)
        p.push_back(fPars[i]->GetUniformRandomValue(R));
    return p;
}

