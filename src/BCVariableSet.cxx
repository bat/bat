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
{
}

// ---------------------------------------------------------
BCVariableSet & BCVariableSet::operator=(const BCVariableSet & rhs) {
	for (unsigned i=0; i<rhs.fPars.size(); ++i)
		Add(new BCVariable(*fPars[i]));
	return *this;
}

// ---------------------------------------------------------
bool BCVariableSet::Add(BCVariable * parameter)
{
   if ( !parameter)
      return false;
   // check if parameter with same name or same safe name exists
   for (unsigned int i = 0; i < fPars.size() ; ++i) {
		 if ( parameter->GetName() == fPars[i]->GetName() ) {
			 BCLog::OutError(TString::Format("BCVariableSet::Add : Variable with name %s exists already.", parameter->GetName().data()));
			 return false;
		 }
		 if ( parameter->GetSafeName() == fPars[i]->GetSafeName() ) {
			 BCLog::OutError(TString::Format("BCVariableSet::Add : Variable with safe name %s exists already.", parameter->GetSafeName().data()));
			 return false;
		 }
	 }
   // add parameter to parameter container
   fPars.push_back(parameter);
	 fMaxNameLength = std::max(fMaxNameLength, (unsigned)parameter->GetName().length());
	 // extend existing entries of fFillH2
	 for (unsigned i=0; i<fFillH2.size(); ++i)
		 fFillH2[i].resize(fPars.size(),false);
	 // extend fFillH2
	 fFillH2.resize(fPars.size(),std::vector<bool>(fPars.size(),false));
   return true;
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
   if (index >= fPars.size()) {
      BCLog::OutError(TString::Format("BCVariableSet::%s : Variable index %u out of range", caller.c_str(), index));
      return false;
   }
   else
      return true;
}

// ---------------------------------------------------------
unsigned BCVariableSet::Index(const std::string & name) const {
	for (unsigned int i=0; i < fPars.size() ; ++i)
		if ( fPars[i]->IsNamed(name) )
			return i;
	
	BCLog::OutWarning(TString::Format("BCVariableSet::Index: no parameter named '%s'", name.c_str()));
	return fPars.size();
}


// ---------------------------------------------------------
void BCVariableSet::SetNBins(unsigned nbins) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> SetNbins(nbins);
}

// ---------------------------------------------------------
void BCVariableSet::SetPrecision(unsigned n) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> SetPrecision(n);
}

// ---------------------------------------------------------
void BCVariableSet::FillH1(bool flag) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> FillH1(flag);
}

// ---------------------------------------------------------
void BCVariableSet::FillH2(bool flag) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> FillH2(flag);
}

// ---------------------------------------------------------
bool BCVariableSet::FillH2(unsigned x, unsigned y) const {
	return x<fFillH2.size() and y<fFillH2[x].size() // x and y within range
		and (fFillH2[x][y]														// explicitely designated to be filled
				 or (fPars[x]->FillH2() and fPars[y]->FillH2())); // or implicitely
}

// ---------------------------------------------------------
double BCVariableSet::Volume() const {
	double volume = 1;
	for (unsigned i=0; i<fPars.size(); ++i)
		volume *= fPars[i] -> GetRangeWidth();
	if (volume<0)
		return 0;
	return volume;
}

// ---------------------------------------------------------
bool BCVariableSet::IsWithinLimits(const std::vector<double> & x) const {
	if (x.size() != fPars.size())
		return false;
	for (unsigned i=0; i<fPars.size(); ++i)
		if (!fPars[i]->IsWithinLimits(x[i]))
			return false;
	return true;
}

// ---------------------------------------------------------
std::vector<double> BCVariableSet::PositionInRange(const std::vector<double> & x) const {
	std::vector<double> p;
	for (unsigned i=0; i<fPars.size(); ++i)
		p.push_back(fPars[i]->PositionInRange(x[i]));
	return p;
}

// ---------------------------------------------------------
void BCVariableSet::ValueFromPositionInRange(std::vector<double> & p) const {
	if ( p.size() != fPars.size() )
		return;
	for (unsigned i=0; i<fPars.size(); ++i)
		p[i] = fPars[i]->ValueFromPositionInRange(p[i]);
}

// ---------------------------------------------------------
std::vector<double> BCVariableSet::GetRangeCenters() const {
	std::vector<double> p;
	for (unsigned i=0; i<fPars.size(); ++i)
		p.push_back(fPars[i]->GetRangeCenter());
	return p;
}

// ---------------------------------------------------------
std::vector<double> BCVariableSet::GetUniformRandomValues(TRandom * const R) const {
	std::vector<double> p;
	for (unsigned i=0; i<fPars.size(); ++i)
		p.push_back(fPars[i]->GetUniformRandomValue(R));
	return p;
}
