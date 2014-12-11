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
void BCVariableSet::FillHistograms(bool flag) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> FillHistograms(flag);
}
