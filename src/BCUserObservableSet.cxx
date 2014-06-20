/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCUserObservableSet.h"

#include "BCLog.h"
#include "BCUserObservable.h"

#include <TString.h>

#include <algorithm>

// ---------------------------------------------------------
BCUserObservableSet::BCUserObservableSet() 
	: fMaxNameLength(0)
{
}

// ---------------------------------------------------------
bool BCUserObservableSet::Add(BCUserObservable * observable)
{
   if ( !observable)
      return false;
   // check if observable with same name exists
   for (unsigned int i = 0; i < fPars.size() ; ++i)
      if ( observable->GetName() == fPars[i]->GetName() ) {
         BCLog::OutError(TString::Format("BCUserObservableSet::Add : Observable with name %s exists already. ",
               observable->GetName().data()));
         return false;
      }

   // add observable to observable container
   fPars.push_back(observable);
	 fMaxNameLength = std::max(fMaxNameLength, (unsigned)observable->GetName().length());
   return true;
}

// ---------------------------------------------------------
void BCUserObservableSet::Clear(bool hard)
{
   if (hard) {
      for (unsigned int i = 0; i < fPars.size() ; ++i)
         delete fPars[i];
   }

   fPars.clear();
}

// ---------------------------------------------------------
bool BCUserObservableSet::ValidIndex(unsigned index, const std::string caller) const
{
   if (index >= fPars.size()) {
      BCLog::OutError(TString::Format("BCUserObservableSet::%s : Observable index %u out of range", caller.c_str(), index));
      return false;
   }
   else
      return true;
}

// ---------------------------------------------------------
unsigned BCUserObservableSet::Index(const std::string & name) const
{
   for (unsigned int i=0; i < fPars.size() ; ++i)
      if (name == fPars[i]->GetName())
         return i;

   BCLog::OutWarning(TString::Format("BCUserObservableSet::Index: no observable named '%s'", name.c_str()));
   return fPars.size();
}


// ---------------------------------------------------------
void BCUserObservableSet::SetNBins(unsigned nbins) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> SetNbins(nbins);
}

// ---------------------------------------------------------
void BCUserObservableSet::SetPrecision(unsigned n) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> SetPrecision(n);
}

// ---------------------------------------------------------
void BCUserObservableSet::FillHistograms(bool flag) {
	for (unsigned i = 0 ; i < fPars.size() ; ++i )
		fPars[i] -> FillHistograms(flag);
}
