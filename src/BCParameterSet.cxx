/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCParameterSet.h"

#include "BCLog.h"
#include "BCParameter.h"

#include <TString.h>

// ---------------------------------------------------------
bool BCParameterSet::Add(BCParameter * parameter)
{
   if ( !parameter)
      return false;
   // check if parameter with same name exists
   for (unsigned int i = 0; i < fPars.size() ; ++i)
      if ( parameter->GetName() == fPars[i]->GetName() ) {
         BCLog::OutError(TString::Format("BCParameterSet::Add : Parameter with name %s exists already. ",
               parameter->GetName().data()));
         return false;
      }

   // add parameter to parameter container
   fPars.push_back(parameter);
   return true;
}

// ---------------------------------------------------------
void BCParameterSet::Clear(bool hard)
{
   if (hard) {
      for (unsigned int i = 0; i < fPars.size() ; ++i)
         delete fPars[i];
   }

   fPars.clear();
}

// ---------------------------------------------------------
bool BCParameterSet::ValidIndex(unsigned index, const std::string caller) const
{
   if (index >= fPars.size()) {
      BCLog::OutError(TString::Format("BCParameterSet::%s : Parameter index %u out of range", caller.c_str(), index));
      return false;
   }
   else
      return true;
}

// ---------------------------------------------------------
unsigned BCParameterSet::Index(const std::string & name) const
{
   for (unsigned int i=0; i < fPars.size() ; ++i)
      if (name == fPars[i]->GetName())
         return i;

   BCLog::OutWarning(TString::Format("BCParameterSet::Index: no parameter named '%s'", name.c_str()));
   return fPars.size();
}
