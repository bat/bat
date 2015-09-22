/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCObservableSet.h"

// ---------------------------------------------------------
BCObservableSet::BCObservableSet()
    : BCVariableSet<BCObservable>()
{
}

// ---------------------------------------------------------
BCObservableSet::BCObservableSet(const BCObservableSet& other, bool share_pointers)
    : BCVariableSet<BCObservable>(other)
{
    if (!share_pointers)
        for (unsigned i = 0; i < fVars.size(); ++i)
            fVars[i].SetValueLocation(new double);
}

// ---------------------------------------------------------
BCObservableSet::~BCObservableSet()
{
}
