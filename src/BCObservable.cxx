/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCObservable.h"

// ---------------------------------------------------------
BCObservable::BCObservable()
    : BCVariable()
    , fObservableValue(0)
{
    fPrefix = "Observable";
}

// ---------------------------------------------------------
BCObservable::BCObservable(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname, const std::string& unitstring)
    :	BCVariable(name, lowerlimit, upperlimit, latexname, unitstring)
    , fObservableValue(0)
{
    fPrefix = "Observable";
}
