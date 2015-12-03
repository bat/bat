/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCObservable.h"

#include <string>


// ---------------------------------------------------------
BCObservable::BCObservable()
    : BCVariable()
    , fObservableValue(new double)
{
    fPrefix = "Observable";
}

// ---------------------------------------------------------
BCObservable::BCObservable(const BCObservable& other, bool share_pointer)
    : BCVariable(other)
    , fObservableValue(share_pointer ? other.fObservableValue : new double)
{
}

// ---------------------------------------------------------
BCObservable::BCObservable(const std::string& name, double lowerlimit, double upperlimit, const std::string& latexname, const std::string& unitstring)
    :	BCVariable(name, lowerlimit, upperlimit, latexname, unitstring)
    , fObservableValue(new double)
{
    fPrefix = "Observable";
}

// ---------------------------------------------------------
BCObservable::~BCObservable()
{
}

// ---------------------------------------------------------
BCObservable& BCObservable::operator=(BCObservable other)
{
    swap(*this, other);
    return *this;
}

// ---------------------------------------------------------
void swap(BCObservable& A, BCObservable& B)
{
    std::swap(static_cast<BCVariable&>(A), static_cast<BCVariable&>(B));
    std::swap(A.fObservableValue, B.fObservableValue);
}
