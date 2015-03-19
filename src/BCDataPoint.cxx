/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCDataPoint.h"

#include <TString.h>

#include <cstdlib>

// ---------------------------------------------------------
BCDataPoint::BCDataPoint(int nvariables)
{
    fData.assign(nvariables, 0.);
}

// ---------------------------------------------------------
BCDataPoint::BCDataPoint(const std::vector<double>& x)
{
    fData = x;
}

// ---------------------------------------------------------
BCDataPoint::BCDataPoint(const BCDataPoint& other)
{
    Copy(other);
}

// ---------------------------------------------------------
double BCDataPoint::GetValue(unsigned index) const
{
    double value;

    // check if index is in range. return value if true ...
    if (index < fData.size())
        value = fData[index];
    // ... or give error if not.
    else {
        // exit on error
        BCLog::OutError(Form("BCDataPoint::GetValue : Index %u out of range (%u to %u).", index, 0, unsigned(fData.size() - 1)));
        exit(1);
    }

    return value;
}

// ---------------------------------------------------------
void BCDataPoint::SetValue(unsigned index, double value)
{
    // check if index is in range. set value if true ...
    if (index < fData.size())
        fData[index] = value;
    // ... or give error if not.
    else {
        // exit on error
        BCLog::OutError(Form("BCDataPoint::SetValue : Index %u out of range (%u to %u).", index, 0 , unsigned(fData.size() - 1)));
        exit(1);
    }
}

// ---------------------------------------------------------
void BCDataPoint::SetValues(const std::vector<double>& values)
{
    // check if sizes are the same. if true, clear the data point and copy from
    // the vector passed to the method ...
    if (values.size() == fData.size()) {
        fData = values;
    }
    // ... or give error if the size if not the same.
    else {
        BCLog::OutError("BCDataPoint::SetValues : Vectors have different ranges.");
        exit(1);
    }
}

// ---------------------------------------------------------
void BCDataPoint::Dump(void (*output)(const char*)) const
{
    for (unsigned i = 0; i < fData.size(); ++i)
        output(Form("%u : %12.5g", i, fData[i]));
}
