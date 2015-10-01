/*
 * Copyright (C) 2007-2015, the BAT core developer team
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
    : fData(nvariables, 0)
{
}

// ---------------------------------------------------------
BCDataPoint::BCDataPoint(const std::vector<double>& x)
    : fData(x)
{
}

// ---------------------------------------------------------
void BCDataPoint::SetValues(const std::vector<double>& values)
{
    // if vectors are same size, copy values
    if (values.size() == fData.size()) {
        fData = values;
    }
    // else give error
    else {
        BCLog::OutError("BCDataPoint::SetValues : Vectors have different ranges.");
        exit(1);
    }
}

// ---------------------------------------------------------
void BCDataPoint::PrintSummary(void (*output)(const std::string&)) const
{
    for (unsigned i = 0; i < fData.size(); ++i)
        output(Form("%u : %12.5g", i, fData[i]));
}
