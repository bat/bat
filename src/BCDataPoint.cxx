/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BAT/BCDataPoint.h"

#include "BAT/BCLog.h"

#include <TString.h>

// ---------------------------------------------------------

BCDataPoint::BCDataPoint(int nvariables)
{
	// assign the specified number of variables to the data
	// point and fill with zero
	fData.assign(nvariables, 0.);
}

// ---------------------------------------------------------

BCDataPoint::BCDataPoint(std::vector<double> x)
{
	// copy all values of x to the data point
	for (std::vector<double>::const_iterator it = x.begin(); it != x.end(); ++it)
		fData.push_back(*it);
}

// ---------------------------------------------------------

BCDataPoint::~BCDataPoint()
{}

// ---------------------------------------------------------

double BCDataPoint::GetValue(int index)
{
// this is not good at all
// -1 can be ok value
// and one can turn off warnings
// so if index is out of range the program should stop
	double value = -1.0;

	// check if index is in range. return value if true ...
	if (index >= 0 && index < int(fData.size()))
		value = fData[index];
	// ... or give warning if not.
	else
      BCLog::Out(BCLog::warning, BCLog::warning, TString::Format(
            "BCDataPoint::GetValue : Index %d out of range (%d to %d).", index,
            0, fData.size()-1));

	return value;
}

// ---------------------------------------------------------

void BCDataPoint::SetValue(int index, double value)
{
// this is not good at all
// -1 can be ok value
// and one can turn off warnings
// so if index is out of range the program should stop

	// check if index is in range. set value if true ...
	if (index >= 0 && index < int(fData.size()))
		fData[index] = value;
	// ... or give warning if not.
	else
		BCLog::Out(BCLog::warning, BCLog::warning,TString::Format(
				"BCDataPoint::SetValue : Index %d out of range (%d to %d).",
				index, 0 ,fData.size()-1));
}

// ---------------------------------------------------------

void BCDataPoint::SetValues(std::vector <double> values)
{
	// check if sizes are the same. if true, clear the data point and copy from
	// the vector passed to the method ...
	if (values.size() == fData.size())
	{
		fData.clear();
		for (std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it)
			fData.push_back(*it);
	}
	// ... or give warning if the size if not the same.
	else
		BCLog::Out(BCLog::warning, BCLog::warning,"BCDataPoint::SetValues : Vectors have different ranges.");
}

// ---------------------------------------------------------


