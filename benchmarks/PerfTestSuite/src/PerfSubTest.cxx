/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <cmath>

#include "include/PerfSubTest.h"

//______________________________________________________________________________
PerfSubTest::PerfSubTest(std::string name) 
	: fTestValue(-1)
	, fStatusRegion(std::vector<double>(3))
	, fStatusUnknown(true)
	, fStatusOff(false)
{
	// define default status regions 
	fStatusRegion[0] =  1.0; // test value range for good
	fStatusRegion[1] =  2.0; // test value range for good
	fStatusRegion[2] =  3.0; // test value range for good

	// set name
	fName = name; 
}
	
//______________________________________________________________________________
PerfSubTest::~PerfSubTest()
{
}
	
//______________________________________________________________________________
std::string PerfSubTest::ToString(PerfSubTest::Status status)
{
	switch (status) 
		{
		case PerfSubTest::kGood : 
			return std::string("good"); 

		case PerfSubTest::kFlawed : 
			return std::string("flawed"); 

		case PerfSubTest::kBad : 
			return std::string("bad"); 

		case PerfSubTest::kFatal : 
			return std::string("fatal"); 

		case PerfSubTest::kUnknown : 
			return std::string("unknown"); 

		case PerfSubTest::kOff : 
			return std::string("off"); 

		default :
			return std::string("-"); 				
		}
}

//______________________________________________________________________________
std::string PerfSubTest::ToStringHTML(PerfSubTest::Status status)
{
	switch (status) 
		{
		case PerfSubTest::kGood : 
			return std::string("<font color=\"#01DF01\">good</font>"); 

		case PerfSubTest::kFlawed : 
			return std::string("<font color=\"#FF00FF\">flawed</font>"); 

		case PerfSubTest::kBad : 
			return std::string("<font color=\"#FF8000\">bad</font>"); 

		case PerfSubTest::kFatal : 
			return std::string("<font color=\"#FF0000\">fatal</font>"); 

		case PerfSubTest::kUnknown : 
			return std::string("<font color=\"#0174DF\">unknown</font>"); 

		case PerfSubTest::kOff : 
			return std::string("<font color=\"#A4A4A4\">off</font>"); 

		default :
			return std::string("-"); 				
		}
}

//______________________________________________________________________________
void PerfSubTest::SetStatusRegion(PerfSubTest::Status status, double delta) 
{
	switch (status)
		{
		case PerfSubTest::kGood : 
			fStatusRegion[0] = delta;
			break; 
		case PerfSubTest::kFlawed : 
			fStatusRegion[1] = delta;
			break; 
		case PerfSubTest::kBad : 
			fStatusRegion[2] = delta;
			break;
		case PerfSubTest::kFatal : 
			break;
		default:
			std::cout << "Status " << ToString(status) << " has no region." << std::endl;
		}

	// set status
	GetStatus();
}

//______________________________________________________________________________
PerfSubTest::Status PerfSubTest::GetStatus() 
{
	// return status off
	if (fStatusOff) {
		return PerfSubTest::kOff; 
	}

	// calculate difference 
	double diff = fabs(fTestValue - fTargetValue); 

	// define status flag
	fStatusUnknown = false;

	// check for "good", "flawed" and "bad" 
	if (diff < fStatusRegion.at(0)) {
		return PerfSubTest::kGood;
	}

	else if (diff < fStatusRegion.at(1)) {
		return PerfSubTest::kFlawed;
	}

	else if (diff < fStatusRegion.at(2)) {
		return PerfSubTest::kBad;
	}

	else 
		return PerfSubTest::kFatal; 

}

//______________________________________________________________________________
double PerfSubTest::GetStatusRegion(PerfSubTest::Status status) 
{
	switch (status)
		{
		case PerfSubTest::kGood : 
			return fStatusRegion[0]; 

		case PerfSubTest::kFlawed : 
			return fStatusRegion[1]; 

		case PerfSubTest::kBad : 
			return fStatusRegion[2]; 
				
		case PerfSubTest::kFatal : 
			return 1e99; 

		default:
			std::cout << "Could not find status region." << std::endl;
			return -1; 
		}	
			
	return -1; 
}

//______________________________________________________________________________
	
