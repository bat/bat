#include "ParameterSummary.h"

#include <BAT/BCH1D.h>

#include <iostream>

// ---------------------------------------------------------
ParameterSummary::ParameterSummary(const char* name) :
	fName(name)
{
}

// ---------------------------------------------------------
ParameterSummary::~ParameterSummary()
{
}

// ---------------------------------------------------------
void ParameterSummary::Summarize(BCH1D* hist)
{
	SetMode( hist->GetMode() ); 
	SetMean( hist->GetMean() ); 
	SetMedian( hist->GetMedian() ); 
	SetQuantile5( hist->GetQuantile(0.05) );
	SetQuantile10( hist->GetQuantile(0.10) );
	SetQuantile16( hist->GetQuantile(0.16) );
	SetQuantile84( hist->GetQuantile(0.84) );
	SetQuantile90( hist->GetQuantile(0.90) );
	SetQuantile95( hist->GetQuantile(0.95) );
	SetRMS( hist->GetRMS() );
}

// ---------------------------------------------------------
