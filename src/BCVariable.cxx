/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include "BCVariable.h"
#include "BCLog.h"

#include <TString.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TRandom.h>

#include <algorithm>
#include <limits>

// ---------------------------------------------------------
BCVariable::BCVariable()
	:	fPrefix("Variable")
	, fLowerLimit(0)
	, fUpperLimit(1)
	, fPrecision(3)
	, fFillH1(true)
	, fFillH2(true)
	, fNbins(100)
{
	SetName("parameter");
}

// ---------------------------------------------------------
BCVariable::BCVariable(const BCVariable & other) 
	: fPrefix(other.fPrefix)
	, fName(other.fName)
	, fSafeName(other.fSafeName)
	, fLowerLimit(other.fLowerLimit)
	, fUpperLimit(other.fUpperLimit)
	, fPrecision(other.fPrecision)
	, fLatexName(other.fLatexName)
	, fFillH1(other.fFillH1)
	, fFillH2(other.fFillH2)
	, fNbins(other.fNbins)
{
}

// ---------------------------------------------------------
BCVariable::BCVariable(const char * name, double lowerlimit, double upperlimit, const char * latexname)
	:	fPrefix("Variable")
	, fLowerLimit(lowerlimit)
	, fUpperLimit(upperlimit)
	, fLatexName(latexname)
	, fFillH1(true)
	, fFillH2(true)
	, fNbins(100)
{
	SetName(name);
	CalculatePrecision();
}

// ---------------------------------------------------------
BCVariable::~BCVariable() {
}

// ---------------------------------------------------------
void BCVariable::SetName(const char * name) {
	fName = name;
	fSafeName = name;
	fSafeName.erase(std::remove_if(fSafeName.begin(),fSafeName.end(),::isspace),fSafeName.end());
}

// ---------------------------------------------------------
bool BCVariable::IsAtLimit(double value) const {
	if ( (value-fLowerLimit)*(value-fLowerLimit)/fLowerLimit/fLowerLimit <= 1e-10
			 or
			 (value-fUpperLimit)*(value-fUpperLimit)/fUpperLimit/fUpperLimit <= 1e-10 )
		return true;
	
	return false;
}

// ---------------------------------------------------------
void BCVariable::PrintSummary() const {
	BCLog::OutSummary(Form("%s summary:",fPrefix.c_str()));
	BCLog::OutSummary(Form("%11s : %s", fPrefix.c_str(),fName.c_str()));
	BCLog::OutSummary(Form("Lower limit : % .*f", fPrecision,fLowerLimit));
	BCLog::OutSummary(Form("Upper limit : % .*f", fPrecision,fUpperLimit));
}

// ---------------------------------------------------------
std::string BCVariable::OneLineSummary() const {
	return std::string(Form("%s \"%s\" : [%.*f, %.*f]",fPrefix.data(),fName.data(),fPrecision,fLowerLimit,fPrecision,fUpperLimit));
}

// ---------------------------------------------------------
TH1 * BCVariable::CreateH1(const char * name) const {
	return new TH1D(name, TString::Format(";%s",GetLatexName().c_str()), fNbins, fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------
TH2 * BCVariable::CreateH2(const char * name, BCVariable * ordinate) const {
	return new TH2D(name, TString::Format(";%s;%s",GetLatexName().c_str(),ordinate->GetLatexName().c_str()),
									fNbins, fLowerLimit, fUpperLimit,
									ordinate->GetNbins(),ordinate->GetLowerLimit(),ordinate->GetUpperLimit());
}

// ---------------------------------------------------------
double BCVariable::GetUniformRandomValue(TRandom * const R) const {
	if (!R)
		return std::numeric_limits<double>::quiet_NaN();
	return ValueFromPositionInRange(R->Rndm());
}
