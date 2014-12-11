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
#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>

// ---------------------------------------------------------

BCVariable::BCVariable():
	fPrefix("Variable"),
	fLowerLimit(0),
	fUpperLimit(1),
	fPrecision(3),
	fFillHistograms(true),
	fNbins(100)
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
	, fFillHistograms(other.fFillHistograms)
	, fNbins(other.fNbins)
{
}

// ---------------------------------------------------------

BCVariable::BCVariable(const char * name, double lowerlimit, double upperlimit, const char * latexname) :
	fPrefix("Variable"),
	fLowerLimit(lowerlimit),
	fUpperLimit(upperlimit),
	fLatexName(latexname),
	fFillHistograms(true),
	fNbins(100)
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

TH1D * BCVariable::CreateH1(const char * name) const {
	return new TH1D(name, TString::Format(";%s",GetLatexName().c_str()), fNbins, fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------

TH2D * BCVariable::CreateH2(const char * name, BCVariable * ordinate) const {
	return new TH2D(name, TString::Format(";%s;%s",GetLatexName().c_str(),ordinate->GetLatexName().c_str()),
									fNbins, fLowerLimit, fUpperLimit,
									ordinate->GetNbins(),ordinate->GetLowerLimit(),ordinate->GetUpperLimit());
}
