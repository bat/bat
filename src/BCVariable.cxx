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

// ---------------------------------------------------------

BCVariable::BCVariable():
	fPrefix("Variable"),
	fName("parameter"),
	fLowerLimit(0),
	fUpperLimit(1),
	fPrecision(3),
	fFillHistograms(true),
	fNbins(100)
{
}

// ---------------------------------------------------------

BCVariable::BCVariable(const char * name, double lowerlimit, double upperlimit, const char * latexname) :
	fPrefix("Variable"),
	fName(name),
	fLowerLimit(lowerlimit),
	fUpperLimit(upperlimit),
	fPrecision(3),
	fLatexName(latexname),
	fFillHistograms(true),
	fNbins(100)
{
}

// ---------------------------------------------------------

BCVariable::~BCVariable() {
}

// ---------------------------------------------------------

void BCVariable::PrintSummary() const {
	BCLog::OutSummary(Form("%s summary:",fPrefix.c_str()));
	BCLog::OutSummary(Form("%11s : %s", fPrefix.c_str(),fName.c_str()));
	BCLog::OutSummary(Form("Lower limit : % .*f", fPrecision,fLowerLimit));
	BCLog::OutSummary(Form("Upper limit : % .*f", fPrecision,fUpperLimit));
}

// ---------------------------------------------------------

TH1D * BCVariable::CreateH1(const char * name) {
	return new TH1D(name, TString::Format(";%s",fLatexName.c_str()), fNbins, fLowerLimit, fUpperLimit);
}

// ---------------------------------------------------------

TH2D * BCVariable::CreateH2(const char * name, BCVariable * ordinate) {
	return new TH2D(name, TString::Format(";%s;%s",fLatexName.c_str(),ordinate->GetLatexName().c_str()),
									fNbins, fLowerLimit, fUpperLimit,
									ordinate->GetNbins(),ordinate->GetLowerLimit(),ordinate->GetUpperLimit());
}
