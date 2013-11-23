/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPPDiagnostics.h"

#include <TH1D.h>
#include <TH2D.h>

#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
BCPPDiagnostics::BCPPDiagnostics() : BCPostProcessor()
{
}

// ---------------------------------------------------------
BCPPDiagnostics::~BCPPDiagnostics()
{
}

// ---------------------------------------------------------
BCH1D* BCPPDiagnostics::PlotLogProbability(std::string options)
{
  return 0;
}

// ---------------------------------------------------------
