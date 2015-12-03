/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMTFProcess.h"

// ---------------------------------------------------------
BCMTFProcess::BCMTFProcess(const std::string& name)
    : fHistogramColor(-1)
    , fHistogramFillStyle(-1)
    , fHistogramLineStyle(-1)
{
    SetName(name);
}

// ---------------------------------------------------------
BCMTFProcess::~BCMTFProcess()
{}

// ---------------------------------------------------------
