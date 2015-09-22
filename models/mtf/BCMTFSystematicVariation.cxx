/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <TH1D.h>

#include "BCMTFSystematicVariation.h"

// ---------------------------------------------------------
BCMTFSystematicVariation::BCMTFSystematicVariation(int nprocesses)
    : fHistogramUpContainer(nprocesses, (TH1D*)0),
      fHistogramDownContainer(nprocesses, (TH1D*)0)
{
}

// ---------------------------------------------------------
BCMTFSystematicVariation::~BCMTFSystematicVariation()
{}
