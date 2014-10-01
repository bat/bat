/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <TH1D.h>

#include "BCMTFSystematicVariation.h"

// ---------------------------------------------------------
BCMTFSystematicVariation::BCMTFSystematicVariation(const char * channelname, const char * systematicname, int nprocesses)
{
   fChannelName = channelname;
   fSystematicName = systematicname;
   for (int i = 0; i < nprocesses; ++i) {
      fHistogramUpContainer.push_back(0);
      fHistogramDownContainer.push_back(0);
   }
}

// ---------------------------------------------------------
BCMTFSystematicVariation::~BCMTFSystematicVariation()
{}

// ---------------------------------------------------------
