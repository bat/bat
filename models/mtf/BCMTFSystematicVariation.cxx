/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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
