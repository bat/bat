/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <TH1D.h>

#include "BCMTFTemplate.h"

// ---------------------------------------------------------
BCMTFTemplate::BCMTFTemplate(const char * channelname, const char * processname)
 : fEfficiency(0)
 , fHistogram(0)
 , fNBins(0)
{
   fChannelName = channelname;
   fProcessName = processname;
   fFunctionContainer = new std::vector<TF1 *>(0);
}

// ---------------------------------------------------------
BCMTFTemplate::~BCMTFTemplate()
{
   // debugKK
   //   if (fHistogram)
   //      delete fHistogram;
}

// ---------------------------------------------------------
void BCMTFTemplate::SetHistogram(TH1D * hist)
{
   fHistogram = hist;
   if (hist)
      fNBins = fHistogram->GetNbinsX();
}

// ---------------------------------------------------------
void BCMTFTemplate::SetFunctionContainer(std::vector<TF1 *> * funccont, int nbins)
{
   fFunctionContainer = funccont;
   fNBins = nbins;
}

// ---------------------------------------------------------
