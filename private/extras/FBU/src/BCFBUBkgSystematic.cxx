/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCFBUBkgSystematic.h"

// ---------------------------------------------------------
BCFBUBkgSystematic::BCFBUBkgSystematic(std::string systname):BCFBUSystematic(systname.c_str())
 { 
   
 }


BCFBUBkgSystematic::BCFBUBkgSystematic(std::string systname, std::string samplename, TH1 *h_up, TH1 *h_down, TH2D *responseup, TH2D *responsedown):
  BCFBUSystematic(systname.c_str())
{ 
  fResponseUp = responseup;
  fResponseDown = responsedown;
  fHistogramUp[samplename]=h_up;
  fHistogramDown[samplename]=h_down;
}



// ---------------------------------------------------------
BCFBUBkgSystematic::~BCFBUBkgSystematic()
{}

// ---------------------------------------------------------
