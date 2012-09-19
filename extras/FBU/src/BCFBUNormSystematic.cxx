/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCFBUNormSystematic.h"

// ---------------------------------------------------------
BCFBUNormSystematic::BCFBUNormSystematic(std::string name, double uncertainty):
  BCFBUSystematic(name.c_str())
{
  fNormUncertainty = uncertainty;
}


// ---------------------------------------------------------
BCFBUNormSystematic::~BCFBUNormSystematic()
{}

// ---------------------------------------------------------
