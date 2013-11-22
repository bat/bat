/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPPMarginalize.h"

#include <TH1D.h>
#include <TH2D.h>

#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
BCPPMarginalize::BCPPMarginalize() : BCPostProcessor()
{
}

// ---------------------------------------------------------
BCPPMarginalize::~BCPPMarginalize()
{
}

// ---------------------------------------------------------
BCH1D* BCPPMarginalize::BuildMarginalized1D(int parindex, int nbins, double parmin, double parmax)
{
  // check if parameter index is within range
  if (parindex < 0 || parindex >= GetNParameters()) {
    BCLog::OutWarning("BCPPMarginalize::BuildMarginalized1D. Index not within range.");
    return 0;
  }

  // define histogram
  TH1D* hist = new TH1D("hist", Form(";p_{%i};p(p_{%i}|data)", parindex, parindex), nbins, parmin, parmax);
  hist->SetStats(kFALSE);

  // loop over all trees
  for (int i = 0; i < fNTrees; ++i) {
    // fill histogram
    fTrees.at(i)->Draw(Form("Parameter%i>>+hist", parindex), "Phase==2");
  }

  // define BAT histogram
  BCH1D* bhist = new BCH1D(hist);

  return bhist;
}

// ---------------------------------------------------------
BCH2D* BCPPMarginalize::BuildMarginalized2D(int parindex1, int nbins1, double parmin1, double parmax1, int parindex2, int nbins2, double parmin2, double parmax2)
{
  // check if parameter index is within range
  if (parindex1 < 0 || parindex1 >= GetNParameters() ||
      parindex2 < 0 || parindex2 >= GetNParameters()) {
    BCLog::OutWarning("BCPPMarginalize::BuildMarginalized1D. Index not within range.");
    return 0;
  }

  // define histogram
  TH2D* hist = new TH2D("hist", Form(";p_{%i};p_{%i}", parindex1, parindex2), nbins1, parmin1, parmax1, nbins2, parmin2, parmax2);
  hist->SetStats(kFALSE);

  // loop over all trees
  for (int i = 0; i < fNTrees; ++i) {
    // fill histogram
    fTrees.at(i)->Draw(Form("Parameter%i:Parameter%i>>+hist",parindex1, parindex2), "Phase==2");
  }

  // define BAT histogram
  BCH2D* bhist = new BCH2D(hist);

  return bhist;
}

// ---------------------------------------------------------
