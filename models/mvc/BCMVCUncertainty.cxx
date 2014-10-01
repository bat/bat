/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMVCUncertainty.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
BCMVCUncertainty::BCMVCUncertainty(std::string name) : fName(name)
						     , fFlagActive(true)
{
}

// ---------------------------------------------------------
BCMVCUncertainty::~BCMVCUncertainty()
{
}

// ---------------------------------------------------------
void BCMVCUncertainty::SetCorrelationMatrix(const TMatrixT<double> &matrix)
{
  fCorrelationMatrix.ResizeTo(matrix);
  fCorrelationMatrix=matrix;
}

// ---------------------------------------------------------
void BCMVCUncertainty::SetCovarianceMatrix(const TMatrixT<double> &matrix)
{
  fCovarianceMatrix.ResizeTo(matrix);
  fCovarianceMatrix = matrix;
  fInvCovarianceMatrix.ResizeTo(matrix);
  fInvCovarianceMatrix = fCovarianceMatrix;
  fInvCovarianceMatrix.Invert();
}

// ---------------------------------------------------------
