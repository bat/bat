#include "MVUncertainty.h"

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
MVUncertainty::MVUncertainty(std::string name) : fName(name)
{
}

// ---------------------------------------------------------
MVUncertainty::~MVUncertainty()
{
}

// ---------------------------------------------------------
void MVUncertainty::SetCorrelationMatrix(const TMatrixT<double> &matrix)
{
	fCorrelationMatrix.ResizeTo(matrix);
	fCorrelationMatrix=matrix; 
};

// ---------------------------------------------------------
void MVUncertainty::SetCovarianceMatrix(const TMatrixT<double> &matrix)
{ 
	fCovarianceMatrix.ResizeTo(matrix);
	fCovarianceMatrix = matrix; 
	fInvCovarianceMatrix.ResizeTo(matrix);
	fInvCovarianceMatrix = fCovarianceMatrix;
	fInvCovarianceMatrix.Invert();
};

// ---------------------------------------------------------
