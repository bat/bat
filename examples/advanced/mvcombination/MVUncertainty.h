#ifndef __MVUNCERTAINTY__H
#define __MVUNCERTAINTY__H

// ---------------------------------------------------------

#include <string>
#include <vector>
#include <iostream>

#include <TMatrixT.h>
#include <TMatrixD.h>

// ---------------------------------------------------------
class MVUncertainty
{
 public:
	
	// constructor
	// name: the name of the uncertainty
	MVUncertainty(std::string name);

	// destructor
	~MVUncertainty();

	// getters

	// return the name of the uncertainty
	std::string GetName()
		{ return fName; };
	
	TMatrixD GetCorrelationMatrix()
		{ return fCorrelationMatrix; }; 

	TMatrixD GetCovarianceMatrix()
		{ return fCovarianceMatrix; }; 

	TMatrixD GetInvCovarianceMatrix()
		{ return fInvCovarianceMatrix; }; 

	// setters
	void SetCorrelationMatrix(const TMatrixD &matrix);

	// setters
	void SetCovarianceMatrix(const TMatrixT<double> &matrix);

 private:

	// the name of the uncertainty
	std::string fName;

	// the symmetric correlation matrix; the number of columns and rows
	// is equal to the number of measurements
	TMatrixD fCorrelationMatrix; 

	// the symmetric covariance matrix; the number of columns and rows
	// is equal to the number of measurements
	TMatrixD fCovarianceMatrix; 

	// the inverse of the covariance matrix
	TMatrixD fInvCovarianceMatrix;

	double* f;

};
// ---------------------------------------------------------

#endif
