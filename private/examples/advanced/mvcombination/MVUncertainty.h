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
	
	// return the correlation matrix
	TMatrixD GetCorrelationMatrix()
		{ return fCorrelationMatrix; }; 

	// return the covariance matrix
	TMatrixD GetCovarianceMatrix()
		{ return fCovarianceMatrix; }; 

	// return the inverse of the covariance matrix
	TMatrixD GetInvCovarianceMatrix()
		{ return fInvCovarianceMatrix; }; 

	// return the flag if the uncertainty is active or not
	bool GetFlagActive()
	{ return fFlagActive; };

	// setters
	void SetCorrelationMatrix(const TMatrixD &matrix);

	// setters
	void SetCovarianceMatrix(const TMatrixT<double> &matrix);

 	// set flag if uncertainty is active for the combination
	void SetFlagActive(bool flag)
	{ fFlagActive = flag; }; 

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

	// debugKK: is that used?
	//	double* f;

	// flag: active in combination (true) or not (false)
	bool fFlagActive; 
};
// ---------------------------------------------------------

#endif
