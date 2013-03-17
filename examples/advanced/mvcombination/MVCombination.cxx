#include "MVCombination.h"

#include <iomanip>

#include <TMath.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MVCombination::MVCombination() : BCModel("MVCombination")
{
}

// ---------------------------------------------------------
MVCombination::~MVCombination()
{
}

// ---------------------------------------------------------
void MVCombination::AddObservable(std::string name, double min, double max)
{
	AddParameter(name.c_str(), min, max); 
}

// ---------------------------------------------------------
void MVCombination::AddUncertainty(std::string name)
{
	MVUncertainty* u = new MVUncertainty(name);

	fUncertainties.push_back(u);
}

// ---------------------------------------------------------
void MVCombination::AddMeasurement(std::string name, std::string observable, double value, std::vector<double> uncertainties)
{
	// get index of the corresponding observable
	int index = GetIndexObservable(observable);
	
	// check if observable exists
	if (index < 0) {
		BCLog::OutWarning("Observable does not exist. Measurement was not added.");
		return;
	}

	MVMeasurement* m = new MVMeasurement(name);
	m->SetObservable(index);
	m->SetCentralValue(value);
	m->SetUncertainties(uncertainties);

	fMeasurements.push_back(m);

	fVectorObservable.push_back(index);

	int n = GetNMeasurements();
	fVectorMeasurements.ResizeTo(n);
	fVectorMeasurements[n-1]=value;
}

// ---------------------------------------------------------
double MVCombination::LogLikelihood(const std::vector<double> &parameters)
{
   double logprob = 0.;

	 int nmeas = GetNMeasurements();

	 // copy parameters into a vector
	 TVectorD observables(nmeas);

	 for (int i = 0; i < nmeas; ++i) {
		 observables[i] = parameters[fVectorObservable[i]];
	 }

	 TVectorD prod1 = observables - fVectorMeasurements;
	 TVectorD prod2 = fInvCovarianceMatrix * prod1;
	 double prod = prod1 * prod2;

	 logprob = -0.5 * prod - log(TMath::Power(2*TMath::Pi(), nmeas/2.) * sqrt(fDetCovariance));

   return logprob;
}

// ---------------------------------------------------------
double MVCombination::LogAPrioriProbability(const std::vector<double> &parameters)
{
   double logprob = 0.;

	 // debugKK
	 //	 if (parameters[0] + parameters[1] > 1)
	 //		 return -1e55;

   return logprob;
}

// ---------------------------------------------------------
void MVCombination::CalculateCorrelationMatrix(int index)
{
	MVUncertainty* u = fUncertainties.at(index);

	int n = GetNMeasurements();

	TMatrixD cov(n, n);

	TMatrixD corr = u->GetCorrelationMatrix();

	PrintMatrix(corr, "corr");

	for (int i = 0; i < n; ++i) 
		for (int j = 0; j < n; ++j) {
			double sigma_i = GetMeasurement(i)->GetUncertainty(index);
			double sigma_j = GetMeasurement(j)->GetUncertainty(index);
			
			double corr_ij = corr[i][j];
			cov[i][j] = corr_ij*sigma_i*sigma_j;
		}
	u->SetCovarianceMatrix(cov);
}

// ---------------------------------------------------------
void MVCombination::CalculateCovarianceMatrix()
{
	int n = GetNUncertainties();

	fCovarianceMatrix.Clear();
	fCovarianceMatrix.ResizeTo(GetNMeasurements(), GetNMeasurements());
	for (int i = 0; i < n; ++i) {
		MVUncertainty* u = fUncertainties.at(i);
		TMatrixD cov = u->GetCovarianceMatrix();
		fCovarianceMatrix += cov;
	}
	fInvCovarianceMatrix.Clear();
	fInvCovarianceMatrix.ResizeTo(fCovarianceMatrix);
	fInvCovarianceMatrix = fCovarianceMatrix;
	fInvCovarianceMatrix.Invert();

	fDetCovariance = fCovarianceMatrix.Determinant();
}

// ---------------------------------------------------------
void MVCombination::CalculateBLUE()
{
	int nmeas = GetNMeasurements();
	int nobs = GetNObservables();
	int nunc = GetNUncertainties();

	// calculate U matrix
	TMatrixD u(nmeas, nobs);

	for (int i = 0; i < nmeas; ++i)
		for (int j = 0; j < nobs; ++j) {
			MVMeasurement* m = GetMeasurement(i);
			if (m->GetObservable() == j)
				u[i][j] = 1;
			else 
				u[i][j] = 0;
		}

	// calculate weight matrix
	TMatrixD ut = u;
	ut.Transpose(ut); 

	TMatrixD m1 = ut * fInvCovarianceMatrix;
	TMatrixD m2 = m1 * u;
	m2.Invert();

	fBLUEWeights.Clear();
	fBLUEWeights.ResizeTo(nobs, nmeas);
	
	fBLUEWeights = m2*m1;

	// calculate central values
	fBLUECentral.Clear();
	fBLUECentral.ResizeTo(nobs);

	fBLUECentral = fBLUEWeights * fVectorMeasurements;

	// calculate covariance matrix
	fBLUECovarianceMatrix.Clear();
	fBLUECovarianceMatrix.ResizeTo(nobs, nobs);

	TMatrixD weightt = fBLUEWeights;
	weightt.Transpose(weightt); 

	fBLUECovarianceMatrix = fBLUEWeights * fCovarianceMatrix * weightt;

	// calculate uncertainties, covariance and correlation matrices for each uncertainty
	for (int i = 0; i < nunc; ++i) {

		// calculate covariance matrix
		MVUncertainty* u = GetUncertainty(i);
		TMatrixD cov = u->GetCovarianceMatrix();
		TMatrixD mat = fBLUEWeights * cov * weightt;
		fBLUECovarianceMatrices.push_back(mat);

		// calculate uncertainties
		TVectorD vec(nobs);
		for (int j = 0; j < nobs; ++j) {
			double sigma_j = sqrt(mat[j][j]);
			vec[j] = sigma_j;
		}

		TMatrixD mat2(nobs, nobs);
		for (int j = 0; j < nobs; ++j) 
			for (int k = 0; k < nobs; ++k) {
				mat2[j][k] = mat[j][k]/vec[j]/vec[k];
		}
		fBLUECorrelationMatrices.push_back(mat2);
	}

	// calculate uncertainties
	fBLUEUncertainties.Clear();
	fBLUEUncertainties.ResizeTo(nobs);

	for (int i = 0; i < nobs; ++i)
		fBLUEUncertainties[i] = sqrt(fBLUECovarianceMatrix[i][i]);
	
	// calculate correlation matrix
	fBLUECorrelationMatrix.Clear();
	fBLUECorrelationMatrix.ResizeTo(nobs, nobs);

	for (int i = 0; i < nobs; ++i)
		for (int j = 0; j < nobs; ++j) 
			fBLUECorrelationMatrix[i][j] = fBLUECovarianceMatrix[i][j] / fBLUEUncertainties[i] / fBLUEUncertainties[j];
}

// ---------------------------------------------------------
void MVCombination::PrintSummary()
{
	int nobs = GetNObservables();
	int nmeas = GetNMeasurements();
	int nunc = GetNUncertainties();

	std::cout << std::endl;
	std::cout << "Summary of combination" << std::endl;
	std::cout << "======================" << std::endl << std::endl;

	std::cout << "* Observables:" << std::endl;
	std::cout << "  Observable (range): " << std::endl;
	for (int i = 0; i < nobs; ++i) 
		std::cout << "  " << std::setiosflags(std::ios::left) << GetParameter(i)->GetName()
							<< " (" << GetParameter(i)->GetLowerLimit() << " - " << GetParameter(i)->GetUpperLimit() << ")" << std::endl;
	std::cout << std::endl;

	std::cout << "* Measurements:" << std::endl;
	std::cout << "  Measurement (observable): central value +- total uncertainty" << std::endl;
	for (int i = 0; i < nmeas; ++i) {
		MVMeasurement* m = GetMeasurement(i);
		std::cout << "  " << std::setiosflags(std::ios::left) << m->GetName() 
							<< std::setiosflags(std::ios::left) << " (" << GetParameter(m->GetObservable())->GetName() << ")" 
							<< ": " << std::setiosflags(std::ios::left) << std::setw(7) << m->GetCentralValue()
							<< " +- " << std::setiosflags(std::ios::left) << std::setw(7) << m->GetTotalUncertainty() << std::endl;
	}
	std::cout << std::endl;
	
	std::cout << "* Uncertainties:" << std::endl;
	for (int i = 0; i < nunc; ++i) {
		MVUncertainty* u = GetUncertainty(i);
		std::cout << "  " << u->GetName() << " " << "(correlation matrix)" << std::endl;
		TMatrixD mat = u->GetCorrelationMatrix();

		for (int j = 0; j < nmeas; ++j) {
			std::cout << "  ";
			for (int k = 0; k < nmeas; ++k) 
				std::cout << std::setw(7) << mat[j][k] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "* BLUE results: " << std::endl;
	std::cout << "  Observable: estimate +- uncertainty" << std::endl;
	for (int i = 0; i < nobs; ++i) {
		if (i < fBLUECentral.GetNoElements())
			std::cout << "  " << GetParameter(i)->GetName() << " " << fBLUECentral[i] << " +- " << fBLUEUncertainties[i] << std::endl;
	}
	std::cout << std::endl;

	std::cout << "  Correlation matrix" << std::endl;
	TMatrixD mat = fBLUECorrelationMatrix;
	for (int j = 0; j < nobs; ++j) {
		std::cout << "  ";
		for (int k = 0; k < nobs; ++k) 
			std::cout << std::setw(7) << mat[j][k] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	std::cout << "  Weights [%]:" <<std::endl;
	for (int j = 0; j < nmeas; ++j) {
		std::cout << "  " << GetMeasurement(j)->GetName() << " : ";
		for (int k = 0; k < nobs; ++k) 
			std::cout << std::setw(7) << fBLUEWeights[k][j]*100. << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;

}

// ---------------------------------------------------------
int MVCombination::GetIndexObservable(std::string name)
{
	int n = GetNObservables();

	// go through the list of parameters and compare strings
	for (int i = 0; i < n; ++i) {
		if (name == std::string(GetParameter(i)->GetName()))
			return i;
	}
	
	// return -1 if not in the list
	return -1;
}

// ---------------------------------------------------------
void MVCombination::PrintMatrix(TMatrixD &matrix, std::string name)
{
	int nrows = matrix.GetNrows();
	int ncols = matrix.GetNcols();

	std::cout << std::endl;
	std::cout << name.c_str() << " (" << nrows << "x" << ncols << "):" << std::endl;

	for (int i = 0; i < nrows; ++i) {
		for (int j = 0; j < ncols; ++j) 
			std::cout << std::setprecision(3) << std::setw(7) << matrix[i][j] << " ";
		std::cout << std::endl;
	}
}

// ---------------------------------------------------------
void MVCombination::PrintVector(TVectorD &vector, std::string name)
{
	int nrows = vector.GetNoElements();

	std::cout << std::endl;
	std::cout << name.c_str() << " (" << nrows << "):" << std::endl;

	for (int i = 0; i < nrows; ++i) {
		std::cout << std::setprecision(3) << std::setw(7) << vector[i] << " ";
		std::cout << std::endl;
	}
}

// ---------------------------------------------------------
