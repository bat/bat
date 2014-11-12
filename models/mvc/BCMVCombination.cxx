/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMVCombination.h"

#include "BCMVCMeasurement.h"
#include "BCMVCObservable.h"
#include "BCMVCUncertainty.h"

#include "../../BAT/BCLog.h"
#include "../../BAT/BCMath.h"
#include "../../BAT/BCParameter.h"

#include <fstream>
#include <iomanip>
#include <iostream>

// ---------------------------------------------------------
BCMVCombination::BCMVCombination()
  : BCModel("BCMVCombination")
  , fDetCovariance(0)
  , fNObservables(0)
  , fNNuisanceCorrelation(0)
{
}

// ---------------------------------------------------------
BCMVCombination::~BCMVCombination()
{
   int nuncertainties = GetNUncertainties();
   int nmeasurements = GetNMeasurements();

   for (int i = 0; i < nuncertainties; ++i) {
      BCMVCUncertainty* u = GetUncertainty(i);
      delete u;
   }
   fUncertainties.clear();

   for (int i = 0; i < nmeasurements; ++i) {
      BCMVCMeasurement* m = GetMeasurement(i);
      delete m;
   }
   fMeasurements.clear();
}

// ---------------------------------------------------------
void BCMVCombination::AddObservable(std::string name, double min, double max)
{
   // check if observable exists already
   int index = GetIndexObservable(name);

   if (index >= 0)
      return;

   BCMVCObservable* obs = new BCMVCObservable();
   obs->SetName(name);
   obs->SetMinMax(min, max);
   fObservables.push_back(obs);

   fNObservables++;

   AddParameter(name.c_str(), min, max);

   SetPriorConstant(name.c_str());
}

// ---------------------------------------------------------
void BCMVCombination::AddUncertainty(std::string name)
{
   BCMVCUncertainty* u = new BCMVCUncertainty(name);

   fUncertainties.push_back(u);
}

// ---------------------------------------------------------
void BCMVCombination::AddMeasurement(std::string name, std::string observable, double value, std::vector<double> uncertainties)
{
   // get index of the corresponding observable
   int index = GetIndexObservable(observable);

   // check if observable exists
   if (index < 0) {
      BCLog::OutWarning(Form("BCMVCombination::AddMeasurement. Observable \"%s\" does not exist. Measurement was not added.", observable.c_str()));
      return;
   }

   BCMVCMeasurement* m = new BCMVCMeasurement(name);
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
double BCMVCombination::LogLikelihood(const std::vector<double> &parameters)
{
   double logprob = 0.;

   if (fNNuisanceCorrelation > 0) {
      CalculateCovarianceMatrix(parameters);
      if (!PositiveDefinite())
         return -1e90;
   }

   int nmeas = GetNActiveMeasurements();

   // copy parameters into a vector
   TVectorD observables(nmeas);

   for (int i = 0; i < nmeas; ++i) {
      observables[i] = parameters[fVectorActiveObservable[i]];
   }

   TVectorD prod1 = observables - fVectorActiveMeasurements;
   TVectorD prod2 = fInvCovarianceMatrix * prod1;
   double prod = prod1 * prod2;

   logprob = -0.5 * prod - log(TMath::Power(2*TMath::Pi(), nmeas/2.) * sqrt(fabs(fDetCovariance)));

   return logprob;
}

// ---------------------------------------------------------
int BCMVCombination::ReadInput(std::string filename)
{
   // open input file
   std::ifstream infile;
   infile.open(filename.c_str(), std::ifstream::in);

   // check if file is open
   if (!infile.is_open()) {
      BCLog::OutWarning(Form("BCMVCombination::ReadInput. Could not open input file %s.", filename.c_str()));
      return 0;
   }

   int nobservables;
   int nmeasurements;
   int nuncertainties;
   int nnuisance;

   infile >> nobservables >> nmeasurements >> nuncertainties >> nnuisance;

   std::vector<std::string> observable_names;

   for (int i = 0; i < nobservables; ++i) {
      std::string name;
      double min;
      double max;
      infile >> name >> min >> max;

      // add observable
      AddObservable(name.c_str(), min, max);
   }

   for (int i = 0; i < nuncertainties; ++i) {
      std::string name;
      infile >> name;

      // add uncertainty
      AddUncertainty(name);
   }

   for (int i = 0; i < nmeasurements; ++i) {
      std::string name;
      std::string observable;
      double central;
      std::vector<double> uncertainties(0);

      infile >> name;
      infile >> observable;
      infile >> central;

      for (int j = 0; j < nuncertainties; ++j) {
         double uncertainty;
         infile >> uncertainty;
         uncertainties.push_back(uncertainty);
      }

      // add measurement
      AddMeasurement(name, observable, central, uncertainties);
   }

   for (int i = 0; i < nuncertainties; ++i) {
      TMatrixD mat(nmeasurements, nmeasurements);

      for (int j = 0; j < nmeasurements; ++j)
         for (int k = 0; k < nmeasurements; ++k) {
            double corr;
            infile >> corr;
            mat[j][k] = corr;
         }

      // set correlation matrix
      GetUncertainty(i)->SetCorrelationMatrix(mat);
   }

   for (int i = 0; i < nnuisance; ++i) {
      std::string uncertainty;
      std::string measurement1;
      std::string observable1;
      std::string measurement2;
      std::string observable2;
      std::string parname;
      double min;
      double max;
      double pre;
      std::string priorshape;

      infile >> uncertainty >> measurement1 >> observable1 >> measurement2 >> observable2 >> parname;

      // check if parameter exists already
      int index = -1;

      for (unsigned int i = 0; i < GetNParameters(); i++)
         if (parname.c_str() == GetParameter(i)->GetName())
            index = i;

      if (index >= 0)
         infile >> pre;

      else {
         // read properties of parameter
         infile >> min >> max >> priorshape;

         // add nuisance parameter
         AddParameter(parname.c_str(), min, max);

         // set pre-factor
         pre = 1;

         // set index
         index = GetNParameters()-1;

         if (priorshape == "flat") {
            SetPriorConstant(parname.c_str());
         }
         else if (priorshape == "gauss") {
            double mean;
            double std;

            infile >> mean >> std;

            SetPriorGauss(parname.c_str(), mean, std);
         }
         else {
            BCLog::OutWarning(Form("BCMVCombination::ReadInput. Unknown prior shape %s.", priorshape.c_str()));
         }
      }

      // increase counter of nuisance parametera
      fNNuisanceCorrelation++;

      NuisanceParameter p;
      p.index_uncertainty  = GetIndexUncertainty(uncertainty);
      p.index_measurement1 = GetIndexMeasurement(measurement1, observable1);
      p.index_measurement2 = GetIndexMeasurement(measurement2, observable2);
      p.index_rhoparameter = index;
      p.pre = pre;

      fNuisanceCorrelation.push_back(p);
   }

   // close input file
   infile.close();

   PrepareAnalysis();

   // no error
   return 1;
}

// ---------------------------------------------------------
void BCMVCombination::PrepareAnalysis()
{
   int nuncertainties = GetNUncertainties();

   // prepare analysis
   for (int i = 0; i < nuncertainties; ++i)
      CalculateCorrelationMatrix(i);

   CalculateCovarianceMatrix();

   if (!PositiveDefinite()) {
      BCLog::OutWarning("BCMVCombination::PrepareAnalysis. Covariance matrix is not positive definite.");
   }

   CalculateHelperVectors();

}

// ---------------------------------------------------------
int BCMVCombination::GetNActiveMeasurements()
{
   int n = GetNMeasurements();

   int counter = 0;

   for (int i = 0; i < n; ++i) {
      BCMVCMeasurement* m = GetMeasurement(i);
      if (m->GetFlagActive())
         counter++;
   }
   return counter;
}

// ---------------------------------------------------------
void BCMVCombination::CalculateHelperVectors()
{
   int nmeasurements  = GetNMeasurements();

   fVectorMeasurements.Clear();
   fVectorMeasurements.ResizeTo(nmeasurements);
   fVectorObservable.clear();

   for (int i = 0; i < nmeasurements; ++i) {
      BCMVCMeasurement* m = GetMeasurement(i);
      fVectorMeasurements[i] =  m->GetCentralValue();
      fVectorObservable.push_back(m->GetObservable());
   }

   int nactive = GetNActiveMeasurements();

   fVectorActiveMeasurements.Clear();
   fVectorActiveMeasurements.ResizeTo(nactive);
   fVectorActiveObservable.clear();

   int counter = 0;
   for (int i = 0; i < nmeasurements; ++i) {
      BCMVCMeasurement* m = GetMeasurement(i);
      if (m->GetFlagActive()) {
         fVectorActiveMeasurements[counter] =  m->GetCentralValue();
         fVectorActiveObservable.push_back(m->GetObservable());
         counter++;
      }
   }

}

// ---------------------------------------------------------
void BCMVCombination::CalculateCorrelationMatrix(int index)
{
   BCMVCUncertainty* u = fUncertainties.at(index);

   int n = GetNMeasurements();
   int nactive = GetNActiveMeasurements();

   TMatrixD cov(nactive, nactive);

   TMatrixD corr = u->GetCorrelationMatrix();

   int counteri = 0;
   for (int i = 0; i < n; ++i) {
      BCMVCMeasurement* mi = GetMeasurement(i);
      double sigma_i = mi->GetUncertainty(index);

      // skip line if not active
      if (!mi->GetFlagActive())
         continue;

      int counterj = 0;
      for (int j = 0; j < n; ++j) {
         BCMVCMeasurement* mj = GetMeasurement(j);
         double sigma_j = mj->GetUncertainty(index);

         if (mj->GetFlagActive()) {
            cov[counteri][counterj] = corr[i][j]*sigma_i*sigma_j;
            counterj++;
         }
      }
      counteri++;
   }

   u->SetCovarianceMatrix(cov);
}

// ---------------------------------------------------------
void BCMVCombination::CalculateCovarianceMatrix(std::vector<double> parameters)
{
   int n = GetNUncertainties();
   int nmeasurements = GetNMeasurements();
   int nnuisance = fNNuisanceCorrelation;

   fCovarianceMatrix.Clear();
   fCovarianceMatrix.ResizeTo(GetNActiveMeasurements(), GetNActiveMeasurements());

   // loop over all uncertainties
   for (int unc_i = 0; unc_i < n; ++unc_i) {
      BCMVCUncertainty* u = fUncertainties.at(unc_i);

      // shrink covariance matrix such that it fits only active measurements
      TMatrixD mat_small = u->GetCovarianceMatrix();

      // loop over all measurements (i)
      int counteri = 0;
      for (int meas_i = 0; meas_i < nmeasurements; ++meas_i) {
         BCMVCMeasurement* mi = GetMeasurement(meas_i);

         // skip line if not active
         if (!mi->GetFlagActive())
            continue;

         // loop over all measurements (j)
         int counterj = 0;
         for (int meas_j = 0; meas_j < nmeasurements; ++meas_j) {
            BCMVCMeasurement* mj = GetMeasurement(meas_j);

            // skip line if not active
            if (!mj->GetFlagActive())
               continue;

            // modify matrix if nuisance parameter present
            if (parameters.size() > 0) {

               // loop over all nuisance parameters
               for (int  nuis_k = 0; nuis_k < nnuisance; ++nuis_k) {

                  // get nuisance parameter
                  NuisanceParameter p = fNuisanceCorrelation.at(nuis_k);

                  // compare
                  if (p.index_uncertainty == unc_i) {
                     double sigma_i = GetMeasurement(p.index_measurement1)->GetUncertainty(unc_i);
                     double sigma_j = GetMeasurement(p.index_measurement2)->GetUncertainty(unc_i);
                     double pre = p.pre;

                     // check if indices agree
                     if (meas_i == p.index_measurement1 && meas_j == p.index_measurement2) {
                        mat_small[counteri][counterj] = pre * parameters.at(p.index_rhoparameter) * sigma_i*sigma_j;
                        mat_small[counterj][counteri] = mat_small[counteri][counterj];
                     }
                  }
               } // end loop: nuisance parameters
            }

            // increase counter of active measurements
            counterj++;
         } // end loop: measurements j

         // increase counter of active measurements
         counteri++;
      } // end loop: measurements i

      // add matrix if active
      if (u->GetFlagActive())
         fCovarianceMatrix += mat_small;
   }
   fInvCovarianceMatrix.Clear();
   fInvCovarianceMatrix.ResizeTo(fCovarianceMatrix);
   fInvCovarianceMatrix = fCovarianceMatrix;
   fInvCovarianceMatrix.Invert();

   fDetCovariance = fCovarianceMatrix.Determinant();
}

// ---------------------------------------------------------
bool BCMVCombination::PositiveDefinite()
{
   TMatrixDEigen m(fCovarianceMatrix);

   // calculate real part of all eigenvalues
   TVectorD eigen_re = m.GetEigenValuesRe();

   int n_eigen = eigen_re.GetNoElements();

   bool flag_ispositive = true;

   // check if eigenvalues are positive
   for (int i = 0; i < n_eigen; ++i) {
      if (eigen_re[i] < 0)
         flag_ispositive = false;
   }

   // true if all eigenvalues are positive
   return flag_ispositive;
}

// ---------------------------------------------------------
void BCMVCombination::CalculateBLUE()
{
   int nmeas = GetNMeasurements();
   int nactivemeas = GetNActiveMeasurements();
   int nobs = GetNObservables();
   int nunc = GetNUncertainties();

   // calculate U matrix
   //  TMatrixD u(nmeas, nobs);
   TMatrixD u(nactivemeas, nobs);

   int counter = 0;
   for (int i = 0; i < nmeas; ++i) {
      BCMVCMeasurement* m = GetMeasurement(i);

      // if measurement is active fill matrix u
      if (m->GetFlagActive()) {
         for (int j = 0; j < nobs; ++j) {
            if (m->GetObservable() == j)
               u[counter][j] = 1;
            else
               u[counter][j] = 0;
         }
         counter++;
      }
   }

   // calculate weight matrix
   TMatrixD ut = u;
   ut.Transpose(ut);

   TMatrixD m1 = ut * fInvCovarianceMatrix;
   TMatrixD m2 = m1 * u;
   m2.Invert();

   fBLUEWeights.Clear();
   //  fBLUEWeights.ResizeTo(nobs, nmeas);
   fBLUEWeights.ResizeTo(nobs, nactivemeas);

   fBLUEWeights = m2*m1;

   // calculate central values
   fBLUECentral.Clear();
   fBLUECentral.ResizeTo(nobs);

   //  fBLUECentral = fBLUEWeights * fVectorMeasurements;
   fBLUECentral = fBLUEWeights * fVectorActiveMeasurements;

   // calculate covariance matrix
   fBLUECovarianceMatrix.Clear();
   fBLUECovarianceMatrix.ResizeTo(nobs, nobs);

   TMatrixD weightt = fBLUEWeights;
   weightt.Transpose(weightt);

   fBLUECovarianceMatrix = fBLUEWeights * fCovarianceMatrix * weightt;

   // calculate uncertainties, covariance and correlation matrices for each uncertainty
   for (int i = 0; i < nunc; ++i) {

      // calculate covariance matrix
      BCMVCUncertainty* u = GetUncertainty(i);
      if (!u->GetFlagActive())
         continue;

      TMatrixD cov = u->GetCovarianceMatrix();
      TMatrixD mat = fBLUEWeights * cov * weightt;
      fBLUECovarianceMatrices.push_back(mat);

      // calculate uncertainties
      TVectorD vec(nobs);
      for (int j = 0; j < nobs; ++j) {
         double sigma_j = sqrt(mat[j][j]);
         vec[j] = sigma_j;
      }
      fBLUEUncertaintiesPerSource.push_back(vec);

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
void BCMVCombination::PrintBLUEResults(std::string filename)
{
   // open file
   std::ofstream ofi(filename.c_str());

   // check if file is open
   if (!ofi.is_open()) {
      std::cerr << "Couldn't open file " << filename << std::endl;
      return;
   }

   int nobs = GetNObservables();
   int nmeas = GetNMeasurements();
   int nunc = GetNUncertainties();

   ofi << std::endl;
   ofi << "Summary of the combination" << std::endl;
   ofi << "==========================" << std::endl << std::endl;

   ofi << "* Observables:" << std::endl;
   ofi << "  Observable (range): " << std::endl;
   for (int i = 0; i < nobs; ++i)
      ofi << "  " << std::setiosflags(std::ios::left) << GetParameter(i)->GetName()
      << " (" << GetParameter(i)->GetLowerLimit() << " - " << GetParameter(i)->GetUpperLimit() << ")" << std::endl;
   ofi << std::endl;

   ofi << "* Measurements:" << std::endl;
   ofi << "  Measurement (observable): central value +- total uncertainty" << std::endl;
   for (int i = 0; i < nmeas; ++i) {
      BCMVCMeasurement* m = GetMeasurement(i);
      if (m->GetFlagActive()) {
         double total2 = 0;
         for (int j = 0; j < nunc; ++j) {
            if (GetUncertainty(j)->GetFlagActive())
               total2+=m->GetUncertainty(j)*m->GetUncertainty(j);
         }
         ofi << "  " << std::setiosflags(std::ios::left) << std::setw(20) << m->GetName()
	        << std::setiosflags(std::ios::left) << " (" << GetParameter(m->GetObservable())->GetName() << ")"
	        << ": " << std::setiosflags(std::ios::left) << std::setw(7) << std::setprecision(4) << m->GetCentralValue()
	        << " +- " << std::setiosflags(std::ios::left) << std::setw(7) << std::setprecision(4) << sqrt(total2) << std::endl;
      }
   }
   ofi << std::endl;

   ofi << "* Uncertainties:" << std::endl;
   ofi << "  Measurement (observable): Uncertainty (";
   for (int j = 0; j < nunc-1; ++j )
      if (GetUncertainty(j)->GetFlagActive())
         ofi << GetUncertainty(j)->GetName() << ", ";
   if (GetUncertainty(nunc-1)->GetFlagActive())
      ofi << GetUncertainty(nunc-1)->GetName() << ")" << std::endl;
   else
      ofi << ")" << std::endl;

   for (int i = 0; i < nmeas; ++i) {
      BCMVCMeasurement* m = GetMeasurement(i);
      if (m->GetFlagActive()) {
         ofi << "  " << std::setiosflags(std::ios::left) << std::setw(20) << m->GetName()
	        << std::setiosflags(std::ios::left) << " (" << GetParameter(m->GetObservable())->GetName() << "): ";
         for (int j = 0; j < nunc; ++j )
            if (GetUncertainty(j)->GetFlagActive())
               ofi << std::setiosflags(std::ios::left) << std::setw(7) << m->GetUncertainty(j);
         ofi << std::endl;
      }
   }
   ofi << std::endl;

   for (int i = 0; i < nunc; ++i) {
      BCMVCUncertainty* u = GetUncertainty(i);

      if (!u->GetFlagActive())
         continue;

      ofi << "  " << u->GetName() << " " << "(correlation matrix)" << std::endl;
      TMatrixD mat = u->GetCorrelationMatrix();

      int counterk = 0;
      for (int k = 0; k < nmeas; ++k) {
         BCMVCMeasurement* mk = GetMeasurement(k);

         // skip line if not active
         if (!mk->GetFlagActive())
            continue;

         ofi << "  ";

         int counterj = 0;
         for (int j = 0; j < nmeas; ++j) {
            BCMVCMeasurement* mj = GetMeasurement(j);

            if (mj->GetFlagActive()) {
               ofi << std::setw(7) << std::showpos << mat[k][j] << " ";
               counterj++;
            }
         }
         ofi << std::noshowpos << std::endl;
         counterk++;
      }
      ofi << std::endl;
   }

   ofi << "* BLUE results: " << std::endl;
   ofi << "  Observable: estimate +- total uncertainty" << std::endl;
   for (int i = 0; i < nobs; ++i) {
      if (i < fBLUECentral.GetNoElements())
         ofi << "  " << GetParameter(i)->GetName() << ": " << fBLUECentral[i] << " +- " << std::setprecision(4) << fBLUEUncertainties[i] << std::endl;
   }
   ofi << std::endl;

   ofi << "  Observable: Uncertainty (";
   for (int j = 0; j < nunc-1; ++j )
      if (GetUncertainty(j)->GetFlagActive())
         ofi << GetUncertainty(j)->GetName() << ", ";
   if (GetUncertainty(nunc-1)->GetFlagActive())
      ofi << GetUncertainty(nunc-1)->GetName() << ")" << std::endl;
   else
      ofi << ")" << std::endl;

   for (int i = 0; i < nobs; ++i) {
      ofi << "  " << std::setiosflags(std::ios::left) << GetParameter(i)->GetName()<< ":";
      int counterj = 0;
      for (int j = 0; j < nunc; ++j )
         if (GetUncertainty(j)->GetFlagActive()) {
            ofi << " " << std::setiosflags(std::ios::left) << std::setw(7) << std::setprecision(4) << GetBLUEUncertainties(counterj)[i];
            counterj++;
         }
      ofi << std::endl;
   }
   ofi << std::endl;

   if (nobs > 1) {
      ofi << "  Individual correlation matrices " << std::endl;
      for (int i = 0; i < nunc; ++i) {
         if (GetUncertainty(i)->GetFlagActive()) {
            TMatrixD mat = GetBLUECorrelationMatrix(i);
            ofi << "  " << GetUncertainty(i)->GetName() << std::endl;
            for (int j = 0; j < nobs; ++j) {
               ofi << "  ";
               for (int k = 0; k < nobs; ++k) {
                  ofi << std::setw(7) << std::setprecision(4) << std::showpos << mat[j][k] << " ";
               }
               ofi << std::noshowpos << std::endl;
            }
            ofi << std::endl;
         }
      }
   }

   if (nobs > 1) {
      ofi << "  Overall correlation matrix" << std::endl;
      TMatrixD mat = fBLUECorrelationMatrix;
      for (int j = 0; j < nobs; ++j) {
         ofi << "  ";
         for (int k = 0; k < nobs; ++k)
            ofi << std::setw(7) << std::setprecision(4) << std::showpos << mat[j][k] << " ";
         ofi << std::noshowpos << std::endl;
      }
      ofi << std::endl;
   }

   ofi << "  Weights [%]:" <<std::endl;
   int counter = 0;
   for (int j = 0; j < nmeas; ++j) {
      if (GetMeasurement(j)->GetFlagActive()) {
         ofi << "  " << std::setw(20) << GetMeasurement(j)->GetName() << " : ";
         for (int k = 0; k < nobs; ++k)
            ofi << std::setw(7) << std::setprecision(4) << std::showpos << fBLUEWeights[k][counter]*100. << " ";
         ofi << std::endl;
         counter++;
      }
   }
   ofi << std::endl;

   // close file
   ofi.close();
}

// ---------------------------------------------------------
void BCMVCombination::PrintMatrix(TMatrixD &matrix, std::string name)
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
void BCMVCombination::PrintVector(TVectorD &vector, std::string name)
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
int BCMVCombination::GetIndexMeasurement(std::string measurement, std::string observable)
{
   int index_observable = GetIndexObservable(observable);

   int nmeasurements = GetNMeasurements();

   for (int i = 0; i < nmeasurements; ++i) {
      if ((measurement == GetMeasurement(i)->GetName()) && (index_observable == GetMeasurement(i)->GetObservable()))
         return i;
   }

   return -1;
}

// ---------------------------------------------------------
int BCMVCombination::GetIndexUncertainty(std::string name)
{
   int nuncertainties = GetNUncertainties();

   int index = -1;

   for (int i = 0; i < nuncertainties; ++i) {
      if (name == GetUncertainty(i)->GetName())
         index = i;
   }

   return index;
}

// ---------------------------------------------------------
int BCMVCombination::GetIndexObservable(std::string name)
{
   int n = GetNObservables();

   // go through the list of parameters and compare strings
   for (int i = 0; i < n; ++i) {
      //    if (name == std::string(GetParameter(i)->GetName()))
      if (name == std::string(fObservables.at(i)->GetName()))
         return i;
   }

   // return -1 if not in the list
   return -1;
}

// ---------------------------------------------------------
