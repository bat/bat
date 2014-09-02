/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModel.h"

#include "BCDataPoint.h"
#include "BCDataSet.h"
#include "BCGoFTest.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"
#include "BCMath.h"
#include "BCModelOutput.h"
#include "BCParameter.h"
#include "BCPriorModel.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TArrow.h>

#include <fstream>
#include <iomanip>
#include <set>

// ---------------------------------------------------------
BCModel::BCModel(const char * name)
	: BCIntegrate(name)
  , fModelAPriori(0)
  , fModelAPosteriori(0)
  , fDataSet(0)
  , fDataPointLowerBoundaries(0)
  , fDataPointUpperBoundaries(0)
  , fPValue(-1)
  , fChi2NDoF(-1)
  , fPValueNDoF(-1)
  , flag_discrete(false)
  , fGoFNIterationsMax(100000)
  , fGoFNIterationsRun(2000)
  , fGoFNChains(5)
  , fPriorConstantAll(false)
	, fPriorModel(0)
{
}

// ---------------------------------------------------------
BCModel::BCModel(const BCModel & bcmodel)
	: BCIntegrate(bcmodel)
	, fPriorModel(0)
{
	Copy(bcmodel);
}

// ---------------------------------------------------------
void BCModel::Copy(const BCModel & bcmodel)
{
   //  called for the second time in copy constructor? do copy-and-swap instead
   //   BCIntegrate::Copy(bcmodel);
   fName                            = bcmodel.fName;
   fModelAPriori                    = bcmodel.fModelAPriori;
   fModelAPosteriori                = bcmodel.fModelAPosteriori;
   if (fDataSet)
      fDataSet = bcmodel.fDataSet;
   else
      fDataSet = 0;

   if (bcmodel.fDataPointLowerBoundaries)
      fDataPointLowerBoundaries = new BCDataPoint(*bcmodel.fDataPointLowerBoundaries);
   else
      fDataPointLowerBoundaries = 0;
   if (bcmodel.fDataPointUpperBoundaries)
      fDataPointUpperBoundaries = new BCDataPoint(*bcmodel.fDataPointUpperBoundaries);
   else
      fDataPointUpperBoundaries = 0;

   fDataFixedValues                 = bcmodel.fDataFixedValues;

   fPValue                          = bcmodel.fPValue;
   fChi2NDoF                        = bcmodel.fChi2NDoF;
   fPValueNDoF                      = bcmodel.fPValueNDoF;
   flag_discrete                    = bcmodel.flag_discrete;
   fGoFNIterationsMax               = bcmodel.fGoFNIterationsMax;
   fGoFNIterationsRun               = bcmodel.fGoFNIterationsRun;
   fGoFNChains                      = bcmodel.fGoFNChains;
   for (int i = 0; i < int(bcmodel.fPriorContainer.size()); ++i) {
      if (bcmodel.fPriorContainer.at(i))
         fPriorContainer.push_back(new TNamed(*bcmodel.fPriorContainer.at(i)));
      else
         fPriorContainer.push_back(0);
   }
   fPriorConstantAll                = bcmodel.fPriorConstantAll;
   fPriorContainerConstant          = bcmodel.fPriorContainerConstant;
   fPriorContainerInterpolate       = bcmodel.fPriorContainerInterpolate;
}

// ---------------------------------------------------------
BCModel::~BCModel()
{
   for (unsigned int i = 0; i < GetNParameters(); ++i)
      delete fPriorContainer[i];
   fPriorContainer.clear();

	 if (fPriorModel)
		 delete fPriorModel;

   delete fDataPointLowerBoundaries;
   delete fDataPointUpperBoundaries;
}

// ---------------------------------------------------------
BCModel & BCModel::operator = (const BCModel & bcmodel)
{
   Copy(bcmodel);

   return *this;
}

// ---------------------------------------------------------
unsigned BCModel::GetNDataPoints() const
{
   if (fDataSet)
      return fDataSet->GetNDataPoints();
   else
      return 0;
   }

// ---------------------------------------------------------
BCDataPoint * BCModel::GetDataPoint(unsigned int index) const
{
   if (fDataSet)
      return fDataSet->GetDataPoint(index);

   BCLog::OutWarning("BCModel::GetDataPoint : No data set defined.");
   return 0;
}

// ---------------------------------------------------------
double BCModel::GetDataPointLowerBoundary(unsigned int index) const
{
    return fDataPointLowerBoundaries -> GetValue(index);
}

// ---------------------------------------------------------
double BCModel::GetDataPointUpperBoundary(unsigned int index) const
{
    return fDataPointUpperBoundaries -> GetValue(index);
}

// ---------------------------------------------------------
bool BCModel::GetFlagBoundaries() const
{
   if (!fDataPointLowerBoundaries)
      return false;

   if (!fDataPointUpperBoundaries)
      return false;

   if (fDataPointLowerBoundaries->GetNValues() != fDataSet->GetDataPoint(0)->GetNValues())
      return false;

   if (fDataPointUpperBoundaries->GetNValues() != fDataSet->GetDataPoint(0)->GetNValues())
      return false;

   return true;
}

// ---------------------------------------------------------
void BCModel::SetSingleDataPoint(BCDataPoint * datapoint)
{
   // create new data set consisting of a single data point
   BCDataSet * dataset = new BCDataSet();

   // add the data point
   dataset->AddDataPoint(datapoint);

   // set this new data set
   SetDataSet(dataset);
}

// ---------------------------------------------------------
void BCModel::SetSingleDataPoint(BCDataSet * dataset, unsigned int index)
{
   if (index > dataset->GetNDataPoints())
      return;

   SetSingleDataPoint(dataset->GetDataPoint(index));
}

// ---------------------------------------------------------
void BCModel::SetDataBoundaries(unsigned int index, double lowerboundary, double upperboundary, bool fixed)
{
   // check if data set exists
   if (!fDataSet) {
      BCLog::OutError("BCModel::SetDataBoundaries : Need to define data set first.");
      return;
   }

   // check if index is within range
   if (index > fDataSet->GetDataPoint(0)->GetNValues()) {
      BCLog::OutError("BCModel::SetDataBoundaries : Index out of range.");
      return;
   }

   // check if boundary data points exist
   if (!fDataPointLowerBoundaries)
      fDataPointLowerBoundaries = new BCDataPoint(fDataSet->GetDataPoint(0)->GetNValues());

   if (!fDataPointUpperBoundaries)
      fDataPointUpperBoundaries = new BCDataPoint(fDataSet->GetDataPoint(0)->GetNValues());

   if (fDataFixedValues.size() == 0)
      fDataFixedValues.assign(fDataSet->GetDataPoint(0)->GetNValues(), false);

   // set boundaries
   fDataPointLowerBoundaries->SetValue(index, lowerboundary);
   fDataPointUpperBoundaries->SetValue(index, upperboundary);
   fDataFixedValues[index] = fixed;
}

// ---------------------------------------------------------
int BCModel::AddParameter(BCParameter * parameter)
{
	if ( !BCEngineMCMC::AddParameter(parameter))
		return 1;

   // add empty object to prior container
   fPriorContainer.push_back(0);

   // don't interpolate the prior histogram by default
   fPriorContainerInterpolate.push_back(false);

   // prior assumed to be non-constant in general case
   fPriorContainerConstant.push_back(false);

   return 0;
}

// ---------------------------------------------------------
double BCModel::ProbabilityNN(const std::vector<double> &params)
{
   return exp(LogProbabilityNN(params) );
}

// ---------------------------------------------------------
double BCModel::Probability(const std::vector<double> &parameter)
{
   return exp(LogProbability(parameter));
}

// ---------------------------------------------------------
double BCModel::LogProbability(const std::vector<double> &parameters)
{
   // check if normalized
  if (GetIntegral() <= 0.) {
      BCLog::OutError("BCModel::LogProbability. Normalization not available or zero.");
      return 0.;
   }

  return LogProbabilityNN(parameters) - log(GetIntegral());
}

// ---------------------------------------------------------
double BCModel::APrioriProbability(const std::vector<double> &parameters)
{
   return exp(this->LogAPrioriProbability(parameters) );
}

// ---------------------------------------------------------
double BCModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
   double logprob = 0;

   // loop over all parameters, assume prior factorizes
   // into n independent parts
   for (unsigned i = 0; i < GetNParameters(); ++i) {
		 BCParameter * par = GetParameter(i);

      // avoid fixed and zero-width parameters
      if (par->Fixed() or not par->GetRangeWidth())
         continue;

      if (fPriorContainerConstant[i]) {
         logprob -= log(par->GetRangeWidth());
         continue;
      }

      if (fPriorContainer[i]) {
         // check what type of object is stored
         TF1 * f = dynamic_cast<TF1*>(fPriorContainer[i]);
         TH1 * h = dynamic_cast<TH1*>(fPriorContainer[i]);

         if (f) // TF1
            logprob += log(f->Eval(parameters[i]));
         else if (h) { // TH1
            if(fPriorContainerInterpolate[i])
               logprob += log(h->Interpolate(parameters[i]));
            else
               logprob += log(h->GetBinContent(h->FindBin(parameters[i])));
         }
         else
            BCLog::OutError(Form(
                  "BCModel::LogAPrioriProbability : Prior for parameter %s "
                  "is defined but not recognized.",
                  par->GetName().c_str())); // this should never happen
      }
      // use constant only if user has defined it
      else {
         BCLog::OutError(Form(
               "BCModel::LogAPrioriProbability: Prior for parameter %s "
               "is undefined. Using constant prior to proceed.",
               par->GetName().c_str()));
         logprob -= log(par->GetRangeWidth());
      }
   }

   return logprob;
}

// ---------------------------------------------------------
double BCModel::Likelihood(const std::vector<double> &params)
{
   return exp(LogLikelihood(params));
}

// ---------------------------------------------------------
double BCModel::Eval(const std::vector<double> &parameters)
{
   return exp(LogEval(parameters));
}

// ---------------------------------------------------------
double BCModel::LogEval(const std::vector<double> &parameters)
{
   return LogProbabilityNN(parameters);
}

// ---------------------------------------------------------
double BCModel::SamplingFunction(const std::vector<double> & /*parameters*/)
{
   double probability = 1;
   for (unsigned i = 0 ; i < GetNParameters() ; ++i)
		 probability *= 1. / GetParameter(i)->GetRangeWidth();
   return probability;
}

// ---------------------------------------------------------
double BCModel::GetPvalueFromChi2(const std::vector<double> &par, int sigma_index)
{
   double ll = LogLikelihood(par);
   int n = GetNDataPoints();

   double sum_sigma = 0;
   for (int i = 0; i < n; i++)
      sum_sigma += log(GetDataPoint(i)->GetValue(sigma_index));

   double chi2 = -2. * (ll + (double) n / 2. * log(2. * M_PI) + sum_sigma);

   fPValue = TMath::Prob(chi2, n);

   return fPValue;
}

// ---------------------------------------------------------
std::vector<double> BCModel::GetChi2Runs(int /*dataIndex*/, int /*sigmaIndex*/)
{
   std::vector<double> x;
   return x;
}

// ---------------------------------------------------------
double BCModel::GetPvalueFromChi2NDoF(std::vector<double> par, int sigma_index)
{
   double ll = LogLikelihood(par);
   int n = GetNDataPoints();
   int npar = GetNParameters();

   double sum_sigma = 0;
   for (int i = 0; i < n; i++)
      sum_sigma += log(GetDataPoint(i)->GetValue(sigma_index));

   double chi2 = -2. * (ll + (double) n / 2. * log(2. * M_PI) + sum_sigma);

   fChi2NDoF = chi2 / double(n - npar);
   fPValueNDoF = TMath::Prob(chi2, n - npar);

   return fPValueNDoF;
}

// ---------------------------------------------------------
double BCModel::GetPvalueFromKolmogorov(const std::vector<double>& par,int index)
{
   if (flag_discrete) {
      BCLog::OutError(Form("BCModel::GetPvalueFromKolmogorov : "
            "test defined only for continuous distributions."));
      return -1.;
   }

   //calculate the ECDF from the 1D data
   std::vector<double> yData = fDataSet->GetDataComponents(index);
   TH1D * ECDF = BCMath::ECDF(yData);

   int N = GetNDataPoints();

   // calculated expected CDF for unique points only
   std::set<double> uniqueObservations;
   for (int i = 0; i < N; i++)
      uniqueObservations.insert(CDF(par, i, false));

   int nUnique = uniqueObservations.size();
   if (nUnique != ECDF->GetNbinsX() + 1) {
      BCLog::OutError(Form("BCModel::GetPvalueFromKolmogorov : "
            "Number of unique data points doesn't match (%d vs %d)", nUnique,
            ECDF->GetNbinsX() + 1));
      return -1.;
   }

   // find maximum distance
   double distMax = 0.;

   // current distance
   double dist = 0.;

   std::set<double>::const_iterator iter = uniqueObservations.begin();
   for (int iBin = 0; iBin < nUnique; ++iBin) {
      // distance between current points
      dist = TMath::Abs(*iter - ECDF->GetBinContent(iBin + 1));
      // update maximum if necessary
      distMax = TMath::Max(dist, distMax);

      // advance to next entry in the set
      ++iter;
   }

   // correct for total #points, not unique #points.
   // would need sqrt(n1*n2/(n1+n2)) if we had two experimental datasets
   double z = distMax * sqrt(N);

   fPValue = TMath::KolmogorovProb(z);

   // clean up
   delete ECDF;

   return fPValue;
}

// ---------------------------------------------------------
BCH1D * BCModel::CalculatePValue(std::vector<double> par, bool flag_histogram)
{
   BCH1D * hist = 0;

   // print log
   BCLog::OutSummary("Do goodness-of-fit-test");

   // create model test
   BCGoFTest * goftest = new BCGoFTest("modeltest");

   // set this model as the model to be tested
   goftest->SetTestModel(this);

   // set the point in parameter space which is tested an initialize
   // the model testing
   if (!goftest->SetTestPoint(par))
      return 0;

   // disable the creation of histograms to save _a lot_ of memory
   // (histograms are not needed for p-value calculation)
   goftest->MCMCSetFlagFillHistograms(false);

   // set parameters of the MCMC for the GoFTest
   goftest->MCMCSetNChains(fGoFNChains);
   goftest->MCMCSetNIterationsPreRunMax(fGoFNIterationsMax);
   goftest->MCMCSetNIterationsRun(fGoFNIterationsRun);

   // get p-value
   fPValue = goftest->GetCalculatedPValue(flag_histogram);

   // get histogram
   if (flag_histogram) {
      hist = new BCH1D();
      hist->SetHistogram(goftest->GetHistogramLogProb());
   }

   // delete model test
   delete goftest;

   // return histogram
   return hist;
}

// ---------------------------------------------------------
void BCModel::CorrelateDataPointValues(std::vector<double> & /*x*/)
{
   // ...
}

// ---------------------------------------------------------
double BCModel::HessianMatrixElement(const BCParameter * par1, const BCParameter * par2, std::vector<double> point)
{
   // check number of entries in vector
   if (point.size() != GetNParameters()) {
      BCLog::OutError("BCModel::HessianMatrixElement : Invalid number of entries in the vector.");
      return -1;
   }

   // define steps
   double nsteps = 1e5;
   double dx1 = par1->GetRangeWidth() / nsteps;
   double dx2 = par2->GetRangeWidth() / nsteps;

   // define points at which to evaluate
   std::vector<double> xpp = point;
   std::vector<double> xpm = point;
   std::vector<double> xmp = point;
   std::vector<double> xmm = point;

   unsigned idx1 = fParameters.Index(par1->GetName());
   unsigned idx2 = fParameters.Index(par2->GetName());

   xpp[idx1] += dx1;
   xpp[idx2] += dx2;

   xpm[idx1] += dx1;
   xpm[idx2] -= dx2;

   xmp[idx1] -= dx1;
   xmp[idx2] += dx2;

   xmm[idx1] -= dx1;
   xmm[idx2] -= dx2;

   // calculate probability at these points
   double ppp = Likelihood(xpp);
   double ppm = Likelihood(xpm);
   double pmp = Likelihood(xmp);
   double pmm = Likelihood(xmm);

   // return derivative
   return (ppp + pmm - ppm - pmp) / (4. * dx1 * dx2);
}

// ---------------------------------------------------------
void BCModel::SetDataPointLowerBoundary(int index, double lowerboundary)
{
   fDataPointLowerBoundaries -> SetValue(index, lowerboundary);
}

// ---------------------------------------------------------
void BCModel::SetDataPointUpperBoundary(int index, double upperboundary)
{
   fDataPointUpperBoundaries -> SetValue(index, upperboundary);
}

// ---------------------------------------------------------
int BCModel::SetPrior(int index, TF1 * f)
{
   // check index range
   if (index < 0 && index >= int(GetNParameters())) {
      BCLog::OutError("BCModel::SetPrior : Index out of range.");
      return 0;
   }

   if (fPriorContainer[index])
      delete fPriorContainer[index];

   // copy function
   fPriorContainer[index] = new TF1(*f);

   fPriorContainerConstant[index] = false;

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCModel::SetPrior(const char * name, TF1 * f)
{
   // find index
   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (name == GetParameter(i)->GetName())
         index = i;
	 if (index < 0) {
		 BCLog::OutError(TString::Format("BCModel::SetPrior : No parameter named %s.",name));
		 return 0;
	 }
		 
   // set prior
   return SetPrior(index, f);
}

// ---------------------------------------------------------
int BCModel::SetPriorDelta(int index, double value)
{
   // set range to value
   GetParameter(index)->Fix(value);

   // set prior
   return 1;
}

// ---------------------------------------------------------
int BCModel::SetPriorDelta(const char* name, double value)
{
   GetParameter(name)->Fix(value);

   return 1;
}

// ---------------------------------------------------------
int BCModel::SetPriorGauss(int index, double mean, double sigma)
{
   // check index range
   if (index < 0 || index >= int(GetNParameters())) {
      BCLog::OutError("BCModel::SetPriorGauss : Index out of range.");
      return 0;
   }

   // create new function
   TF1 * f = new TF1(Form("prior_%s", GetParameter(index)->GetName().c_str()),
         "1./sqrt(2.*TMath::Pi())/[1] * exp(- (x-[0])*(x-[0])/2./[1]/[1])",
         GetParameter(index)->GetLowerLimit(),
         GetParameter(index)->GetUpperLimit());
   f->SetParameter(0, mean);
   f->SetParameter(1, sigma);

   // set prior
   return SetPrior(index, f);
}

// ---------------------------------------------------------
int BCModel::SetPriorGauss(const char* name, double mean, double sigma)
{
   // find index
   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (name == GetParameter(i)->GetName())
         index = i;
	 if (index < 0) {
		 BCLog::OutError(TString::Format("BCModel::SetPriorGauss : No parameter named %s.",name));
		 return 0;
	 }

   // set prior
   return SetPriorGauss(index, mean, sigma);
}

// ---------------------------------------------------------
int BCModel::SetPriorGauss(int index, double mean, double sigmadown, double sigmaup)
{
   // check index range
   if (index < 0 || index >= int(GetNParameters())) {
      BCLog::OutError("BCModel::SetPriorGauss : Index out of range.");
      return 0;
   }

   // create new function
   TF1 * f = new TF1(Form("prior_%s", GetParameter(index)->GetName().c_str()),
         BCMath::SplitGaussian,
         GetParameter(index)->GetLowerLimit(),
         GetParameter(index)->GetUpperLimit(),
         3);
   f->SetParameter(0, mean);
   f->SetParameter(1, sigmadown);
   f->SetParameter(2, sigmaup);

   // set prior
   return SetPrior(index, f);
}

// ---------------------------------------------------------
int BCModel::SetPriorGauss(const char * name, double mean, double sigmadown, double sigmaup)
{
   // find index
   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (name == GetParameter(i)->GetName())
         index = i;
	 if (index < 0) {
		 BCLog::OutError(TString::Format("BCModel::SetPriorGauss : No parameter named %s.",name));
		 return 0;
	 }

   // set prior
   return SetPriorGauss(index, mean, sigmadown, sigmaup);
}

// ---------------------------------------------------------
int BCModel::SetPrior(int index, TH1 * h, bool interpolate)
{
   // check index range
   if (index < 0 && index >= int(GetNParameters())) {
      BCLog::OutError("BCModel::SetPrior : Index out of range.");
      return 0;
   }

   // if the histogram exists
   if(h) {

      // check if histogram is 1d
      if (h->GetDimension() != 1) {
         BCLog::OutError(Form("BCModel::SetPrior : Histogram given for parameter %d is not 1D.",index));
         return 0;
      }

      // normalize the histogram
      h->Scale(1./h->Integral("width"));

      if(fPriorContainer[index])
         delete fPriorContainer[index];

      // set function
      fPriorContainer[index] = (TNamed*) h->Clone();

      if (interpolate)
         fPriorContainerInterpolate[index] = true;

      fPriorContainerConstant[index] = false;
   }

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCModel::SetPrior(const char * name, TH1 * h, bool interpolate)
{
   // find index
   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (name == GetParameter(i)->GetName())
         index = i;
	 if (index < 0) {
		 BCLog::OutError(TString::Format("BCModel::SetPrior : No parameter named %s.",name));
		 return 0;
	 }

   // set prior
   return SetPrior(index, h, interpolate);
}

// ---------------------------------------------------------
int BCModel::SetPriorConstant(int index)
{
   // check index range
   if (index < 0 && index >= int(GetNParameters())) {
      BCLog::OutError("BCModel::SetPriorConstant : Index out of range.");
      return 0;
   }

   if(fPriorContainer[index]) {
      delete fPriorContainer[index];
      fPriorContainer[index] = 0;
   }

   // set prior to a constant
   fPriorContainerConstant[index] = true;

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCModel::SetPriorConstantAll()
{
	if ( fParameters.Empty() )
      BCLog::OutWarning("BCModel::SetPriorConstantAll : No parameters defined.");

   // loop over all 1-d priors
   for (unsigned i = 0; i < fParameters.Size(); ++i) {
      if (fPriorContainer[i]) {
         delete fPriorContainer[i];
         fPriorContainer[i]=0;
      }
      fPriorContainerConstant[i] = true;
   }

   // no error
   return 1;
}

// ---------------------------------------------------------
void BCModel::PrintSummary()
{
	BCIntegrate::PrintSummary();

   // model testing
   if (fPValue >= 0) {
      double likelihood = Likelihood(GetBestFitParameters());
      BCLog::OutSummary(" Model testing:");
      BCLog::OutSummary(Form("  p(data|lambda*) = %f", likelihood));
      BCLog::OutSummary(Form("  p-value         = %f", fPValue));
   }
}

// ---------------------------------------------------------
void BCModel::PrintMarginalizationToStream(std::ofstream & ofi) {
	if (fPValue >= 0) {
		ofi << " Results of the model test" << std::endl
				<< " =========================" << std::endl
				<< " p-value: " << fPValue << std::endl;
		if (fPValueNDoF >= 0)
			ofi << " p-value corrected for degrees of freedom: " << fPValueNDoF << std::endl;
		ofi << std::endl;
	}
	BCIntegrate::PrintMarginalizationToStream(ofi);
}

// ---------------------------------------------------------
void BCModel::PrintShortFitSummary(int chi2flag)
{
   BCLog::OutSummary("---------------------------------------------------");
   BCLog::OutSummary(Form("Fit summary for model \'%s\':", GetName().data()));
   BCLog::OutSummary(Form("   Number of parameters:  Npar  = %i", GetNParameters()));
   if (GetNDataPoints()) {
      BCLog::OutSummary(Form("   Number of data points: Ndata = %i", GetNDataPoints()));
      BCLog::OutSummary("   Number of degrees of freedom:");
      BCLog::OutSummary(Form("      NDoF = Ndata - Npar = %i", GetNDataPoints() - GetNParameters()));
   }
   if (!GetBestFitParameters().empty())
      BCLog::OutSummary("   Best fit parameters (global):");
	 PrintParameters(GetBestFitParameters(),BCLog::OutSummary);

   if (GetPValue() >= 0) {
      BCLog::OutSummary("   Goodness-of-fit test:");
      BCLog::OutSummary(Form("      p-value = %.3g", GetPValue()));
      if (chi2flag) {
         BCLog::OutSummary(Form("      p-value corrected for NDoF = %.3g", GetPValueNDoF()));
         BCLog::OutSummary(Form("      chi2 / NDoF = %.3g", GetChi2NDoF()));
      }
   }
   BCLog::OutSummary("---------------------------------------------------");
}

// ---------------------------------------------------------
void BCModel::PrintHessianMatrix(std::vector<double> parameters)
{
   // check number of entries in vector
   if (parameters.size() != GetNParameters()) {
      BCLog::OutError("BCModel::PrintHessianMatrix : Invalid number of entries in the vector");
      return;
   }

   // print to screen
   BCLog::OutSummary("Hessian matrix elements: ");
   BCLog::OutSummary("Parameter values:");

   for (int i = 0; i < int(parameters.size()); i++)
      BCLog::OutSummary(Form("Parameter %d : %f", i, parameters.at(i)));

   BCLog::OutSummary("Hessian matrix:");
   // loop over all parameter pairs
   for (unsigned int i = 0; i < GetNParameters(); i++)
      for (unsigned int j = 0; j < i; j++) {
				// calculate Hessian matrix element
				double hessianmatrixelement = HessianMatrixElement(GetParameter(i),GetParameter(j),parameters);

				// print to screen
				BCLog::OutSummary(Form("%d %d : %f", i, j, hessianmatrixelement));
      }
}

// ---------------------------------------------------------
BCPriorModel * BCModel::GetPriorModel(bool prepare) {
	if (!fPriorModel)
		fPriorModel = new BCPriorModel(this);
	else if (prepare)
		fPriorModel -> PreparePriorModel();
	return fPriorModel;
}

// ---------------------------------------------------------
int BCModel::DrawKnowledgeUpdatePlot1D(unsigned index, std::string options_post, std::string options_prior) {
	// option flags
	bool flag_slice_post  = (options_post.find("slice") < options_post.size());
	bool flag_slice_prior = (options_prior.find("slice") < options_prior.size());

	// Get Prior
	BCH1D * bch1d_prior = 0;
	TLine * const_prior = (IsPriorConstant(index)) ? new TLine() : 0;
	TF1   * f1_prior    = (const_prior) ? 0 : dynamic_cast<TF1*> (PriorContainer(index));
	TH1   * h1_prior    = (const_prior) ? 0 : dynamic_cast<TH1*> (PriorContainer(index));
	
	double max_prior = 0;

	if (const_prior) {
		max_prior = 1./GetVariable(index)->GetRangeWidth();
		const_prior -> SetLineColor(kRed);
	}
	else if (f1_prior) {
		max_prior = f1_prior -> GetMaximum(GetVariable(index)->GetLowerLimit(),GetVariable(index)->GetUpperLimit());
		f1_prior -> SetLineColor(kRed);
		f1_prior -> SetLineWidth(1);
	}
	else if (h1_prior) {
		max_prior = h1_prior -> GetMaximum();
		h1_prior -> SetLineColor(kRed);
		h1_prior -> SetStats(false);
		h1_prior -> GetXaxis() -> SetNdivisions(508);
	}
	else {
		if (flag_slice_prior and index<fPriorModel->GetNParameters()) {
			if (fPriorModel->GetNParameters()==2) {
				TH1D * hist = (index==0) ? fPriorModel->GetSlice(0,1)->ProjectionX(Form("projx_%i",BCLog::GetHIndex())) : fPriorModel->GetSlice(0,1)->ProjectionY(Form("projy_%i",BCLog::GetHIndex()));
				bch1d_prior = new BCH1D(hist);
			} else if (fPriorModel->GetNParameters()==1)
				bch1d_prior = new BCH1D(fPriorModel->GetSlice(index));
		}
		if (!bch1d_prior and fPriorModel->MarginalizedHistogramExists(index))
			bch1d_prior = fPriorModel->GetMarginalized(index);
		if (!bch1d_prior)
			return 0;
		max_prior = bch1d_prior->GetHistogram()->GetMaximum();
		bch1d_prior -> GetHistogram() -> SetStats(false);
		bch1d_prior -> GetHistogram() -> SetLineColor(kRed);
		bch1d_prior -> GetHistogram() -> GetXaxis() -> SetNdivisions(508);
	}
	
	// if prior doesn't exist, exit
	if (!const_prior and !f1_prior and !h1_prior and !bch1d_prior)
		return 0;

	// Get Posterior
	BCH1D* bch1d_posterior = 0;
	if (flag_slice_post and index<GetNParameters()) {
		if (GetNParameters()==2) {
			TH1D * hist = (index==0) ? GetSlice(0,1)->ProjectionX(Form("projx_%i",BCLog::GetHIndex()))
				: GetSlice(0,1)->ProjectionY(Form("projy_%i",BCLog::GetHIndex()));
			hist -> Scale(1./hist->Integral("width"));
			bch1d_posterior = new BCH1D(hist);
		} else if (GetNParameters()==1)
			bch1d_posterior = new BCH1D(GetSlice(index));
	} else if (MarginalizedHistogramExists(index))
		bch1d_posterior = GetMarginalized(index);
	
	// if marginal doesn't exist, exit
	if (!bch1d_posterior)
		return 0;
	
	bch1d_posterior -> GetHistogram() -> Scale(1./bch1d_posterior->GetHistogram()->Integral("width"));
	bch1d_posterior -> GetHistogram() -> SetStats(kFALSE);
	
	// get maximum
	double maxy = 1.1 * TMath::Max(max_prior, bch1d_posterior->GetHistogram()->GetMaximum());

	// prepare legend
	TLegend* legend = new TLegend();
	legend->SetBorderSize(0);
	legend->SetFillColor(kWhite);
	legend->SetTextAlign(12);
	legend->SetTextFont(62);
	legend->SetTextSize(0.03);

	// draw axes
	TH2D * h2_axes = new TH2D(TString::Format("h2_axes_%s_knowledge_update_%d",GetName().data(),index), TString::Format(";%s;P(%s|Data)",GetVariable(index)->GetLatexName().data(),GetVariable(index)->GetLatexName().data()),
														10, GetVariable(index)->GetLowerLimit(), GetVariable(index)->GetUpperLimit(),
														10, 0, maxy);
	h2_axes -> SetStats(false);
	h2_axes -> GetXaxis() -> SetNdivisions(508);
	h2_axes -> Draw();
	
	// draw prior
	if (const_prior) {
		legend -> AddEntry(const_prior,"prior","L");
		const_prior -> DrawLine(GetVariable(index)->GetLowerLimit(),max_prior,GetVariable(index)->GetUpperLimit(),max_prior);
	} else if (f1_prior) {
		legend -> AddEntry(f1_prior,"prior","L");
		f1_prior -> Draw("same");
	} else if (h1_prior) {
		legend -> AddEntry(h1_prior,"prior","L");
		h1_prior -> Draw("same");
	} else if (bch1d_prior) {
		legend -> AddEntry(bch1d_prior->GetHistogram(), "prior", "L");
		bch1d_prior -> Draw(std::string(options_prior+"same"));
	}
	
	// draw posterior
	legend -> AddEntry(bch1d_posterior->GetHistogram(), "posterior", "L");
	bch1d_posterior -> Draw(std::string(options_post+"same"));

	gPad->SetTopMargin(0.02);
   
	// Draw legend on top of histogram
	legend -> SetX1NDC(gPad->GetLeftMargin() + 0.10 * (1.0 - gPad->GetRightMargin() - gPad->GetLeftMargin()));
	legend -> SetX2NDC(1. - gPad->GetRightMargin());
	double y1 = gPad->GetTopMargin() + legend->GetTextSize()*legend->GetNRows();
	legend -> SetY1NDC(1-y1);
	legend -> SetY2NDC(1. - gPad->GetTopMargin());
	legend -> Draw();

	// rescale top margin
	gPad -> SetTopMargin(y1+0.01);

	gPad -> RedrawAxis();

	return 1;
}

// ---------------------------------------------------------
int BCModel::DrawKnowledgeUpdatePlot2D(unsigned index1, unsigned index2, bool flag_slice, double interval_content) {
	
	if (index1 == index2)
		return 0;
	if (index1 > index2)
		return DrawKnowledgeUpdatePlot2D(index2,index1,flag_slice);

	// Get Posterior
	TH2D * h2d_2dposterior = 0;
	if (flag_slice and GetNParameters()==2 and index1<GetNParameters() and index2<GetNParameters())
		h2d_2dposterior = GetSlice(index1,index2);
	else if (MarginalizedHistogramExists(index1,index2))
		h2d_2dposterior = GetMarginalizedHistogram(index1,index2);

	if (!h2d_2dposterior)
		return 0;

	// Get Prior
	bool const_prior1 = IsPriorConstant(index1);
	TF1 * f1_prior1   = (const_prior1) ? 0 : dynamic_cast<TF1*> (PriorContainer(index1));
	TH1 * h1_prior1   = (const_prior1) ? 0 : dynamic_cast<TH1*> (PriorContainer(index1));
	bool auto_prior1 = const_prior1 or f1_prior1 or h1_prior1;

	bool const_prior2 = IsPriorConstant(index2);
	TF1 * f1_prior2   = (const_prior2) ? 0 : dynamic_cast<TF1*> (PriorContainer(index2));
	TH1 * h1_prior2   = (const_prior2) ? 0 : dynamic_cast<TH1*> (PriorContainer(index2));
	bool auto_prior2 = const_prior2 or f1_prior2 or h1_prior2;

	TH2D * h2d_2dprior = 0;

	if (!auto_prior1 or !auto_prior2) { // one or both prior pre-defined
		if (flag_slice and GetNParameters()==2 and index1<GetNParameters() and index2<GetNParameters())
			h2d_2dprior = fPriorModel -> GetSlice(index1,index2);
		else if (fPriorModel->MarginalizedHistogramExists(index1, index2))
			h2d_2dprior = fPriorModel -> GetMarginalizedHistogram(index1,index2);
	}

	// if not predefined, use the projection of the marginalization
	if (!auto_prior1 and h2d_2dprior)
		h1_prior1 = h2d_2dprior -> ProjectionX(TString::Format("h1_prior1_%s_%d",GetName().data(),index1));
	if (!auto_prior2 and h2d_2dprior)
		h1_prior2 = h2d_2dprior -> ProjectionY(TString::Format("h1_prior2_%s_%d",GetName().data(),index2));

	if (!h2d_2dprior)
		h2d_2dprior = fPriorModel->GetVariable(index1) -> CreateH2(TString::Format("h2d_2dprior_%s_%d_%d",GetName().data(),index1,index2).Data(),fPriorModel->GetVariable(index2));

	// Set 2D-prior histogram binning to match 1D binning, for if prior was defined by histogram
	if (h2d_2dprior) {
		TAxis * xaxis = (h1_prior1) ? h1_prior1->GetXaxis() : h2d_2dprior->GetXaxis();
		TAxis * yaxis = (h1_prior2) ? h1_prior2->GetXaxis() : h2d_2dprior->GetYaxis();

		int n_xbins = xaxis->GetNbins();
		double xbins[n_xbins+1];
		xaxis -> GetLowEdge(xbins);
		xbins[n_xbins] = xaxis->GetXmax();

		int n_ybins = yaxis->GetNbins();
		double ybins[n_ybins+1];
		yaxis -> GetLowEdge(ybins);
		ybins[n_ybins] = yaxis->GetXmax();

		h2d_2dprior -> SetBins(n_xbins,xbins,n_ybins,ybins);
	}
	
	for (int i = 1; i <= h2d_2dprior->GetNbinsX(); ++i) {
		// x prior
		double x = 1;
		if (f1_prior1)
			x = f1_prior1 -> Eval(h2d_2dprior->GetXaxis()->GetBinCenter(i));
		else if (h1_prior1)
			x = h1_prior1 -> GetBinContent(h1_prior1->FindFixBin(h2d_2dprior->GetXaxis()->GetBinCenter(i)));
		
		for (int j = 1; j <= h2d_2dprior->GetNbinsY(); ++j) {
			// y prior
			double y = 1;
			if (f1_prior2)
				y = f1_prior2 -> Eval(h2d_2dprior->GetYaxis()->GetBinCenter(j));
			else if (h1_prior2)
				y = h1_prior2 -> GetBinContent(h1_prior2->FindFixBin(h2d_2dprior->GetYaxis()->GetBinCenter(j)));
			
			h2d_2dprior -> SetBinContent(i,j,x*y);
		}
	}

	if (!h2d_2dprior)
		return 0;


	// TH2D drawing options
	h2d_2dprior -> SetLineColor(kRed);
	h2d_2dprior -> SetStats(false);
	h2d_2dposterior -> SetStats(false);

	// Create BCH2D's (these normalize the TH2D's)
	BCH2D * bch2d_2dprior     = new BCH2D(h2d_2dprior);
	BCH2D * bch2d_2dposterior = new BCH2D(h2d_2dposterior);

	// Calculate integrated histograms for getting contour line values
	bch2d_2dprior     -> CalculateIntegratedHistogram();
	bch2d_2dposterior -> CalculateIntegratedHistogram();

	// Set contour levels 
	if (interval_content <= 0 or interval_content >= 1)
		interval_content = 68e-2;
	double level[1] = {bch2d_2dprior -> GetLevel(1-interval_content)};
	h2d_2dprior -> SetContour(1, level);
	h2d_2dprior -> Draw("CONT3");
	level[0] = bch2d_2dposterior -> GetLevel(1-interval_content);
	h2d_2dposterior -> SetContour(1, level);
	h2d_2dposterior -> Draw("CONT3 SAME");

	// create legend
	TLegend * legend2d = new TLegend();
	legend2d -> SetBorderSize(0);
	legend2d -> SetFillColor(kWhite);
	legend2d -> SetTextAlign(12);
	legend2d -> SetTextFont(62);
	legend2d -> SetTextSize(0.03);
	
  // create markers and arrows
	TMarker * marker_prior = new TMarker();
	marker_prior -> SetMarkerStyle(20);
	marker_prior -> SetMarkerColor(kRed);
	marker_prior -> SetMarkerSize(1.5*gPad->GetWNDC());
	
	TMarker * marker_posterior = new TMarker();
	marker_posterior -> SetMarkerStyle(20);
	marker_posterior -> SetMarkerColor(kBlue);
	marker_posterior -> SetMarkerSize(1.5*gPad->GetWNDC());
	
	TArrow * arrow = new TArrow();
	arrow -> SetArrowSize(0.02*gPad->GetWNDC());
	arrow -> SetLineColor(kBlack);

	double prior_mode_X = fPriorModel -> GetBestFitParameter(index1);
	double prior_mode_Y = fPriorModel -> GetBestFitParameter(index2);

	std::string marker_prior_text = "";
	if (const_prior1) {
		prior_mode_X = fPriorModel -> GetParameter(index1) -> GetRangeCenter();
		marker_prior_text += Form("prior(%s) constant",fPriorModel->GetParameter(index1)->GetLatexName().data());
	}
	if (const_prior2) {
		prior_mode_Y = fPriorModel -> GetParameter(index2) -> GetRangeCenter();
		if (!marker_prior_text.empty())
			marker_prior_text += ", ";
		marker_prior_text += Form("prior(%s) constant",fPriorModel->GetParameter(index2)->GetLatexName().data());
	}
	marker_prior_text = (marker_prior_text.empty()) ? "prior mode" : "prior mode* [" + marker_prior_text + "]";

	marker_prior     -> DrawMarker(prior_mode_X,prior_mode_Y);
	marker_posterior -> DrawMarker(GetBestFitParameter(index1),GetBestFitParameter(index2));
	arrow            -> DrawArrow(prior_mode_X,prior_mode_Y, GetBestFitParameter(index1), GetBestFitParameter(index2));
	
	legend2d->AddEntry(h2d_2dprior,      TString::Format("smallest %.0f%% interval(s) of prior",     100*interval_content), "L");
	legend2d->AddEntry(h2d_2dposterior,  TString::Format("smallest %.0f%% interval(s) of posterior", 100*interval_content), "L");
	legend2d->AddEntry(marker_prior,     marker_prior_text.data(), "P");
	legend2d->AddEntry(marker_posterior, "posterior mode", "P");
	legend2d->AddEntry(arrow,            "change in mode", "L");
	
	gPad->SetTopMargin(0.02);

	// place legend on top of histogram
	legend2d->SetX1NDC(gPad->GetLeftMargin());
	legend2d->SetX2NDC(1. - gPad->GetRightMargin());
	double y1 = gPad->GetTopMargin() + legend2d->GetTextSize()*legend2d->GetNRows();
	legend2d->SetY1NDC(1.-y1);
	legend2d->SetY2NDC(1. - gPad->GetTopMargin());

	legend2d->Draw();
	
	gPad->SetTopMargin(y1+0.01);
	
	gPad->RedrawAxis();
	return 1;
}

// ---------------------------------------------------------
int BCModel::PrintKnowledgeUpdatePlots(const char * filename, unsigned hdiv, unsigned vdiv, std::string options, double interval_content)
{
	// prepare prior
	GetPriorModel(true);
	// return 0 if failed
	if ( fPriorModel->GetNParameters() == 0 )
		return 0;
	fPriorModel -> MarginalizeAll();
	fPriorModel -> FindMode();

   // option flags
   bool flag_slice = false;

   // check content of options string
   if (options.find("slice") < options.size())
      flag_slice = true;

   std::string file(filename);

   // if file extension is neither .pdf nor .ps, force to .pdf
	 if ( file.rfind(".pdf") != file.size()-4 and file.rfind(".ps") != file.size()-3 )
		 file += ".pdf";

   // create canvas and prepare postscript
   TCanvas * c = new TCanvas(TString::Format("c_%s_update",fPriorModel->GetName().data()));
   c->cd();
   c->Print(std::string(file + "[").c_str());

	 if (hdiv<1) hdiv = 1;
	 if (vdiv<1) vdiv = 1;
	 int npads = hdiv * vdiv;

	 c -> Divide(hdiv,vdiv);

   // loop over all parameters and draw 1D plots
	 int ndrawn = 0;
	 int nprinted = -1;
	 c -> cd(1);
   for (unsigned i = 0; i < GetNVariables(); ++i)
		 if(DrawKnowledgeUpdatePlot1D(i, options, options)) {
			 ++ndrawn;
			 if (ndrawn!=0 and ndrawn%npads==0) {
				 c -> Print(file.c_str());
				 nprinted = ndrawn;
				 c -> Clear("D");
			 }
			 c -> cd(ndrawn%npads+1);
		 }
	 if (nprinted<ndrawn)
		 c -> Print(file.c_str());

	 c -> Clear("D");

   // loop over all parameter pairs
	 ndrawn = 0;
	 nprinted = -1;
	 c -> cd(1);
   for (unsigned i = 0; i < GetNVariables(); ++i)
		 for (unsigned j = i+1; j < GetNVariables(); ++j)
			 if (DrawKnowledgeUpdatePlot2D(i,j,flag_slice,interval_content)) {
				 ++ndrawn;
				 if (ndrawn!=0 and ndrawn%npads==0) {
					 c -> Print(file.c_str());
					 nprinted = ndrawn;
					 c -> Clear("D");
				 }
				 c -> cd(ndrawn%npads+1);
			 }
	 if (nprinted<ndrawn)
		 c -> Print(file.c_str());

   // close output
   c->Print(std::string(file + "]").c_str());
   c->Update();

   // no error
   return 1;
}
