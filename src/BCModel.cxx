/*
 * Copyright (C) 2007-2014, the BAT core developer team
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

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTree.h>

#include <fstream>
#include <iomanip>
#include <set>

// ---------------------------------------------------------
BCModel::BCModel(const char * name) :
   BCIntegrate(),
   fName((char *) name),
   fModelAPriori(0),
   fModelAPosteriori(0),
   fDataSet(0),
   fDataPointLowerBoundaries(0),
   fDataPointUpperBoundaries(0),
   fPValue(-1),
   fChi2NDoF(-1),
   fPValueNDoF(-1),
   flag_discrete(false),
   fGoFNIterationsMax(100000),
   fGoFNIterationsRun(2000),
   fGoFNChains(5),
   fPriorConstantAll(false)
{
}

// ---------------------------------------------------------
BCModel::BCModel(const BCModel & bcmodel) : BCIntegrate(bcmodel)
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
   for (unsigned i = 0; i < fParameters.Size(); ++i) {
      BCParameter * par = fParameters[i];

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
      probability *= 1. / fParameters[i]->GetRangeWidth();
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
   goftest->MCMCSetNIterationsMax(fGoFNIterationsMax);
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
   fParameters.Get(name)->Fix(value);

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

   return 0;
}

// ---------------------------------------------------------
int BCModel::SetPriorGauss(const char * name, double mean, double sigmadown, double sigmaup)
{
   // find index
   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (name == GetParameter(i)->GetName())
         index = i;

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
   if ( !fParameters.Size())
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
   // model summary
   BCLog::OutSummary(Form("Model : %s", fName.data()));
   BCLog::OutSummary(Form("Number of parameters : %u", GetNParameters()));
   BCLog::OutSummary("Parameters:");

   // parameter summary
   for (unsigned i = 0; i < GetNParameters(); i++)
      fParameters[i]->PrintSummary();

   // best fit parameters
   if ( !GetBestFitParameters().empty()) {
     BCLog::OutSummary(Form("Log of the maximum posterior: %f", GetLogMaximum()));
      BCLog::OutSummary("Best fit parameters:");

      for (unsigned i = 0; i < GetNParameters(); i++) {
        if ( fParameters[i]->Fixed() )
          BCLog::OutSummary(Form(" %s = %f (fixed)", fParameters[i]->GetName().data(), GetBestFitParameter(i)));
        else
          BCLog::OutSummary(Form(" %s = %f (global)", fParameters[i]->GetName().data(), GetBestFitParameter(i)));

        if ( fMarginalModes.size() == GetNParameters())
          BCLog::OutSummary(Form(" %s = %f (marginalized)", fParameters[i]->GetName().data(), GetBestFitParametersMarginalized()[i]));

      }
   }

   // model testing
   if (fPValue >= 0) {
      double likelihood = Likelihood(GetBestFitParameters());
      BCLog::OutSummary(" Model testing:");
      BCLog::OutSummary(Form("  p(data|lambda*) = %f", likelihood));
      BCLog::OutSummary(Form("  p-value         = %f", fPValue));
   }

   // normalization
   if (GetIntegral() > 0) {
      BCLog::OutSummary(" Evidence:");
      BCLog::OutSummary(Form(" - evidence : %f", GetIntegral()));
   }
}

// ---------------------------------------------------------
void BCModel::PrintResults(const char * file)
{
   // print summary of Markov Chain Monte Carlo

   // open file
   std::ofstream ofi(file);

   // check if file is open
   if (!ofi.is_open()) {
      std::cerr << "Couldn't open file " << file << std::endl;
      return;
   }

   // number of parameters and chains
   unsigned npar = GetNParameters();
   unsigned nchains = MCMCGetNChains();

   // check convergence
   bool flag_conv = MCMCGetNIterationsConvergenceGlobal() > 0;

   ofi << std::endl
         << " -----------------------------------------------------" << std::endl
         << " Summary" << std::endl
         << " -----------------------------------------------------" << std::endl
         << std::endl;

   ofi << " Model summary" << std::endl << " =============" << std::endl
         << " Model: " << fName.data() << std::endl
         << " Number of parameters: " << npar << std::endl
         << " List of Parameters and ranges:" << std::endl;
   for (unsigned i = 0; i < npar; ++i) {
     ofi << "  (" << i << ") Parameter \""
         << fParameters[i]->GetName() << "\"" << ": "
         << "[" << fParameters[i]->GetLowerLimit() << ", "
         << fParameters[i]->GetUpperLimit() << "]";
     if (fParameters[i]->Fixed()) {
       ofi << " (fixed)";
     }
     ofi << std::endl;
   }
   ofi << std::endl;

   ofi << " Results of the optimization" << std::endl
         << " ===========================" << std::endl
         << " Optimization algorithm used: "
         << DumpUsedOptimizationMethod()<< std::endl;

   if ( ! GetBestFitParameters().empty()) {
     ofi << " Log of the maximum posterior: " << GetLogMaximum() << std::endl;
      ofi << " List of parameters and global mode:" << std::endl;
      for (unsigned i = 0; i < npar; ++i) {
        ofi << "  (" << i << ") Parameter \""
            << fParameters[i]->GetName() << "\": "
            << GetBestFitParameter(i);
        if (fParameters[i]->Fixed()) {
          ofi << " (fixed)";
        }
        else if (GetBestFitParameterErrors().size() == npar) {
          if(GetBestFitParameterError(i) >= 0.)
            ofi << " +- " << GetBestFitParameterError(i);
          else
            ofi << " (no error estimate available) ";
        }
        else {
          ofi << " (no error estimate available) ";
        }
        ofi << std::endl;
      }
      ofi << std::endl;
   }
   else {
      ofi << " No best fit information available." << std::endl;
      ofi << std::endl;
   }

   if (fPValue >= 0.) {
      ofi << " Results of the model test" << std::endl
            << " =========================" << std::endl
            << " p-value: " << fPValue << std::endl;
      if (fPValueNDoF >= 0)
         ofi << " p-value corrected for degrees of freedom: " << fPValueNDoF << std::endl;

      ofi << std::endl;
   }

   if (GetIntegral() >= 0.) {
      ofi << " Results of the integration" << std::endl
            << " ============================" << std::endl
            << " Integration method used: "
            << DumpUsedIntegrationMethod() << std::endl;
      ofi << " Evidence: " << GetIntegral();
      if (GetError() >= 0)
        ofi << " +- " << GetError() << std::endl;
      else
        ofi << " (no error estimate available) " << std::endl;
      ofi << std::endl;
   }

   // give warning if MCMC did not converge
   if (!flag_conv && fMCMCFlagRun)
      ofi << " WARNING: the Markov Chain did not converge!" << std::endl
      << " Be cautious using the following results!" << std::endl
      << std::endl;

   // print results of marginalization (if MCMC was run)
   if (fFlagMarginalized) {
      ofi << " Results of the marginalization" << std::endl
          << " ==============================" << std::endl
          << " Marginalization algorithm used: "
          << DumpUsedMarginalizationMethod() << std::endl
          << " List of parameters and properties of the marginalized"
          << std::endl << " distributions:" << std::endl;
      for (unsigned i = 0; i < npar; ++i) {
         if ( ! fParameters[i]->FillHistograms())
            continue;

         // get marginalized histogram
         BCH1D * bch1d = GetMarginalized(fParameters[i]);

         ofi << "  (" << i << ") Parameter \""
             << fParameters[i]->GetName() << "\":";

         if (!bch1d) {
            ofi << " fixed (or histogram does not exist) " << std::endl;
            continue;
         }
         else
            ofi << std::endl;

         ofi << "      Mean +- sqrt(V):                " << std::setprecision(4)
             << bch1d->GetMean() << " +- " << std::setprecision(4)
             << bch1d->GetRMS() << std::endl

             << "      Median +- central 68% interval: "
             << std::setprecision(4) << bch1d->GetMedian() << " +  "
             << std::setprecision(4) << bch1d->GetQuantile(0.84) - bch1d->GetMedian()
             << " - " << std::setprecision(4)
             << bch1d->GetMedian() - bch1d->GetQuantile(0.16) << std::endl

             << "      (Marginalized) mode:            " << bch1d->GetMode() << std::endl;

         ofi << "       5% quantile:                   " << std::setprecision(4)
             << bch1d->GetQuantile(0.05) << std::endl
             << "      10% quantile:                   " << std::setprecision(4)
             << bch1d->GetQuantile(0.10) << std::endl
             << "      16% quantile:                   " << std::setprecision(4)
             << bch1d->GetQuantile(0.16) << std::endl
             << "      84% quantile:                   " << std::setprecision(4)
             << bch1d->GetQuantile(0.85) << std::endl
             << "      90% quantile:                   " << std::setprecision(4)
             << bch1d->GetQuantile(0.90) << std::endl
             << "      95% quantile:                   " << std::setprecision(4)
             << bch1d->GetQuantile(0.95) << std::endl;

         std::vector<double> v;
         v = bch1d->GetSmallestIntervals(0.68);
         ofi << "      Smallest interval(s) containing at least 68% and local mode(s):"
             << std::endl;
         for (unsigned j = 0; j < v.size(); j += 5)
            ofi << "       (" << v[j] << ", " << v[j + 1]
                << ") (local mode at " << v[j + 3] << " with rel. height "
                << v[j + 2] << "; rel. area " << v[j + 4] << ")"
                << std::endl;
         ofi << std::endl;
      }
   }
   if (fMCMCFlagRun) {
      ofi << " Status of the MCMC" << std::endl << " ==================" << std::endl
          << " Convergence reached:                    " << (flag_conv ? "yes" : "no")
          << std::endl;

      if (flag_conv)
         ofi << " Number of iterations until convergence: "
             << MCMCGetNIterationsConvergenceGlobal() << std::endl;
      else
         ofi << " WARNING: the Markov Chain did not converge! Be\n"
             << " cautious using the following results!" << std::endl
             << std::endl;
      ofi << " Number of chains:                       " << MCMCGetNChains() << std::endl
          << " Number of iterations per chain:         " << MCMCGetNIterationsRun() << std::endl
          << " Average pre-run efficiencies:" << std::endl;

      std::vector<double> efficiencies;
      efficiencies.assign(npar, 0.);

      for (unsigned ipar = 0; ipar < npar; ++ipar)
        for (unsigned ichain = 0; ichain < nchains; ++ichain) {
          efficiencies[ipar] +=
            fMCMCEfficiencies[ichain*npar+ipar] / double(nchains) * 100.;
         }

      for (unsigned ipar = 0; ipar < npar; ++ipar)
         ofi << "  (" << ipar << ") Parameter \""
             << fParameters[ipar]->GetName().data() << "\": "
             << efficiencies.at(ipar) << "%" << std::endl;
      ofi << std::endl;
   }

   ofi << " -----------------------------------------------------" << std::endl
       << " Notation:" << std::endl
       << " Mean        : mean value of the marg. pdf" << std::endl
       << " Median      : median of the marg. pdf" << std::endl
       << " Marg. mode  : most probable value of the marg. pdf" << std::endl
       << " V           : Variance of the marg. pdf" << std::endl
       << " Quantiles   : most commonly used quantiles" <<std::endl
       << " -----------------------------------------------------" << std::endl
       << std::endl;

   // close file
   ofi.close();
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
   for (unsigned int i = 0; i < GetNParameters(); ++i)
      BCLog::OutSummary(Form("      %s = %.3g", GetParameter(i)->GetName().data(), GetBestFitParameter(i)));

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
         double hessianmatrixelement = HessianMatrixElement(
               fParameters[i], fParameters[j], parameters);

         // print to screen
         BCLog::OutSummary(Form("%d %d : %f", i, j, hessianmatrixelement));
      }
}

// ---------------------------------------------------------
BCDataPoint * BCModel::VectorToDataPoint(const std::vector<double> &data)
{
   BCDataPoint * datapoint = new BCDataPoint(data.size());
   datapoint->SetValues(data);
   return datapoint;
}

// ---------------------------------------------------------
int BCModel::CompareStrings(const char * string1, const char * string2)
{
   int flag_same = 0;

   if (strlen(string1) != strlen(string2))
      return -1;

   for (int i = 0; i < int(strlen(string1)); i++)
      if (string1[i] != string2[i])
         flag_same = -1;

   return flag_same;
}
