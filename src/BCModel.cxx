/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCModel.h"

#include "BCDataPoint.h"
#include "BCDataSet.h"
#include "BCParameter.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCGoFTest.h"
#include "BCLog.h"
#include "BCErrorCodes.h"
#include "BCMath.h"
#include "BCModelOutput.h"

#include <TROOT.h>
#include <TNamed.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TMath.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TRandom3.h>
#include <TF1.h>
#include <Math/QuantFuncMathCore.h>

#include <algorithm>
#include <set>

#include <fstream>
#include <iomanip>

// ---------------------------------------------------------
BCModel::BCModel(const char * name)
   : BCIntegrate()
{
   fNormalization = -1.;
   fDataSet = 0;
   fParameterSet = new BCParameterSet;

   fIndex = -1;
   fPValue = -1;
   fPValueNDoF = -1;
   fChi2NDoF = -1;

   fName = (char *) name;

   fDataPointUpperBoundaries = 0;
   fDataPointLowerBoundaries = 0;

   fErrorBandXY = 0;

   fGoFNChains = 5;
   fGoFNIterationsMax = 100000;
   fGoFNIterationsRun = 2000;

   flag_discrete = false;

   fPriorConstantAll = false;
   fPriorConstantValue = 0;
}

// ---------------------------------------------------------
BCModel::BCModel()
   : BCIntegrate()
{
   fNormalization = -1.;
   fDataSet = 0;
   fParameterSet = new BCParameterSet();

   fIndex = -1;
   fPValue = -1;
   fPValueNDoF = -1;
   fChi2NDoF = -1;

   fName = "model";
   fDataPointUpperBoundaries = 0;
   fDataPointLowerBoundaries = 0;

   fGoFNChains = 5;
   fGoFNIterationsMax = 100000;
   fGoFNIterationsRun = 2000;

   flag_discrete = false;

   fPriorConstantAll = false;
   fPriorConstantValue = 0;
}

// ---------------------------------------------------------
BCModel::BCModel(const BCModel & bcmodel) : BCIntegrate(bcmodel)
{
   fIndex                           = bcmodel.fIndex;
   fName                            = bcmodel.fName;
   fModelAPriori                    = bcmodel.fModelAPriori;
   fModelAPosteriori                = bcmodel.fModelAPosteriori;
   for (int i = 0; i < int(bcmodel.fParameterSet->size()); ++i) {
      if (bcmodel.fParameterSet->at(i)) {
         fParameterSet->push_back(new BCParameter(*(bcmodel.fParameterSet->at(i))));

      }
      else
         fParameterSet->push_back(0);
   }
   if (fDataSet)
      fDataSet = bcmodel.fDataSet;
   else
      fDataSet = 0;
   if (bcmodel.fNDataPointsMinimum)
      fNDataPointsMinimum = bcmodel.fNDataPointsMinimum;
   else
      fNDataPointsMinimum = 0;
   if (bcmodel.fNDataPointsMaximum)
      fNDataPointsMaximum = bcmodel.fNDataPointsMaximum;
   else
      fNDataPointsMaximum = 0;
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
   fPriorConstantValue              = bcmodel.fPriorConstantValue;
   fPriorContainerConstant          = bcmodel.fPriorContainerConstant;
   fPriorContainerInterpolate       = bcmodel.fPriorContainerInterpolate;
   fNormalization                   = bcmodel.fNormalization;
}

// ---------------------------------------------------------
BCModel::~BCModel()
{
   for (unsigned int i = 0; i < GetNParameters(); ++i)
      delete fPriorContainer[i];
   fPriorContainer.clear();

   delete fParameterSet;

   delete fDataPointLowerBoundaries;
   delete fDataPointUpperBoundaries;
}

// ---------------------------------------------------------
BCModel & BCModel::operator = (const BCModel & bcmodel)
{
   BCIntegrate::operator=(bcmodel);
   fIndex                           = bcmodel.fIndex;
   fName                            = bcmodel.fName;
   fModelAPriori                    = bcmodel.fModelAPriori;
   fModelAPosteriori                = bcmodel.fModelAPosteriori;
   for (int i = 0; i < int(bcmodel.fParameterSet->size()); ++i) {
      if (bcmodel.fParameterSet->at(i)) {
         fParameterSet->push_back(new BCParameter(*(bcmodel.fParameterSet->at(i))));

      }
      else
         fParameterSet->push_back(0);
   }
   if (fDataSet)
      fDataSet = bcmodel.fDataSet;
   else
      fDataSet = 0;
   if (bcmodel.fNDataPointsMinimum)
      fNDataPointsMinimum = bcmodel.fNDataPointsMinimum;
   else
      fNDataPointsMinimum = 0;
   if (bcmodel.fNDataPointsMaximum)
      fNDataPointsMaximum = bcmodel.fNDataPointsMaximum;
   else
      fNDataPointsMaximum = 0;
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
   fPriorConstantValue              = bcmodel.fPriorConstantValue;
   fPriorContainerConstant          = bcmodel.fPriorContainerConstant;
   fPriorContainerInterpolate       = bcmodel.fPriorContainerInterpolate;
   fNormalization                   = bcmodel.fNormalization;

   // return this
   return *this;
}

// ---------------------------------------------------------
int BCModel::GetNDataPoints()
{
   int npoints = 0;
   if (fDataSet)
      npoints = fDataSet->GetNDataPoints();
   else {
      BCLog::OutWarning("BCModel::GetNDataPoints() : No data set defined.");
      return ERROR_NOEVENTS;
   }

   return npoints;
}

// ---------------------------------------------------------
BCDataPoint * BCModel::GetDataPoint(unsigned int index)
{
   if (fDataSet)
      return fDataSet->GetDataPoint(index);

   BCLog::OutWarning("BCModel::GetDataPoint : No data set defined.");
   return 0;
}

// ---------------------------------------------------------
BCParameter * BCModel::GetParameter(int index)
{
   if (!fParameterSet)
      return 0;

   if (index < 0 || index >= (int)GetNParameters()) {
      BCLog::OutWarning(
            Form("BCModel::GetParameter : Parameter index %d not within range.",index));
      return 0;
   }

   return fParameterSet->at(index);
}

// ---------------------------------------------------------
BCParameter * BCModel::GetParameter(const char * name)
{
   if (!fParameterSet)
      return 0;

   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (name == GetParameter(i)->GetName())
         index = i;

   if (index < 0) {
      BCLog::OutWarning(
            Form("BCModel::GetParameter : Model %s has no parameter named '%s'",
                 GetName().data(), name));
      return 0;
   }

   return GetParameter(index);
}

// ---------------------------------------------------------
double BCModel::GetBestFitParameter(unsigned int index)
{
   if(index >= GetNParameters()) {
      BCLog::OutError("BCModel::GetBestFitParameter : Parameter index out of range, returning -1e+111.");
      return -1e+111;
   }

   if(fBestFitParameters.size()==0) {
      BCLog::OutError("BCModel::GetBestFitParameter : Mode finding not yet run, returning center of the range.");
      return (GetParameter(index)->GetUpperLimit() + GetParameter(index)->GetLowerLimit() ) / 2.;
   }

   return fBestFitParameters[index];
}

// ---------------------------------------------------------
double BCModel::GetBestFitParameterError(unsigned int index)
{
   if(index >= GetNParameters()) {
      BCLog::OutError("BCModel::GetBestFitParameterError : Parameter index out of range, returning -1.");
      return -1;
   }

   if(fBestFitParameterErrors.size()==0) {
      BCLog::OutError("BCModel::GetBestFitParameterError : Mode finding not yet run, returning -1.");
      return -1.;
   }

   if(fBestFitParameterErrors[index]<0.)
      BCLog::OutWarning("BCModel::GetBestFitParameterError : Parameter error not available, returning -1.");

   return fBestFitParameterErrors[index];
}

// ---------------------------------------------------------
double BCModel::GetBestFitParameterMarginalized(unsigned int index)
{
   if(index >= GetNParameters()) {
      BCLog::OutError("BCModel::GetBestFitParameterMarginalized : Parameter index out of range, returning -1e+111.");
      return -1e+111;
   }

   if(fBestFitParametersMarginalized.size()==0) {
      BCLog::OutError("BCModel::GetBestFitParameterMarginalized : MCMC not yet run, returning center of the range.");
      return (GetParameter(index)->GetUpperLimit() + GetParameter(index)->GetLowerLimit() ) / 2.;
   }

   return fBestFitParametersMarginalized[index];
}

// ---------------------------------------------------------
void BCModel::SetNbins(const char * parname, int nbins)
{
   BCParameter * p = GetParameter(parname);
   if (!p) {
      BCLog::OutWarning(
            Form("BCModel::SetNbins : Parameter '%s' not found so Nbins not set",parname));
      return;
   }

   BCIntegrate::SetNbins(nbins, p->GetIndex());
}

// ---------------------------------------------------------
std::vector<double> BCModel::GetErrorBand(double level)
{
   std::vector<double> errorband;

   if (!fErrorBandXY)
      return errorband;

   int nx = fErrorBandXY->GetNbinsX();
   errorband.assign(nx, 0.);

   // loop over x and y bins
   for (int ix = 1; ix <= nx; ix++) {
      TH1D * temphist = fErrorBandXY->ProjectionY("temphist", ix, ix);

      int nprobSum = 1;
      double q[1];
      double probSum[1];
      probSum[0] = level;

      temphist->GetQuantiles(nprobSum, q, probSum);

      errorband[ix - 1] = q[0];
   }

   return errorband;
}

// ---------------------------------------------------------
TGraph * BCModel::GetErrorBandGraph(double level1, double level2)
{
   if (!fErrorBandXY)
      return 0;

   // define new graph
   int nx = fErrorBandXY->GetNbinsX();

   TGraph * graph = new TGraph(2 * nx);
   graph->SetFillStyle(1001);
   graph->SetFillColor(kYellow);

   // get error bands
   std::vector<double> ymin = GetErrorBand(level1);
   std::vector<double> ymax = GetErrorBand(level2);

   for (int i = 0; i < nx; i++) {
      graph->SetPoint(i, fErrorBandXY->GetXaxis()->GetBinCenter(i + 1), ymin[i]);
      graph->SetPoint(nx + i, fErrorBandXY->GetXaxis()->GetBinCenter(nx - i), ymax[nx - i - 1]);
   }

   return graph;
}

// ---------------------------------------------------------
TH2D * BCModel::GetErrorBandXY_yellow(double level, int nsmooth)
{
   if (!fErrorBandXY)
      return 0;

   int nx = fErrorBandXY->GetNbinsX();
   int ny = fErrorBandXY->GetNbinsY();

   // copy existing histogram
   TH2D * hist_tempxy = (TH2D*) fErrorBandXY->Clone(
         TString::Format("%s_sub_%f.2", fErrorBandXY->GetName(), level));
   hist_tempxy->Reset();
   hist_tempxy->SetFillColor(kYellow);

   // loop over x bins
   for (int ix = 1; ix < nx; ix++) {
      BCH1D * hist_temp = new BCH1D();

      TH1D * hproj = fErrorBandXY->ProjectionY("temphist", ix, ix);
      if (nsmooth > 0)
         hproj->Smooth(nsmooth);

      hist_temp->SetHistogram(hproj);

      TH1D * hist_temp_yellow = hist_temp->GetSmallestIntervalHistogram(level);

      for (int iy = 1; iy <= ny; ++iy)
         hist_tempxy->SetBinContent(ix, iy, hist_temp_yellow->GetBinContent(iy));

      delete hist_temp_yellow;
      delete hist_temp;
   }

   return hist_tempxy;
}

// ---------------------------------------------------------
TGraph * BCModel::GetFitFunctionGraph(const std::vector<double> &parameters)
{
   if (!fErrorBandXY)
      return 0;

   // define new graph
   int nx = fErrorBandXY->GetNbinsX();
   TGraph * graph = new TGraph(nx);

   // loop over x values
   for (int i = 0; i < nx; i++) {
      double x = fErrorBandXY->GetXaxis()->GetBinCenter(i + 1);

      std::vector<double> xvec;
      xvec.push_back(x);
      double y = FitFunction(xvec, parameters);
      xvec.clear();

      graph->SetPoint(i, x, y);
   }

   return graph;
}

// ---------------------------------------------------------
TGraph * BCModel::GetFitFunctionGraph(const std::vector<double> &parameters, double xmin, double xmax, int n)
{
   // define new graph
   TGraph * graph = new TGraph(n + 1);

   double dx = (xmax - xmin) / (double) n;

   // loop over x values
   for (int i = 0; i <= n; i++) {
      double x = (double) i * dx + xmin;
      std::vector<double> xvec;
      xvec.push_back(x);
      double y = FitFunction(xvec, parameters);

      xvec.clear();

      graph->SetPoint(i, x, y);
   }

   return graph;
}

// ---------------------------------------------------------
bool BCModel::GetFlagBoundaries()
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
void BCModel::SetErrorBandContinuous(bool flag)
{
   fErrorBandContinuous = flag;

   if (flag)
      return;

   // clear x-values
   fErrorBandX.clear();

   // copy data x-values
   for (unsigned int i = 0; i < fDataSet->GetNDataPoints(); ++i)
      fErrorBandX.push_back(fDataSet->GetDataPoint(i)->GetValue(fFitFunctionIndexX));
}

// ---------------------------------------------------------
int BCModel::AddParameter(const char * name, double lowerlimit, double upperlimit)
{
   // create new parameter
   BCParameter * parameter = new BCParameter(name, lowerlimit, upperlimit);

   int flag_ok = AddParameter(parameter);
   if (flag_ok)
      delete parameter;

   return flag_ok;
}

// ---------------------------------------------------------
int BCModel::AddParameter(BCParameter * parameter)
{
   // check if parameter set exists
   if (!fParameterSet) {
      BCLog::OutError("BCModel::AddParameter : Parameter set does not exist");
      return ERROR_PARAMETERSETDOESNOTEXIST;
   }

   // check if parameter with same name exists
   int flag_exists = 0;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (CompareStrings(parameter->GetName().data(), GetParameter(i)->GetName().data()) == 0)
         flag_exists = -1;

   if (flag_exists < 0) {
      BCLog::OutError(Form(
            "BCModel::AddParameter : Parameter with name %s exists already. ",
            parameter->GetName().data()));
      return ERROR_PARAMETEREXISTSALREADY;
   }

   // define index of new parameter
   unsigned int index = fParameterSet->size();
   parameter->SetIndex(index);

   // add parameter to parameter container
   fParameterSet->push_back(parameter);

   // add empty object to prior container
   fPriorContainer.push_back(0);

   // don't interpolate the prior histogram by default
   fPriorContainerInterpolate.push_back(false);

   // prior assumed to be non-constant in general case
   fPriorContainerConstant.push_back(false);

   // add parameters to integration methods
   SetParameters(fParameterSet);

   // reset results
   ResetResults();

   return 0;
}

// ---------------------------------------------------------
double BCModel::LogProbabilityNN(const std::vector<double> &parameters)
{
   // add log of likelihood
   double logprob = LogLikelihood(parameters);

   // add log of prior probability
   logprob += LogAPrioriProbability(parameters);

   return logprob;
}

// ---------------------------------------------------------
double BCModel::LogProbability(const std::vector<double> &parameters)
{
   // check if normalized
   if (fNormalization < 0. || fNormalization == 0.) {
      BCLog::OutError("BCModel::LogProbability. Normalization not available or zero.");
      return 0.;
   }

   return LogProbabilityNN(parameters) - log(fNormalization);
}

// ---------------------------------------------------------
double BCModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
   double logprob = 0;

   if(!fPriorConstantAll) {

      // get number of parameters
      int npar = GetNParameters();

      // loop over all 1-d priors
      for (int i = 0; i < npar; ++i) {
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
                  "is defined but not recodnized.",
                  GetParameter(i)->GetName().c_str())); // this should never happen
         }
         // use constant only if user has defined it
         else if (!fPriorContainerConstant[i]) {
            BCLog::OutWarning(Form(
                  "BCModel::LogAPrioriProbability: Prior for parameter %s "
                  "is undefined. Using constant prior to proceed.",
                  GetParameter(i)->GetName().c_str()));
            logprob -= log(GetParameter(i)->GetRangeWidth());
         }
      }
   }

   // add the contribution from constant priors in one step
   logprob += fPriorConstantValue;

   return logprob;
}

// ---------------------------------------------------------
double BCModel::LogEval(const std::vector<double> &parameters)
{
   return LogProbabilityNN(parameters);
}

// ---------------------------------------------------------
double BCModel::EvalSampling(const std::vector<double> &parameters)
{
   return SamplingFunction(parameters);
}

// ---------------------------------------------------------
double BCModel::SamplingFunction(const std::vector<double> & /*parameters*/)
{
   double probability = 1.;
   for (std::vector<BCParameter *>::const_iterator it = fParameterSet->begin(); it != fParameterSet->end(); ++it)
      probability *= 1. / ((*it)->GetUpperLimit() - (*it)->GetLowerLimit());
   return probability;
}

// ---------------------------------------------------------
double BCModel::Normalize()
{
   if(fParameterSet->size()<1) {
      BCLog::OutError(Form("Normalize : No parameters defined in model \'%s\'. Aborting.",GetName().data()));
      return -1.;
   }

   BCLog::OutSummary(Form("Model \'%s\': Normalizing probability",GetName().data()));

   unsigned int n = GetNvar();

   // initialize BCIntegrate if not done already
   if (n == 0) {
      SetParameters(fParameterSet);
      n = GetNvar();
   }

   // integrate and get best fit parameters
   // maybe we have to remove the mode finding from here in the future
   fNormalization = Integrate();

   BCLog::OutDetail(Form(" --> Normalization factor : %.6g", fNormalization));

   return fNormalization;
}

// ---------------------------------------------------------
int BCModel::CheckParameters(const std::vector<double> &parameters)
{
   // check if vectors are of equal size
   if (!parameters.size() == fParameterSet->size())
      return ERROR_INVALIDNUMBEROFPARAMETERS;

   // check if parameters are within limits
   for (unsigned int i = 0; i < fParameterSet->size(); i++) {
      BCParameter * modelparameter = fParameterSet->at(i);

      if (modelparameter->GetLowerLimit() > parameters.at(i)
          || modelparameter->GetUpperLimit() < parameters.at(i)) {
         BCLog::OutError(Form(
               "BCModel::CheckParameters : Parameter %s not within limits.",
               fParameterSet-> at(i)-> GetName().data()));
         return ERROR_PARAMETERNOTWITHINRANGE;
      }
   }

   return 0;
}

// ---------------------------------------------------------
void BCModel::FindMode(std::vector<double> start)
{
   if(fParameterSet->size()<1) {
      BCLog::OutError(Form("FindMode : No parameters defined in model \'%s\'. Aborting.",GetName().data()));
      return;
   }

   // this implementation is CLEARLY not good we have to work on this.

   BCLog::OutSummary(Form("Model \'%s\': Finding mode using %s",GetName().data(), DumpOptimizationMethod().c_str()));

   // synchronize parameters in BCIntegrate
   SetParameters(fParameterSet);

   switch (GetOptimizationMethod()) {
      case BCIntegrate::kOptSA:
         FindModeSA(start);
         return;

      case BCIntegrate::kOptMinuit: {
         int printlevel = -1;
         if (BCLog::GetLogLevelScreen() <= BCLog::detail)
            printlevel = 0;
         if (BCLog::GetLogLevelScreen() <= BCLog::debug)
            printlevel = 1;
         BCIntegrate::FindModeMinuit(start, printlevel);
         return;
      }

      case BCIntegrate::kOptMetropolis:
          MarginalizeAll();
         return;

      default:
         BCLog::OutError(Form("BCModel::FindMode : Invalid mode finding method: %d",GetOptimizationMethod()));
         break;
   }

   return;
}

// ---------------------------------------------------------
void BCModel::FindModeMinuit(std::vector<double> start, int printlevel)
{
   if(fParameterSet->size()<1) {
      BCLog::OutError(Form("FindModeMinuit : No parameters defined in model \'%s\'. Aborting.",GetName().data()));
      return;
   }

   // synchronize parameters in BCIntegrate
   SetParameters(fParameterSet);
   BCIntegrate::FindModeMinuit(start, printlevel);
}

// ---------------------------------------------------------
void BCModel::WriteMode(const char * file)
{
	std::ofstream ofi(file);
   if (!ofi.is_open()) {
      std::cerr << "Couldn't open file " << file << std::endl;
      return;
   }

   int npar = fParameterSet->size();
   for (int i = 0; i < npar; i++)
      ofi << fBestFitParameters.at(i) << std::endl;

   ofi << std::endl;
   ofi << "#######################################################################"
       << std::endl;
   ofi << "#" << std::endl;
   ofi << "#  This file was created automatically by BCModel::WriteMode() call."
       << std::endl;
   ofi << "#  It can be read in by call to BCModel::ReadMode()." << std::endl;
   ofi << "#  Do not modify it unless you know what you're doing." << std::endl;
   ofi << "#" << std::endl;
   ofi << "#######################################################################"
       << std::endl;
   ofi << "#" << std::endl;
   ofi << "#  Best fit parameters (mode) for model:" << std::endl;
   ofi << "#  \'" << fName.data() << "\'" << std::endl;
   ofi << "#" << std::endl;
   ofi << "#  Number of parameters: " << npar << std::endl;
   ofi << "#  Parameters ordered as above:" << std::endl;

   for (int i = 0; i < npar; i++) {
      ofi << "#     " << i << ": ";
      ofi << fParameterSet->at(i)->GetName().data() << " = ";
      ofi << fBestFitParameters.at(i) << std::endl;
   }

   ofi << "#" << std::endl;
   ofi << "########################################################################"
       << std::endl;
}

// ---------------------------------------------------------
int BCModel::ReadMode(const char * file)
{
	std::ifstream ifi(file);
   if (!ifi.is_open()) {
      BCLog::OutError(Form("BCModel::ReadMode : Couldn't open file %s.", file));
      return 0;
   }

   int npar = fParameterSet->size();
   std::vector<double> mode;

   int i = 0;
   while (i < npar && !ifi.eof()) {
      double a;
      ifi >> a;
      mode.push_back(a);
      i++;
   }

   if (i < npar) {
      BCLog::OutError(Form("BCModel::ReadMode : Couldn't read mode from file %s.", file));
      BCLog::OutError(Form("BCModel::ReadMode : Expected %d parameters, found %d.", npar, i));
      return 0;
   }

   BCLog::OutSummary(
         Form("#  Read in best fit parameters (mode) for model \'%s\' from file %s:",
               fName.data(), file));
   SetMode(mode);
   for (int j = 0; j < npar; j++)
      BCLog::OutSummary(Form("#   ->Parameter %d : %s = %e", j,
            fParameterSet->at(j)->GetName().data(), fBestFitParameters[j]));

   BCLog::OutWarning("#  ! Best fit values obtained before this call will be overwritten !");

   return npar;
}

// ---------------------------------------------------------
int BCModel::MarginalizeAll()
{
   if(fParameterSet->size()<1) {
      BCLog::OutError(Form("MarginalizeAll : No parameters defined in model \'%s\'. Aborting.",GetName().data()));
      return 0;
   }

   BCLog::OutSummary(Form("Running MCMC for model \'%s\'",GetName().data()));

   // prepare function fitting
   double dx = 0.;
   double dy = 0.;

   if (fFitFunctionIndexX >= 0) {
      dx = (fDataPointUpperBoundaries->GetValue(fFitFunctionIndexX)
            - fDataPointLowerBoundaries->GetValue(fFitFunctionIndexX))
            / (double) fErrorBandNbinsX;

      dy = (fDataPointUpperBoundaries->GetValue(fFitFunctionIndexY)
            - fDataPointLowerBoundaries->GetValue(fFitFunctionIndexY))
            / (double) fErrorBandNbinsY;

      fErrorBandXY
            = new TH2D(TString::Format("errorbandxy_%d", BCLog::GetHIndex()), "",
                  fErrorBandNbinsX,
                  fDataPointLowerBoundaries->GetValue(fFitFunctionIndexX) - .5 * dx,
                  fDataPointUpperBoundaries->GetValue(fFitFunctionIndexX) + .5 * dx,
                  fErrorBandNbinsY,
                  fDataPointLowerBoundaries->GetValue(fFitFunctionIndexY) - .5 * dy,
                  fDataPointUpperBoundaries->GetValue(fFitFunctionIndexY) + .5 * dy);
      fErrorBandXY->SetStats(kFALSE);

      for (int ix = 1; ix <= fErrorBandNbinsX; ++ix)
         for (int iy = 1; iy <= fErrorBandNbinsX; ++iy)
            fErrorBandXY->SetBinContent(ix, iy, 0.);
   }

   // run the Markov chains
   MCMCMetropolis();

   // store the mode found by the chains
   StoreMode();

   return 1;
}

// ---------------------------------------------------------
BCH1D * BCModel::GetMarginalized(BCParameter * parameter)
{
   if (!parameter) {
// don't print any error message, should be done upstream
//      BCLog::OutError("BCModel::GetMarginalized : Parameter does not exist.");
      return 0;
   }

   int index = parameter->GetIndex();
   if (fMCMCFlagsFillHistograms[index] == false) {
// don't print any error message, should be done upstream
//      BCLog::OutError(Form(
//            "BCModel::GetMarginalized : Distribuion for '%s' not filled.",
//            parameter->GetName().data()));
      return 0;
   }

   // get histogram
   TH1D * hist = MCMCGetH1Marginalized(index);
   if (!hist)
      return 0;

   // set axis labels
   hist->SetName(Form("hist_%s_%s", GetName().data(), parameter->GetName().data()));
   hist->SetXTitle(parameter->GetName().data());
   hist->SetYTitle(Form("p(%s|data)", parameter->GetName().data()));
   hist->SetStats(kFALSE);

   // set histogram
   BCH1D * hprob = new BCH1D();
   hprob->SetHistogram(hist);

   // set best fit parameter
   double bestfit = hprob->GetMode();

   if (fBestFitParametersMarginalized.empty())
       fBestFitParametersMarginalized.assign(GetNParameters(), 0.0);

   fBestFitParametersMarginalized[index] = bestfit;

   hprob->SetGlobalMode(fBestFitParameters.at(index));

   return hprob;
}

// ---------------------------------------------------------
int BCModel::ReadMarginalizedFromFile(const char * file)
{
   TFile * froot = new TFile(file);
   if (!froot->IsOpen()) {
      BCLog::OutError(Form("BCModel::ReadMarginalizedFromFile : Couldn't open file %s.", file));
      return 0;
   }

   // We reset the MCMCEngine here for the moment.
   // In the future maybe we only want to do this if the engine
   // wans't initialized at all or when there were some changes
   // in the model.
   // But maybe we want reset everything since we're overwriting
   // the marginalized distributions anyway.
   MCMCInitialize();

   int k = 0;
   int n = GetNParameters();
   for (int i = 0; i < n; i++) {
      BCParameter * a = GetParameter(i);
      TKey * key = froot->GetKey(TString::Format("hist_%s_%s",GetName().data(), a->GetName().data()));
      if (key) {
         TH1D * h1 = (TH1D*) key->ReadObjectAny(TH1D::Class());
         h1->SetDirectory(0);
         if (SetMarginalized(i, h1))
            k++;
      }
      else
         BCLog::OutWarning(
               Form("BCModel::ReadMarginalizedFromFile : Couldn't read histogram \"hist_%s_%s\" from file %s.",
                     GetName().data(), a->GetName().data(), file));
   }

   for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
         BCParameter * a = GetParameter(i);
         BCParameter * b = GetParameter(j);
         TKey * key = froot->GetKey(
               TString::Format("hist_%s_%s_%s",GetName().data(), a->GetName().data(),b->GetName().data()));
         if (key) {
            TH2D * h2 = (TH2D*) key->ReadObjectAny(TH2D::Class());
            h2->SetDirectory(0);
            if (SetMarginalized(i, j, h2))
               k++;
         }
         else
            BCLog::OutWarning(
                  Form("BCModel::ReadMarginalizedFromFile : Couldn't read histogram \"hist_%s_%s_%s\" from file %s.",
                        GetName().data(), a->GetName().data(), b->GetName().data(), file));
      }
   }

   froot->Close();

   return k;
}

// ---------------------------------------------------------
int BCModel::ReadErrorBandFromFile(const char * file)
{
   TFile * froot = new TFile(file);
   if (!froot->IsOpen()) {
      BCLog::OutError(Form("BCModel::ReadErrorBandFromFile. Couldn't open file %s.", file));
      return 0;
   }

   int r = 0;

   TH2D * h2 = (TH2D*) froot->Get("errorbandxy");
   if (h2) {
      h2->SetDirectory(0);
      h2->SetName(TString::Format("errorbandxy_%d", BCLog::GetHIndex()));
      SetErrorBandHisto(h2);
      r = 1;
   }
   else
      BCLog::OutWarning(
            Form("BCModel::ReadErrorBandFromFile : Couldn't read histogram \"errorbandxy\" from file %s.",file));

   froot->Close();

   return r;
}

// ---------------------------------------------------------
int BCModel::PrintAllMarginalized1D(const char * filebase)
{
   if (fMCMCH1Marginalized.size() == 0) {
      BCLog::OutError("BCModel::PrintAllMarginalized : Marginalized distributions not available.");
      return 0;
   }

   int n = GetNParameters();
   for (int i = 0; i < n; i++) {
      BCParameter * a = GetParameter(i);
      if (GetMarginalized(a))
         GetMarginalized(a)->Print(Form("%s_1D_%s.ps", filebase, a->GetName().data()));
   }

   return n;
}

// ---------------------------------------------------------
int BCModel::PrintAllMarginalized2D(const char * filebase)
{
   if (fMCMCH2Marginalized.size() == 0) {
      BCLog::OutError("BCModel::PrintAllMarginalized : Marginalized distributions not available.");
      return 0;
   }

   int k = 0;
   int n = GetNParameters();
   for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
         BCParameter * a = GetParameter(i);
         BCParameter * b = GetParameter(j);

         double meana = (a->GetLowerLimit() + a->GetUpperLimit()) / 2.;
         double deltaa = (a->GetUpperLimit() - a->GetLowerLimit());
         if (deltaa <= 1e-7 * meana)
            continue;

         double meanb = (b->GetLowerLimit() + b->GetUpperLimit()) / 2.;
         double deltab = (b->GetUpperLimit() - b->GetLowerLimit());
         if (deltab <= 1e-7 * meanb)
            continue;

         if (GetMarginalized(a, b))
            GetMarginalized(a, b)->Print(
                  Form("%s_2D_%s_%s.ps",filebase, a->GetName().data(), b->GetName().data()));
         k++;
      }
   }

   return k;
}

// ---------------------------------------------------------
int BCModel::PrintAllMarginalized(const char * file, unsigned int hdiv, unsigned int vdiv)
{
   if (!fMCMCFlagFillHistograms) {
      BCLog::OutError("BCModel::PrintAllMarginalized : Marginalized distributions not filled.");
      return 0;
   }

   int npar = GetNParameters();

   if (fMCMCH1Marginalized.size() == 0 || (fMCMCH2Marginalized.size() == 0
         && npar > 1)) {
      BCLog::OutError(
            "BCModel::PrintAllMarginalized : Marginalized distributions not available.");
      return 0;
   }

   // if there's only one parameter, we just want to call Print()
   if (fMCMCH1Marginalized.size() == 1 && fMCMCH2Marginalized.size() == 0) {
      BCParameter * a = GetParameter(0);
      if (GetMarginalized(a))
         GetMarginalized(a)->Print(file);
      return 1;
   }

   int c_width = 600; // default canvas width
   int c_height = 600; // default canvas height

   int type = 112; // landscape

   if (hdiv > vdiv) {
      if (hdiv > 3) {
         c_width = 1000;
         c_height = 700;
      }
      else {
         c_width = 800;
         c_height = 600;
      }
   }
   else if (hdiv < vdiv) {
      if (hdiv > 3) {
         c_height = 1000;
         c_width = 700;
      }
      else {
         c_height = 800;
         c_width = 600;
      }
      type = 111;
   }

   // calculate number of plots
   int nplots2d = npar * (npar - 1) / 2;
   int nplots = npar + nplots2d;

   // give out warning if too many plots
   BCLog::OutSummary(
         Form(
               "Printing all marginalized distributions (%d x 1D + %d x 2D = %d) into file %s",
               npar, nplots2d, nplots, file));
   if (nplots > 100)
      BCLog::OutDetail("This can take a while...");

   // setup the canvas and postscript file
   TCanvas * c = new TCanvas("c", "canvas", c_width, c_height);

   TPostScript * ps = new TPostScript(file, type);

   if (type == 112)
      ps->Range(24, 16);
   else
      ps->Range(16, 24);

   // draw all 1D distributions
   ps->NewPage();
   c->cd();
   c->Clear();
   c->Divide(hdiv, vdiv);

   int n = 0;
   for (int i = 0; i < npar; i++) {
      // get corresponding parameter
      BCParameter * a = GetParameter(i);

      // check if histogram exists
      if (!GetMarginalized(a))
         continue;

      // check if histogram is filled
      if (GetMarginalized(a)->GetHistogram()->Integral() <= 0)
         continue;

      // if current page is full, switch to new page
      if (i != 0 && i % (hdiv * vdiv) == 0) {
         c->Update();
         ps->NewPage();
         c->cd();
         c->Clear();
         c->Divide(hdiv, vdiv);
      }

      // go to next pad
      c->cd(i % (hdiv * vdiv) + 1);

      // just draw a line for a delta prior
      if (a->GetRangeWidth() == 0)
          GetMarginalized(a)->Draw(4, a->GetLowerLimit());
      else
          GetMarginalized(a)->Draw();

      n++;

      if (n % 100 == 0)
         BCLog::OutDetail(Form(" --> %d plots done", n));
   }

   if (n > 0) {
      c->Update();
   }

   // check how many 2D plots are actually drawn, despite no histogram filling or delta prior
   int k = 0;
   for (int i = 0; i < npar - 1; i++) {
      for (int j = i + 1; j < npar; j++) {
         // get corresponding parameters
         BCParameter * a = GetParameter(i);
         BCParameter * b = GetParameter(j);

         // check if histogram exists, or skip if one par has a delta prior
         if (!GetMarginalized(a, b))
            continue;

         // check if histogram is filled
         if (GetMarginalized(a, b)->GetHistogram()->Integral() <= 0)
            continue;

         // if current page is full, switch to new page, but only if there is data to plot
         if ((k != 0 && k % (hdiv * vdiv) == 0) || k == 0) {
            c->Update();
            ps->NewPage();
            c->cd();
            c->Clear();
            c->Divide(hdiv, vdiv);
         }

         // go to next pad
         c->cd(k % (hdiv * vdiv) + 1);

         double meana = (a->GetLowerLimit() + a->GetUpperLimit()) / 2.;
         double deltaa = (a->GetUpperLimit() - a->GetLowerLimit());
         if (deltaa <= 1e-7 * meana)
            continue;

         double meanb = (b->GetLowerLimit() + b->GetUpperLimit()) / 2.;
         double deltab = (b->GetUpperLimit() - b->GetLowerLimit());
         if (deltab <= 1e-7 * meanb)
            continue;

         GetMarginalized(a, b)->Draw(52);
         k++;

         if ((n + k) % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n + k));
      }
   }

   if ((n + k) > 100 && (n + k) % 100 != 0)
      BCLog::OutDetail(Form(" --> %d plots done", n + k));

   c->Update();
   ps->Close();

   delete c;
   delete ps;

   // return total number of drawn histograms
   return n + k;
}

// ---------------------------------------------------------
BCH2D * BCModel::GetMarginalized(BCParameter * par1, BCParameter * par2)
{
   int index1 = par1->GetIndex();
   int index2 = par2->GetIndex();

   if (fMCMCFlagsFillHistograms[index1] == false || fMCMCFlagsFillHistograms[index2] == false) {
// don't print any error message, should be done upstream
//      BCLog::OutError(
//            Form("BCModel::GetMarginalized : Distribuion for '%s' and/or '%s' not filled.",
//                  par1->GetName().data(), par2->GetName().data()));
      return 0;
   }

   if (index1 == index2) {
// don't print any error message, should be done upstream
//      BCLog::OutError("BCModel::GetMarginalized : Provided parameters are identical. Distribution not available.");
      return 0;
   }

   if (par1->GetRangeWidth() == 0 || par2->GetRangeWidth() == 0){
       return 0;
   }

   BCParameter * npar1 = par1;
   BCParameter * npar2 = par2;

   if (index1 > index2) {
      npar1 = par2;
      npar2 = par1;

      int itmp = index1;
      index1 = index2;
      index2 = itmp;
   }

   // get histogram
   TH2D * hist = MCMCGetH2Marginalized(index1, index2);

   if (hist == 0)
      return 0;

   BCH2D * hprob = new BCH2D();

   // set axis labels
   hist->SetName(Form("hist_%s_%s_%s", GetName().data(), npar1->GetName().data(), npar2->GetName().data()));
   hist->SetXTitle(Form("%s", npar1->GetName().data()));
   hist->SetYTitle(Form("%s", npar2->GetName().data()));
   hist->SetStats(kFALSE);

   double gmode[] = { fBestFitParameters.at(index1), fBestFitParameters.at(index2) };
   hprob->SetGlobalMode(gmode);

   // set histogram
   hprob->SetHistogram(hist);

   return hprob;
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
double BCModel::GetPvalueFromChi2Johnson(std::vector<double> par)
{
   double chi2 = GetChi2Johnson(par);
   // look up corresponding p value
   fPValue = TMath::Prob(chi2, NumberBins() - 1);
   return fPValue;
}

// ---------------------------------------------------------
double BCModel::GetChi2Johnson(std::vector<double> par, int nBins)
{
   typedef unsigned int uint;

   // number of observations
   int n = GetNDataPoints();

   if (nBins < 0)
      nBins = NumberBins();

   // fixed width quantiles, including final point!
   std::vector<double> a;
   for (int i = 0; i <= nBins; i++)
      a.push_back(i / double(nBins));

   // determine the bin counts and fill the histogram with data using the CDF
   TH1I * hist = new TH1I("h1", "h1 title", nBins, 0., 1.);

   // discrete case requires randomization to allocate counts of bins that cover more
   // than one quantile
   if (flag_discrete) {
      // loop over observations, each may have different likelihood and CDF
      for (int j = 0; j < n; j++) {
         // actual value
            double CDFval = CDF(par, j, false);
         // for the bin just before
         double CDFlower = CDF(par, j, true);

         // what quantiles q are covered, count from q_1 to q_{nBins}
         int qMax = 1;
         for (int i = 0; i < nBins; i++) {
            if (CDFval > a[i])
               qMax = i + 1;
            else
               break;
         }
         int qMin = 1;
         for (int i = 0; i < nBins; i++) {
            if (CDFlower > a[i])
               qMin = i + 1;
            else
               break;
         }

         // simplest case: observation bin entirely contained in one quantile
         if (qMin == qMax) {
            // watch out for overflow because CDF exactly 1
            if (CDFval < 1)
               hist->Fill(CDFval);
            else
               hist->AddBinContent(qMax);

            continue; // this observation finished
         }

         // if more than quantile is covered need more work:
         // determine probabilities of this observation to go for each quantile covered
         // as follows: If each quantile has size 0.25 and the CDF(integral of likelihood)
         // for current observation gives gives 0.27, but for observation-1 we would have
         // 0.20, then 5/7 of the 7% go for first quantile and 2/7 for the second.
         // This extend to bins covering more than two quantiles
         std::vector<double> prob;
         // normalization
         double norm = 1 / double(CDFval - CDFlower);

         for (int i = 0; i < (qMax - qMin + 1); i++) {
            if (i == 0) {
               prob.push_back(norm * (a[qMin] - CDFlower));
               continue;
            }
            if (i == (qMax - qMin)) {
               prob.push_back(norm * (CDFval - a[qMax - 1]));
               continue;
            }
            // default case
            prob.push_back(norm * (a[i] - a[i - 1]));
         }
         // have distribution, use inverse-transform method
         double U = fRandom->Rndm();
         // build up the integral (CDF)
         for (uint i = 1; i < prob.size(); i++)
            prob[i] += prob[i - 1];
         // and search with linear comput. complexity
         for (uint i = 0; i < prob.size(); i++) {
            // we finally allocate the count, as center of quantile
            if (U <= prob[i]) {
               hist->Fill((a[qMin + i - 1] + a[qMin + i]) / 2.);
               break;
            }
         }
      }
   }
   else { //continuous case is simple
      for (int j = 0; j < n; j++)
            hist->Fill(CDF(par, j, false));
   }

   // calculate chi^2
   double chi2 = 0.;
   double mk, pk;
   double N = double(n);
   for (int i = 1; i <= nBins; i++) {
      mk = hist->GetBinContent(i);
      pk = a[i] - a[i - 1];
      chi2 += (mk - N * pk) * (mk - N * pk) / (N * pk);
   }

   delete hist;

   return chi2;
}

// ---------------------------------------------------------
double BCModel::GetAvalueFromChi2Johnson(TTree * tree, TH1D * distribution)
{
   // model parameters
   int nPar = (int)GetNParameters();
   std::vector<double> param(nPar);

   //parameters saved in branches should be the same
   int nParBranches = -1;
   tree->SetBranchAddress("fNParameters", &nParBranches);

   //assume all events have same number of parameters, so check only first
   tree->GetEntry(0);
   if (nParBranches != nPar) {
      BCLog::OutError(Form(
            "Cannot compute A value: number of parameters in tree (%d)"
               "doesn't match  number of parameters in model (%d)",
            nParBranches, nPar));
      return -1.;
   }

   // buffer to construct correct branchnames for parameters, e.g. "fParameter3"
   char * branchName = new char[10 + nPar];
   // set up variables filled for each sample of parameters
   // assume same order as in model
   for (int i = 0; i < (int) nPar; i++) {
      sprintf(branchName, "fParameter%d", i);
      tree->SetBranchAddress(branchName, &param[i]);
   }

   // get the p value from Johson's definition for each param from posterior
   long nEntries = tree->GetEntries();

   // RN ~ chi2 with K-1 DoF needed for comparison
   std::vector<double> randoms(nEntries);
   int K = NumberBins();
   BCMath::RandomChi2(randoms, K - 1);

   // number of Johnson chi2 values bigger than reference
   int nBigger = 0;
   for (int i = 0; i < nEntries; i++) {
      tree->GetEntry(i);
      double chi2 = GetChi2Johnson(param);
      if (distribution != 0)
         distribution->Fill(chi2);
      // compare to set of chi2 variables
      if (chi2 > randoms[i])
         nBigger++;
   }

   return nBigger / double(nEntries);
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

//      BCLog::OutDebug(Form("BCModel::GetPvalueFromKolmogorov : "
//         "expected vs empirical (%f vs %f)", *iter, ECDF->GetBinContent(iBin
//            + 1)));

      // advance to next entry in the set
      ++iter;
   }

   // correct for total #points, not unique #points.
   // would need sqrt(n1*n2/(n1+n2)) if we had two experimental datasets
   double z = distMax * sqrt(N);

   fPValue = TMath::KolmogorovProb(z);

//   BCLog::OutDebug(Form("BCModel::GetPvalueFromKolmogorov : "
//      "max distance vs corrected (%f vs %f)", distMax, z));

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
double BCModel::HessianMatrixElement(BCParameter * par1, BCParameter * par2, std::vector<double> point)
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

   int idx1 = par1->GetIndex();
   int idx2 = par2->GetIndex();

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
void BCModel::FixDataAxis(unsigned int index, bool fixed)
{
   // check if index is within range
   if (index > fDataSet->GetDataPoint(0)->GetNValues()) {
      BCLog::OutError("BCModel::FixDataAxis : Index out of range.");
      return;
   }

   if (fDataFixedValues.size() == 0)
      fDataFixedValues.assign(fDataSet->GetDataPoint(0)->GetNValues(),
            false);

   fDataFixedValues[index] = fixed;
}

// ---------------------------------------------------------
bool BCModel::GetFixedDataAxis(unsigned int index)
{
   // check if index is within range
   if (index > fDataSet->GetDataPoint(0)->GetNValues()) {
      BCLog::OutError("BCModel::GetFixedDataAxis : Index out of range.");
      return false;
   }

   return fDataFixedValues[index];
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

   RecalculatePriorConstant();

   // reset all results
   ResetResults();

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
   SetParameterRange(index, value, value);

   // set prior
   return SetPriorConstant(index);
}

// ---------------------------------------------------------
int BCModel::SetPriorDelta(const char* name, double value)
{
  // find index
   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++)
      if (name == GetParameter(i)->GetName())
        index = i;

   // set prior
    return SetPriorDelta(index, value);
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
      fPriorContainer[index] = new TH1(*h);

      if (interpolate)
         fPriorContainerInterpolate[index] = true;

      fPriorContainerConstant[index] = false;
   }

   RecalculatePriorConstant();

   // reset all results
   ResetResults();

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

   RecalculatePriorConstant();

   // reset all results
   ResetResults();

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCModel::SetPriorConstant(const char* name)
{
   // find index
   int index = -1;
   for (unsigned int i = 0; i < GetNParameters(); i++) {
      if (name == GetParameter(i)->GetName()) {
         index = i;
         break;
      }
   }

   if (index == -1) {
      BCLog::OutError(Form(
            "BCModel::SetPriorConstant : parameter '%s' doesn't exist.", name));
      return 0;
   }

   return SetPriorConstant(index);
}

// ---------------------------------------------------------
int BCModel::SetPriorConstantAll()
{
   // get number of parameters
   int nPar = GetNParameters();

   if (nPar == 0)
      BCLog::OutWarning("BCModel::SetPriorConstantAll : No parameters defined.");

   // loop over all 1-d priors
   for (int i = 0; i < nPar; ++i) {
      if (fPriorContainer[i]) {
         delete fPriorContainer[i];
         fPriorContainer[i]=0;
      }
      fPriorContainerConstant[i] = true;
   }

   RecalculatePriorConstant();

   // reset all results
   ResetResults();

   // no error
   return 1;
}

// ---------------------------------------------------------
void BCModel::RecalculatePriorConstant()
{
   fPriorConstantValue = 0.;

   // get number of parameters
   int npar = GetNParameters();

   int nconstant = 0;

   for (int i=0; i<npar; ++i)
      if (fPriorContainerConstant[i]) {
         // default case
         if (GetParameter(i)->GetRangeWidth() > 0)
             fPriorConstantValue -= log(GetParameter(i)->GetRangeWidth());
         // do not add infinity due to zero width for delta prior
         ++nconstant;
      }

   if (nconstant == npar)
      fPriorConstantAll = true;
   else
      fPriorConstantAll = false;
}

// ---------------------------------------------------------
int BCModel::SetParameterRange(int index, double parmin, double parmax)
{
   // check index
   if (index < 0 || index >= int(GetNParameters())) {
      BCLog::OutError("BCModel::SetParameterRange : Index out of range.");
      return 0;
   }

   // set parameter ranges in BAT
   GetParameter(index)->SetLowerLimit(parmin);
   GetParameter(index)->SetUpperLimit(parmax);
   fMCMCBoundaryMin[index] = parmin;
   fMCMCBoundaryMax[index] = parmax;

   // reset results
   ResetResults();

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCModel::ResetResults()
{
   BCIntegrate::IntegrateResetResults();

   BCEngineMCMC::MCMCResetResults();

   // no error
   return 1;
}

// ---------------------------------------------------------
void BCModel::PrintSummary()
{
   int nparameters = GetNParameters();

   // model summary
	 BCLog::OutSummary(Form("Model : %s", fName.data()));
	 BCLog::OutSummary(Form("Index                : %d", fIndex));
	 BCLog::OutSummary(Form("Number of parameters : %d", nparameters));
	 BCLog::OutSummary("Parameters:");

   // parameter summary
   for (int i = 0; i < nparameters; i++)
      fParameterSet->at(i)->PrintSummary();

   // best fit parameters
   if (GetBestFitParameters().size() > 0) {
		 BCLog::OutSummary("Best fit parameters:");
		 
		 for (int i = 0; i < nparameters; i++) {
			 BCLog::OutSummary(Form(" %s = %f (global)", fParameterSet->at(i)->GetName().data(), GetBestFitParameter(i)));
			 
			 if ((int) fBestFitParametersMarginalized.size() == nparameters)
				 BCLog::OutSummary(Form(" %s = %f (marginalized)", fParameterSet->at(i)->GetName().data(), GetBestFitParameterMarginalized(i)));
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
   if (fNormalization > 0) {
		 BCLog::OutSummary(" Normalization:");
		 BCLog::OutSummary(Form(" - normalization : %f", fNormalization));
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
   int npar = MCMCGetNParameters();
   int nchains = MCMCGetNChains();

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
   for (int i = 0; i < npar; ++i)
      ofi << "  (" << i << ") Parameter \""
          << fParameterSet->at(i)->GetName().data() << "\"" << ": "
          << "(" << fParameterSet->at(i)->GetLowerLimit() << ", "
          << fParameterSet->at(i)->GetUpperLimit() << ")" << std::endl;
   ofi << std::endl;

   ofi << " Results of the optimization" << std::endl
       << " ===========================" << std::endl
       << " Optimization algorithm used: "
       << DumpUsedOptimizationMethod()<< std::endl;

   if (int(fBestFitParameters.size()) > 0) {
      ofi << " List of parameters and global mode:" << std::endl;
      for (int i = 0; i < npar; ++i) {
         ofi << "  (" << i << ") Parameter \""
             << fParameterSet->at(i)->GetName().data() << "\": "
             << fBestFitParameters[i];
             if (int(fBestFitParameterErrors.size()) == npar)
                if(fBestFitParameterErrors[i]>=0.)
                   ofi << " +- " << fBestFitParameterErrors[i];
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
          << " p-value at global mode: " << fPValue << std::endl << std::endl;
   }

   if (fNormalization >= 0.) {
       ofi << " Results of the normalization" << std::endl
             << " ============================" << std::endl
             << " Integration method used:"
             << DumpIntegrationMethod() << std::endl;
       ofi << " Normalization factor: " << fNormalization << std::endl << std::endl;
   }

   // give warning if MCMC did not converge
   if (!flag_conv && fMCMCFlagRun)
      ofi << " WARNING: the Markov Chain did not converge!" << std::endl
          << " Be cautious using the following results!" << std::endl
          << std::endl;

   // print results of marginalization (if MCMC was run)
   if (fMCMCFlagRun && fMCMCFlagFillHistograms) {
      ofi << " Results of the marginalization" << std::endl
          << " ==============================" << std::endl
          << " List of parameters and properties of the marginalized"
          << std::endl << " distributions:" << std::endl;
      for (int i = 0; i < npar; ++i) {
         if (!fMCMCFlagsFillHistograms[i])
            continue;

         BCH1D * bch1d = GetMarginalized(fParameterSet->at(i));

         ofi << "  (" << i << ") Parameter \""
             << fParameterSet->at(i)->GetName().data() << "\":" << std::endl

             << "      Mean +- sqrt(V):                " << std::setprecision(4)
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
         int ninter = int(v.size());
         ofi << "      Smallest interval(s) containing at least 68% and local mode(s):"
             << std::endl;
         for (int j = 0; j < ninter; j += 5)
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
          << " Average efficiencies:" << std::endl;

      std::vector<double> efficiencies;
      efficiencies.assign(npar, 0.);

      for (int ipar = 0; ipar < npar; ++ipar)
         for (int ichain = 0; ichain < nchains; ++ichain) {
            int index = ichain * npar + ipar;
            efficiencies[ipar] +=
                  double(MCMCGetNTrialsTrue().at(index)) / double(MCMCGetNTrialsTrue().at(index)
                  + MCMCGetNTrialsFalse().at(index)) / double(nchains) * 100.;
         }

      for (int ipar = 0; ipar < npar; ++ipar)
         ofi << "  (" << ipar << ") Parameter \""
             << fParameterSet->at(ipar)->GetName().data() << "\": "
             << efficiencies.at(ipar) << "%" << std::endl;
      ofi << std::endl;
   }

    ofi << " -----------------------------------------------------" << std::endl
          << " Notation:" << std::endl
          << " Mean        : mean value of the marg. pdf" << std::endl
          << " Median      : maximum of the marg. pdf" << std::endl
          << " Marg. mode  : most probable value of the marg. pdf" << std::endl
          << " V           : Variance of the marg. pdf" << std::endl
          << " Quantiles   : most commonly used quantiles" <<std::endl
          << " -----------------------------------------------------" << std::endl
          << std::endl;

   // close file
   //   ofi.close;
}

// ---------------------------------------------------------
void BCModel::PrintShortFitSummary(int chi2flag)
{
   BCLog::OutSummary("---------------------------------------------------");
   BCLog::OutSummary(Form("Fit summary for model \'%s\':", GetName().data()));
   BCLog::OutSummary(Form("   Number of parameters:  Npar  = %i", GetNParameters()));
   BCLog::OutSummary(Form("   Number of data points: Ndata = %i", GetNDataPoints()));
   BCLog::OutSummary("   Number of degrees of freedom:");
   BCLog::OutSummary(Form("      NDoF = Ndata - Npar = %i", GetNDataPoints() - GetNParameters()));

   BCLog::OutSummary("   Best fit parameters (global):");
   for (unsigned int i = 0; i < GetNParameters(); ++i)
      BCLog::OutSummary(Form("      %s = %.3g", GetParameter(i)->GetName().data(), GetBestFitParameter(i)));

   BCLog::OutSummary("   Goodness-of-fit test:");
   BCLog::OutSummary(Form("      p-value = %.3g", GetPValue()));
   if (chi2flag) {
      BCLog::OutSummary(Form("      p-value corrected for NDoF = %.3g", GetPValueNDoF()));
      BCLog::OutSummary(Form("      chi2 / NDoF = %.3g", GetChi2NDoF()));
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
               fParameterSet->at(i), fParameterSet->at(j), parameters);

         // print to screen
				 BCLog::OutSummary(Form("%d %d : %f", i, j, hessianmatrixelement));
      }
}

// ---------------------------------------------------------
BCDataPoint * BCModel::VectorToDataPoint(const std::vector<double> &data)
{
   int sizeofvector = int(data.size());
   BCDataPoint * datapoint = new BCDataPoint(sizeofvector);
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

void BCModel::StoreMode()
{
   //start with max. of posterior at first chain
   double probmax = fMCMCprobMax.at(0);
   int probmaxindex = 0;

   // loop over all chains and find the maximum point
   for (int i = 1; i < fMCMCNChains; ++i) {
      if (fMCMCprobMax.at(i) > probmax) {
         probmax = fMCMCprobMax.at(i);
         probmaxindex = i;
      }
   }

   // save if improved the log posterior
   if (fBestFitParameters.empty() || probmax > LogEval(fBestFitParameters)) {
      fBestFitParameters.assign(fMCMCNParameters, 0.0);
      for (int i = 0; i < fMCMCNParameters; ++i)
         fBestFitParameters[i] = fMCMCxMax[probmaxindex * fMCMCNParameters + i];

      SetOptimizationMethodMode(BCIntegrate::kOptMetropolis);
   }
}

// ---------------------------------------------------------

