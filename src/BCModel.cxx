/*
 * Copyright (C) 2008-2012, Daniel Kollar, Kevin Kroeninger, and Daniel Greenwald.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCModel.h"

#include "BCDataPoint.h"
#include "BCDataSet.h"
#include "BCErrorCodes.h" // todo remove
#include "BCGoFTest.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"
#include "BCMath.h"
#include "BCModelOutput.h"
#include "BCParameter.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TMath.h>
#include <TKey.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TTree.h>

#include <fstream>
#include <iomanip>
#include <set>

// ---------------------------------------------------------
BCModel::BCModel(const char * name)
   : BCIntegrate()
{
   fNormalization = -1.;
   fDataSet = 0;
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
   if (bcmodel.fNDataPointsMinimum)
      fNDataPointsMinimum = bcmodel.fNDataPointsMinimum;
   else
      fNDataPointsMinimum = 0;
   if (bcmodel.fNDataPointsMaximum)
      fNDataPointsMaximum = bcmodel.fNDataPointsMaximum;
   else
      fNDataPointsMaximum = 0;

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
const std::string & BCModel::Get1DDefaultPlotOptions()
{
   static const std::string opt = "BTciB1CS1D0pdf0Lmeanmode";
   return opt;
}

// ---------------------------------------------------------
const std::string & BCModel::Get2DDefaultPlotOptions()
{
   static const std::string opt = "BTfB1CS1meangmodelmode";
   return opt;
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
std::vector<double> BCModel::GetErrorBand(double level) const
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
TGraph * BCModel::GetErrorBandGraph(double level1, double level2) const
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
TH2D * BCModel::GetErrorBandXY_yellow(double level, int nsmooth) const
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
   // todo still needed?
   // reset results
   ResetResults();

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
   if (fNormalization <= 0.) {
      BCLog::OutError("BCModel::LogProbability. Normalization not available or zero.");
      return 0.;
   }

   return LogProbabilityNN(parameters) - log(fNormalization);
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
                  "is defined but not recognized.",
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
double BCModel::Normalize()
{
   if (fParameters.Size() < 1) {
      BCLog::OutError(Form("Normalize : No parameters defined in model \'%s\'. Aborting.",GetName().data()));
      return -1.;
   }

   BCLog::OutSummary(Form("Model \'%s\': Normalizing probability",GetName().data()));

   // integrate and get best fit parameters
   // maybe we have to remove the mode finding from here in the future
   fNormalization = Integrate();

   BCLog::OutDetail(Form(" --> Normalization factor : %.6g", fNormalization));

   return fNormalization;
}

// ---------------------------------------------------------
int BCModel::MarginalizeAll()
{
   if (fParameters.Size() < 1) {
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

      for (unsigned ix = 1; ix <= fErrorBandNbinsX; ++ix)
         for (unsigned iy = 1; iy <= fErrorBandNbinsX; ++iy)
            fErrorBandXY->SetBinContent(ix, iy, 0.);
   }

   // run the Markov chains
   MCMCMetropolis();

   return 1;
}

// ---------------------------------------------------------
BCH1D * BCModel::GetMarginalized(const BCParameter * parameter)
{
   if ( !parameter)
      return 0;
   return GetMarginalized(fParameters.Index(parameter->GetName()));
}

// ---------------------------------------------------------
BCH1D * BCModel::GetMarginalized(unsigned index)
{
   // get histogram
   BCH1D * hist = MCMCGetH1Marginalized(index);
   if (!hist)
      return 0;

   // set axis labels
   hist->GetHistogram()->SetName(Form("hist_%s_%s", GetName().data(), fParameters[index]->GetName().data()));
   hist->GetHistogram()->SetYTitle(Form("p(%s|data)", fParameters[index]->GetLatexName().data()));

   return hist;
}

// ---------------------------------------------------------
BCH1D * BCModel::GetSlice(const BCParameter* parameter, const std::vector<double> parameters, int nbins)
{
	// check if parameter exists
	if (!parameter) {
		BCLog::OutError("BCModel::GetSlice : Parameter does not exist.");
		return 0;
	}

	// create local copy of parameter set
	std::vector<double> parameters_temp;
	parameters_temp = parameters;

	// normalization flag: if true, normalize slice histogram to unity
	bool flag_norm = false;

	// check if parameter set if defined
	if (parameters_temp.size()==0 && GetNParameters()==1) {
		parameters_temp.push_back(0);
		flag_norm = true; // slice is the 1D pdf, so normalize it to unity
	}
	else if (parameters_temp.size()==0 && GetNParameters()!=1) {
		BCLog::OutError("BCModel::GetSlice : No parameters defined.");
		return 0;
	}

	// calculate number of bins
	if (nbins <= 0)
		nbins = parameter->GetNbins();

	// create histogram
	TH1D * hist = new TH1D("", "", nbins, parameter->GetLowerLimit(), parameter->GetUpperLimit());

	// set axis labels
	hist->SetName(Form("hist_%s_%s", GetName().data(), parameter->GetName().data()));
	hist->SetXTitle(parameter->GetLatexName().data());
	if (GetNParameters() == 1)
		hist->SetYTitle(Form("p(%s|data)", parameter->GetLatexName().data()));
	else
		hist->SetYTitle(Form("p(%s|data, all other parameters fixed)", parameter->GetLatexName().data()));
	hist->SetStats(kFALSE);

	// fill histogram
	for (int i = 1; i <= nbins; ++i) {
		double par_temp = hist->GetBinCenter(i);
		parameters_temp[fParameters.Index(parameter->GetName())] = par_temp;
		double prob = Eval(parameters_temp);
		hist->SetBinContent(i, prob);
	}

	// normalize
	if (flag_norm)
		hist->Scale(1.0/hist->Integral());

	// set histogram
	BCH1D * hprob = new BCH1D();
	hprob->SetHistogram(hist);

	return hprob;
}
// ---------------------------------------------------------
BCH2D* BCModel::GetSlice(const BCParameter* parameter1, const BCParameter* parameter2, const std::vector<double> parameters, int nbins)
{
   return GetSlice(parameter1->GetName().c_str(), parameter2->GetName().c_str(), parameters, nbins);
}

// ---------------------------------------------------------
BCH2D* BCModel::GetSlice(const char* name1, const char* name2, const std::vector<double> parameters, int nbins)
{
   return GetSlice(fParameters.Index(name1), fParameters.Index(name2), parameters, nbins);
}

// ---------------------------------------------------------
BCH2D* BCModel::GetSlice(unsigned index1, unsigned index2, const std::vector<double> parameters, int nbins)
{
	// check if parameter exists
	if (!fParameters.ValidIndex(index1) || !fParameters.ValidIndex(index2)) {
		BCLog::OutError("BCModel::GetSlice : Parameter does not exist.");
		return 0;
	}

	// create local copy of parameter set
	std::vector<double> parameters_temp;
	parameters_temp = parameters;

	// normalization flag: if true, normalize slice histogram to unity
	bool flag_norm = false;

	// check number of dimensions
	if (GetNParameters() < 2) {
		BCLog::OutError("BCModel::GetSlice : Number of parameters need to be at least 2.");
	}

	// check if parameter set if defined
	if (parameters_temp.size()==0 && GetNParameters()==2) {
		parameters_temp.push_back(0);
		parameters_temp.push_back(0);
		flag_norm = true; // slice is the 1D pdf, so normalize it to unity
	}
	else if (parameters_temp.size()==0 && GetNParameters()>2) {
		BCLog::OutError("BCModel::GetSlice : No parameters defined.");
		return 0;
	}

	// calculate number of bins
	const BCParameter * p1 = fParameters.Get(index1);
	const BCParameter * p2 = fParameters.Get(index2);
	unsigned nbins1, nbins2;
	if (nbins <= 0) {
	   nbins1 = p1->GetNbins();
	   nbins2 = p2->GetNbins();
	} else {
	   nbins1 = nbins2 = nbins;
	}

	// create histogram
	TH2D * hist = new TH2D("", "", nbins1, p1->GetLowerLimit(), p1->GetUpperLimit(),
											 nbins2, p2->GetLowerLimit(), p2->GetUpperLimit());

	// set axis labels
	hist->SetName(Form("hist_%s_%s_%s", GetName().c_str(), p1->GetName().c_str(), p2->GetName().c_str()));
	hist->SetXTitle(Form("%s", p1->GetLatexName().data()));
	hist->SetYTitle(Form("%s", p2->GetLatexName().data()));
	hist->SetStats(kFALSE);

	// fill histogram
	for (int ix = 1; ix <= nbins; ++ix) {
		for (int iy = 1; iy <= nbins; ++iy) {
		// debugKK: I am here
			double par_temp1 = hist->GetXaxis()->GetBinCenter(ix);
			double par_temp2 = hist->GetYaxis()->GetBinCenter(iy);

		parameters_temp[index1] = par_temp1;
		parameters_temp[index2] = par_temp2;

		double prob = Eval(parameters_temp);
		hist->SetBinContent(ix, iy, prob);
		}
	}

	// normalize
	if (flag_norm)
		hist->Scale(1.0/hist->Integral());

	// set histogram
	BCH2D * hprob = new BCH2D();
	hprob->SetHistogram(hist);

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
      const BCParameter * a = GetParameter(i);
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
         const BCParameter * a = GetParameter(i);
         const BCParameter * b = GetParameter(j);
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
      const BCParameter * a = GetParameter(i);
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
         const BCParameter * a = GetParameter(i);
         const BCParameter * b = GetParameter(j);

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
int BCModel::PrintAllMarginalized(const char * file, std::string options1d, std::string options2d, unsigned int hdiv, unsigned int vdiv)
{

   if (fMCMCH1Marginalized.size() == 0 and fMCMCH2Marginalized.size() == 0) {
      BCLog::OutError("BCModel::PrintAllMarginalized : Marginalized distributions not available.");
      return 0;
   }

   // find all valid (non NULL) histograms
   std::vector<TH1 *> validH1;

   for (unsigned i = 0 ; i < fMCMCH1Marginalized.size() ; ++i)
   {
      if (TH1D * h = fMCMCH1Marginalized[i])
         validH1.push_back(h);
   }

   std::vector<TH2D *> validH2;
   for (unsigned i = 0 ; i < fMCMCH2Marginalized.size() ; ++i)
   {
      if (TH2D * h = fMCMCH2Marginalized[i])
         validH2.push_back(h);
   }

   std::string filename(file);

   // check if file extension does not exist or is not pdf or ps
   if ( (filename.find_last_of(".") == std::string::npos) or
         ((filename.substr(filename.find_last_of(".")+1) != "pdf") and
               (filename.substr(filename.find_last_of(".")+1) != "ps"))) {
      // make it a PDF file
      filename += ".pdf";
   }

   // todo do we really need this?
   // if there's only one parameter, we just want to call Print()
   if (fMCMCH1Marginalized.size() == 1 && fMCMCH2Marginalized.size() == 0) {
      if (BCH1D * m = GetMarginalized(0u)) {
         m->Print(filename.c_str());
         delete m;
      }
      return 1;
   }

   int c_width  = gStyle->GetCanvasDefW(); // default canvas width
   int c_height = gStyle->GetCanvasDefH(); // default canvas height

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
   }

   const unsigned nplots = validH1.size() + validH2.size();

   // give out warning if too many plots
   BCLog::OutSummary(Form("Printing all marginalized distributions (%d x 1D + %d x 2D = %d) into file %s",
         validH1.size(), validH2.size(), nplots, filename.c_str()));
   if (nplots > 100)
      BCLog::OutDetail("This can take a while...");

   // setup the canvas and file
   TCanvas c("c", "canvas", c_width, c_height);
   c.Divide(hdiv, vdiv);

   // count plots
   unsigned n = 0;
   for (unsigned i = 0; i < fParameters.Size(); ++i) {
      BCH1D * h = GetMarginalized(i);

      // check if histogram exists
      if ( !h)
         continue;

      // if current page is full, switch to new page
      if (n != 0 && n % (hdiv * vdiv) == 0) {
         if ( n <= (hdiv * vdiv)) {
            c.Print(std::string( filename + "(").c_str());
         }
         else {
            c.Print(filename.c_str());
         }
      }

      // go to next pad
      c.cd(n % (hdiv * vdiv) + 1);

      h->Draw(options1d);
      delete h;
      {
      }

      if (++n % 100 == 0)
         BCLog::OutDetail(Form(" --> %d plots done", n));
   }

   // for clean up later
   std::vector<BCH2D *> h2;

   // check how many 2D plots are actually drawn, despite no histogram filling or delta prior
   unsigned k = 0;
   for (unsigned i = 0; i < fParameters.Size() - 1; ++i) {
      for (unsigned j = i + 1; j < fParameters.Size(); ++j) {
         // check if histogram exists, or skip if one par has a delta prior
         h2.push_back(GetMarginalized(i, j));
         if ( !h2.back())
            continue;

         // get corresponding parameters
         BCParameter * a = GetParameter(i);
         BCParameter * b = GetParameter(j);

         // if current page is full, switch to new page, but only if there is data to plot
         if ((k != 0 && k % (hdiv * vdiv) == 0) || k == 0) {
            c.Print(filename.c_str());
         }

         // go to next pad
         c.cd(k % (hdiv * vdiv) + 1);

         const double meana = (a->GetLowerLimit() + a->GetUpperLimit()) / 2.;
         if (a->GetRangeWidth() <= 1e-7 * meana)
            continue;

         const double meanb = (b->GetLowerLimit() + b->GetUpperLimit()) / 2.;
         if (b->GetRangeWidth() <= 1e-7 * meanb)
            continue;

         h2.back()->Draw(options2d);
         k++;

         if ((n + k) % 100 == 0)
            BCLog::OutDetail(Form(" --> %d plots done", n + k));
      }
   }

   if ((n + k) > 100 && (n + k) % 100 != 0)
      BCLog::OutDetail(Form(" --> %d plots done", n + k));

   c.Print(std::string( filename + ")").c_str());

   // clean up
   for (unsigned i = 0; i < h2.size() ; ++i)
      delete h2[i];

   // return total number of drawn histograms
   return n + k;

   return 0;
}

// ---------------------------------------------------------
BCH2D * BCModel::GetMarginalized(const BCParameter * par1, const BCParameter * par2)
{
   if ( !par1 or !par2 or (par1 == par2) )
      return 0;

  return GetMarginalized(fParameters.Index(par1->GetName()), fParameters.Index(par2->GetName()));
}

// ---------------------------------------------------------
BCH2D * BCModel::GetMarginalized(unsigned index1, unsigned index2)
{
   BCH2D * h = MCMCGetH2Marginalized(index1, index2);
   if ( !h)
      return 0;

   h->GetHistogram()->SetName(Form("hist_%s_%s_%s", GetName().data(),
                                    fParameters[index1]->GetName().data(),
                                    fParameters[index2]->GetName().data()));

   return h;
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
bool BCModel::GetFixedDataAxis(unsigned int index) const
{
   // check if index is within range
   if (index > fDataSet->GetDataPoint(0)->GetNValues()) {
      BCLog::OutError("BCModel::GetFixedDataAxis : Index out of range.");
      return false;
   }

   return fDataFixedValues[index];
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
// todo remove prior delta, or jus call Fix() for back ward compatability? How to deal with in MCMC?
// ---------------------------------------------------------
int BCModel::SetPriorDelta(int index, double value)
{
   // set range to value
//   SetParameterRange(index, value, value);
   GetParameter(index)->Fix(value);

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
      BCLog::OutSummary(Form("Log of the maximum posterior: %f", fLogMaximum));
      BCLog::OutSummary("Best fit parameters:");

      for (unsigned i = 0; i < GetNParameters(); i++) {
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
   for (unsigned i = 0; i < npar; ++i)
      ofi << "  (" << i << ") Parameter \""
          << fParameters[i]->GetName() << "\"" << ": "
          << "(" << fParameters[i]->GetLowerLimit() << ", "
          << fParameters[i]->GetUpperLimit() << ")" << std::endl;
   ofi << std::endl;

   ofi << " Results of the optimization" << std::endl
         << " ===========================" << std::endl
         << " Optimization algorithm used: "
         << DumpUsedOptimizationMethod()<< std::endl;

   if ( ! fBestFitParameters.empty()) {
      ofi << "Log of the maximum posterior: " << fLogMaximum << std::endl;
      ofi << " List of parameters and global mode:" << std::endl;
      for (unsigned i = 0; i < npar; ++i) {
         ofi << "  (" << i << ") Parameter \""
             << fParameters[i]->GetName() << "\": "
               << fBestFitParameters[i];
             if (fBestFitParameterErrors.size() == npar)
                if(fBestFitParameterErrors[i] >= 0.)
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
            << " p-value: " << fPValue << std::endl;
      if (fPValueNDoF >= 0)
         ofi << " p-value corrected for degrees of freedom: " << fPValueNDoF << std::endl;

      ofi << std::endl;
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
   if (fMCMCFlagRun) {
      ofi << " Results of the marginalization" << std::endl
            << " ==============================" << std::endl
            << " List of parameters and properties of the marginalized"
            << std::endl << " distributions:" << std::endl;
      for (unsigned i = 0; i < npar; ++i) {
         if ( ! fParameters[i]->FillHistograms())
            continue;

         BCH1D * bch1d = GetMarginalized(fParameters[i]);

         ofi << "  (" << i << ") Parameter \""
             << fParameters[i]->GetName() << "\":" << std::endl

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
         ofi << "      Smallest interval(s) containing at least 68% and local mode(s):"
               << std::endl;
         for (unsigned j = 0; j < v.size(); j += 5)
            ofi << "       (" << v[j] << ", " << v[j + 1]
                                                   << ") (local mode at " << v[j + 3] << " with rel. height "
                                                   << v[j + 2] << "; rel. area " << v[j + 4] << ")"
                                                   << std::endl;
         ofi << std::endl;
         delete bch1d;
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

      for (unsigned ipar = 0; ipar < npar; ++ipar)
         for (unsigned ichain = 0; ichain < nchains; ++ichain) {
            unsigned index = ichain * npar + ipar;
            efficiencies[ipar] +=
                  double(MCMCGetNTrialsTrue().at(index)) / double(MCMCGetNTrialsTrue().at(index)
                        + MCMCGetNTrialsFalse().at(index)) / double(nchains) * 100.;
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
