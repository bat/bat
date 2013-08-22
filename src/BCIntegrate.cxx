/*
 * Copyright (C) 2008-2013, Daniel Kollar, Kevin Kroeninger, Daniel Greenwald, and Frederik Beaujean.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------
#include <config.h>

#include "BCIntegrate.h"
#include "BCErrorCodes.h"
#include "BCLog.h"
#include "BCMath.h"
#include "BCParameter.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>

#ifdef HAVE_CUBA_H
#include <cuba.h>
#endif

#include <math.h>
#include <limits>

namespace
{
    /**
     * Hold an instance of BCIntegrate to emulate a global variable for use with Minuit
     */
    class BCIntegrateHolder
    {
    private:
       BCIntegrate * global_this;

    public:

        BCIntegrateHolder() :
           global_this(NULL)
        {
        }

       /**
        * Set and/or retrieve the static BCIntegrate object
        */
       static BCIntegrate * instance(BCIntegrate * obj = NULL)
       {
          static BCIntegrateHolder result;
          if (obj)
             result.global_this = obj;

          return result.global_this;
       }
    };
}

// ---------------------------------------------------------
BCIntegrate::BCIntegrate() : BCEngineMCMC(),
    fMinuit(0),
    fMinuitErrorFlag(0),
    fFlagIgnorePrevOptimization(false),
    fSAT0(100),
    fSATmin(0.1),
    fTreeSA(0),
    fFlagWriteSAToFile(false),
    fSANIterations(0),
    fSATemperature(0),
    fSALogProb(0),
#ifdef HAVE_CUBA_H
    fIntegrationMethod(BCIntegrate::kIntCuba),
#else
    fIntegrationMethod(BCIntegrate::kIntMonteCarlo),
#endif
    fMarginalizationMethod(fIntegrationMethod),
    fSASchedule(BCIntegrate::kSACauchy),
    fNIterationsMin(0),
    fNIterationsMax(1000000),
    fNIterationsPrecisionCheck(1000),
    fNIterationsOutput(0),
    fNIterations(0),
    fRelativePrecision(1e-2),
    fAbsolutePrecision(1e-6),
    fError(-999.),
    fCubaIntegrationMethod(BCIntegrate::kCubaVegas)
{
	fMinuitArglist[0] = 20000;
	fMinuitArglist[1] = 0.01;
}

// ---------------------------------------------------------
BCIntegrate::BCIntegrate(const BCIntegrate & other) : BCEngineMCMC(other)
{
	Copy(other);
}

// ---------------------------------------------------------
BCIntegrate & BCIntegrate::operator = (const BCIntegrate & other)
{
	Copy(other);
	return *this;
}

// ---------------------------------------------------------
void BCIntegrate::Copy(const BCIntegrate & other)
{
	BCEngineMCMC::Copy(other);

	fBestFitParameterErrors   = other.fBestFitParameterErrors;
	fFillErrorBand            = other.fFillErrorBand;
	fFitFunctionIndexX        = other.fFitFunctionIndexX;
	fFitFunctionIndexY        = other.fFitFunctionIndexY;
	fErrorBandX               = other.fErrorBandX;
	if (other.fErrorBandXY)
		fErrorBandXY = new TH2D(*(other.fErrorBandXY));
	else
		fErrorBandXY = 0;
	fErrorBandNbinsX          = other.fErrorBandNbinsX;
	fErrorBandNbinsY          = other.fErrorBandNbinsY;
	fMinuit                   = new TMinuit();
	fMinuitArglist[0]         = other.fMinuitArglist[0];
	fMinuitArglist[1]         = other.fMinuitArglist[1];
	fMinuitErrorFlag          = other.fMinuitErrorFlag;
	fFlagIgnorePrevOptimization = other.fFlagIgnorePrevOptimization;
	fSAT0                     = other.fSAT0;
	fSATmin                   = other.fSATmin;
	fTreeSA = 0;
	fFlagWriteSAToFile        = other.fFlagWriteSAToFile;
	fSANIterations            = other.fSANIterations;
	fSATemperature            = other.fSATemperature;
	fSALogProb                = other.fSALogProb;
	fSAx                      = other.fSAx;
	fIntegrationMethod        = other.fIntegrationMethod;
	fMarginalizationMethod    = other.fMarginalizationMethod;
	fSASchedule               = other.fSASchedule;
	fNIterationsMin           = other.fNIterationsMin;
	fNIterationsMax           = other.fNIterationsMax;
	fNIterationsPrecisionCheck= other.fNIterationsPrecisionCheck;
	fNIterationsOutput        = other.fNIterationsOutput;
	fNIterations              = other.fNIterations;
	fRelativePrecision        = other.fRelativePrecision;
	fAbsolutePrecision        = other.fAbsolutePrecision;
	fError                    = other.fError;
	fCubaIntegrationMethod    = other.fCubaIntegrationMethod;
	fCubaVegasOptions         = other.fCubaVegasOptions;
	fCubaSuaveOptions         = other.fCubaSuaveOptions;
	fCubaDivonneOptions       = other.fCubaDivonneOptions;
	fCubaCuhreOptions         = other.fCubaCuhreOptions;
}

// ---------------------------------------------------------
BCIntegrate::~BCIntegrate()
{
   delete fMinuit;
}

// ---------------------------------------------------------
double BCIntegrate::GetBestFitParameter(unsigned index) const
{
   if( ! fParameters.ValidIndex(index))
      return -1e+111;

   if(fBestFitParameters.size()==0) {
      BCLog::OutError("BCIntegrate::GetBestFitParameter : Mode finding not yet run, returning center of the range.");
      return (fParameters[index]->GetUpperLimit() + fParameters[index]->GetLowerLimit() ) / 2.;
   }

   return fBestFitParameters[index];
}

// ---------------------------------------------------------
double BCIntegrate::GetBestFitParameterError(unsigned index) const
{
   if( ! fParameters.ValidIndex(index)) {
      BCLog::OutError("BCIntegrate::GetBestFitParameterError : Parameter index out of range, returning -1.");
      return -1;
   }

   if(fBestFitParameterErrors.size()==0) {
      BCLog::OutError("BCIntegrate::GetBestFitParameterError : Mode finding not yet run, returning -1.");
      return -1.;
   }

   if(fBestFitParameterErrors[index]<0.)
      BCLog::OutWarning("BCIntegrate::GetBestFitParameterError : Parameter error not available, returning -1.");

   return fBestFitParameterErrors[index];
}

// ---------------------------------------------------------
double BCIntegrate::Eval(const std::vector<double> & /*x*/)
{
   BCLog::OutWarning( "BCIntegrate::Eval. No function. Function needs to be overloaded.");
   return 0;
}

// ---------------------------------------------------------
double BCIntegrate::LogEval(const std::vector<double> &x)
{
   // this method should better also be overloaded
   return log(Eval(x));
}

// ---------------------------------------------------------
double BCIntegrate::EvalSampling(const std::vector<double> &)
{
   BCLog::OutError("BCIntegrate::EvalSampling. The importance sampling function needs to be implemented by the user.");
   return 0;
}

// ---------------------------------------------------------
void BCIntegrate::ResetResults()
{
   BCEngineMCMC::ResetResults();
   fBestFitParameterErrors.clear();
}

// ---------------------------------------------------------
unsigned BCIntegrate::GetNIntegrationVariables() {
    unsigned n = 0;
    for(unsigned i = 0; i < fParameters.Size(); ++i)
        if ( ! fParameters[i]->Fixed())
            ++n;
    return n;
}

// ---------------------------------------------------------
double BCIntegrate::CalculateIntegrationVolume() {
   double integrationVolume = -1.;

   for(unsigned i = 0; i < fParameters.Size(); i++)
      if ( ! fParameters[i]->Fixed()) {
         if (integrationVolume<0)
            integrationVolume = 1;
         integrationVolume *= fParameters[i]->GetRangeWidth();
      }

   if (integrationVolume<0)
      integrationVolume = 0;

   return integrationVolume;
}

// ---------------------------------------------------------
double BCIntegrate::Integrate(BCIntegrationMethod type)
{
   switch(type)
   {

   // Monte Carlo Integration
   case BCIntegrate::kIntMonteCarlo:
   {
      std::vector<double> sums (2,0.0);
      sums.push_back(CalculateIntegrationVolume());
      return Integrate(kIntMonteCarlo,
            &BCIntegrate::GetRandomVectorInParameterSpace,
            &BCIntegrate::EvaluatorMC,
            &IntegralUpdaterMC,
            sums);
   }

   // Importance Sampling Integration
   case BCIntegrate::kIntImportance:
   {
      std::vector<double> sums (2,0.0);
      return Integrate(kIntImportance,
            &BCIntegrate::GetRandomVectorInParameterSpace,
            &BCIntegrate::EvaluatorImportance,	 // use same evaluator as for metropolis
            &IntegralUpdaterImportance,	 // use same updater as for metropolis
            sums);
   }
   case BCIntegrate::kIntCuba:
      return IntegrateCuba();
   default:
      BCLog::OutError(TString::Format("BCIntegrate::Integrate : Invalid integration method: %d", fIntegrationMethod));
      break;
   }

   return 0;
}

// ---------------------------------------------------------
void BCIntegrate::SetBestFitParameters(const std::vector<double> &x, const double &new_value, double &old_value) {
	if (new_value < old_value)
		return;
	old_value = new_value;
	SetBestFitParameters(x);
}

// ---------------------------------------------------------
unsigned BCIntegrate::IntegrationOutputFrequency() const
{
	if (fNIterationsOutput > 0)
		return fNIterationsOutput;
	const unsigned nwrite = fNIterationsMax / 10;
	if(nwrite < 10000)
		return 1000;
	if(nwrite < 100000)
		return 10000;
	if(nwrite < 1000000)
		return 100000;
	return 1000000;
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtStartOfIntegration(BCIntegrationMethod type, BCCubaMethod cubatype)
{
   const unsigned NVarNow = GetNIntegrationVariables();

   BCLog::LogLevel level = BCLog::summary;

   if(fParameters.Size() != NVarNow) {

      level = BCLog::detail;
      bool printed = false;

      if (type==kIntCuba)
      {
         BCLog::OutDetail(TString::Format("Running %s (%s) integration over %i dimensions out of %i.",
               DumpIntegrationMethod(type).c_str(),
               DumpCubaIntegrationMethod(cubatype).c_str(),
               NVarNow, fParameters.Size()));
         printed = true;
      }

      if (not printed)
         BCLog::OutDetail(TString::Format("Running %s integration over %i dimensions out of %i.",
               DumpIntegrationMethod(type).c_str(),
               NVarNow, fParameters.Size()));

      BCLog::OutDetail(" --> Fixed parameters:");
      for(unsigned i = 0; i < fParameters.Size(); i++)
         if(fParameters[i]->Fixed())
            BCLog::OutDetail(TString::Format("      %3i :  %g", i, fParameters[i]->GetFixedValue()));
   }
   else {
      bool printed = false;

      if (type==kIntCuba) {
         BCLog::OutDetail(TString::Format("Running %s (%s) integration over %i dimensions.",
               DumpIntegrationMethod(type).c_str(),
               DumpCubaIntegrationMethod(cubatype).c_str(),
               NVarNow));
         printed = true;
      }

      if (not printed)
         BCLog::OutDetail(TString::Format("Running %s integration over %i dimensions.",
               DumpIntegrationMethod(type).c_str(),
               NVarNow));
   }

   BCLog::Out(level, TString::Format(" --> Minimum number of iterations: %i", GetNIterationsMin()));
   BCLog::Out(level, TString::Format(" --> Maximum number of iterations: %i", GetNIterationsMax()));
   BCLog::Out(level, TString::Format(" --> Target relative precision:    %e", GetRelativePrecision()));
   BCLog::Out(level, TString::Format(" --> Target absolute precision:    %e", GetAbsolutePrecision()));
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtEndOfIntegration(double integral, double absprecision, double relprecision, int nIterations)
{
   BCLog::OutSummary(TString::Format(" --> Result of integration:        %e +- %e   in %i iterations.", integral, absprecision, nIterations));
   BCLog::OutSummary(TString::Format(" --> Obtained relative precision:  %e. ", relprecision));
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtIntegrationStatusUpdate(BCIntegrationMethod type, double integral, double absprecision, int nIterations)
{
	BCLog::OutDetail(TString::Format("%s. Iteration %i, integral: %e +- %e.", DumpIntegrationMethod(type).c_str(), nIterations, integral, absprecision));
}

// ---------------------------------------------------------
double BCIntegrate::Integrate(BCIntegrationMethod type, tRandomizer randomizer, tEvaluator evaluator, tIntegralUpdater updater,
      std::vector<double> &sums)
{
	LogOutputAtStartOfIntegration(type, NCubaMethods);

	// reset variables
	double pmax = 0.;
	double integral = 0.;
	double absprecision = 2.*fAbsolutePrecision;
	double relprecision = 2.*fRelativePrecision;

	std::vector<double> randx (fParameters.Size(), 0.);

	// how often to print out the info line to screen
	int nwrite = IntegrationOutputFrequency();

	// reset number of iterations
	fNIterations = 0;
	bool accepted;

	// iterate while number of iterations is lower than minimum number of iterations
	// or precision is not reached and the number of iterations is lower than maximum number of iterations
	while ((GetRelativePrecision() < relprecision and GetAbsolutePrecision() < absprecision and GetNIterations() < GetNIterationsMax())
	      or GetNIterations() < GetNIterationsMin())
		{

			// get random numbers
			(this->*randomizer)(randx);

			// evaluate function at sampled point
			// updating sums & checking for maximum probability

			SetBestFitParameters(randx, (this->*evaluator)(sums,randx,accepted), pmax);

			// increase number of iterations
			if (accepted)
			   fNIterations++;

			// update precisions
			if (fNIterations % fNIterationsPrecisionCheck == 0) {
				(*updater)(sums, fNIterations, integral, absprecision);
				relprecision = absprecision / integral;
			}

			// write status
			if (fNIterations % nwrite == 0) {
				double temp_integral;
				double temp_absprecision;
				(*updater)(sums, fNIterations, temp_integral, temp_absprecision);
				LogOutputAtIntegrationStatusUpdate(type, temp_integral, temp_absprecision, fNIterations);
			}
		}

	// calculate integral
	(*updater)(sums, fNIterations, integral, absprecision);
	relprecision = absprecision / integral;

	if (unsigned(fNIterations) >= fMCMCNIterationsMax)
	   BCLog::OutWarning("BCIntegrate::Integrate: Did not converge within maximum number of iterations");

	// print to log
	LogOutputAtEndOfIntegration(integral, absprecision, relprecision, fNIterations);

	fError = absprecision;
	return integral;
}

// ---------------------------------------------------------
double BCIntegrate::EvaluatorMC(std::vector<double> &sums, const std::vector<double> &point, bool &accepted) {
	const double value = Eval(point);

	// all samples accepted in sample-mean integration
	accepted = true;

	// add value to sum and sum of squares
	sums[0] += value;
	sums[1] += value * value;

	return value;
}

// ---------------------------------------------------------
void BCIntegrate::IntegralUpdaterMC(const std::vector<double> &sums, const int &nIterations, double &integral, double &absprecision) {
	// sample mean including the volume of the parameter space
	integral = sums[2] * sums[0] / nIterations;

	// unbiased estimate
	absprecision = sqrt((1.0 / (nIterations - 1)) * (sums[2] * sums[2] * sums[1] / double(nIterations) - integral * integral));
}

// ---------------------------------------------------------
double BCIntegrate::EvaluatorImportance(std::vector<double> &sums, const std::vector<double> &point, bool &accepted)
{
	// calculate value at random point
	double val_f = Eval(point);

	// calculate sampling distributions at that point
	double val_g = EvalSampling(point);

	// add ratio to sum and sum of squares
	if (val_g > 0 && val_f/val_g > fRandom->Rndm()) {
	   accepted = true;
	   sums[0] += val_f / val_g;
	   sums[1] += val_f * val_f / val_g / val_g;
	}
	else
	   accepted = false;

	return val_f;
}

// ---------------------------------------------------------
void BCIntegrate::IntegralUpdaterImportance(const std::vector<double> &sums, const int &nIterations, double &integral, double &absprecision) {
	integral = sums[0] / nIterations;
	absprecision = sqrt((1.0 / double(nIterations)) * (sums[1] / double(nIterations) - integral * integral));
}

// ---------------------------------------------------------
bool BCIntegrate::CheckMarginalizationAvailability(BCIntegrationMethod type) {
   switch(type)
   {
   case BCIntegrate::kIntMonteCarlo:
   case BCIntegrate::kIntImportance:
      return true;
   case BCIntegrate::kIntCuba:
#ifdef HAVE_CUBA_H
      return true;
#else
      return false;
#endif
   default:
      BCLog::OutError(TString::Format("BCIntegrate::Marginalize. Invalid marginalization method: %d.", type));
      break;
   }
   return false;
}

// ---------------------------------------------------------
TH1D* BCIntegrate::Marginalize(BCIntegrationMethod type, unsigned index)
{
    BCParameter * par = fParameters.Get(index);
    if ( ! par)
        return NULL;

    // create histogram
    TH1D * hist = new TH1D(TString::Format("h1_%s", par->GetName().data()),
                           TString::Format(";%s;(unit %s)^{-1};", par->GetName().data(), par->GetName().data()),
                           par->GetNbins(), par->GetLowerLimit(), par->GetUpperLimit());

    // fill histogram
    std::vector<unsigned> indices(1, index);
    Marginalize(hist, type, indices);

    return hist;
}

// ---------------------------------------------------------
TH2D* BCIntegrate::Marginalize(BCIntegrationMethod type, unsigned index1, unsigned index2)
{
   BCParameter * par1 = fParameters.Get(index1);
   BCParameter * par2 = fParameters.Get(index2);
   if ( ! par1 || par2)
      return NULL;

	// create histogram
	TH2D * hist = new TH2D(TString::Format("h2_%s_%s", par1->GetName().data(), par2->GetName().data()),
												 TString::Format(";%s;%s;(unit %s)^{-1}(unit %s)^{-1}",
															par1->GetName().data(), par2->GetName().data(),
															par1->GetName().data(), par2->GetName().data()),
															par1->GetNbins(), par1->GetLowerLimit(), par1->GetUpperLimit(),
															par2->GetNbins(), par2->GetLowerLimit(), par2->GetUpperLimit());

	// fill histogram
	std::vector<unsigned> indices;
	indices.push_back(index1);
	indices.push_back(index2);
	Marginalize(hist,type,indices);

	return hist;
}

// ---------------------------------------------------------
bool BCIntegrate::CheckMarginalizationIndices(TH1* hist, const std::vector<unsigned> &index) {
	if (index.size()==0) {
		BCLog::OutError("BCIntegrate::Marginalize : No marginalization parameters chosen.");
		return false;
	}

	if (index.size() >= 4 or index.size() > fParameters.Size()) {
		BCLog::OutError("BCIntegrate::Marginalize : Too many marginalization parameters.");
		return false;
	}

	if ((int)index.size()<hist->GetDimension()) {
		BCLog::OutError(TString::Format("BCIntegrate::Marginalize : Too few (%d) indices supplied for histogram dimension (%d)",(int)index.size(),hist->GetDimension()));
		return false;
	}

	for (unsigned i=0; i<index.size(); i++) {
		// check if indices are in bounds
		if ( ! fParameters.ValidIndex(index[i])) {
			BCLog::OutError(TString::Format("BCIntegrate::Marginalize : Parameter index (%d) out of bound.",index[i]));
			return false;
		}
		// check for duplicate indices
		for (unsigned j=0; j<index.size(); j++)
			if (i!=j and index[i]==index[j]) {
				BCLog::OutError(TString::Format("BCIntegrate::Marginalize : Parameter index (%d) appears more than once",index[i]));
				return false;
			}
	}
	return true;
}

// ---------------------------------------------------------
bool BCIntegrate::Marginalize(TH1* hist, BCIntegrationMethod type, const std::vector<unsigned> &index)
{
	if (!CheckMarginalizationAvailability(type))
		return false;

	if (!CheckMarginalizationIndices(hist,index))
		return false;

	// todo is this true?
	// Currently this does not seem to work for Cuba integration beyond 1 dimension
	if (type==kIntCuba and hist->GetDimension()>1) {
		BCLog::OutError("BCIntegrate::Marginalize : Cuba integration currently only functions for 1D marginalization.");
		return false;
	}

	// generate string output
	char * parnames = new char;
	for (unsigned i=0; i<index.size(); i++) {
		parnames = Form("%s %d (%s),",parnames,index[i],fParameters[i]->GetName().data());
	}
	if (index.size()==1)
		BCLog::OutDebug(TString::Format(" --> Marginalizing model w/r/t parameter%s using %s.",parnames,DumpMarginalizationMethod().c_str()));
	else
		BCLog::OutDebug(TString::Format(" --> Marginalizing model w/r/t parameters%s using %s.",parnames,DumpMarginalizationMethod().c_str()));

	// Store original integration limits and fixed flags
	// and unfix marginalization parameters;
	std::vector<double> origMins, origMaxs;
	std::vector<bool> origFix;
	for (unsigned i=0; i<index.size(); i++) {
	   BCParameter * par = fParameters[index[i]];
		origMins.push_back(par->GetLowerLimit());
		origMaxs.push_back(par->GetUpperLimit());
		origFix.push_back(par->Fixed());
		par->Fix(false);
	}

	// Set histogram title to indicate fixed variables
   std::string title;
   for (unsigned i=0; i < fParameters.Size(); ++i)
      if (fParameters[i]->Fixed()) {
         title += TString::Format(" (%s=%e)", fParameters[i]->GetName().data(), fParameters[i]->GetFixedValue());
      }
   if ( ! title.empty())
   {
      hist -> SetTitle(("Fixed: "+title).data());
   }

	// get bin count
	int N = hist->GetNbinsX() * hist->GetNbinsY() * hist->GetNbinsZ();

	// Store fNIterationsMax and change
	int nIterationsMax = fNIterationsMax;
	fNIterationsMax = fNIterationsMax / N;
	if (fNIterationsMax<fNIterationsMin)
		fNIterationsMax = fNIterationsMin;

	// todo not thread safe!
	// integrate each bin
	for (int i=1; i<=hist->GetNbinsX(); i++) {
		fParameters[index[0]]->SetLowerLimit(hist->GetXaxis()->GetBinLowEdge(i));
		fParameters[index[0]]->SetUpperLimit(hist->GetXaxis()->GetBinLowEdge(i+1));
		double binwidth1d = hist -> GetXaxis() -> GetBinWidth(i);
		if (hist->GetDimension()>1)
			for (int j=1; j<=hist->GetNbinsY(); j++) {
				fParameters[index[1]]->SetLowerLimit(hist -> GetYaxis() -> GetBinLowEdge(j));
				fParameters[index[1]]->SetUpperLimit(hist -> GetYaxis() -> GetBinLowEdge(j+1));
				double binwidth2d = binwidth1d * hist->GetYaxis()->GetBinWidth(j);
				if (hist->GetDimension()>2)
					for (int k=1; k<=hist->GetNbinsZ(); k++) {
						fParameters[index[2]]->SetLowerLimit(hist -> GetZaxis() -> GetBinLowEdge(k));
						fParameters[index[2]]->SetUpperLimit(hist -> GetZaxis() -> GetBinLowEdge(k+1));
						double binwidth3d = binwidth2d * hist->GetZaxis()->GetBinWidth(k);
						hist -> SetBinContent(i, j, k, Integrate(type)/binwidth3d);
						hist -> SetBinError  (i, j, k, GetError()/binwidth3d);
					}
				else {
					hist -> SetBinContent(i, j, Integrate(type)/binwidth2d);
					hist -> SetBinError  (i, j, GetError()/binwidth2d);
				}
			}
		else {
			hist -> SetBinContent(i, Integrate(type)/binwidth1d);
			hist -> SetBinError  (i, GetError()/binwidth1d);
		}
	}

	// restore original integration limits and fixed flags
	for (unsigned i=0; i<index.size(); i++) {
		fParameters[index[i]]->SetLowerLimit(origMins[i]);
		fParameters[index[i]]->SetUpperLimit(origMaxs[i]);
		fParameters[index[i]]->Fix(origFix[i]);
	}

	// restore original fNIterationsMax
	fNIterationsMax = nIterationsMax;

	return true;
}

// ---------------------------------------------------------
double BCIntegrate::GetRandomPoint(std::vector<double> &x)
{
   GetRandomVectorInParameterSpace(x);
   return Eval(x);
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomVectorUnitHypercube(std::vector<double> &x) const
{
   fRandom->RndmArray(x.size(), &x.front());
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomVectorInParameterSpace(std::vector<double> &x) const
{
	 GetRandomVectorUnitHypercube(x);

	 for (unsigned i = 0; i < fParameters.Size(); i++){
		 if (fParameters[i]->Fixed())
			 x[i] = fParameters[i]->GetFixedValue();
		 else
			 x[i] = fParameters[i]->GetLowerLimit() + x[i] * fParameters[i]->GetRangeWidth();
	 }
}
// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindMode(std::vector<double> start)
{
   if (fParameters.Size() < 1) {
      BCLog::OutError("FindMode : No parameters defined. Aborting.");
      return std::vector<double>();
   }

   BCLog::OutSummary(Form("Finding mode using %s", DumpOptimizationMethod().c_str()));

   switch (fOptimizationMethod) {
      case BCIntegrate::kOptSA:
         return FindModeSA(start);

      case BCIntegrate::kOptMinuit: {
         int printlevel = -1;
         if (BCLog::GetLogLevelScreen() <= BCLog::detail)
            printlevel = 0;
         if (BCLog::GetLogLevelScreen() <= BCLog::debug)
            printlevel = 1;
         return BCIntegrate::FindModeMinuit(start, printlevel);
      }

      case BCIntegrate::kOptMetropolis:
         return FindModeMCMC();

      default:
         BCLog::OutError(Form("BCModel::FindMode : Invalid mode finding method: %d", GetOptimizationMethod()));
         return std::vector<double>();
   }
}

// ---------------------------------------------------------
TMinuit * BCIntegrate::GetMinuit()
{
   if (!fMinuit)
      fMinuit = new TMinuit();

   return fMinuit;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindModeMinuit(std::vector<double> start, int printlevel)
{
   if (fParameters.Size() < 1) {
      BCLog::OutError("BCIntegrate::FindModeMinuit : No parameters defined. Aborting.");
      return std::vector<double>();
   }

   bool have_start = true;

   // check start values
   if (start.size() != fParameters.Size())
      have_start = false;

   // set global this
   ::BCIntegrateHolder::instance(this);

   // define minuit
   if (fMinuit)
      delete fMinuit;
   fMinuit = new TMinuit(fParameters.Size());

   // set print level
   fMinuit->SetPrintLevel(printlevel);

   // set function
   fMinuit->SetFCN(&BCIntegrate::FCNLikelihood);

   // set UP for likelihood
   fMinuit->SetErrorDef(0.5);

   // set parameters
   int flag;
   for (unsigned i = 0; i < fParameters.Size(); i++) {
      double starting_point = (fParameters[i]->GetUpperLimit() + fParameters[i]->GetLowerLimit()) / 2.;
      if(have_start)
         starting_point = start[i];
      if (fParameters[i]->Fixed())
         starting_point = fParameters[i]->GetFixedValue();
      fMinuit->mnparm(i,
                      fParameters[i]->GetName().data(),
                      starting_point,
                      fParameters[i]->GetRangeWidth() / 100.,
                      fParameters[i]->GetLowerLimit(),
                      fParameters[i]->GetUpperLimit(),
                      flag);
   }

   for (unsigned i = 0; i < fParameters.Size(); i++)
      if (fParameters[i]->Fixed())
         fMinuit->FixParameter(i);

   // do mcmc minimization
   //   fMinuit->mnseek();

   // do minimization
   fMinuit->mnexcm("MIGRAD", fMinuitArglist, 2, flag);

   // improve search for local minimum
   //   fMinuit->mnimpr();

   // copy flag and result
   fMinuitErrorFlag = flag;
   std::vector<double> localMode(fParameters.Size(), 0);
   std::vector<double> errors(fParameters.Size(), 0);
   for (unsigned i = 0; i < fParameters.Size(); i++) {
      fMinuit->GetParameter(i, localMode[i], errors[i]);
   }

   // set best fit parameters
   if ( (-fMinuit->fAmin > fLogMaximum and !fFlagIgnorePrevOptimization) or fFlagIgnorePrevOptimization) {
      fBestFitParameters = localMode;
      fBestFitParameterErrors = errors;
      fLogMaximum = -fMinuit->fAmin;

      // set optimization method used to find the mode
      fOptimizationMethodMode = BCIntegrate::kOptMinuit;
   }

   // delete minuit
   delete fMinuit;
   fMinuit = 0;

   return localMode;
}

// ---------------------------------------------------------
void BCIntegrate::InitializeSATree()
{
  if (fTreeSA)
    delete fTreeSA;
  fTreeSA = new TTree("SATree", "SATree");

   fTreeSA->Branch("Iteration",      &fSANIterations,   "iteration/I");
   // todo make consistent with examples and test
   // remove because would require fixed number of parameters, and just superfluous info
//   fTreeSA->Branch("NParameters",    &fNvar,            "parameters/I");
   fTreeSA->Branch("Temperature",    &fSATemperature,   "temperature/D");
   fTreeSA->Branch("LogProbability", &fSALogProb,       "log(probability)/D");

   for (unsigned i = 0; i < fParameters.Size(); ++i)
      fTreeSA->Branch(TString::Format("Parameter%i", i), &fSAx[i], TString::Format("parameter %i/D", i));
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindModeSA(std::vector<double> start)
{
   // note: if f(x) is the function to be minimized, then
   // f(x) := - LogEval(parameters)

   bool have_start = true;
   std::vector<double> x, y, best_fit; // vectors for current state, new proposed state and best fit up to now
   double fval_x, fval_y, fval_best_fit; // function values at points x, y and best_fit (we save them rather than to re-calculate them every time)
   int t = 1; // time iterator

   // check start values
   if (start.size() != fParameters.Size())
      have_start = false;

   // if no starting point is given, set to center of parameter space
   if ( !have_start ) {
      start.clear();
      for (unsigned i=0; i<fParameters.Size(); i++)
         if (fParameters[i]->Fixed())
            start.push_back(fParameters[i]->GetFixedValue());
         else
            start.push_back((fParameters[i]->GetLowerLimit() + fParameters[i]->GetUpperLimit()) / 2.);
   }

   // set current state and best fit to starting point
   x = start;
   best_fit = start;

   // calculate function value at starting point
   fval_x = fval_best_fit = LogEval(x);

   // run while still "hot enough"
   while ( SATemperature(t) > fSATmin ) {
      // generate new state
      y = GetProposalPointSA(x, t);

      // check if the proposed point is inside the phase space
      // if not, reject it
      bool in_range = true;

      for (unsigned i = 0; i < fParameters.Size(); i++)
         if ( !fParameters[i]->IsValid(y[i]))
            in_range = false;

      if ( in_range ){
         // calculate function value at new point
         fval_y = LogEval(y);

         // is it better than the last one?
         // if so, update state and chef if it is the new best fit...
         if (fval_y >= fval_x) {
            x = y;

            fval_x = fval_y;

            if (fval_y > fval_best_fit) {
               best_fit = y;
               fval_best_fit = fval_y;
            }
         }
         // ...else, only accept new state w/ certain probability
         else {
            if (fRandom->Rndm() <= exp( (fval_y - fval_x) / SATemperature(t) )) {
               x = y;
               fval_x = fval_y;
            }
         }
      }

      // update tree variables
      fSANIterations = t;
      fSATemperature = SATemperature(t);
      fSALogProb = fval_x;
      fSAx = x;

      // fill tree
      if (fFlagWriteSAToFile)
         fTreeSA->Fill();

      // increate t
      t++;
   }

   if ( fval_best_fit > fLogMaximum or fFlagIgnorePrevOptimization) {
      // set best fit parameters
      fBestFitParameters = best_fit;
      fBestFitParameterErrors.assign(fParameters.Size(),-1);
      fLogMaximum = fval_best_fit;

      // set optimization moethod used to find the mode
      fOptimizationMethodMode = BCIntegrate::kOptSA;
   }

   return best_fit;
}

// ---------------------------------------------------------
double BCIntegrate::SATemperature(double t)
{
   // do we have Cauchy (default), Boltzmann or custom annealing schedule?
   if (fSASchedule == BCIntegrate::kSABoltzmann)
      return SATemperatureBoltzmann(t);
   else if (fSASchedule == BCIntegrate::kSACauchy)
      return SATemperatureCauchy(t);
   else
      return SATemperatureCustom(t);
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureBoltzmann(double t)
{
   return fSAT0 / log((double)(t + 1));
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureCauchy(double t)
{
   return fSAT0 / (double)t;
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureCustom(double /*t*/)
{
   BCLog::OutError("BCIntegrate::SATemperatureCustom : No custom temperature schedule defined");
   return 0.;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSA(const std::vector<double> &x, int t)
{
   // do we have Cauchy (default), Boltzmann or custom annealing schedule?
   if (fSASchedule == BCIntegrate::kSABoltzmann)
      return GetProposalPointSABoltzmann(x, t);
   else if (fSASchedule == BCIntegrate::kSACauchy)
      return GetProposalPointSACauchy(x, t);
   else
      return GetProposalPointSACustom(x, t);
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSABoltzmann(const std::vector<double> &x, int t)
{
   std::vector<double> y;
   y.clear();
   double new_val, norm;

   for (unsigned i = 0; i < fParameters.Size(); i++) {
      if (fParameters[i]->Fixed()) {
         y.push_back(fParameters[i]->GetFixedValue());
      }
      else {
         norm = fParameters[i]->GetRangeWidth() * SATemperature(t) / 2.;
         new_val = x[i] + norm * fRandom->Gaus();
         y.push_back(new_val);
      }
   }

   return y;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSACauchy(const std::vector<double> &x, int t)
{
   std::vector<double> y;
   y.clear();

   if (fParameters.Size() == 1) {
      double cauchy, new_val, norm;

      if (fParameters[0]->Fixed()) {
         y.push_back(fParameters[0]->GetFixedValue());
      }
      else {
         norm = fParameters[0]->GetRangeWidth() * SATemperature(t) / 2.;
         cauchy = tan(3.14159 * (fRandom->Rndm() - 0.5));
         new_val = x[0] + norm * cauchy;
         y.push_back(new_val);
      }
   }
   else {
      // use sampling to get radial n-dim Cauchy distribution

      // first generate a random point uniformly distributed on a
      // fNvar-dimensional hypersphere
      y = SAHelperGetRandomPointOnHypersphere();

      // scale the vector by a random factor determined by the radial
      // part of the fNvar-dimensional Cauchy distribution
      double radial = SATemperature(t) * SAHelperGetRadialCauchy();

      // scale y by radial part and the size of dimension i in phase space
      // afterwards, move by x
      for (unsigned i = 0; i < fParameters.Size(); i++) {
         if (fParameters[i]->Fixed()) {
            y[i] = fParameters[i]->GetFixedValue(); }
         else {
            y[i] = fParameters[i]->GetRangeWidth() * y[i] * radial / 2. + x[i];
         }
      }
   }

   return y;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSACustom(const std::vector<double> & /*x*/, int /*t*/)
{
   BCLog::OutError("BCIntegrate::GetProposalPointSACustom : No custom proposal function defined");
   return std::vector<double>(fParameters.Size());
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::SAHelperGetRandomPointOnHypersphere()
{
   std::vector<double> rand_point(fParameters.Size());

   // This method can only be called with fNvar >= 2 since the 1-dim case
   // is already hard wired into the Cauchy annealing proposal function.
   // To speed things up, hard-code fast method for 2 and dimensions.
   // The algorithm for 2D can be found at
   // http://mathworld.wolfram.com/CirclePointPicking.html
   // For 3D just using ROOT's algorithm.
   if (fParameters.Size() == 2) {
      double x1, x2, s;
      do {
         x1 = fRandom->Rndm() * 2. - 1.;
         x2 = fRandom->Rndm() * 2. - 1.;
         s = x1*x1 + x2*x2;
      }
      while (s >= 1);

      rand_point[0] = (x1*x1 - x2*x2) / s;
      rand_point[1] = (2.*x1*x2) / s;
   }
   else if (fParameters.Size() == 3) {
      fRandom->Sphere(rand_point[0], rand_point[1], rand_point[2], 1.0);
   }
   else {
      double s = 0.,
         gauss_num;

      for (unsigned i = 0; i < fParameters.Size(); i++) {
         gauss_num = fRandom->Gaus();
         rand_point[i] = gauss_num;
         s += gauss_num * gauss_num;
      }
      s = sqrt(s);

      for (unsigned i = 0; i < fParameters.Size(); i++)
         rand_point[i] = rand_point[i] / s;
   }

   return rand_point;
}

// ---------------------------------------------------------
double BCIntegrate::SAHelperGetRadialCauchy()
{
   // theta is sampled from a rather complicated distribution,
   // so first we create a lookup table with 10000 random numbers
   // once and then, each time we need a new random number,
   // we just look it up in the table.
   double theta;

   // static vectors for theta-sampling-map
   static double map_u[10001];
   static double map_theta[10001];
   static bool initialized = false;
   static unsigned map_dimension = 0;

   // is the lookup-table already initialized? if not, do it!
   if (!initialized or map_dimension != fParameters.Size()) {
      double init_theta;
      double beta = SAHelperSinusToNIntegral(fParameters.Size() - 1, 1.57079632679);

      for (int i = 0; i <= 10000; i++) {
         init_theta = 3.14159265 * (double)i / 5000.;
         map_theta[i] = init_theta;

         map_u[i] = SAHelperSinusToNIntegral(fParameters.Size() - 1, init_theta) / beta;
      }

      map_dimension = fParameters.Size();
      initialized = true;
   } // initializing is done.

   // generate uniform random number for sampling
   double u = fRandom->Rndm();

   // Find the two elements just greater than and less than u
   // using a binary search (O(log(N))).
   int lo = 0;
   int up = 10000;
   int mid;

   while (up != lo) {
      mid = ((up - lo + 1) / 2) + lo;

      if (u >= map_u[mid])
         lo = mid;
      else
         up = mid - 1;
   }
   up++;

   // perform linear interpolation:
   theta = map_theta[lo] + (u - map_u[lo]) / (map_u[up] - map_u[lo]) * (map_theta[up] - map_theta[lo]);

   return tan(theta);
}

// ---------------------------------------------------------
double BCIntegrate::SAHelperSinusToNIntegral(int dim, double theta)
{
   if (dim < 1)
      return theta;
   else if (dim == 1)
      return (1. - cos(theta));
   else if (dim == 2)
      return 0.5 * (theta - sin(theta) * cos(theta));
   else if (dim == 3)
      return (2. - sin(theta) * sin(theta) * cos(theta) - 2. * cos(theta)) / 3.;
   else
      return - pow(sin(theta), (double)(dim - 1)) * cos(theta) / (double)dim
         + (double)(dim - 1) / (double)dim
         * SAHelperSinusToNIntegral(dim - 2, theta);
}

// ---------------------------------------------------------
void BCIntegrate::FCNLikelihood(int & /*npar*/, double * /*grad*/, double &fval, double * par, int /*flag*/)
{
   // copy parameters
   static std::vector<double> parameters;

   // calculate number of active + fixed parameters
   // remember: npar is just the number of _active_ parameters while
   // par is a vector of _all_ parameters
   int nparameters = ::BCIntegrateHolder::instance()->GetNParameters();

   // adjust size if needed
   parameters.resize(nparameters, 0.0);

   // copy values
   std::copy(par, par + nparameters, parameters.begin());

   // evaluate, for efficiency don't check if npar matches
   fval = - ::BCIntegrateHolder::instance()->LogEval(parameters);
}

// ---------------------------------------------------------

//void BCIntegrate::FindModeMCMC(int flag_run)
std::vector<double> BCIntegrate::FindModeMCMC()
{
   const std::vector<double> tempMode = fBestFitParameters;
   const std::vector<double> tempErrors= fBestFitParameterErrors;
   const double tempMaximum = fLogMaximum;

   // call PreRun
   MCMCMetropolisPreRun();

   // loop over all chains and find the maximum point
   double logProb = -std::numeric_limits<double>::max();
   unsigned index = 0;
   for (unsigned i = 0; i < fMCMCNChains; ++i) {
      if (MCMCGetMaximumLogProb().at(i) > logProb){
         index = i;
         logProb = MCMCGetMaximumLogProb().at(i);
      }
   }
   std::vector<double> localMode = MCMCGetMaximumPoint(index);

   if ( logProb > tempMaximum  or fFlagIgnorePrevOptimization) {
      fBestFitParameters = localMode;
      fBestFitParameterErrors.assign(fParameters.Size(),-1.);   // error undefined
      fLogMaximum = logProb;

      // set optimization method used to find the mode
      fOptimizationMethodMode = BCIntegrate::kOptMetropolis;
   }
   else if (!fFlagIgnorePrevOptimization) {
      fBestFitParameters = tempMode;
      fBestFitParameterErrors = tempErrors;
      fLogMaximum = tempMaximum;
   }

   return MCMCGetMaximumPoint(index);
}

// ---------------------------------------------------------
void BCIntegrate::SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method)
{
   if (method >= BCIntegrate::NIntMethods) {
      BCLog::OutError(Form("BCIntegrate::SetIntegrationMethod: Invalid method '%d' ", method));
      return;
   }
   if (method == BCIntegrate::kIntCuba) {
#ifndef HAVE_CUBA_H
      BCLog::OutError("BCIntegrate::SetIntegrationMethod: Cuba not enabled during configure");
#endif
   }
   fIntegrationMethod = method;
}

// ---------------------------------------------------------
void BCIntegrate::SetCubaIntegrationMethod(BCIntegrate::BCCubaMethod type)
{
#ifdef HAVE_CUBA_H
   switch(type) {
      case BCIntegrate::kCubaVegas:
      case BCIntegrate::kCubaSuave:
      case BCIntegrate::kCubaDivonne:
      case BCIntegrate::kCubaCuhre:
         fCubaIntegrationMethod = type;
         return;
      default:
         BCLog::OutError(TString::Format("Integration method of type %d is not defined for Cuba",type));
         return;
   }
#else
   BCLog::OutError("SetCubaIntegrationMethod: Cuba not enabled during configure");
#endif
}

// ---------------------------------------------------------
int BCIntegrate::CubaIntegrand(const int * ndim, const double xx[],
      const int * /*ncomp*/, double ff[], void * userdata)
{
   BCIntegrate * local_this = static_cast<BCIntegrate *>(userdata);

   // scale variables
   double jacobian = 1.0;

   // create local parameter vector
   // important for thread safety, though not super efficient
   std::vector<double> scaled_parameters(local_this->fParameters.Size());

   // stay in sync with the possible lower number of parameters
   // that cuba sees due to fixing in BAT
   unsigned cubaIndex = 0;
   unsigned batIndex = 0;
   for (batIndex = 0; batIndex < local_this->fParameters.Size(); ++batIndex) {
  	 const BCParameter * p = local_this->fParameters[batIndex];

		 // get the scaled parameter value
		 if (p->Fixed())
			 scaled_parameters[batIndex] = p->GetFixedValue();
		 else {
			 // convert from unit hypercube to actual parameter hyperrectangle
			 scaled_parameters[batIndex] = p->GetLowerLimit() + xx[cubaIndex] * p->GetRangeWidth();

			 // multiply range to jacobian
			 jacobian *= p->GetRangeWidth();

			 // one more parameter that cuba varies
			 ++cubaIndex;
		 }
   }

   if (cubaIndex != unsigned(*ndim))
  	 BCLog::OutError(Form("BCIntegrate::CubaIntegrand: mismatch between variable parameters"
  			 	 	 	 	 	 	 	 	  "in BAT (%d) and Cuba(%d)", batIndex, cubaIndex));

   // call function to integrate
   ff[0] = local_this->Eval(scaled_parameters);

   // multiply jacobian
   ff[0] *= jacobian;

   return 0;
}

// ---------------------------------------------------------
double BCIntegrate::IntegrateCuba(BCCubaMethod cubatype) {
#if HAVE_CUBA_H
   LogOutputAtStartOfIntegration(kIntCuba, cubatype);

   // integrand has only one component
   static const int ncomp     = 1;

   // don't store integration in a slot for reuse
   static const int gridno    = -1;

   // we don't want dumps of internal state
   static const char * statefile = "";

   // cuba returns info in these variables
   int fail = 0;
   int nregions = 0;
   std::vector<double> integral(ncomp, -1);
   std::vector<double> error(ncomp, -1);
   std::vector<double> prob(ncomp, -1);

   // reset number of iterations
   fNIterations = 0;

   // Cuba needs int variable
   int nIntegrationVariables = GetNIntegrationVariables();

   switch (cubatype) {

   case BCIntegrate::kCubaVegas:
      Vegas(nIntegrationVariables, ncomp,
            &BCIntegrate::CubaIntegrand, static_cast<void *>(this),
            fRelativePrecision, fAbsolutePrecision,
            fCubaVegasOptions.flags, fRandom->GetSeed(),
            fNIterationsMin, fNIterationsMax,
            fCubaVegasOptions.nstart, fCubaVegasOptions.nincrease, fCubaVegasOptions.nbatch,
            gridno, statefile,
            &fNIterations, &fail,
            &integral[0], &error[0], &prob[0]);
      break;

   case BCIntegrate::kCubaSuave:
      Suave(nIntegrationVariables, ncomp,
            &BCIntegrate::CubaIntegrand, static_cast<void *>(this),
            fRelativePrecision, fAbsolutePrecision,
            fCubaSuaveOptions.flags, fRandom->GetSeed(),
            fNIterationsMin, fNIterationsMax,
            fCubaSuaveOptions.nnew, fCubaSuaveOptions.flatness,
            statefile,
            &nregions, &fNIterations, &fail,
            &integral[0], &error[0], &prob[0]);
      break;

   case BCIntegrate::kCubaDivonne:
      if (nIntegrationVariables < 2 or nIntegrationVariables > 33)
         BCLog::OutError("BCIntegrate::IntegrateCuba(Divonne): Divonne only works in 1 < d < 34");
      else {
         // no extra info supported
         static const int ngiven = 0;
         static const int nextra = ngiven;
         Divonne(nIntegrationVariables, ncomp,
               &BCIntegrate::CubaIntegrand, static_cast<void *>(this),
               fRelativePrecision, fAbsolutePrecision,
               fCubaDivonneOptions.flags, fRandom->GetSeed(),
               fNIterationsMin, fNIterationsMax,
               fCubaDivonneOptions.key1, fCubaDivonneOptions.key2, fCubaDivonneOptions.key3,
               fCubaDivonneOptions.maxpass, fCubaDivonneOptions.border,
               fCubaDivonneOptions.maxchisq, fCubaDivonneOptions.mindeviation,
               ngiven, nIntegrationVariables /*ldxgiven*/, NULL, nextra, NULL,
               statefile,
               &nregions, &fNIterations, &fail,
               &integral[0], &error[0], &prob[0]);
      }
      break;

   case BCIntegrate::kCubaCuhre:
      if (nIntegrationVariables < 2)
         BCLog::OutError("BCIntegrate::IntegrateCuba(Cuhre): Cuhre(cubature) only works in d > 1");

      Cuhre(nIntegrationVariables, ncomp,
            &BCIntegrate::CubaIntegrand, static_cast<void *>(this),
            fRelativePrecision, fAbsolutePrecision,
            fCubaCuhreOptions.flags, fNIterationsMin, fNIterationsMax,
            fCubaCuhreOptions.key,
            statefile,
            &nregions, &fNIterations, &fail,
            &integral[0], &error[0], &prob[0]);
      break;

   case BCIntegrate::NCubaMethods:
   default:
      BCLog::OutError("Cuba integration method not set.");
      error[0] = -1;
      integral[0] = -1;
      break;
   }

   fError = error[0];
   double result = integral[0];

   if (fail != 0) {
      BCLog::OutWarning(" Warning, integral did not converge with the given set of parameters. ");
      BCLog::OutWarning(TString::Format(" neval    = %d", fNIterations));
      BCLog::OutWarning(TString::Format(" fail     = %d", fail));
      BCLog::OutWarning(TString::Format(" integral = %e", result));
      BCLog::OutWarning(TString::Format(" error    = %e", fError));
      BCLog::OutWarning(TString::Format(" prob     = %e", prob[0]));

      // handle cases in which cuba does not alter integral
      if (fNIterations == 0)
      {
         fError = -1;
         result = -1;
      }

   } else
      LogOutputAtEndOfIntegration(result,fError,fError/result,fNIterations);

   return result;
#else
   BCLog::OutError("IntegrateCuba: Cuba not enabled during configure");
   return -1;
#endif
}

// ---------------------------------------------------------
void BCIntegrate::SAInitialize()
{
   fSAx.clear();
   fSAx.assign(fParameters.Size(), 0.0);
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpIntegrationMethod(BCIntegrate::BCIntegrationMethod type)
{
   switch(type) {
      case BCIntegrate::kIntMonteCarlo:
         return "Sampled Mean Monte Carlo";
      case BCIntegrate::kIntImportance:
         return "Importance Sampling";
      case BCIntegrate::kIntCuba:
         return "Cuba";
      default:
         return "Undefined";
   }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpOptimizationMethod(BCIntegrate::BCOptimizationMethod type)
{
   switch(type) {
      case BCIntegrate::kOptSA:
         return "Simulated Annealing";
      case BCIntegrate::kOptMetropolis:
         return "Metropolis MCMC";
      case BCIntegrate::kOptMinuit:
         return "Minuit";
      default:
         return "Undefined";
   }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpCubaIntegrationMethod(BCIntegrate::BCCubaMethod type)
{
   switch(type) {
      case BCIntegrate::kCubaVegas:
         return "Vegas";
      case BCIntegrate::kCubaSuave:
         return "Suave";
      case BCIntegrate::kCubaDivonne:
         return "Divonne";
      case BCIntegrate::kCubaCuhre:
         return "Cuhre";
      default:
         return "Undefined";
   }
}

namespace BCCubaOptions
{
General::General() :
      ncomp(1),
      flags(0),
      nregions(0),
      neval(0),
      fail(0),
      error(0),
      prob(0)
{}

General::~General()
{}

/*
 * copy values from demo-c.c shipping with cuba 3.2
 * for three-dimensionsal examples.
 */

Vegas::Vegas() :
      General(),
      nstart(1000),
      nincrease(500),
      nbatch(1000),
      gridno(0)
{}

Suave::Suave() :
      General(),
      nnew(1000),
      flatness(25)
{}

Divonne::Divonne() :
      General(),
      key1(47),
      key2(1),
      key3(1),
      maxpass(5),
      border(0),
      maxchisq(10),
      mindeviation(0.25)
{}

Cuhre::Cuhre() :
      General(),
      key(0) // let cuba choose default cubature rule
{}
}
