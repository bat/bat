/*
 * Copyright (C) 2008-2010, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BAT/BCModelManager.h"

#include "BAT/BCModel.h"
#include "BAT/BCDataSet.h"
#include "BAT/BCDataPoint.h"
#include "BAT/BCLog.h"
#include "BAT/BCErrorCodes.h"

#include <fstream>
#include <iostream>

#include <TString.h>

// ---------------------------------------------------------

BCModelManager::BCModelManager()
{
   fModelContainer = new BCModelContainer();
   fDataSet = 0;
}

// ---------------------------------------------------------

BCModelManager::~BCModelManager()
{
   delete fModelContainer;

   if (fDataSet)
      delete fDataSet;
}

// ---------------------------------------------------------

BCModelManager::BCModelManager(const BCModelManager & modelmanager)
{
   modelmanager.Copy(*this);
}

// ---------------------------------------------------------

BCModelManager & BCModelManager::operator = (const BCModelManager & modelmanager)
{
   if (this != &modelmanager)
      modelmanager.Copy(* this);

   return * this;
}

// ---------------------------------------------------------

void BCModelManager::SetDataSet(BCDataSet * dataset)
{
   // set data set
   fDataSet = dataset;

   // set data set of all models in the manager
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetDataSet(fDataSet);
}

// ---------------------------------------------------------

void BCModelManager::SetSingleDataPoint(BCDataPoint * datapoint)
{
   // create new data set consisting of a single data point
   BCDataSet * dataset = new BCDataSet();

   // add the data point
   dataset->AddDataPoint(datapoint);

   // set this new data set
   SetDataSet(dataset);
}

// ---------------------------------------------------------

void BCModelManager::SetSingleDataPoint(BCDataSet * dataset, unsigned int index)
{
   if (index < 0 || index > dataset->GetNDataPoints())
      return;

   SetSingleDataPoint(dataset->GetDataPoint(index));
}

// ---------------------------------------------------------

void BCModelManager::AddModel(BCModel * model, double probability)
{
   // create index
   unsigned int index = fModelContainer->size();

   // set index of new model
   model->SetIndex(index);

   // set a priori probability of new model
   model->SetModelAPrioriProbability(probability);

   // set data set
   model->SetDataSet(fDataSet);

   // fill model into container
   fModelContainer->push_back(model);
}

// ---------------------------------------------------------
// DEBUG DELETE?
/*
void BCModelManager::SetNIterationsMax(int niterations)
{
   // set maximum number of iterations of all models in the manager

   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetNIterationsMax(niterations);
}

// ---------------------------------------------------------
*/

void BCModelManager::SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method)
{
   // set integration method for all models registered

   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetIntegrationMethod(method);
}

// ---------------------------------------------------------

void BCModelManager::SetMarginalizationMethod(BCIntegrate::BCMarginalizationMethod method)
{
   // set marginalization method for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetMarginalizationMethod(method);
};

// ---------------------------------------------------------

void BCModelManager::SetOptimizationMethod(BCIntegrate::BCOptimizationMethod method)
{
   // set mode finding method for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetOptimizationMethod(method);
}

// ---------------------------------------------------------

void BCModelManager::SetNiterationsPerDimension(unsigned int niterations)
{
   // set number of iterations per dimension for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetNiterationsPerDimension(niterations);
}

// ---------------------------------------------------------

void BCModelManager::SetNSamplesPer2DBin(unsigned int n)
{
   // set samples per 2d bin for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetNSamplesPer2DBin(n);
}

// ---------------------------------------------------------

void BCModelManager::SetRelativePrecision(double relprecision)
{
   // set relative precision for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetRelativePrecision(relprecision);
}

// ---------------------------------------------------------

void BCModelManager::SetNbins(unsigned int n)
{
   // set number of bins for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->BCIntegrate::SetNbins(n);
}

// ---------------------------------------------------------

void BCModelManager::SetFillErrorBand(bool flag)
{
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetFillErrorBand(flag);
}

// ---------------------------------------------------------

void BCModelManager::SetFitFunctionIndexX(int index)
{
   // set fit function x index for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetFitFunctionIndexX(index);
}

// ---------------------------------------------------------

void BCModelManager::SetFitFunctionIndexY(int index)
{
   // set  fit function y index for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetFitFunctionIndexY(index);
}

// ---------------------------------------------------------

void BCModelManager::SetFitFunctionIndices(int indexx, int indexy)
{
   // set fit function indices for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetFitFunctionIndices(indexx, indexy);
}

// ---------------------------------------------------------

void BCModelManager::SetDataPointLowerBoundaries(BCDataPoint * datasetlowerboundaries)
{
   // set lower boundary point for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetDataPointLowerBoundaries(datasetlowerboundaries);
}

// ---------------------------------------------------------

void BCModelManager::SetDataPointUpperBoundaries(BCDataPoint * datasetupperboundaries)
{
   // set upper boundary point for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetDataPointUpperBoundaries(datasetupperboundaries);
}

// ---------------------------------------------------------

void BCModelManager::SetDataPointLowerBoundary(int index, double lowerboundary)
{
   // set lower bounday values for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetDataPointLowerBoundary(index, lowerboundary);
}

// ---------------------------------------------------------

void BCModelManager::SetDataPointUpperBoundary(int index, double upperboundary)
{
   // set upper boundary values for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetDataPointUpperBoundary(index, upperboundary);
}

// ---------------------------------------------------------

void BCModelManager::SetDataBoundaries(int index, double lowerboundary, double upperboundary)
{
   // set lower and upper boundary values for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->SetDataBoundaries(index, lowerboundary, upperboundary);
}

// ---------------------------------------------------------

void BCModelManager::FixDataAxis(int index, bool fixed)
{
   // fix axis for all models
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->FixDataAxis(index, fixed);
}

// ---------------------------------------------------------

void BCModelManager::SetNChains(unsigned int n)
{
   // set number of Markov chains for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->MCMCSetNChains(n);
}

// ---------------------------------------------------------

int BCModelManager::ReadDataFromFileTree(const char * filename, const char * treename, const char * branchnames)
{
   if (fModelContainer->size() < 0) {
      BCLog::Out(BCLog::warning, BCLog::warning, "BCModelManager::ReadDataFromFileTree : No model defined.");
      return ERROR_NOMODELS;
   }

   // create data set
   if (!fDataSet)
      fDataSet = new BCDataSet();
   else
      fDataSet->Reset();

   // read data from tree
   int read_file = fDataSet->ReadDataFromFileTree(filename, treename, branchnames);

   if (read_file >=0) {
      SetDataSet(fDataSet);

      for (unsigned int i = 0; i < GetNModels(); i++)
         fModelContainer->at(i)->SetDataSet(fDataSet);
   }
   else if (read_file == ERROR_FILENOTFOUND) {
      delete fDataSet;
      return ERROR_FILENOTFOUND;
   }

   return 0;
}

// ---------------------------------------------------------

int BCModelManager::ReadDataFromFileTxt(const char * filename, int nbranches)
{
   if (fModelContainer->size() < 0) {
      BCLog::Out(BCLog::warning, BCLog::warning, "BCModelManager::ReadDataFromFileTree. No model defined.");
      return ERROR_NOMODELS;
   }

   // create data set
   if (!fDataSet)
      fDataSet = new BCDataSet();
   else
      fDataSet->Reset();

   // read data from txt file
   int read_file = fDataSet->ReadDataFromFileTxt(filename, nbranches);

   if (read_file >=0) {
      SetDataSet(fDataSet);

      for (unsigned int i = 0; i < GetNModels(); i++)
         fModelContainer->at(i)->SetDataSet(fDataSet);
   }
   else {
      delete fDataSet;
      return ERROR_FILENOTFOUND;
   }

   return 0;
}

// ---------------------------------------------------------

void BCModelManager::Normalize()
{
   // initialize likelihood norm
   double normalization = 0.0;

//   BCLog::Out(BCLog::summary, BCLog::summary, "Running normalization of all models.");
   BCLog::OutSummary("Running normalization of all models.");

   for (unsigned int i = 0; i < GetNModels(); i++) {
      fModelContainer->at(i)->Normalize();

      // add to total normalization
      normalization += (fModelContainer->at(i)->GetNormalization() * fModelContainer->at(i)->GetModelAPrioriProbability());
   }

   // set model a posteriori probabilities
   for (unsigned int i = 0; i < GetNModels(); i++)
      fModelContainer->at(i)->SetModelAPosterioriProbability(
            (fModelContainer->at(i)->GetNormalization() * fModelContainer->at(i)->GetModelAPrioriProbability()) /
            normalization);
}

// ---------------------------------------------------------

double BCModelManager::BayesFactor(const unsigned int imodel1, const unsigned int imodel2)
{
   // Bayes Factors are the likelihoods integrated over the parameters
   // Is this equal to the posteriors?
   //    NOOOO
   // But it is equal to normalization factors.

   const double norm1 = fModelContainer->at(imodel1)->GetNormalization();
   const double norm2 = fModelContainer->at(imodel2)->GetNormalization();

   // check model 1
   if(norm1<0.) {
      BCLog::OutError(
         Form("Model %s (index %d) not normalized. Cannot calculate Bayes factor.",
              fModelContainer->at(imodel1)->GetName().data(),imodel1));
      return -1.;
   }

   // check model 2
   if(norm2<0.) {
      BCLog::OutError(
         Form("Model %s (index %d) not normalized. Cannot calculate Bayes factor.",
              fModelContainer->at(imodel2)->GetName().data(),imodel2));
      return -1.;
   }

   // denominator cannot be zero
   if(norm2==0. && norm1!=0.) {// not good since norm2 is double !!!
      BCLog::OutError(
         Form("Model %s (index %d) has ZERO probability. Bayes factor is infinite.",
              fModelContainer->at(imodel2)->GetName().data(),imodel2));
      return -1.;
   }

   // denominator cannot be zero unless also numerator is zero
   if(norm2==0. && norm1==0.) {// not good since norm2 and norm1 are both double !!!
      BCLog::OutWarning(
         Form("Models %s and %s have ZERO probability. Bayes factor is unknown. Returning 1.",
              fModelContainer->at(imodel2)->GetName().data(),fModelContainer->at(imodel1)->GetName().data()));
      return 1.;
   }

   // now calculate the factor
   return norm1/norm2;
}

// ---------------------------------------------------------

void BCModelManager::FindMode()
{
   // finds mode for all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->FindMode();
}

// ---------------------------------------------------------

void BCModelManager::MarginalizeAll()
{
   // marginalizes all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->MarginalizeAll();
}

// ---------------------------------------------------------

void BCModelManager::WriteMarkovChain(bool flag)
{
   // marginalizes all models registered
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->WriteMarkovChain(flag);
}

// ---------------------------------------------------------

void BCModelManager::CalculatePValue(bool flag_histogram)
{
   // calculate p-value for all models
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->CalculatePValue(GetModel(i)->GetBestFitParameters(), flag_histogram);
}

// ---------------------------------------------------------

void BCModelManager::PrintSummary(const char * file)
{
   ofstream out;
   std::streambuf * old_buffer = 0;

   if(file) {
      out.open(file);
      if (!out.is_open()) {
         std::cerr<<"Couldn't open file "<<file<<std::endl;
         return;
      }
      old_buffer = std::cout.rdbuf(out.rdbuf());
   }

   // model summary
   int nmodels = fModelContainer->size();
   std::cout<<std::endl
            <<"======================================"<<std::endl
            <<" Summary"<<std::endl
            <<"======================================"<<std::endl
            <<std::endl
            <<" Number of models               : "<<nmodels<<std::endl
            <<std::endl
            <<" - Models:"<<std::endl;

   for (int i = 0; i < nmodels; i++)
      fModelContainer->at(i)->PrintSummary();

   // data summary
   std::cout<<" - Data:"<<std::endl
            <<std::endl
            <<"     Number of entries: "<<fDataSet->GetNDataPoints()<<std::endl
            <<std::endl;

   std::cout<<"======================================"<<std::endl
            <<" Model comparison"<<std::endl
            <<std::endl;

   // probability summary
   std::cout<<" - A priori probabilities:"<<std::endl<<std::endl;
  
   for (int i=0; i<nmodels; i++)
      std::cout<<"     p("<< fModelContainer->at(i)->GetName()
               <<") = "<< fModelContainer->at(i)->GetModelAPrioriProbability()
               <<std::endl;
   std::cout<<std::endl;

   std::cout<<" - A posteriori probabilities:"<<std::endl<<std::endl;

   for (int i = 0; i < nmodels; i++)
      std::cout<<"     p("<< fModelContainer->at(i)->GetName()
               <<" | data) = "<< fModelContainer->at(i)->GetModelAPosterioriProbability()
               <<std::endl;
   std::cout<<std::endl;

   std::cout<<"======================================"<<std::endl<<std::endl;

   if (file)
      std::cout.rdbuf(old_buffer);

}

// ---------------------------------------------------------

void BCModelManager::PrintModelComparisonSummary(const char * file)
{
   ofstream out;
   std::streambuf * old_buffer = 0;

   if(file) {
      out.open(file);
      if (!out.is_open()) {
         std::cerr<<"Couldn't open file "<<file<<std::endl;
         return;
      }
      old_buffer = std::cout.rdbuf(out.rdbuf());
   }

   // model summary
   int nmodels = fModelContainer->size();
   std::cout<<std::endl
            <<"==========================================="<<std::endl
            <<" Model Comparison Summary"<<std::endl
            <<"==========================================="<<std::endl
            <<std::endl
            <<" Number of models               : "<<nmodels<<std::endl
            <<std::endl;

   // probability summary
   std::cout<<" - A priori probabilities:"<<std::endl<<std::endl;
  
   for (int i=0; i<nmodels; i++)
      std::cout<<"     p("<< fModelContainer->at(i)->GetName()
               <<") = "<< fModelContainer->at(i)->GetModelAPrioriProbability()
               <<std::endl;
   std::cout<<std::endl;

   std::cout<<" - A posteriori probabilities:"<<std::endl<<std::endl;

   for (int i = 0; i < nmodels; i++)
      std::cout<<"     p("<< fModelContainer->at(i)->GetName()
               <<" | data) = "<< fModelContainer->at(i)->GetModelAPosterioriProbability()
               <<std::endl;
   std::cout<<std::endl;

   // Bayes factors summary
   std::cout<<" - Bayes factors:"<<std::endl<<std::endl;
   for (int i = 0; i < nmodels-1; i++)
      for (int j = i+1; j < nmodels; j++)
         std::cout<<"     K = p(data | "<<fModelContainer->at(i)->GetName()<<") / "
                  <<"p(data | "<<fModelContainer->at(j)->GetName()<<") = "
                  <<BayesFactor(i,j)<<std::endl;
   std::cout<<std::endl;

   // p-values summary
   std::cout
      <<" - p-values:"<<std::endl
      <<std::endl;

   for (int i = 0; i < nmodels; i++)
   {
      double p = fModelContainer->at(i)->GetPValue();
      std::cout <<"     "<< fModelContainer->at(i)->GetName();
      if(p>=0.)
         std::cout<<":  p-value = "<< p;
      else
         std::cout<<":  p-value not calculated";
      std::cout<<std::endl;
   }
   std::cout<<std::endl;

   std::cout<<"==========================================="<<std::endl<<std::endl;

   if (file)
      std::cout.rdbuf(old_buffer);

}

// ---------------------------------------------------------

void BCModelManager::PrintResults()
{
   // print summary of all models
   for (unsigned int i = 0; i < GetNModels(); i++)
      GetModel(i)->PrintResults(Form("%s.txt", GetModel(i)->GetName().data()));
}

// ---------------------------------------------------------

void BCModelManager::Copy(BCModelManager & modelmanager) const
{
   // don't copy the content only the pointers
   modelmanager.fModelContainer = fModelContainer;
   modelmanager.fDataSet        = fDataSet;
}

// ---------------------------------------------------------
