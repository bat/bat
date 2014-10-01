/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModelOutput.h"

#include "BCModel.h"
#include "BCParameter.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>

#include <iostream>

// ---------------------------------------------------------
BCModelOutput::BCModelOutput()
{
   Init();
}

// ---------------------------------------------------------
BCModelOutput::BCModelOutput(BCModel * model, const char * filename)
{
   Init();
   SetModel(model);
   SetFile(filename);
}

// ---------------------------------------------------------
BCModelOutput::~BCModelOutput()
{
  if (fOutputFile) {
    Close();
    delete fOutputFile;
  }
}

// ---------------------------------------------------------
BCModelOutput::BCModelOutput(const BCModelOutput & modeloutput)
{
   modeloutput.Copy(* this);
}

// ---------------------------------------------------------
BCModelOutput & BCModelOutput::operator = (const BCModelOutput & modeloutput)
{
   if (this != &modeloutput)
      modeloutput.Copy(* this);

   return * this;
}

// ---------------------------------------------------------
void BCModelOutput::Init()
{
   fIndex = 0;
   fOutputFile = 0;
   fAnalysisTree = 0;
   fTreeSA = 0;
   fModel = 0;
   fFileName = 0;
}

// ---------------------------------------------------------
void BCModelOutput::SetModel(BCModel * model)
{
   fModel = model;
   fModel->MCMCInitialize();
   fModel->SAInitialize();
}

// ---------------------------------------------------------
void BCModelOutput::SetFile(const char * filename)
{
   if(!fModel) {
      BCLog::OutError("BCModelOutput::SetFile : Cannot set file if model is not set.");
      return;
   }

   // delete the old file
   if (fOutputFile) {
      fOutputFile->Close();
      delete fOutputFile;
   }

   // remember current directory
   TDirectory * dir = gDirectory;

   // create a new file
   fFileName = const_cast<char *>(filename);
   fOutputFile = TFile::Open(fFileName, "RECREATE");

   // initialize trees
   InitializeAnalysisTree();
   fModel->MCMCInitializeMarkovChainTrees();
   fModel->InitializeSATree();

   // change back to the old directory
   gDirectory = dir;
}

// ---------------------------------------------------------
void BCModelOutput::WriteMarkovChain(bool flag)
{
   if (fModel)
      fModel->WriteMarkovChain(flag);
}

// ---------------------------------------------------------
void BCModelOutput::FillAnalysisTree()
{
   if(!fOutputFile) {
      BCLog::OutError("BCModelOutput::FillAnalysisTree : No file to write to.");
      return;
   }

   // get output values from model
   fNParameters = fModel->GetNParameters();
   fProbability_apriori   = fModel->GetModelAPrioriProbability();
   fProbability_aposteriori = fModel->GetModelAPosterioriProbability();

   fMode_global.clear();
   fMode_marginalized.clear();
   fMean_marginalized.clear();
   fMedian_marginalized.clear();
   fQuantile_05.clear();
   fQuantile_10.clear();
   fQuantile_16.clear();
   fQuantile_84.clear();
   fQuantile_90.clear();
   fQuantile_95.clear();

   // loop over parameters
   int nparameters = fModel->GetNParameters();
   for (int i = 0; i < nparameters; ++i) {
     const BCParameter * parameter = fModel->GetParameter(i);
     if (fModel->GetBestFitParameters().size() > 0)
       fMode_global.push_back(fModel->GetBestFitParameters().at(i));
     if (fModel->GetMarginalized(parameter->GetName().data())) {
       fMode_marginalized.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetMode());
       fMean_marginalized.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetMean());
       fMedian_marginalized.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetMedian());
       fQuantile_05.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.05));
       fQuantile_10.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.10));
       fQuantile_16.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.16));
       fQuantile_84.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.84));
       fQuantile_90.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.90));
       fQuantile_95.push_back(fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.95));
     }
   }

   // fill tree
   fAnalysisTree->Fill();

   // increase index
   fIndex++;
}

// ---------------------------------------------------------
void BCModelOutput::WriteMarginalizedDistributions()
{
   if(!fOutputFile) {
      BCLog::OutError("BCModelOutput::WriteMarginalizedDistributions : No file to write to.");
      return;
   }

   if (!fOutputFile->IsOpen())
     return;

   // remember current directory
   TDirectory * dir = gDirectory;

   // change to file
   fOutputFile->cd();

   int nparameters = fModel->GetNParameters();
   for (int i = 0; i < nparameters; ++i) {
      BCH1D* bchist = fModel->GetMarginalized(fModel->GetParameter(i));
      if (bchist) {
         TH1D* hist = bchist->GetHistogram();
         if (hist)
            hist->Write();
      }
   }

   if (nparameters > 1)
      for (int i = 0; i < nparameters - 1; ++i) {
         for (int j = i + 1; j < nparameters; ++j) {
            BCH2D* bchist = fModel->GetMarginalized(fModel->GetParameter(i),
                                                    fModel->GetParameter(j));
            if (bchist) {
               TH2D* hist = bchist->GetHistogram();
               if (hist)
                  hist->Write();
            }
         }
      }

   // return to old directory
   gDirectory = dir;
}

// ---------------------------------------------------------
void BCModelOutput::Write(TObject * o)
{
   if(!fOutputFile) {
      BCLog::OutError("BCModelOutput::Write : No file to write to.");
      return;
   }

   // remember current directory
   TDirectory * dir = gDirectory;

   // change to file
   fOutputFile->cd();

   o->Write();

   // return to old directory
   gDirectory = dir;
}

// ---------------------------------------------------------
void BCModelOutput::Close()
{
  if (!fOutputFile)
    return;

  if (!fOutputFile->IsOpen())
    return;

   // remember current directory
   TDirectory * dir = gDirectory;

   // change to file
   fOutputFile->cd();

   // write analysis tree to file
   if (fAnalysisTree)
     if (fAnalysisTree->GetEntries() > 0)
       fAnalysisTree->Write();

   // write markov chain tree to file
   if (fModel) {
     for (unsigned i = 0; i < fModel->MCMCGetNChains(); ++i) {
       if (fModel->MCMCGetMarkovChainTree(i)){
         if (fModel->MCMCGetMarkovChainTree(i)->GetEntries() > 0) {
           fModel->MCMCGetMarkovChainTree(i)->Write();
         }
       }
     }
   }

   // write SA tree to file
   if (fModel)
     if (fModel->GetSATree()->GetEntries() > 0)
       fModel->GetSATree()->Write();

   // close file
   fOutputFile->Close();

   // return to old directory
   gDirectory = dir;
}

// ---------------------------------------------------------
void BCModelOutput::InitializeAnalysisTree()
{
   // create new tree
   fAnalysisTree = new TTree("AnalysisTree", "AnalysisTree");

   // set branch addresses
   fAnalysisTree->Branch("fIndex",                   &fIndex,                   "index/I");
   fAnalysisTree->Branch("fNParameters",             &fNParameters,             "parameters/I");
   fAnalysisTree->Branch("fProbability_apriori" ,    &fProbability_apriori,     "apriori probability/D");
   fAnalysisTree->Branch("fProbability_aposteriori", &fProbability_aposteriori, "aposteriori probability/D");
   fAnalysisTree->Branch("fMode_global",             &fMode_global);
   fAnalysisTree->Branch("fMode_marginalized",       &fMode_marginalized);
   fAnalysisTree->Branch("fMean_marginalized",       &fMean_marginalized);
   fAnalysisTree->Branch("fMedian_marginalized",     &fMedian_marginalized);
   fAnalysisTree->Branch("fQuantile_5" ,             &fQuantile_05);
   fAnalysisTree->Branch("fQuantile_10" ,            &fQuantile_10);
   fAnalysisTree->Branch("fQuantile_16" ,            &fQuantile_16);
   fAnalysisTree->Branch("fQuantile_84" ,            &fQuantile_84);
   fAnalysisTree->Branch("fQuantile_90" ,            &fQuantile_90);
   fAnalysisTree->Branch("fQuantile_95" ,            &fQuantile_95);
}

// ---------------------------------------------------------
void BCModelOutput::Copy(BCModelOutput & modeloutput) const
{
   // don't copy the content
   modeloutput.fModel            = fModel;
   modeloutput.fAnalysisTree     = fAnalysisTree;
}
