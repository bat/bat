/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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
      fOutputFile->Close();
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
   fOutputFile = new TFile(fFileName, "RECREATE");

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

   // loop over parameters
   int nparameters = fModel->GetNParameters();
   for (int i = 0; i < nparameters; ++i) {
      BCParameter * parameter = fModel->GetParameter(i);
      if (fModel->GetBestFitParameters().size() > 0)
         fMode_global[i] = fModel->GetBestFitParameters().at(i);
      if (fModel->GetMarginalized(parameter->GetName().data())) {
         fMode_marginalized[i] = fModel->GetMarginalized(parameter->GetName().data())->GetMode();
         fMean_marginalized[i] = fModel->GetMarginalized(parameter->GetName().data())->GetMean();
         fMedian_marginalized[i] = fModel->GetMarginalized(parameter->GetName().data())->GetMedian();
         fQuantile_05[i] = fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.05);
         fQuantile_10[i] = fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.10);
         fQuantile_16[i] = fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.16);
         fQuantile_84[i] = fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.84);
         fQuantile_90[i] = fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.90);
         fQuantile_95[i] = fModel->GetMarginalized(parameter->GetName().data())->GetQuantile(0.95);
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

   // remember current directory
   TDirectory * dir = gDirectory;

   // change to file
   fOutputFile->cd();

   int nparameters = fModel->GetNParameters();
   for (int i = 0; i < nparameters; ++i)
      fModel->GetMarginalized(fModel->GetParameter(i))->GetHistogram()->Write();

   if (nparameters > 1)
      for (int i = 0; i < nparameters - 1; ++i)
         for (int j = i + 1; j < nparameters; ++j)
            fModel->GetMarginalized(fModel->GetParameter(i),
                  fModel->GetParameter(j))->GetHistogram()->Write();

   // return to old directory
   gDirectory = dir;
}

// ---------------------------------------------------------
void BCModelOutput::WriteErrorBand()
{
   if(!fOutputFile) {
      BCLog::OutError("BCModelOutput::WriteErrorBand : No file to write to.");
      return;
   }

   // remember current directory
   TDirectory * dir = gDirectory;

   // change to file
   fOutputFile->cd();

   TH2D * h0 = fModel->GetErrorBandXY();
   if (h0) {
      TH2D * h1 = (TH2D*)h0->Clone("errorbandxy");
      h1->Write();

      double levels[] = { .68, .90, .95 };
      int nlevels = sizeof(levels)/sizeof(double);
      for (int i=0;i<nlevels;i++) {
         TH2D * htmp = fModel->GetErrorBandXY_yellow(levels[i]);
         htmp->SetName(TString::Format("%s_sub_%f.2",h1->GetName(),levels[i]));
         htmp->Write();
         delete htmp;
      }

      delete h1;
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
   // remember current directory
   TDirectory * dir = gDirectory;

   // change to file
   fOutputFile->cd();

   // write analysis tree to file
   if (fAnalysisTree->GetEntries() > 0)
      fAnalysisTree->Write();

   // write markov chain tree to file
   for (int i = 0; i < fModel->MCMCGetNChains(); ++i)
      if (fModel->MCMCGetMarkovChainTree(i)->GetEntries() > 0)
         fModel->MCMCGetMarkovChainTree(i)->Write();

   // write SA tree to file
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
   fAnalysisTree->Branch("fMode_global",              fMode_global,             "mode (global) [parameters]/D");
   fAnalysisTree->Branch("fMode_marginalized",        fMode_marginalized,       "mode (marginalized) [parameters]/D");
   fAnalysisTree->Branch("fMean_marginalized",        fMean_marginalized,       "mean (marginalized)[parameters]/D");
   fAnalysisTree->Branch("fMedian_marginalized",      fMedian_marginalized,     "median (marginalized)[parameters]/D");
   fAnalysisTree->Branch("fQuantile_05" ,             fQuantile_05,             "quantile 5% [parameters]/D");
   fAnalysisTree->Branch("fQuantile_10" ,             fQuantile_10,             "quantile 10% [parameters]/D");
   fAnalysisTree->Branch("fQuantile_16" ,             fQuantile_16,             "quantile 16% [parameters]/D");
   fAnalysisTree->Branch("fQuantile_84" ,             fQuantile_84,             "quantile 84% [parameters]/D");
   fAnalysisTree->Branch("fQuantile_90" ,             fQuantile_90,             "quantile 90% [parameters]/D");
   fAnalysisTree->Branch("fQuantile_95" ,             fQuantile_95,             "quantile 95% [parameters]/D");
}

// ---------------------------------------------------------
void BCModelOutput::Copy(BCModelOutput & modeloutput) const
{
   // don't copy the content
   modeloutput.fModel            = fModel;
   modeloutput.fAnalysisTree     = fAnalysisTree;
//   modeloutput.fMarkovChainTrees = fMarkovChainTrees;
}

// ---------------------------------------------------------
