/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TPostScript.h>
#include <TCanvas.h>

#include <math.h>
#include <iostream>

#include "../../BAT/BCLog.h"
#include "../../BAT/BCH1D.h"

#include "BCTemplateFitter.h"
#include "BCTemplateEnsembleTest.h"

// ---------------------------------------------------------
BCTemplateEnsembleTest::BCTemplateEnsembleTest()
   : fTemplateFitter(0)
   , fFile(0)
   , fTree(0)
   , fEnsembleCounter(0)
   , fEnsembleExpectation(0)
   , fNEnsembles(0)
   , fFlagMCMC(false)
{
   fRandom = new TRandom3(0);

   BCLog::OutWarning("Class BCTemplateFitter is depreceted and it will be removed");
   BCLog::OutWarning("in the future version of BAT. The same is true for the class");
   BCLog::OutWarning("BCTemplateEnsembleTest. Use the new class BCMTF and its ensemble");
   BCLog::OutWarning("testing capabilities instead.");
}

// ---------------------------------------------------------
BCTemplateEnsembleTest::~BCTemplateEnsembleTest()
{
   if (fRandom)
      delete fRandom;
}

// ---------------------------------------------------------
int BCTemplateEnsembleTest::SetEnsembleTemplate(TH1D hist)
{
   // calculate integral
   double integral = hist.Integral();

   // check if integral is ok
   if (integral <= 0) {
      std::cout << "Template not valid. Integral is lower or equal to 0." << std::endl;
      return 0;
   }

   // set template
   fEnsembleTemplate = hist;

   // scale template
   fEnsembleTemplate.Scale(1.0/integral);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCTemplateEnsembleTest::PerformEnsembleTest(TTree* tree)
{
   // set log level to nothing
   BCLog::LogLevel ll = BCLog::GetLogLevelScreen();
   BCLog::SetLogLevel(BCLog::nothing);

   // initialize template fitter
   fTemplateFitter->Initialize();

   // Prepare tree
   PrepareTree();

   // get number of parameters
   int npar = fTemplateFitter->GetNParameters();

   // connect tree
   if (tree) {
      fTemplateParameters = std::vector<double>(npar);
      for (int i = 0; i < npar; ++i) {
         tree->SetBranchAddress(Form("Parameter%i", i), &(fTemplateParameters[i]));
      }
   }

   // loop over ensembles
   for(int j = 0; j < fNEnsembles; j++){

      // print status
      if ((j+1) % 100 == 0 && j > 0)
         std::cout << "Fraction of ensembles analyzed: " << double(j+1) / double(fNEnsembles) * 100 << "%" << std::endl;

      // get parameters from tree
      if (tree) {
        int index = (int) fRandom->Uniform(tree->GetEntries());
        tree->GetEntry(index);
      }

      // create new ensemble
      TH1D* ensemble = BuildEnsemble();

      // set ensemble as new data set
      fTemplateFitter->SetData(*ensemble);

      // find mode
      fTemplateFitter->FindMode();

      // perform MCMC
      if(fFlagMCMC) {
         fTemplateFitter->MCMCInitialize();

         fTemplateFitter->MarginalizeAll();

         // find mode with MCMC best fit
         fTemplateFitter->FindMode(fTemplateFitter->GetBestFitParameters());

         // loop over parameters and set tree variables
         for (int i = 0; i < npar; ++i) {
            BCH1D* hist = fTemplateFitter->GetMarginalized(fTemplateFitter->GetParameter(i));
            fOutParModeMarg[i]       = hist->GetMode();
            fOutParMedianMarg[i]     = hist->GetMedian();
            fOutParMeanMarg[i]       = hist->GetMean();
            fOutParRMSMarg[i]        = hist->GetRMS();
            fOutParErrorUpMarg[i]    = hist->GetQuantile(0.84)-hist->GetMode();
            fOutParErrorDownMarg[i]  = hist->GetMode()-hist->GetQuantile(0.16);
            fOutParQuantile5Marg[i]  = hist->GetQuantile(0.05);
            fOutParQuantile10Marg[i] = hist->GetQuantile(0.10);
            fOutParQuantile90Marg[i] = hist->GetQuantile(0.90);
            fOutParQuantile95Marg[i] = hist->GetQuantile(0.95);
            fOutParPullMarg[i]       = -1;
            double error = 0.5 * (fOutParErrorUpMarg.at(i) + fOutParErrorDownMarg.at(i));
            if (error > 0)
               fOutParPullMarg[i]     = (fOutParModeMarg.at(i) - fTemplateParameters.at(i)) / ( error );
         }
         fOutChi2Marg             = fTemplateFitter->CalculateChi2( fTemplateFitter->GetBestFitParametersMarginalized() );
         fOutChi2ProbMarg         = fTemplateFitter->CalculateChi2Prob(fTemplateFitter->GetBestFitParametersMarginalized());
      }

      if (fFlagMCMC) {
         int nratios = fTemplateFitter->GetNRatios();
         for (int i = 0; i < nratios; ++i) {
            TH1D histtemp = fTemplateFitter->GetHistRatio1D(i);
             BCH1D * hist = new BCH1D( &histtemp );
             fOutRatioModeMarg[i]       = hist->GetMode();
             fOutRatioMedianMarg[i]     = hist->GetMedian();
             fOutRatioMeanMarg[i]       = hist->GetMean();
             fOutRatioRMSMarg[i]        = hist->GetRMS();
             fOutRatioErrorUpMarg[i]    = hist->GetQuantile(0.84)-hist->GetMode();
            fOutRatioErrorDownMarg[i]  = hist->GetMode()-hist->GetQuantile(0.16);
             fOutRatioQuantile5Marg[i]  = hist->GetQuantile(0.05);
             fOutRatioQuantile10Marg[i] = hist->GetQuantile(0.10);
             fOutRatioQuantile90Marg[i] = hist->GetQuantile(0.90);
             fOutRatioQuantile95Marg[i] = hist->GetQuantile(0.95);
         }
      }

      // set tree variables
      fOutParValue           = fTemplateParameters;
      fOutParModeGlobal      = fTemplateFitter->GetBestFitParameters();
      fOutParErrorUpGlobal   = fTemplateFitter->GetBestFitParameterErrors();
      fOutParErrorDownGlobal = fTemplateFitter->GetBestFitParameterErrors();
      fOutChi2Global         = fTemplateFitter->CalculateChi2( fTemplateFitter->GetBestFitParameters() );
      fOutNDF                = fTemplateFitter->GetNDF();
      fOutChi2ProbGlobal     = fTemplateFitter->CalculateChi2Prob(fTemplateFitter->GetBestFitParameters());
      fOutKSProb             = fTemplateFitter->CalculateKSProb();
      fOutPValue             = fTemplateFitter->CalculatePValue();
      fOutNEvents            = int(fTemplateFitter->GetData().Integral());

      for (int i = 0; i < npar; ++i) {
         fOutParPullGlobal[i] = -1;
         double error = 0.5 * (fOutParErrorUpGlobal.at(i) + fOutParErrorDownGlobal.at(i));
         if (error > 0)
            fOutParPullGlobal[i] = (fOutParModeGlobal.at(i) - fTemplateParameters.at(i)) / ( error );
      }

      // fill the tree
      fTree->Fill();
   }

   // reset log level
   BCLog::SetLogLevel(ll);

   // no error
   return 1;
}

//---------------------------------------------------------------------------------------------------------
TH1D* BCTemplateEnsembleTest::BuildEnsemble()
{
   // get histogram parameters
  int nbins   = fTemplateFitter->GetData().GetNbinsX();

   // create new ensemble
   TH1D* ensemble = new TH1D(fTemplateFitter->GetData());

   // increase ensemble counter
   fEnsembleCounter++;

   // get new parameter if needed
   std::vector<double> parameters = fTemplateParameters;

   // loop over bins and fill them
  for(int ibin = 1; ibin <= nbins; ++ibin){
      double nexp = fTemplateFitter->Expectation(ibin, parameters);
      double nobs = gRandom->Poisson(nexp);

      // set the bin content
    ensemble->SetBinContent(ibin, nobs);
  }

   // return the ensemble histogram
  return ensemble;
}

//---------------------------------------------------------------------------------------------------------
int BCTemplateEnsembleTest::Write(const char * filename)
{
   // open file
   fFile = new TFile(filename, "RECREATE");

   // write tree
   fTree->Write();

   // close file
   fFile->Close();

   // free memory
   delete fFile;

   // no error
   return 1;
}

//---------------------------------------------------------------------------------------------------------
int BCTemplateEnsembleTest::PrepareTree()
{
   // delete old tree if necessary
   if (fTree)
      delete fTree;

   // create new tree
   fTree = new TTree("fTree", "fTree");

   // get number of parameters and ratios
   int npar = fTemplateFitter->GetNParameters();
   int nratios = fTemplateFitter->GetNRatios();

   // initialize variables
   fOutParValue.assign(npar, 0);
   fOutParModeGlobal.assign(npar, 0);
   fOutParErrorUpGlobal.assign(npar, 0);
   fOutParErrorDownGlobal.assign(npar, 0);
   fOutParPullGlobal.assign(npar, 0);
   fOutParModeMarg.assign(npar, 0);
   fOutParMeanMarg.assign(npar, 0);
   fOutParMedianMarg.assign(npar, 0);
   fOutParRMSMarg.assign(npar, 0);
   fOutParErrorUpMarg.assign(npar, 0);
   fOutParErrorDownMarg.assign(npar, 0);
   fOutParPullMarg.assign(npar, 0);
   fOutParQuantile5Marg.assign(npar, 0);
   fOutParQuantile10Marg.assign(npar, 0);
   fOutParQuantile90Marg.assign(npar, 0);
   fOutParQuantile95Marg.assign(npar, 0);
   fOutRatioModeMarg.assign(nratios, 0);
   fOutRatioMeanMarg.assign(nratios, 0);
   fOutRatioMedianMarg.assign(nratios, 0);
   fOutRatioRMSMarg.assign(nratios, 0);
   fOutRatioErrorUpMarg.assign(nratios, 0);
   fOutRatioErrorDownMarg.assign(nratios, 0);
   fOutRatioQuantile5Marg.assign(nratios, 0);
   fOutRatioQuantile10Marg.assign(nratios, 0);
   fOutRatioQuantile90Marg.assign(nratios, 0);
   fOutRatioQuantile95Marg.assign(nratios, 0);

   fTree->Branch("chi2_global", &fOutChi2Global,     "chi2 (global)/D");
   fTree->Branch("ndf",      &fOutNDF,      "ndf/I");
   fTree->Branch("chi2prob_global", &fOutChi2ProbGlobal, "chi2 prob probability (global)/D");
   fTree->Branch("KSprob",   &fOutKSProb,   "KS probability/D");
   fTree->Branch("pvalue",   &fOutPValue,   "p-value/D");
   fTree->Branch("nevents",  &fOutNEvents,  "n events/I");

   for (int i = 0; i < npar; ++i) {
      // add branches
      fTree->Branch(Form("parameter_%i", i),                 &fOutParValue[i],      Form("parameter_%i/D", i));
      fTree->Branch(Form("par_global_mode_par_%i", i),       &fOutParModeGlobal[i],      Form("par_global_Mode_par_%i/D", i));
      fTree->Branch(Form("par_global_error_up_par_%i", i),   &fOutParErrorUpGlobal[i],   Form("par_global_error_up_par_%i/D", i));
      fTree->Branch(Form("par_global_error_down_par_%i", i), &fOutParErrorDownGlobal[i], Form("par_global_error_down_par_%i/D", i));
      fTree->Branch(Form("par_global_pull_par_%i", i),       &fOutParPullGlobal[i],      Form("par_global_pull_par_%i/D", i));

      if(fFlagMCMC) {
         fTree->Branch(Form("par_marg_mode_par_%i", i),       &fOutParModeMarg[i],       Form("par_marg_mode_par_%i/D", i));
         fTree->Branch(Form("par_marg_mean_par_%i", i),       &fOutParMeanMarg[i],       Form("par_marg_mean_par_%i/D", i));
         fTree->Branch(Form("par_marg_median_par_%i", i),     &fOutParMedianMarg[i],     Form("par_marg_median_par_%i/D", i));
         fTree->Branch(Form("par_marg_rms_par_%i", i),        &fOutParRMSMarg[i],        Form("par_marg_rms_par_%i/D", i));
         fTree->Branch(Form("par_marg_error_up_par_%i", i),   &fOutParErrorUpMarg[i],    Form("par_marg_ErrorUp_par_%i/D", i));
         fTree->Branch(Form("par_marg_error_down_par_%i", i), &fOutParErrorDownMarg[i],  Form("par_marg_error_down_par_%i/D", i));
         fTree->Branch(Form("par_marg_pull_par_%i", i),       &fOutParPullMarg[i],       Form("par_marg_pull_par_%i/D", i));
         fTree->Branch(Form("par_marg_quantile5_par_%i", i),  &fOutParQuantile5Marg[i],  Form("par_marg_Quantile5_par_%i/D", i));
         fTree->Branch(Form("par_marg_quantile10_par_%i", i), &fOutParQuantile10Marg[i], Form("par_marg_Quantile10_par_%i/D", i));
         fTree->Branch(Form("par_marg_quantile90_par_%i", i), &fOutParQuantile90Marg[i], Form("par_marg_Quantile90_par_%i/D", i));
         fTree->Branch(Form("par_marg_quantile95_par_%i", i), &fOutParQuantile95Marg[i], Form("par_marg_Quantile95_par_%i/D", i));
      }
   }

   if (fFlagMCMC) {
      fTree->Branch("chi2_marg", &fOutChi2Marg,     "chi2 (marginalized)/D");
      fTree->Branch("chi2prob_marg", &fOutChi2ProbMarg, "chi2 prob probability (marginalized)/D");

      for (int i = 0; i < nratios; ++i) {
         fTree->Branch(Form("ratio_marg_mode_ratio_%i", i),       &fOutRatioModeMarg[i],       Form("ratio_marg_mode_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_mean_ratio_%i", i),       &fOutRatioMeanMarg[i],       Form("ratio_marg_mean_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_median_ratio_%i", i),     &fOutRatioMedianMarg[i],     Form("ratio_marg_median_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_rms_ratio_%i", i),        &fOutRatioRMSMarg[i],        Form("ratio_marg_rms_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_error_up_ratio_%i", i),   &fOutRatioErrorUpMarg[i],    Form("ratio_marg_ErrorUp_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_error_down_ratio_%i", i), &fOutRatioErrorDownMarg[i],  Form("ratio_marg_error_down_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_quantile5_ratio_%i", i),  &fOutRatioQuantile5Marg[i],  Form("ratio_marg_Quantile5_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_quantile10_ratio_%i", i), &fOutRatioQuantile10Marg[i], Form("ratio_marg_Quantile10_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_quantile90_ratio_%i", i), &fOutRatioQuantile90Marg[i], Form("ratio_marg_Quantile90_ratio_%i/D", i));
         fTree->Branch(Form("ratio_marg_quantile95_ratio_%i", i), &fOutRatioQuantile95Marg[i], Form("ratio_marg_Quantile95_ratio_%i/D", i));
      }
   }

   // no error
   return 1;
}

//---------------------------------------------------------------------------------------------------------
void BCTemplateEnsembleTest::PrintPulls(const char* filename)
{
   // create postscript
   TPostScript * ps = new TPostScript(filename);

   // create canvas and prepare postscript
   TCanvas * canvas = new TCanvas();

   canvas->Update();
   ps->NewPage();
   canvas->cd();

   // get number of parameters
   int npar = fTemplateFitter->GetNParameters();

   // create histogram container
   std::vector<TH1D*> histcontainer = std::vector<TH1D*>(0);

   // loop over all parameters
   for (int j = 0; j < npar; ++j) {
      TH1D* hist = new TH1D("", "Pull;N", 11, -5.5, 5.5);
      histcontainer.push_back(hist);
   }

   // loop over all ensembles
   for (int i = 0; i < fNEnsembles; ++i) {

      // get ensemble
      fTree->GetEntry(i);

      // loop over all parameters
      for (int j = 0; j < npar; ++j) {
         histcontainer[j]->Fill(fOutParPullGlobal[j]);
      }
   }

   // loop over all parameters
   for (int j = 0; j < npar; ++j) {
      // update post script
      canvas->Update();
      ps->NewPage();
      canvas->cd();

      canvas->cd();
      histcontainer.at(j)->Draw();
   }
   // close ps
   canvas->Update();
   ps->Close();

   // free memory
   delete ps;
   delete canvas;

}

//---------------------------------------------------------------------------------------------------------
