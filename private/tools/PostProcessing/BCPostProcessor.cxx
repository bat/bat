/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPostProcessor.h"

#include <BAT/BCLog.h>

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

// ---------------------------------------------------------
BCPostProcessor::BCPostProcessor() : fFile(0)
                                   , fTrees(std::vector<TTree*>(0))
                                   , fNTrees(0)
                                   , fNParameters(0)
                                   , fNSamplesPreRun(0)
                                   , fNSamplesMainRun(0)
                                   , fParameters(std::vector<double>(0))
{
}

// ---------------------------------------------------------
BCPostProcessor::~BCPostProcessor()
{
  CloseRootFile();
}

// ---------------------------------------------------------
int BCPostProcessor::OpenRootFile(std::string filename)
{
  // try to open file
  fFile = new TFile(filename.c_str());

  // check if file is open
  if (!fFile->IsOpen()) {
    BCLog::OutWarning("BCPostProcessor::OpenRootFile. Could not open ROOT file.");

    return 0;
  }

  // clear vector of tress
  fTrees.clear();

  // reset number of trees
  fNTrees = 0;

  // get trees from file
  TTree* tree;
  while ( ( tree = (TTree*) fFile->Get(Form("MarkovChainTree_%i", fNTrees)) ) ) {

    // add tree to vector and increase number of trees
    fTrees.push_back(tree);
    fNTrees++;
  }

  // check if there is at least on chain
  if (fNTrees < 1) {
    BCLog::OutWarning("BCPostProcessor::OpenFile. No chain.");

    return 0;
  }

  // calculate number of parameters
  // reset number of parameters
  fNParameters = 0;

  while ( ( fTrees.at(0)->FindBranch(Form("Parameter%i", fNParameters))) ) {
    fNParameters++;
  }

  // calculate number of samples
  fNSamplesPreRun = fTrees.at(0)->GetEntries("Phase==1");
  fNSamplesMainRun = fTrees.at(0)->GetEntries("Phase==2");

  CalculateMinMax();

  // no error
  return 1;
}

// ---------------------------------------------------------
void BCPostProcessor::CloseRootFile()
{
  fFile->Close();
}

// ---------------------------------------------------------
void BCPostProcessor::PrintInfo()
{
  std::cout << "Number of parameters              : " << fNParameters << std::endl;
  std::cout << "Number of trees                   : " << fNTrees << std::endl;
  std::cout << "Number of samples in the pre-run  : " << fNSamplesPreRun << std::endl;
  std::cout << "Number of samples in the main run : " << fNSamplesMainRun << std::endl;
  std::cout << std::endl;
  std::cout << "Minimum and maximum values of the parameters: " << std::endl;
  for (int i = 0; i < fNParameters; ++i) {
    std::cout << " Parameter " << i << " : " << fParametersMin[i] << " - " << fParametersMax[i] << std::endl;
  }
  std::cout << "Minimum and maximum values of the log probability: " << std::endl;
  std::cout << " Log(probability) : " << fLogProbabilityMin << " - " << fLogProbabilityMax << std::endl;

}

// ---------------------------------------------------------
double BCPostProcessor::GetValue(std::string branchname, int chainindex, int entry, bool prerun)
{
  //debugKK
  //  std::cout << "GetValue. " << branchname.c_str() << "  " << chainindex << " " << entry << " " << prerun << std::endl;

  // check chain index
  if (chainindex < 0 || chainindex >= fNTrees) {
    BCLog::OutWarning("BCPostProcessor::GetValue. Chain index not within range.");
    return -1;
  }

  // check entry
  if (entry < 0 || (entry >= fNSamplesPreRun && prerun) || (entry >= fNSamplesMainRun && !prerun)) {
    BCLog::OutWarning("BCPostProcessor::GetValue. Entry not within range.");
    return -1;
  }

  // move to main run if not pre run
  if (!prerun)
    entry += fNSamplesPreRun;

  TTree* tree = fTrees.at(chainindex);

  double value = 0;

  tree->SetBranchAddress(branchname.c_str(), &value);

  tree->GetEntry(entry);

  return value;
}

// ---------------------------------------------------------
double BCPostProcessor::GetParameterValue(int parindex, int chainindex, int entry, bool prerun)
{
  // debugKK
  //  std::cout << "GetParameterValue. " << parindex << " " << fNParameters << std::endl;

  // check parameter index
  if (parindex < 0 || parindex >= fNParameters) {
    BCLog::OutWarning("BCPostProcessor::GetValue. Parameter index not within range.");
    return -1;
  }

  return GetValue(Form("Parameter%i", parindex), chainindex, entry, prerun);
}

// ---------------------------------------------------------
double BCPostProcessor::GetLogProbabilityValue(int chainindex, int entry, bool prerun)
{
  return GetValue("LogProbability", chainindex, entry, prerun);
}

// ---------------------------------------------------------
void BCPostProcessor::CalculateMinMax()
{
  // reset minimum and maximum values
  fParametersMin.clear();
  fParametersMax.clear();
  fLogProbabilityMin = 0;
  fLogProbabilityMax = 0;

  // find minimum and maximum values for parameters
  for (int i = 0; i < fNParameters; ++i) {
    for (int j = 0; j < fNTrees; ++j) {
      fTrees.at(0)->Draw(Form("Parameter%i>>temphist", i), "Phase==2");
      TH1D* temphist = (TH1D*) gDirectory->Get("temphist");
      double xmin = temphist->GetXaxis()->GetXmin();
      double xmax = temphist->GetXaxis()->GetXmax();
      delete temphist;

      if (j == 0) {
        fParametersMin.push_back(xmin);
        fParametersMax.push_back(xmax);
      }
      if (xmin < fParametersMin[i]) {
        fParametersMin[i] = xmin;
      }
      if (xmax > fParametersMax[i]) {
        fParametersMax[i] = xmax;
      }
    }
  }

  // find minimum and maximum values for log probability
  for (int j = 0; j < fNTrees; ++j) {
    fTrees.at(0)->Draw("LogProbability>>temphist", "Phase==2");
    TH1D* temphist = (TH1D*) gDirectory->Get("temphist");
    double xmin = temphist->GetXaxis()->GetXmin();
    double xmax = temphist->GetXaxis()->GetXmax();
    delete temphist;

    if (j == 0) {
      fLogProbabilityMin = xmin;
      fLogProbabilityMax = xmax;
    }
    if (xmin <fLogProbabilityMin ) {
      fLogProbabilityMin = xmin;
    }
    if (xmax > fLogProbabilityMin) {
        fLogProbabilityMax = xmax;
    }
  }

  return;
}


// ---------------------------------------------------------
