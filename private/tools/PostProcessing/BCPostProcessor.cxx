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
}

// ---------------------------------------------------------
