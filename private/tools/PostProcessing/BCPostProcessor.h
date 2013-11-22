#ifndef __BCPOSTPROCESSOR__H
#define __BCPOSTPROCESSOR__H

/*!
 * \class BCPostProcessor
 * \brief A base class for post processing of MCMC samples
 * \author Kevin Kr&ouml;ninger
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <string>

#include <TTree.h>

class TFile;

// ---------------------------------------------------------

class BCPostProcessor
{
 public:

  /** \name Constructors and destructors */
  /** @{ */

  /**
   * The default constructor. */
  BCPostProcessor();

  /**
   * The default destructor. */
  ~BCPostProcessor();

  /** @} */

  /** \name Member functions (get) */
  /** @{ */

  /**
   * Return the number of chains. */
  int GetNChains()
  { return fNTrees; };

  /**
   * Return the number of parameters. */
  int GetNParameters()
  { return fNParameters; };

  /**
   * Return the number of samples in the pre-run. */
  int GetNSamplesPreRun()
  { return fNSamplesPreRun; };

  /**
   * Return the number of samples in the main run. */
  int GetNSamplesMainRun()
  { return fNSamplesMainRun; };

  /** @} */

  /** \name Member functions (set) */
  /** @{ */

  /** @} */

  /** \name Member functions (misc) */
  /** @{ */

  /**
   * Open the ROOT file containing the ROOT truees with the MCMC
   * samples.
   * @param filename The filename
   * @return An error code.
   */
  int OpenRootFile(std::string filename);

  /**
   * Close the ROOT file. */
  void CloseRootFile();

  /**
   * Print info about file on screen. */
  void PrintInfo();

  /** @} */

 protected:

  /**
   * The ROOT containing the ROOT tree with the MCMC samples. */
  TFile* fFile;

  /**
   * The vector of ROOT trees containing the MCMC samples. */
  std::vector<TTree*> fTrees;

  /**
   * Number of MCMC chains. */
  int fNTrees;

  /**
   * Number of parameters. */
  int fNParameters;

  /**
   * Number of samples in the pre-run. */
  int fNSamplesPreRun;

  /**
   * Number of samples in the main run. */
  int fNSamplesMainRun;

};
// ---------------------------------------------------------

#endif

