#ifndef __BCPPDIAGNOSTICS__H
#define __BCPPDIAGNOSTICS__H

/*!
 * \class BCPPDiagnostics
 * \brief A post-processing class for diagnostics of MCMC chains.
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

#include <BCPostProcessor.h>

class BCH1D;
class BCH2D;

#include <string>

// ---------------------------------------------------------

class BCPPDiagnostics : public BCPostProcessor
{
 public:

  /** \name Constructors and destructors */
  /** @{ */

  /**
   * The default constructor. */
  BCPPDiagnostics();

  /**
   * The default destructor. */
  ~BCPPDiagnostics();

  /** @} */

  /** \name Member functions (get) */
  /** @{ */

  /**
   * Calculates the number of batches which have the full batch
   * length. */
  int GetNBatches();

  /** @} */

  /** \name Member functions (set) */
  /** @{ */

  /**
   * Set the batch length.
   * @param length The batch length. */
  void SetBatchLength(int length)
  { fBatchLength = length; };

  /** @} */

  /** \name Member functions (misc) */
  /** @{ */

  /**
   * Draw the log probability
   * @param filename The filename (a pdf file)
   * @param options: \n
   * all : draw it for all chains \n
   * sum : draw it for the sum of all chains */
  void PrintLogProbability(std::string filename, std::string options);

  /**
   * Print parameter-quantities for batches
   * @param filename The filename (a pdf file) */
  void PrintBatchQuantities(std::string filename);

  /**
   * Calculate parameter-quantities for batches. */
  void CalculateBatchQuantities();

  /** @} */

 private:

  /**
   * Length of a batch for calculating mean and variance. */
  int fBatchLength;

  std::vector< std::vector<double> > fParameterMean;
  std::vector< std::vector<double> > fParameterVariance;
  std::vector< std::vector<double> > fParameterStd;

  std::vector< double > fParameterBatchMean;
  std::vector< double > fParameterBatchVariance;
  std::vector< double > fParameterBatchStdDev;


};
// ---------------------------------------------------------

#endif

