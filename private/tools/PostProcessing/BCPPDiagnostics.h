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

  /** \name Member functions (misc) */
  /** @{ */

  /**
   * Create a BCH1D of the log probability
   * @param options: \n
   * all : draw it for all chains \n
   * sum : draw it for the sum of all chains
   * @return The BAT histogram */
  BCH1D* PlotLogProbability(std::string options);

  /** @} */

 private:


};
// ---------------------------------------------------------

#endif

