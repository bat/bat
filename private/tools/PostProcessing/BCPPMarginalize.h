#ifndef __BCPPMARGINALIZE__H
#define __BCPPMARGINALIZE__H

/*!
 * \class BCPPMarginalize
 * \brief A post-processing class for plotting marginalized
 * distributions.
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

// ---------------------------------------------------------

class BCPPMarginalize : public BCPostProcessor
{
 public:

  /** \name Constructors and destructors */
  /** @{ */

  /**
   * The default constructor. */
  BCPPMarginalize();

  /**
   * The default destructor. */
  ~BCPPMarginalize();

  /** @} */

  /** \name Member functions (misc) */
  /** @{ */

  /**
   * Build a marginalized distribution.
   * @param parindex The parameter index
   * @param nbins The number of bins
   * @param parmin The interval minimum
   * @param parmax The interval maximum
   * @return The BCH1D histogram
   */
  BCH1D* BuildMarginalized1D(int parindex, int nbins, double parmin, double parmax);

  /**
   * Build a marginalized distribution.
   * @param parindex1 The parameter index
   * @param nbins1 The number of bins
   * @param parmin1 The interval minimum
   * @param parmax1 The interval maximum
   * @param parindex2 The parameter index
   * @param nbins2 The number of bins
   * @param parmin2 The interval minimum
   * @param parmax2 The interval maximum
   * @return The BCH1D histogram
   */
  BCH2D* BuildMarginalized2D(int parindex1, int nbins1, double parmin1, double parmax1, int parindex2, int nbins2, double parmin2, double parmax2);

  /** @} */

 private:


};
// ---------------------------------------------------------

#endif

