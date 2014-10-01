#ifndef __BCGOFTEST__H
#define __BCGOFTEST__H

/*!
 * \class BCGoFTest
 * \brief The class for testing model hypotheses
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class is used for calculating the p-value of a model.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModel.h"

// ROOT classes
class TH1D;

// BAT classes
class BCDataSet;

// ---------------------------------------------------------

class BCGoFTest : public BCModel
{
   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * Default constructor.
       */
      BCGoFTest(const char * name);

      /**
       * Default destructor. */
      ~BCGoFTest();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * Calculated the p-value.
       * @param flag_histogram A histogram is either filled or not.
       * @return p-value */
      double GetCalculatedPValue(bool flag_histogram = false);

      /**
       * @return distribution of log(likelihood) */
      TH1D * GetHistogramLogProb()
         { return fHistogramLogProb; };

      /**
       * @return pointer to the tested model */
      BCModel * GetTestModel()
         { return fTestModel; };

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * Set the model to be tested.
       * @param testmodel pointer to the model to be tested */
      void SetTestModel(BCModel * testmodel)
         { fTestModel = testmodel; };

      /**
       * Sets the set of parameters which the p-values is calculated for.
       * @param parameters parameters
       * @return error code */
      int SetTestPoint(std::vector<double> parameters);

      /** @} */
      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      double LogLikelihood(const std::vector<double> &parameters);

      double LogAPrioriProbability(const std::vector<double> & /*parameters*/)
         { return 0; };

      void MCMCUserIterationInterface();

      /** @} */

   private:

      /**
       * A map of data points and data values. */
      std::vector<int> fMapDataPoint;
      std::vector<int> fMapDataValue;

      /**
       * Counter for the evaluation of the p-value. */
      int fPValueBelow;
      int fPValueAbove;

      /**
       * A pointer to the model which is tested. */
      BCModel * fTestModel;

      /**
       * A data set used for temporary storage. */
      BCDataSet * fTemporaryDataSet;

      /**
       * The log(likelihood) and its range. */
      double fLogLikelihood;
      double fLogLikelihoodMin;
      double fLogLikelihoodMax;

      /**
       * The distribution of log(likelihood). */
      TH1D * fHistogramLogProb;
};

// ---------------------------------------------------------

#endif
