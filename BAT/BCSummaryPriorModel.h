#ifndef __BCSUMMARYPRIORMODEL__H
#define __BCSUMMARYPRIORMODEL__H

/*!
 * \class BCSummaryPriorModel

 * A helper class for the BCSummaryTool.
 * \brief A helper class for the BCSummaryTool.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0.0
 * \date 15.02.2010
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

// ---------------------------------------------------------

class BCSummaryPriorModel : public BCModel
{
 public:

   // Constructors and destructor

   /**
    * The default constructor. */
   BCSummaryPriorModel();

   /**
    * A constructor.
    * @param name The name of the model. */
   BCSummaryPriorModel(const char * name);

   /**
    * The default destructor. */
   ~BCSummaryPriorModel();

   /**
    * Set a pointer to the model under study.
    * @param model The model under study. */
   void SetModel(BCModel * model);

   /**
    * Calculates and returns the log of the prior probability at a
    * given point in parameter space.
    * @param parameters A vector of coordinates in the parameter space.
    * @return The prior probability. */
   double LogAPrioriProbability(const std::vector<double> &parameters);

   /**
    * Calculates and returns the log of the Likelihood at a given point
    * in parameter space.
    * @param parameters A vector of coordinates in the parameter space.
    * @return The log likelihood. */
   double LogLikelihood(const std::vector<double> &parameters);

 private:

   /**
    * A pointer to the model under study. */
   BCModel * fTestModel;

};
// ---------------------------------------------------------

#endif

