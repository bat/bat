#ifndef __BCSUMMARYTOOL__H
#define __BCSUMMARYTOOL__H

/*!
 * \class BCSummaryTool

 * This class can be used to summarize the results of an analysis. The
 * prior and posterior probabilities are compared.
 * \brief A class for summarizing the results of an analysis.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0.0
 * \date 15.02.2010
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
#include <vector>
#include <limits>

#include "BCModel.h"

// ---------------------------------------------------------

class BCSummaryTool : private BCModel
{
 public:

   /** \name Constructors and destructors */
   /** @{ */

   /**
    * constructor. */
   BCSummaryTool(BCModel * model = 0);

   /**
    * destructor. */
   ~BCSummaryTool();

   /** @} */
   /** \name Member functions (set) */
   /** @{ */

   /**
    * Set the model to be summarized.
    * @param model The BCModel to be summarized.*/
	 void SetModel(BCModel * model);

   /** @} */
   /** \name Member functions (misc) */
   /** @{ */

   /**
    * Calculates and returns the log of the prior probability at a
    * given point in parameter space.
    * @param parameters A vector of coordinates in the parameter space.
    * @return The prior probability. */
	double LogAPrioriProbability(const std::vector<double> &parameters)
	   { return 0; }

   /**
    * Calculates and returns the log of the Likelihood at a given point
    * in parameter space.
    * @param parameters A vector of coordinates in the parameter space.
    * @return The log likelihood. */
   double LogLikelihood(const std::vector<double> &parameters)
	    { return (fModel) ? fModel->LogAPrioriProbability(parameters) : -std::numeric_limits<double>::infinity(); }

	 /**
		* Calculates user observables according to the model. */
	 void CalculateObservables(const std::vector<double> &parameters)
	    { if (fModel) fModel->CalculateObservables(parameters); }

   /**
    * Calculate the marginalized distributions using the prior
    * knowledge alone.
    * @return An error flag.
    */
   int CalculatePriorModel();

   /**
    * Draw a comparison of the prior knowledge to the posterior
    * knowledge for each parameter.
    * @return An error flag. */
   int DrawKnowledgeUpdatePlot1D(unsigned index, std::string options_post="", std::string options_prior="");

   /**
    * Print a comparison of the prior knowledge to the posterior
    * knowledge for each parameter.
    * @return An error flag. */
   int PrintKnowledgeUpdatePlot1D(int index, const char * filename, std::string options_post="-", std::string options_prior="-");

   /**
    * Draw a comparison of the prior knowledge to the posterior.
    * @return An error flag. */
	int DrawKnowledgeUpdatePlot2D(unsigned index1, unsigned index2, bool flag_slice=false, double interval_content=68e-2);

   /**
    * Print a comparison of the prior knowledge to the posterior
    * knowledge for each parameter.
    * @return An error flag. */
	int PrintKnowledgeUpdatePlots(const char * filename = "update.pdf", unsigned hdiv=1, unsigned vdiv=1, std::string options="-", double interval_content=68e-2);

   /** @} */

 private:

   /**
    * The model whose results are summarized */
   BCModel * fModel;

};
// ---------------------------------------------------------

#endif

