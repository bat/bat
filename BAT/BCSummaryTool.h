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
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <string>
#include <vector>

class BCModel;
class BCSummaryPriorModel;

// ---------------------------------------------------------

class BCSummaryTool
{
 public:

   /** \name Constructors and destructors */
   /** @{ */

   /**
    * The default constructor. */
   BCSummaryTool();

   /**
    * A constructor. */
   BCSummaryTool(BCModel * model);

   /**
    * The default destructor. */
   ~BCSummaryTool();

   /** @} */
   /** \name Member functions (get) */
   /** @{ */

   /**
    * Retrieve pointer to the Prior model to allow for its detailed setup */
   BCSummaryPriorModel * GetPriorModel()
      { return fPriorModel; }

   /** @} */
   /** \name Member functions (set) */
   /** @{ */

   /**
    * Set the model to be summarized.
    * @param model The BCModel to be summarized.*/
   void SetModel(BCModel * model)
      { fModel = model; };

   /** @} */
   /** \name Member functions (misc) */
   /** @{ */

   /**
    * Calculate the marginalized distributions using the prior
    * knowledge alone.
    * @return An error flag.
    */
   int CalculatePriorModel();

   /**
    * Copy the summary information from the model.
    * @return An error flag. */
   int CopySummaryData();

   /**
    * Print a summary plot for the parameters.
    * @return An error flag. */
   int PrintParameterPlot(const char * filename = "parameters.pdf");

   /**
    * Print a correlation matrix for the parameters.
    * @return An error flag. */
   int PrintCorrelationMatrix(const char * filename = "matrix.pdf");

   /**
    * Print a correlation plot for the parameters.
    * @return An error flag. */
   int PrintCorrelationPlot(const char * filename = "correlation.pdf");

   /**
    * Draw a comparison of the prior knowledge to the posterior
    * knowledge for each parameter.
    * @return An error flag. */
   int DrawKnowledgeUpdatePlot1D(int index, std::string options_post="", std::string options_prior="");

   /**
    * Print a comparison of the prior knowledge to the posterior
    * knowledge for each parameter.
    * @return An error flag. */
   int PrintKnowledgeUpdatePlot1D(int index, const char * filename, std::string options_post="", std::string options_prior="");

   /**
    * Print a comparison of the prior knowledge to the posterior
    * knowledge for each parameter.
    * @return An error flag. */
   int PrintKnowledgeUpdatePlots(const char * filename = "update.pdf", std::string options="");

   /**
    * Print a Latex table of the parameters.
    * @return An error flag. */
   int PrintParameterLatex(const char * filename);

   /** @} */

 private:

   /** Helper method to get an unique number to be used in histogram naming */
   static unsigned int getNextIndex()
      { return ++fHCounter; }

   /** helper variable to get an unique number to be used in histogram naming */
   static unsigned int fHCounter;

   /**
    * The model whose results are summarized */
   BCModel * fModel;

   /**
    * parameter names */
   std::vector<std::string> fParName;

   /**
    * parameter minima */
   std::vector<double> fParMin;

   /**
    * Parameter maxima */
   std::vector<double> fParMax;

   /**
    * Correlation coefficients.
    * Length of vector equals number of parameters * number of parameters. */
   std::vector<double> fCorrCoeff;

   /**
    * Marginalized modes.\n
    * Length of vector equals number of parameters. */
   std::vector<double> fMargMode;

   /**
    * Mean values.\n
    * Length of vector equals number of parameters. */
   std::vector<double> fMean;

   /**
    * Global modes.\n
    * Length of vector equals number of parameters. */
   std::vector<double> fGlobalMode;

   /**
    * Quantiles.\n
    * The following quantiles are stored: 0.05, 0.10, 0.16, 0.5, 0.84, 0.90, 0.95.\n
    * Length of vector equals number of parameters * number of quantiles. */
   std::vector<double> fQuantiles;

   /**
    * Smallest intervals.\n
    * For each parameter a set of the smallest intervals is recorded.\n
    * Structure: number of intervals n + n * (start, stop, local max, local max pos, integral)
    * Length of vector equals number of parameters * number of quantiles. */
   std::vector<double> fSmallInt;

   /**
    * RMS values.\n
    * Length of vector equals number of parameters. */
   std::vector<double> fRMS;

   /**
    * Sum of probabilities for quantiles */
   std::vector<double> fSumProb;

   /**
    * A model for calculating the marginalized distributions for the
    * prior probabilities. */
   BCSummaryPriorModel * fPriorModel;

   /**
    * A flag: check if marginalized information is present */
   bool fFlagInfoMarg;

   /**
    * A flag: check if optimization information is present */
   bool fFlagInfoOpt;

};
// ---------------------------------------------------------

#endif

