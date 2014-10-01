#ifndef __BCMODELOUTPUT__H
#define __BCMODELOUTPUT__H

/*!
 * \class BCModelOutput
 * \brief A class for creating an (ROOT) output file.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class defines an output interface for the analysis. It
 * creates a ROOT file which can contain summary information,
 * histograms and Markov chains.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>

// BAT classes
class BCModel;

// ROOT classes
class TTree;
class TFile;
class TObject;

// ---------------------------------------------------------

class BCModelOutput
{
   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      BCModelOutput();

      /**
       * A constructor.
       * @param model The model to which this output class is assigned.
       * @param filename Name of the output file. */
      BCModelOutput(BCModel * model, const char * filenname);

      /**
       * The default copy constructor. */
      BCModelOutput(const BCModelOutput & modeloutput);

      /**
       * The default destructor. */
      virtual ~BCModelOutput();

      /** @} */

      /** \name Assignment operators */
      /** @{ */

      /**
       * The defaut assignment operator */
      BCModelOutput & operator = (const BCModelOutput & modeloutput);

      /** @} */

      /** \name Getters */
      /** @{ */

      /**
       * Returns the output TTree tree.
       * @return The pointer to the output TTree */
      TTree * GetAnalysisTree()
         { return fAnalysisTree; };

      /**
       * Returns the output TFile.
       * @return The pointer to the output TFile */
      TFile * GetFile()
         { return fOutputFile; };

      /** @} */

      /** \name Setters */
      /** @{ */

      /**
       * Assign a BCModel to this output class.
       * @param model A pointer to the BCModel */
      void SetModel(BCModel * model);


      /**
       * Sets the output filename.
       * @param filename The filename */
      void SetFile(const char * filename);

      /** @} */

      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Flag for writing Markov chain to file
       * @param flag Writes (true) or does not write (false) the Markov chain */
      void WriteMarkovChain(bool flag = true);

      /**
       * Fill the output TTree with the current information. */
      void FillAnalysisTree();

      /**
       * Writes the marginalized histograms to the TFile. */
      void WriteMarginalizedDistributions();

      /**
       * Writes any object derived from TObject to TFile. */
      void Write(TObject * o);

      /**
       * Closes the TFile. */
      void Close();

      /** @} */

   private:

      /**
       * Initialize the variables */
      void Init();

      /**
       * Copies this BCModelOutput into another one */
      void Copy(BCModelOutput & modeloutput) const;

      /**
       * Initialize the output TTree. */
      void InitializeAnalysisTree();

      /**
       * Initialize SA TTree. */
      void InitializeSATree();

      /**
       * Pointer to the TTree containing the summary output information. */
      TTree * fAnalysisTree;

      /**
       * The trees containing the Markov chains. The length of the vector
       * is fMCMCNChains. */
      std::vector<TTree *> fMarkovChainTrees;

      /**
       * The tree for the simulated annealing. */
      TTree * fTreeSA;

      /**
       * The output filename */
      char * fFileName;

      /**
       * Pointer to the output TFile. */
      TFile * fOutputFile;

      /**
       * Pointer to the model this output class is assigned to */
      BCModel * fModel;

      /**
       * The analysis tree variables */
      int fIndex;
      unsigned int fNParameters;
      double fProbability_apriori;
      double fProbability_aposteriori;
      std::vector<double> fMode_global;
      std::vector<double> fMode_marginalized;
      std::vector<double> fMean_marginalized;
      std::vector<double> fMedian_marginalized;
      std::vector<double> fQuantile_05;
      std::vector<double> fQuantile_10;
      std::vector<double> fQuantile_16;
      std::vector<double> fQuantile_84;
      std::vector<double> fQuantile_90;
      std::vector<double> fQuantile_95;

};

// ---------------------------------------------------------

#endif
