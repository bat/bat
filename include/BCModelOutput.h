/**  
 * \class BCModelOutput
 * \brief A class for creating an (ROOT) output file. 
 * \author D. Kollar  
 * \author K. Kr&ouml;ninger  
 * \version 1.0  
 * \date 12.11.2007  
 *  
 * This class defines an output interface for the analysis. It creates
 * a ROOT file which can contain summary information, histograms and
 * Markov chains. 
 *  
 * Copyright (C) 2007, D. Kollar, K. Kr&ouml;ninger  
 */  

// --------------------------------------------------------- 

#ifndef __BCMODELOUTPUT__H
#define __BCMODELOUTPUT__H

const int MAXNPARAMETERS = 20;

#include "BCModel.h"

#include <TTree.h>
#include <TFile.h>

// --------------------------------------------------------- 

class BCModelOutput 
{
  
 public:
  
	/** \name Constructors and destructors */ 
	/* @{ */

	/** 
	 * The default constructor. 
	 */ 
	BCModelOutput(); 

	/** 
	 * A constructor. 
	 * @param model The model to which this output class is assigned. 
	 * @param filename Name of the output file. 
	 */
	BCModelOutput(BCModel * model, const char * filenname); 

	/** 
	 * The default copy constructor. 
	 */ 
	BCModelOutput(const BCModelOutput & modeloutput); 

	/** 
	 * The default destructor. 
	 */ 
	virtual ~BCModelOutput(); 

	/* @} */ 

	/** \name Assignment operators */ 
	/* @{ */ 

	/**
	 * The defaut assignment operator 
	 */ 
	BCModelOutput & operator = (const BCModelOutput & modeloutput); 

	/* @} */ 

	/** \name Getters */ 
	/* @{ */ 

	/**
	 * Returns the output TTree tree. 
	 * @return The pointer to the output TTree. 
	 */ 
	TTree * GetAnalysisTree()
		{ return fAnalysisTree; }; 

	/**
	 * Returns the output TFile. 
	 * @return The pointer to the output TFile. 
	 */
	TFile * GetFile()
		{ return fOutputFile; }; 

	/* @} */ 

	/** \name Setters */ 
	/* @{ */ 

	/**
	 * Assign a BCModel to this output class. 
	 * @param model A pointer to the BCModel. 
	 */ 
	void SetModel(BCModel * model) 
	{ fModel = model; }; 

	/**
	 * Sets the output filename. 
	 * @param filename The filename. 
	 */ 
	void SetFile(const char * filename); 

	/**
	 * Flag for writing Markov chain to file 
	 * @param flag Writes (true) or does not write (false) the Markov chain 
	 */ 
	void WriteMarkovChain(bool flag); 

	/* @} */ 

	/** \name Miscellaneous methods */ 
	/* @{ */ 	

	/**
	 * Fill the output TTree with the current information. 
	 */ 
	void FillAnalysisTree(); 

	/**
	 * Writes the marginalized histograms to the TFile. 
	 */ 
	void WriteMarginalizedDistributions(); 

	/**
	 * Writes the error band histogram into the TFile. 
	 */ 
	void WriteErrorBand(); 

	/**
	 * Closes the TFile. 
	 */ 
	void Close(); 

	/* @} */ 

 private:

	/* 
	 * Copies this BCModelOutput into another one 
	 */ 
	void Copy(BCModelOutput & modeloutput) const; 

	/**
	 * Initialize the output TTree. 
	 */
	void InitializeAnalysisTree(); 

	/**
	 * Initialize the Markov Chain TTree. 
	 */
	void InitializeMarkovChainTrees(); 

	/**
	 * Pointer to the TTree containing the summary output information.
	 */
	TTree * fAnalysisTree; 

	/**
	 * Pointer to the TTree containing the Markov chain. 
	 */ 
	TTree * fMarkovChainTree; 

	/*
	 * The trees containing the Markov chains. The length of the vector
	 * is fMCMCNChains. 
	 */ 
	std::vector<TTree *> fMarkovChainTrees; 

	/** 
	 * The output filename
	 */ 
	const char * fFilename; 

	/**
	 * Pointer to the output TFile. 
	 */ 
	TFile * fOutputFile; 

	/**
	 * Pointer to the model this output class is assigned to
	 */ 
	BCModel * fModel; 

	/**
	 * The analysis tree variables 
	 */
	int fIndex; 
	int fNParameters; 
	double fProbability_apriori; 
	double fProbability_aposteriori; 
	double fMode_global[MAXNPARAMETERS]; 
	double fMode_marginalized[MAXNPARAMETERS]; 
	double fMean_marginalized[MAXNPARAMETERS]; 
	double fMedian_marginalized[MAXNPARAMETERS]; 
	double fQuantile_05[MAXNPARAMETERS]; 
	double fQuantile_10[MAXNPARAMETERS]; 
	double fQuantile_16[MAXNPARAMETERS]; 
	double fQuantile_84[MAXNPARAMETERS]; 
	double fQuantile_90[MAXNPARAMETERS]; 
	double fQuantile_95[MAXNPARAMETERS]; 

	/**
	 * The markov chain tree variables 
	 */ 
	std::vector<double> * fParameters; 
	std::vector<double> * fLogLikelihood; 
	std::vector <int> * fIteration; 

}; 

// --------------------------------------------------------- 

#endif 
