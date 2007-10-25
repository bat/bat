/*! \class BCModelOutput
 *  \brief Case which does ROOT based I/O 
 *
 * This class saves the results of an analysis with a certain model
 * into a ROOT tree. 
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, 
 *          kevin.kroeninger *at* phys *dot* uni *minus* goettingen *dot* de 
 *
 * CREATED: 23.10.2007 
 * 
 * REVISION: 
 *
 * 25.10.2007 Kevin, added margnialized histograms to output 
 *
 * --------------------------------------------------------- 
 *
 */ 

// --------------------------------------------------------- 

#ifndef __BCMODELOUTPUT__H
#define __BCMODELOUTPUT__H

const int MAXNPARAMETERS = 20; 

#include "BCModel.h"

#include "TTree.h"
#include "TFile.h" 

// --------------------------------------------------------- 

class BCModelOutput 
{
  
 public:
  
	// constructors and destructor 

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
	 * The default destructor. 
	 */ 
	virtual ~BCModelOutput(); 

	// methods (set) 

	/**
	 * Sets the model this output class is assigned to. 
	 */ 
	void SetModel(BCModel * model) 
	{ fModel = model; }; 

	/**
	 * Sets the output filename 
	 */ 
	void SetFile(const char * filename); 

	// methods (get) 

	/**
	 * Returns the ROOT tree 
	 */ 
	TTree * GetAnalysisTree()
		{ return fAnalysisTree; }; 

	/**
	 * Returns the ROOT file 
	 */
	TFile * GetFile()
		{ return fOutputFile; }; 

	// methods 

	/**
	 * Fill the output tree with the current information
	 */ 
	void FillAnalysisTree(); 

	/**
	 * Writes the marginalized histograms to the ROOT file 
	 */ 
	void WriteMarginalizedDistributions(); 

	/**
	 * Closes the file 
	 */ 
	void Close(); 

 private:

	/**
	 * Initializes the output tree 
	 */
	void InitializeAnalysisTree(); 

	/**
	 * The output tree 
	 */
	TTree * fAnalysisTree; 

	/** 
	 * The output filename
	 */ 
	const char * fFilename; 

	/**
	 * The output file 
	 */ 
	TFile * fOutputFile; 

	/**
	 * The model this output class is assigned to 
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

}; 

// --------------------------------------------------------- 

#endif 
