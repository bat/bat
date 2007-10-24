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

	void SetModel(BCModel * model) 
	{ fModel = model; }; 

	void SetFile(const char * filename); 

	// methods 

	void Fill(); 

	void Close(); 

 private:

	/**
	 * Initializes the output tree 
	 */
	void InitializeTree(); 

	/**
	 * The output tree 
	 */
	TTree * fOutputTree; 

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
	 * The tree variables 
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
