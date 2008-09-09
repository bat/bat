/*!
 * \class BCModelTest
 * \brief The class for testing model hypotheses
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class is used for calculating the p-value of a model.
 * 
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger. 
 * All rights reserved. 
 */ 

// --------------------------------------------------------- 

#ifndef __BCMODELTEST__H
#define __BCMODELTEST__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelTest : public BCModel 
{

 public: 

	/** \name Constructors and destructors */ 
	/* @{ */ 

	/** 
	 * The default constructor. 
	 */ 	BCModelTest(const char* name); 

	/** 
	 * The default destructor. 
	 */ 
	~BCModelTest(); 

	/* @} */ 

	/** \name Member functions (get) */ 
	/* @{ */ 

	/*
	 * Calculated the p-value. 
	 * @param flag_histogram A histogram is either filled or not. 
	 * @return The p-value. 
	 */ 
	double GetCalculatedPValue(bool flag_histogram = false); 

	/*
	 * @return The distribution of log(likelihood). 
	 */ 
	TH1D * GetHistogramLogProb()
	{ return fHistogramLogProb; }; 

	/* @} */ 

	/** \name Member functions (set) */ 
	/* @{ */ 

	/*
	 * Set the model to be tested. 
	 * @param testmodel A pointer to the model to be tested. 
	 */ 
	void SetTestModel(BCModel * testmodel)
	{ fTestModel = testmodel; }; 

	/*
	 * Sets the set of parameters which the p-values is calculated for.
	 * @param parameters The parameters.
	 * @return An error code. 
	 */ 
	int SetTestPoint(std::vector<double> parameters); 

	/* @} */ 

	/** \name Member functions (miscellaneous methods) */ 
	/* @{ */ 

	double LogLikelihood(std::vector <double> parameters); 

	double LogAPrioriProbability(std::vector <double> parameters)
	{ return 0; }; 

	void MCMCUserInterface();

	/* @} */ 

 private: 

	/*
	 * A map of data points and data values. 
	 */ 
	std::vector<int> fMapDataPoint; 
	std::vector<int> fMapDataValue; 

	/*
	 * Counter for the evaluation of the p-value. 
	 */ 
	int fPValueBelow; 
	int fPValueAbove; 

	/*
	 * A pointer to the model which is tested. 
	 */ 
	BCModel * fTestModel; 

	/*
	 * A data set used for temporary storage. 
	 */ 
	BCDataSet * fTemporaryDataSet; 

	/*
	 * The log(likelihood) and its ranges. 
	 */ 
	double fLogLikelihood; 	
	double fLogLikelihoodMin; 	
	double fLogLikelihoodMax; 	

	/*
	 * The distribution of log(likelihood). 
	 */ 
	TH1D * fHistogramLogProb; 
	
}; 

// --------------------------------------------------------- 

#endif 
