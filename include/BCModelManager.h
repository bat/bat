/*! \class BCModelManager
 *  \brief A class managing models
 *
 * A class which manages a set of models and a data set. It performs
 * similar operations on all models and handles model comparison.
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, 
 *          kevin.kroeninger *at* phys *dot* uni *minus* goettingen *dot* de 
 *
 * CREATED: 14.03.2007 by Kevin 
 * 
 * REVISION: 
 *
 * 16.03.2007 Kevin,  modified ROOT interface\n
 * 12.06.2007 Kevin,  renamed to BCModelManager \n
 * 30.08.2007 Dano,   added method to set max number of iterations in all
 *                    attached models\n
 *
 * --------------------------------------------------------- 
 *
 */ 

// --------------------------------------------------------- 

#ifndef __BCMODELMANAGER__H
#define __BCMODELMANAGER__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelManager
{

 public:

	// constructor and destructor 

	/**
	 * A constructor. 
	 */ 
	BCModelManager(); 

	/**
	 * A destructor. 
	 */ 
	virtual ~BCModelManager(); 

	// methods (get) 

	/** 
	 * @return The number of models. 
	 */ 
	int GetNModels() 
	{ return int(fModelContainer -> size()); }; 

	/** 
	 * @param index The index of the model 
	 * @return The BCModel. 
	 */ 
	BCModel * GetModel(int index) 
	{ return fModelContainer -> at(index); }; 

	/**
	 * @return The number of entries in the data container
	 */ 
	int GetNDataPoints() 
	{ return fDataSet -> GetNDataPoints(); }; 

	/**
	 * @param index The index of the data point
	 * @return The data point
	 */ 
	BCDataPoint * GetDataPoint(int index) 
		{ return fDataSet -> GetDataPoint(index); }; 

	/** 
	 * @return The data set 
	 */ 
	BCDataSet * GetDataSet()
		{ return fDataSet; }; 

	// methods (set) 

	/**
	 * @param dataset A data set 
	 */ 
	void SetDataSet(BCDataSet * dataset); 

	/**
	 * Set maximum number of iterations for all models added to the manager.
	 * Only works on models already added to the manager.
	 * @param niterations
	 */ 
	void SetNIterationsMax(int niterations);

	/** 
	 * @param method The marginalization method 
	 */ 
	void SetMarginalizationMethod(BCIntegrate::BCMarginalizationType method); 

	/** 
	 * @param method The integration method
	 */ 
	void SetIntegrationMethod(BCIntegrate::BCIntegrationType method); 

	/** 
	 * @param method The mode finding method 
	 */ 
	void SetModeFindingMethod(BCIntegrate::BCModeFindingType method); 

	/**
	 * @param niterations Number of iterations per dimension for Monte
	 * Carlo integration.
	 */
	void SetNiterationsPerDimension(int niterations); 
	/**
	 * @param n Number of samples per 2D bin per variable in the
	 * Metropolis marginalization.  Default is 100.
	 */
	void SetNSamplesPer2DBin(int n); 

	/** 
	 * @param relprecision The relative precision envisioned for Monte
	 * Carlo integration
	 */ 
	void SetRelativePrecision(double relprecision); 

	/**
	 * @param n Number of bins per dimension for the marginalized
	 * distributions.  Default is 100. Minimum number allowed is 2.
	 */
	void SetNbins(int n);

	/**
	 * Sets index of the x values in function fits. 
	 * @param index Index of the x values 
	 */  
	void SetFitFunctionIndexX(int index); 
	 
	/**
	 * Sets index of the y values in function fits. 
	 * @param index Index of the y values 
	 */  
	void SetFitFunctionIndexY(int index); 

	void SetFitFunctionIndices(int indexx, int indexy); 
	
	/**
	 * Sets the data point containing the lower boundaries of possible
	 * data values
	 */ 
	void SetDataPointLowerBoundaries(BCDataPoint* datasetlowerboundaries); 

	/**
	 * Sets the data point containing the upper boundaries of possible
	 * data values
	 */ 
	void SetDataPointUpperBoundaries(BCDataPoint* datasetupperboundaries); 

	/**
	 * Sets the lower boundary of possible data values for a particular
	 * variable
	 */ 
	void SetDataPointLowerBoundary(int index, double lowerboundary); 

	/**
	 * Sets the upper boundary of possible data values for a particular
	 * variable
	 */ 
	void SetDataPointUpperBoundary(int index, double upperboundary); 

	/**
	 * Set the lower and upper boundaries for possible data values for a
	 * particular variable
	 */ 
	void SetDataBoundaries(int index, double lowerboundary, double upperboundary); 


	// methods 

	/** 
	 * Adds a model to the container
	 * @param model The model
	 * @param probability The a priori probability 
	 * @see AddModel(BCModel * model)
	 * @see SetModelPrior(BCModel * model, double probability)
	 */ 
	void AddModel(BCModel * model, double probability=0.); 

	/**
	 * Adds a data point to the data container. 
	 * @param datapoint The data point
	 */ 
	void AddDataPoint(BCDataPoint * datapoint) 
	{ fDataSet -> AddDataPoint(datapoint); }; 

	/** 
	 * Reads tree data from a ROOT file. 
	 * Opens a ROOT file and gets a ROOT tree. It creates data set
	 * containing the values read from the file. 
	 * @param filename The filename of the ROOT file 
	 * @param treename The name of the ROOT tree 
	 * @param branchnames A vector of the names of the branches 
	 * @see ReadDataFromFileHist(char * filename, char * histname, const char*  branchnames); 
	 * @see ReadDataFromFileTxt(char * filename, int nbranches); 
	 */ 
	int ReadDataFromFileTree(const char * filename, const char * treename, const char * branchnames); 

	/** 
	 * Reads data from a txt file. 
	 * Opens a txt file and creates data set
	 * containing the values read from the file. 
	 * @param filename The filename of the ROOT file 
	 * @param nbranches The number of variables 
	 * @see ReadDataFromFileTree(char * filename, char * treename, std::vector<char*> branchnames)
	 * @see ReadDataFromFileHist(char * filename, char * histname, const char * branchnames); 
	 */
	int ReadDataFromFileTxt(const char * filename, int nbranches); 

	/** 
	 * Reads data from a txt file (user specifies). 
	 * Opens a file and creates data set
	 * containing the values read from the file (user specifies). 
	 * @param filename The filename of the ROOT file 
	 * @param nbranches The number of variables 
	 * @see ReadDataFromFileTree(char * filename, char * treename, std::vector<char*> branchnames)
	 * @see ReadDataFromFileHist(char * filename, char * histname, const char * branchnames); 
	 * @see ReadDataFromFileTxt(char * filename, int nbranches); 
	 */
	int ReadDataFromFileUser(const char *  filename, std::vector<int> options_int, std::vector<double> options_double);

	/** 
	 * Calculates the normalization of the likelihood for each model in
	 * the container.
	 */   
	void Normalize(); 

	/**
	 * Does the mode finding 
	 */ 
	void FindMode(); 

	/**
	 * Marginalize all probabilities wrt. single parameters and all
	 * combinations of two parameters for all models. 
	 */ 
	void MarginalizeAll(); 

	/**
	 * Resets the data set 
	 */ 
	void ResetDataSet()
	{ fDataSet -> Reset(); }; 

	/**
	 * Prints a summary into a file. If filename is omitted the summary will
	 * be printed onto the screen
	 * @param filename name of the file to write into.
	 */
	void PrintSummary(const char * filename=0);

 private: 

	/**
	 * A model container. 
	 */ 
	BCModelContainer * fModelContainer; 

	/**
	 * A data set
	 */
	BCDataSet * fDataSet;

}; 

// ---------------------------------------------------------

#endif
