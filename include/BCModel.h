/*! \class BCModel
 *  \brief Base class for all user models
 *
 * This class defines a model. It contains the a container of
 * parameters, their prior distributions and the conditional
 * probabilities given those parameters.  The methods which implement
 * the prior and conditional probabilities have to be overloaded by
 * the user in the user defined model class which will inherit from
 * this class.
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, 
 *          kevin.kroeninger *at* phys *dot* uni *minus* goettingen *dot* de 
 *
 * CREATED: 02.03.2007 
 * 
 * REVISION: 
 *
 * 02.03.2007 Kevin, added comments and header, added marginalized probability class. \n
 * 14.03.2007 Kevin, added a priori and a posteriori probabilities to the class, 
 *                   some renaming of functions. \n
 * 17.04.2007 Kevin, added limits for the variables for goodness-of-fit test \n
 * 30.04.2007 Kevin, added creation of data for goodness-of-fit test \n
 * 15.05.2007 Kevin, added summary print-out and plots \n
 * 23.05.2007 Kevin, added virtual functions for importance sampling integration (with and without Markov chains)\n
 * 12.06.2007 Kevin, some renaming work \n
 * 01.08.2007 Dano,  added MarginalizeAll method and respective GetMarginalized methods,
 *                   added some warnings\n
 * 06.08.2007 Dano,  changed FindMode to use Simulated Annealing (SA) mode finding\n
 * 07.08.2007 Dano,  changed the model to use Log of everything, all standard calls are 
 *                   just the inline calls to TMath::Exp() calls of the Log versions\n
 * 01.10.2007 Kevin, added the error band calculation
 *
 * --------------------------------------------------------- 
 *
 */ 

// --------------------------------------------------------- 

#ifndef __BCMODEL__H
#define __BCMODEL__H

#include <vector.h>
#include <string.h> 

#include <TROOT.h>
#include <TGraph.h> 

#include "BCDataPoint.h" 
#include "BCDataSet.h" 
#include "BCParameter.h" 
#include "BCH1D.h" 
#include "BCH2D.h" 
#include "BCIntegrate.h"

const int MAXNDATAPOINTVALUES = 20; 

// --------------------------------------------------------- 

class BCModel : public BCIntegrate
{
  
 public:
  
	// constructors and destructor 

	/** 
	 * The default constructor. 
	 */ 
	BCModel(); 

	/** 
	 * A constructor. 
	 * @param name The name of the model 
	 */
	BCModel(const char * name); 

	/** 
	 * The default destructor. 
	 */ 
	virtual ~BCModel(); 

	// methods (get) 

	/** 
	 * @return The name of the model
	 */
	std::string GetName()
	{ return fName; }; 

	/** 
	 * @return The index of the model  
	 */ 
	int GetIndex()
	{ return fIndex; }; 

	/** 
	 * @return The a priori probability
	 */ 
	double GetModelAPrioriProbability()
	{ return fModelAPriori; }; 

	/** 
	 * @param model The model. 
	 * @return The a posteriori probability
	 */ 
	double GetModelAPosterioriProbability() 
	{ return fModelAPosteriori; }; 

	/** 
	 * @return The normalization of the likelihood 
	 */ 
	double GetNormalization()
	{ return fNormalization; };

	/** 
	 * @return The data set 
	 */ 
	BCDataSet* GetDataSet()
		{ return fDataSet; }; 

	/**
	 * @return The lower boundaries of possible data values 
	 */ 
	BCDataPoint* GetDataPointLowerBoundaries() 
		{ return fDataPointLowerBoundaries; }; 

	/**
	 * @return The upper boundaries of possible data values 
	 */ 
	BCDataPoint* GetDataPointUpperBoundaries() 
		{ return fDataPointUpperBoundaries; }; 

	/**
	 * @param index The index of the variable 
	 * @return The lower boundary of possible data values for a particular variable. 
	 */ 
	double GetDataPointLowerBoundary(int index)
	{ return fDataPointLowerBoundaries -> GetValue(index); }; 

	/**
	 * @param index The index of the variable 
	 * @return The upper boundary of possible data values for a particular variable. 
	 */ 
	double GetDataPointUpperBoundary(int index)
	{ return fDataPointUpperBoundaries -> GetValue(index); }; 

	/** 
	 * @return The number of data points in the current data set 
	 */ 
	int GetNDataPoints(); 

	/** 
	 * @param index The index of the data point 
	 * @return The data point in the current data set at index 
	 */ 
	BCDataPoint* GetDataPoint(int index); 

	/** 
	 * @return The minimum number of data points. 
	 */ 
	int GetNDataPointsMinimum()
	{ return fNDataPointsMinimum; }; 

	/** 
	 * @return The maximum number of data points. 
	 */ 
	int GetNDataPointsMaximum()
	{ return fNDataPointsMaximum; }; 

	/** 
	 * @return The number of parameters of the model  
	 */
	int GetNParameters()
	{ return fParameterSet -> size(); }; 

	/**
	 * @param index The index of the parameter in the parameter set 
	 * @return The parameter 
	 */ 
	BCParameter* GetParameter(int index); 

	/**
	 * @param name The name of the parameter in the parameter set 
	 * @return The parameter 
	 */ 
	BCParameter* GetParameter(const char * name); 

	/** 
	 * Returns the best fit parameters of the global probability. 
	 * @param index The index of the parameter 
	 * @return The best fit parameter 
	 */
	double GetBestFitParameter(int index) 
	{ return fBestFitParameters.at(index); }; 

	/** 
	 * Returns the best fit parameters of the global probability. 
	 * @return The best fit parameters 
	 */
	std::vector <double> GetBestFitParameters() 
		{ return fBestFitParameters; }; 

	/** 
	 * Returns the best fit parameters of the marginalized probabilities. 
	 * @param index The index of the parameter 
	 * @return The best fit parameter 
	 */
	double GetBestFitParameterMarginalized(int index) 
	{ return fBestFitParametersMarginalized.at(index); }; 

	/** 
	 * Returns the best fit parameters of the marginalized probabilities. 
	 * @return The best fit parameters 
	 */
	std::vector <double> GetBestFitParametersMarginalized() 
		{ return fBestFitParametersMarginalized; }; 

	/**
	 * Returns the 2-d histogram for the error band 
	 * @return The 2-d histogram 
	 */ 
	TH2D * GetErrorBandXY()
		{ return fErrorBandXY; }; 

	/**
	 * Returns a vector of y-values at a certain probability level. 
	 * @param level The level of probability
	 * @return The vector of y-values
	 */ 
	std::vector <double> GetErrorBand(double level); 

	TGraph * GetErrorBandGraph(double level1, double level2); 

	TGraph * GetFitFunctionGraph(std::vector <double> parameters);

	TGraph * GetFitFunctionGraph()
		{ return this -> GetFitFunctionGraph(this -> GetBestFitParameters()); }; 

	// methods (set) 

	/**
	 * Sets the name of the model. 
	 * @param name Name of the model 
	 */ 
	void SetName(const char * name) 
	{ fName = name; }; 

	/** 
	 * Sets the index of the model within the BCModelManager. 
	 * @param index The index of the model 
	 */ 
	void SetIndex(int index) 
	{ fIndex = index; }; 

	/** 
	 * Sets the a priori probability for a model. 
	 * @param model The model
	 * @param probability The a priori probability 
	 */ 
	void SetModelAPrioriProbability(double probability) 
	{ fModelAPriori = probability; }; 

	/** 
	 * Sets the a posteriori probability for a model. 
	 * @param model The model
	 * @param probability The a posteriori probability 
	 */ 
	void SetModelAPosterioriProbability(double probability) 
	{ fModelAPosteriori = probability; }; 

	/** 
	 * Sets the normalization of the likelihood. 
	 * The normalization is the integral of likelihood over all parameters. 
	 * @param norm The normalization of the likelihood 
	 */ 
	void SetNormalization(double norm)
	{ fNormalization = norm; };

	/**
	 * Sets the data set. 
	 * @param dataset A data set 
	 */ 
	void SetDataSet(BCDataSet* dataset)
	{ fDataSet = dataset; fNormalization = -1.0; };

	/**
	 * Sets a single data point as data set. 
	 * @param datapoint A data point 
	 */ 
	void SetSingleDataPoint(BCDataPoint * datapoint); 

	void SetSingleDataPoint(BCDataSet * dataset, int index); 

	/**
	 * Sets the minimum number of data points. 
	 */ 
	void SetNDataPointsMinimum(int minimum)
	{ fNDataPointsMinimum = minimum; }; 

	/**
	 * Sets the maximum number of data points. 
	 */ 
	void SetNDataPointsMaximum(int maximum)
	{ fNDataPointsMaximum = maximum; }; 

	void SetDataBoundaries(int index, double lowerboundary, double upperboundary); 

	/**
	 * Sets the error band flag to continuous function 
	 */ 
	void SetErrorBandContinuous(bool flag); 

	// methods 

	/** 
	 * Adds a parameter to the parameter set 
	 * @param name The name of the parameter 
	 * @param lowerlimit The lower limit of the parameter values
	 * @param upperlimit The upper limit of the parameter values 
	 * @see AddParameter(BCParameter* parameter); 
	 */ 
	int AddParameter(const char * name, double lowerlimit, double upperlimit); 

	/** 
	 * Adds a parameter to the model. 
	 * @param parameter A model parameter 
	 * @see AddParameter(const char * name, double lowerlimit, double upperlimit); 
	 */ 
	int AddParameter(BCParameter* parameter); 

	/** 
	 * Defines the parameters of the model 
	 * Method needs to be overloaded by the user. 
	 */ 
	virtual void DefineParameters() 
	{ ;}; 

	/** 
	 * Returns the prior probability. 
	 * @param parameters A set of parameter values 
	 * @return The prior probability p(parameters) 
	 * @see GetPrior(std::vector <double> parameters) 
	 */ 
	double APrioriProbability(std::vector <double> parameters) 
	{ return exp( this->LogAPrioriProbability(parameters) ); };

	/** 
	 * Returns natural logarithm of the prior probability. 
	 * Method needs to be overloaded by the user. 
	 * @param parameters A set of parameter values 
	 * @return The prior probability p(parameters) 
	 * @see GetPrior(std::vector <double> parameters) 
	 */ 
	virtual double LogAPrioriProbability(std::vector <double> parameters) 
	{ return 1.0; };

	/** 
	 * Returns the likelihood 
	 * @param parameters A set of parameter values 
	 * @return The likelihood
	 */ 
	double Likelihood(std::vector <double> parameter)
	{ return exp( this->LogLikelihood(parameter) ); };

	/** 
	 * Calculates natural logarithm of the likelihood.
	 * Method needs to be overloaded by the user.
	 * @param parameters A set of parameter values 
	 * @return Natural logarithm of the likelihood
	 */ 
	virtual double LogLikelihood(std::vector <double> parameter); 

	/** 
	 * Returns the likelihood times prior probability given a set of parameter values 
	 * @param parameters A set of parameter values 
	 * @return The likelihood times prior probability 
	 */ 
	double ProbabilityNN(std::vector <double> parameter)
	{ return exp( this->LogProbabilityNN(parameter) ); };

	/** 
	 * Returns the natural logarithm of likelihood times prior probability given
	 * a set of parameter values 
	 * @param parameters A set of parameter values 
	 * @return The likelihood times prior probability 
	 */ 
	double LogProbabilityNN(std::vector <double> parameter); 

	/** 
	 * Returns the a posteriori probability given a set of parameter values 
	 * @param parameters A set of parameter values 
	 * @return The a posteriori probability 
	 */ 
	double Probability(std::vector <double> parameter)
	{ return exp( this->LogProbability(parameter) ); };

	/** 
	 * Returns natural logarithm of the  a posteriori probability given a set of parameter values 
	 * @param parameters A set of parameter values 
	 * @return The a posteriori probability 
	 */
	double LogProbability(std::vector <double> parameter); 

	/**
	 * Returns a conditional probability. 
	 * Method needs to be overloaded by the user. 
	 * @param datapoint A data point 
	 * @param parameters A set of parameter values
	 * @return The conditional probability p(datapoint|parameters) 
	 * @see GetConditionalEntry(BCDataPoint* datapoint, std::vector <double> parameters) 
	 */
	double ConditionalProbabilityEntry(BCDataPoint * datapoint, std::vector <double> parameters)
	{ return exp( this->LogConditionalProbabilityEntry(datapoint, parameters) ); };

	/**
	 * Returns a natural logarithm of conditional probability. 
	 * Method needs to be overloaded by the user. 
	 * @param datapoint A data point 
	 * @param parameters A set of parameter values
	 * @return The conditional probability p(datapoint|parameters) 
	 * @see GetConditionalEntry(BCDataPoint* datapoint, std::vector <double> parameters) 
	 */
	virtual double LogConditionalProbabilityEntry(BCDataPoint * datapoint, std::vector <double> parameters)
	{ flag_ConditionalProbabilityEntry = false; return 0.; };

	/**
	 * Returns a probability for the data set container.
	 * @param parameters A set of parameter values
	 * @return The probability for the data set container
	 * @see GetPoisson(std::vector <double> parameters)
	 */
	double PoissonProbability(int nentries, std::vector <double> parameters) 
	{ return exp( this->LogPoissonProbability(nentries, parameters) ); }; 

	/** 
	 * Returns a probability for the data set container. 
	 * Method to be overloaded by the user. 
	 * @param parameters A set of parameter values 
	 * @return The probability for the data set container
	 * @see GetPoisson(std::vector <double> parameters) 
	 */ 
	virtual double LogPoissonProbability(int nentries, std::vector <double> parameters) 
	{ return .0; }; 

	/**
	 * Returns upper and lower value for a given point x on the error band 
	 * param x A vector of x-values 
	 * param ymin The minimum value 
	 * param ymax The maximum value 
	 **/
	void CalculateErrorBandXY(int nx, double xmin, double xmax, int ny, double ymin, double ymax, int niter); 

	/**
	 * Sampling function used for importance sampling. 
	 * Method needs to be overloaded by the user. 
	 * @param parameters A set of parameter values 
	 * @return The probability density at the parameter values 
	 */
	virtual double SamplingFunction(std::vector <double> parameters); 

	/** 
	 * Overloaded function to evaluate integral. 
	 */ 
	double Eval(std::vector <double> parameters)
	{ return exp( this->LogEval(parameters) ); };

	/** 
	 * Overloaded function to evaluate integral. 
	 */ 
  double LogEval(std::vector <double> parameters);

	/**
	 * Overloaded function to evaluate integral.
	 */
	double EvalSampling(std::vector <double> parameters);

	/**
	 * Integrates over the un-normalized probability and updates fNormalization.
	 */
	double Normalize();

	/**
	 * Checks if a set of parameters values is within the given range.
	 * @param parameters A set of parameter values
	 * @return Error code (0: OK, -1 length of parameters not correct, -2 values not within range)
	 */
	int CheckParameters(std::vector <double> parameters);

	/**
	 * Does the mode finding
	 */
	void FindMode();

	/**
	 * Marginalizes the a posteriori probability with respect to a parameter.
	 * @param parameter A model parameter
	 * @return 1D marginalized probability
	 */
	BCH1D * MarginalizeProbability(BCParameter* parameter);

	BCH1D * MarginalizeProbability(const char * name) 
		{ return this -> MarginalizeProbability(this -> GetParameter(name)); }; 

	/** 
	 * Marginalizes the a posteriori probability with respect to two parameters. 
	 * @param parameter1 First parameter
	 * @param parameter2 Second parameter
	 * @return 2D marginalized probability
	 */ 
	BCH2D * MarginalizeProbability(BCParameter * parameter1, BCParameter * parameter2);

	BCH2D * MarginalizeProbability(const char * name1, const char * name2)
		{ return this -> MarginalizeProbability(this -> GetParameter(name1), this -> GetParameter(name2)); };

	/**
	 * Marginalize all probabilities wrt. single parameters and all combinations
	 * of two parameters. The individual distributions can be retrieved using
	 * the GetMarginalized method.
	 * @return Total number of marginalized distributions
	 */
	int MarginalizeAll();

	/**
	 * If MarginalizeAll method was used, the individual marginalized distributions
	 * with respect to one parameter can be retrieved using this method.
	 * @param parameter Model parameter
	 * @return 1D marginalized probability
	 */
	BCH1D * GetMarginalized(BCParameter * parameter);

	BCH1D * GetMarginalized(const char * name)
		{ return this -> GetMarginalized(this -> GetParameter(name)); };

	/**
	 * If MarginalizeAll method was used, the individual marginalized distributions
	 * with respect to otwo parameters can be retrieved using this method.
	 * @param parameter1 First parameter
	 * @param parameter2 Second parameter
	 * @return 2D marginalized probability
	 */
	BCH2D * GetMarginalized(BCParameter * parameter1, BCParameter * parameter2);

	BCH2D * GetMarginalized(const char * name1, const char * name2)
		{ return this -> GetMarginalized(this -> GetParameter(name1), this -> GetParameter(name2)); };

	/**
	 * Creates data sets given a set of parameters. 
	 * @param ndatasets The number of data sets to be created 
	 * @param parameters A set of parameter values 
	 */ 
	void CreateData(int ndatasets, std::vector <double> parameters); 

	/**
	 * Creates data sets in a grid given a set of parameters. 
	 * @param ndatasets The number of data sets to be created 
	 * @param parameters A set of parameter values 
	 * @param grid Boolean for random (false) or grid values (true) 
	 * @param limits Limits for each data value 
	 */ 
	void CreateDataGrid(int ndatasets, std::vector <double> parameters, std::vector <bool> grid, std::vector <double> limits); 


	/**
	 * Creates data sets in a grid given a set of parameters. 
	 * @param ndatasets The number of data sets to be created 
	 * @param parameters A set of parameter values 
	 * @param grid Boolean for random (false) or grid values (true) 
	 * @param limits Limits for each data value 
	 */ 
	void CreateDataGridROOT(int ndatasets, std::vector <double> parameters, std::vector <bool> grid, std::vector <double> limits); 

	/**
	 * Constrains a data point
	 * @param x A vector of double
	 */
	virtual void CorrelateDataPointValues(vector<double> &x);

	/** 
	 * Goodness-of-fit test. 
	 * Assuming a certain set of parameters this function reads data from a file (which was created under 
	 * the assumption of the model and the same set of parameters) and calculates the frequency distribution 
	 * of the probability p(data|parameters). 
	 * @param filename A file which contains a list of files with ensembles 
	 * @param parameters The parameter values for which the probability is calculated
	 * @return A frequency distribution of the conditional probability for the data sets given the parameters 
	 */ 
	BCH1D * GoodnessOfFitTest(const char * filenname, std::vector <double> parameters);

	/** 
	 * Do goodness-of-fit test. 
	 * Creates data sets and performs a goodness-of-fit test. 
	 * @param ndatasets The number of data sets to be created 
	 * @param parameters A set of parameter values 
	 * @param grid Boolean for random (false) or grid values (true) 
	 * @param limits Limits for each data value 
	 * @see CreateData(int ndatasets, std::vector <double> parameters)
	 * @see CreateDataGrid(int ndatasets, std::vector <double> parameters, std::vector <bool> grid, std::vector <double> limits)
	 * @see GoodnessOfFitTest(const char * filenname, std::vector <double> parameters)
 */ 
	BCH1D * DoGoodnessOfFitTest(int ndatasets, std::vector<double> parameters, std::vector <bool> grid, std::vector <double> limits);
	BCH1D * DoGoodnessOfFitTest(int ndatasets, std::vector<double> parameters);
	BCH1D * DoGoodnessOfFitTest(int ndatasets);
	BCH1D * DoGoodnessOfFitTest(const char * filename);
	BCH1D * DoGoodnessOfFitTest(const char * filename, std::vector<double> parameters); 
	BCH1D * DoGoodnessOfFitTestROOT(int ndatasets, std::vector<double> parameters, std::vector <bool> grid, std::vector <double> limits);
	BCH1D * GoodnessOfFitTestROOT(int ndatasets, const char * filename, std::vector <double> parameters); 

	/*
	 * @return The p-value 
	 */ 
	double GetPValue() 
	{ return fPValue; }; 

	/**
	 * Calculates the matrix element of the Hessian matrix
	 * @param parameter1 The parameter for the first derivative
	 * @param parameter2 The parameter for the first derivative
	 * @return The matrix element of the Hessian matrix 
	 */ 
	double HessianMatrixElement(BCParameter * parameter1, BCParameter * parameter2, std::vector<double> point); 

	/** 
	 * Prints a summary on the screen. 
	 */ 
	void PrintSummary(); 

	/** 
	 * Prints matrix elements of the Hessian matrix 
	 * @param parameters The parameter values at which point to evaluate the matrix 
	 */
	void PrintHessianMatrix(std::vector<double> parameters); 

 protected:

	/** 
	 * Index of the model.
	 */ 
	int fIndex; 

	/** 
	 * Name of the model. 
	 */
	std::string fName; 

	/**
	 * The model prior probability. 
	 */
	double fModelAPriori;

	/** 
	 * The model a posteriori probability. 
	 */ 
	double fModelAPosteriori; 

	/** 
	 * A model parameter container. 
	 */ 
	BCParameterSet * fParameterSet;

	/** 
	 * A data set 
	 */ 
	BCDataSet * fDataSet;
  
	/** 
	 * Minimum number of data points 
	 */ 
	int fNDataPointsMinimum; 

	/**
	 * Maximum number of data points
	 */ 
	int fNDataPointsMaximum; 

	/*
	 * A flag for overloading ConditionalProbabilityEntry
	 */ 
	bool flag_ConditionalProbabilityEntry;

 private:

	/** 
	 * Converts a vector of doubles into a BCDataPoint 
	 */
	BCDataPoint * VectorToDataPoint(std::vector<double> data); 

	/** 
	 * Compares to strings
	 */ 
	int CompareStrings(const char * string1, const char * string2); 

	/** 
	 * The Likelihood normalization. 
	 */ 
	double fNormalization; 

	/**
	 * The p-value 
	 */ 
	double fPValue; 

}; 

// --------------------------------------------------------- 

typedef std::vector<BCModel*> BCModelContainer; 

// --------------------------------------------------------- 

#endif 
