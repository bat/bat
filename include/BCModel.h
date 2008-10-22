/*!
 * \class BCModel
 * \brief The base class for all user-defined models.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a model. It contains a container of
 * parameters, their prior distributions and the conditional
 * probabilities given those parameters.  The methods which implement
 * the prior and conditional probabilities have to be overloaded by
 * the user in the user defined model class which will inherit from
 * this class.
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#ifndef __BCMODEL__H
#define __BCMODEL__H

#include <vector>
#include <string>

#include <TROOT.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPostScript.h>

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

		/** \name Constructors and destructors */
		/* @{ */

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

		/* @} */

		/** \name Member functions (get) */
		/* @{ */

		/**
		 * @return The name of the model.
		 */
		std::string GetName()
			{ return fName; };

		/**
		 * @return The index of the model.
		 */
		int GetIndex()
			{ return fIndex; };

		/**
		 * @return The a priori probability.
		 */
		double GetModelAPrioriProbability()
			{ return fModelAPriori; };

		/**
		 * @return The a posteriori probability.
		 */
		double GetModelAPosterioriProbability()
			{ return fModelAPosteriori; };

		/**
		 * @return The normalization factor of the probability
		 */
		double GetNormalization()
			{ return fNormalization; };

		/**
		 * @return The data set.
		 */
		BCDataSet* GetDataSet()
			{ return fDataSet; };

		/**
		 * @return The lower boundaries of possible data values.
		 */
		BCDataPoint* GetDataPointLowerBoundaries()
			{ return fDataPointLowerBoundaries; };

		/**
		 * @return The upper boundaries of possible data values.
		 */
		BCDataPoint* GetDataPointUpperBoundaries()
			{ return fDataPointUpperBoundaries; };

		/**
		 * @param index The index of the variable.
		 * @return The lower boundary of possible data values for a particular variable.
		 */
		double GetDataPointLowerBoundary(int index)
			{ return fDataPointLowerBoundaries -> GetValue(index); };

		/**
		 * @param index The index of the variable.
		 * @return The upper boundary of possible data values for a particular variable.
		 */
		double GetDataPointUpperBoundary(int index)
			{ return fDataPointUpperBoundaries -> GetValue(index); };

		/*
		 * Checks if the boundaries have been defined
		 * @return true, if the boundaries have been set, false otherwise
		 */
		bool GetFlagBoundaries();

		/**
		 * @return The number of data points in the current data set.
		 */
		int GetNDataPoints();

		/**
		 * @param index The index of the data point.
		 * @return The data point in the current data set at index
		 */
		BCDataPoint * GetDataPoint(int index);

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
		 * @return The number of parameters of the model.
		 */
		int GetNParameters()
			{ return fParameterSet -> size(); };

		/**
		 * @param index The index of the parameter in the parameter set.
		 * @return The parameter.
		 */
		BCParameter* GetParameter(int index);

		/**
		 * @param name The name of the parameter in the parameter set.
		 * @return The parameter.
		 */
		BCParameter* GetParameter(const char * name);

		/**
		 * Returns the value of a particular parameter (defined by index) at
		 * the global mode of the posterior pdf.
		 * @param index The index of the parameter.
		 * @return The best fit parameter.
		 */
		double GetBestFitParameter(int index)
			{ return fBestFitParameters.at(index); };

		/**
		 * Returns the set of values of the parameters at the global mode of
		 * the posterior pdf.
		 * @return The best fit parameters
		 */
		std::vector <double> GetBestFitParameters()
			{ return fBestFitParameters; };

		/**
		 * Returns the value of a particular parameter (defined by index) at
		 * the modes of the marginalized posterior pdfs.
		 * @param index The index of the parameter.
		 * @return The best fit parameter
		 */
		double GetBestFitParameterMarginalized(int index)
			{ return fBestFitParametersMarginalized.at(index); };

		/**
		 * Returns the set of values of the parameters at the modes of the
		 * marginalized posterior pdfs.
		 * @return The best fit parameters.
		 */
		std::vector <double> GetBestFitParametersMarginalized()
			{ return fBestFitParametersMarginalized; };

		/**
		 * @return The 2-d histogram of the error band.
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

		TGraph * GetFitFunctionGraph(std::vector <double> parameters, double xmin, double xmax, int n=1000);

		bool GetFixedDataAxis(int index);

		/* @} */

		/** \name Member functions (set) */
		/* @{ */

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

		void SetDataBoundaries(int index, double lowerboundary, double upperboundary, bool fixed=false);

		/**
		 * Sets the error band flag to continuous function
		 */
		void SetErrorBandContinuous(bool flag);

		/* @} */

		/** \name Member functions (miscellaneous methods) */
		/* @{ */

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
		 * Write mode into file
		 */
		void WriteMode(const char * file);

		/**
		 * Read mode from file created by WriteMode() call
		 */
		int ReadMode(const char * file);

		/**
		 * Read
		 */
		int ReadMarginalizedFromFile(const char * file);

		/**
		 * Read
		 */
		int ReadErrorBandFromFile(const char * file);

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
		 *
		 */
		int PrintAllMarginalized1D(const char * filebase);
		int PrintAllMarginalized2D(const char * filebase);
		int PrintAllMarginalized(const char * file, int hdiv=1, int ndiv=1);

		/**
		 * Constrains a data point
		 * @param x A vector of double
		 */
		virtual void CorrelateDataPointValues(std::vector<double> &x);

		/**
		 * Calculate p-value from Chi2 distribution for Gaussian problems
		 * @param par Parameter set for the calculation of the likelihood
		 * @param sigma_index Index of the sigma/uncertainty for the data points
		 *        (for data in format "x y erry" the index would be 2)
		 */
		double GetPvalueFromChi2(std::vector<double> par, int sigma_index);

		BCH1D * CalculatePValue(std::vector<double> par, bool flag_histogram=false);

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
		 * Prints a summary of the Markov Chain Monte Carlo to a file.
		 */
		void PrintResults(const char * file);

		/**
		 * Prints matrix elements of the Hessian matrix
		 * @param parameters The parameter values at which point to evaluate the matrix
		 */
		void PrintHessianMatrix(std::vector<double> parameters);

		void FixDataAxis(int index, bool fixed);

		/* @} */

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

		/**
		 * A flag for overloading ConditionalProbabilityEntry
		 */
		bool flag_ConditionalProbabilityEntry;

		/**
		 * The p-value
		 */
		double fPValue;

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

};

// ---------------------------------------------------------

typedef std::vector<BCModel*> BCModelContainer;

// ---------------------------------------------------------

#endif
