#ifndef __BCMODEL__H
#define __BCMODEL__H

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
 * Copyright (C) 2008-2010, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>
#include <string>

#include "BAT/BCIntegrate.h"

// ROOT classes
class TH2D;
class TGraph;
class TCanvas;
class TPostscript;
class TF1;

//BAT classes
class BCDataPoint;
class BCDataSet;
class BCParameter;
class BCH1D;
class BCH2D;

const int MAXNDATAPOINTVALUES = 20;

// ---------------------------------------------------------

class BCModel : public BCIntegrate
{

	public:

		/** \name Constructors and destructors */
		/* @{ */

		/**
		 * The default constructor. */
		BCModel();

		/**
		 * A constructor.
		 * @param name The name of the model */
		BCModel(const char * name);

		/**
		 * The default destructor. */
		virtual ~BCModel();

		/* @} */

		/** \name Member functions (get) */
		/* @{ */

		/**
		 * @return The name of the model. */
		std::string GetName()
			{ return fName; };

		/**
		 * @return The index of the model. */
		int GetIndex()
			{ return fIndex; };

		/**
		 * @return The a priori probability. */
		double GetModelAPrioriProbability()
			{ return fModelAPriori; };

		/**
		 * @return The a posteriori probability. */
		double GetModelAPosterioriProbability()
			{ return fModelAPosteriori; };

		/**
		 * @return The normalization factor of the probability */
		double GetNormalization()
			{ return fNormalization; };

		/**
		 * @return The data set. */
		BCDataSet* GetDataSet()
			{ return fDataSet; };

		/**
		 * @return The lower boundaries of possible data values. */
		BCDataPoint* GetDataPointLowerBoundaries()
			{ return fDataPointLowerBoundaries; };

		/**
		 * @return The upper boundaries of possible data values. */
		BCDataPoint* GetDataPointUpperBoundaries()
			{ return fDataPointUpperBoundaries; };

		/**
		 * @param index The index of the variable.
		 * @return The lower boundary of possible data values for a particular variable. */
		double GetDataPointLowerBoundary(unsigned int index)
			{ return fDataPointLowerBoundaries -> GetValue(index); };

		/**
		 * @param index The index of the variable.
		 * @return The upper boundary of possible data values for a particular variable. */
		double GetDataPointUpperBoundary(unsigned int index)
			{ return fDataPointUpperBoundaries -> GetValue(index); };

		/*
		 * Checks if the boundaries have been defined
		 * @return true, if the boundaries have been set, false otherwise */
		bool GetFlagBoundaries();

		/**
		 * @return The number of data points in the current data set. */
		int GetNDataPoints();

		/**
		 * @param index The index of the data point.
		 * @return The data point in the current data set at index */
		BCDataPoint * GetDataPoint(unsigned int index);

		/**
		 * @return The minimum number of data points. */
		unsigned int GetNDataPointsMinimum()
			{ return fNDataPointsMinimum; };

		/**
		 * @return The maximum number of data points. */
		unsigned int GetNDataPointsMaximum()
			{ return fNDataPointsMaximum; };

		/**
		 * @return The number of parameters of the model. */
		unsigned int GetNParameters()
			{ return fParameterSet ? fParameterSet -> size() : 0; };

		/**
		 * @param index The index of the parameter in the parameter set.
		 * @return The parameter. */
		BCParameter * GetParameter(int index);

		/**
		 * @param name The name of the parameter in the parameter set.
		 * @return The parameter. */
		BCParameter * GetParameter(const char * name);

		/**
		 * @return parameter set */
		BCParameterSet * GetParameterSet()
			{ return fParameterSet; };

		/**
		 * Returns the value of a parameter (defined by index) at
		 * the global mode of the posterior pdf.
		 * @param index index of the parameter.
		 * @return best fit value of the parameter or -1e+111 on error or center of the range if mode finding not yer run */
		double GetBestFitParameter(unsigned int index);

		/**
		 * Returns the error on the value of a parameter (defined by index) at
		 * the global mode of the posterior pdf.
		 * @param index index of the parameter.
		 * @return error on the best fit value of the parameter or -1 if undefined */
		double GetBestFitParameterError(unsigned int index);

		/**
		 * Returns the set of values of the parameters at the global mode of
		 * the posterior pdf.
		 * @return The best fit parameters */
		std::vector <double> GetBestFitParameters()
			{ return fBestFitParameters; };
		
		std::vector <double> GetBestFitParameterErrors() 
			{ return fBestFitParameterErrors; }; 

		/**
		 * Returns the value of a particular parameter (defined by index) at
		 * the modes of the marginalized posterior pdfs.
		 * @param index index of the parameter.
		 * @return best fit parameter or -1e+111 on error or center of the range if marginalization not yet run */
		double GetBestFitParameterMarginalized(unsigned int index);

		/**
		 * Returns the set of values of the parameters at the modes of the
		 * marginalized posterior pdfs.
		 * @return best fit parameters */
		std::vector <double> GetBestFitParametersMarginalized()
			{ return fBestFitParametersMarginalized; };                                                          

		/**
		 * @return The 2-d histogram of the error band. */
		TH2D * GetErrorBandXY()
			{ return fErrorBandXY; };

		TH2D * GetErrorBandXY_yellow(double level=.68, int nsmooth=0);

		/**
		 * Returns a vector of y-values at a certain probability level.
		 * @param level The level of probability
		 * @return vector of y-values */
		std::vector <double> GetErrorBand(double level);

		TGraph * GetErrorBandGraph(double level1, double level2);

		TGraph * GetFitFunctionGraph(std::vector <double> parameters);

		TGraph * GetFitFunctionGraph()
			{ return this -> GetFitFunctionGraph(this -> GetBestFitParameters()); };

		TGraph * GetFitFunctionGraph(std::vector <double> parameters, double xmin, double xmax, int n=1000);

		bool GetFixedDataAxis(unsigned int index);

		/* @} */

		/** \name Member functions (set) */
		/* @{ */

		/**
		 * Sets the name of the model.
		 * @param name Name of the model */
		void SetName(const char * name)
			{ fName = name; };

		/**
		 * Sets the index of the model within the BCModelManager.
		 * @param index The index of the model */
		void SetIndex(int index)
			{ fIndex = index; };

		/**
		 * Set all parameters of the model using a BCParameterSet container.
		 * @par pointer to parameter set */
		void SetParameterSet( BCParameterSet * parset )
			{ fParameterSet = parset; };

		/** 
		 * Set the range of a parameter
		 * @param index The parameter index
		 * @param parmin The parameter minimum
		 * @param parmax The parameter maximum
		 * @return An error code. */
		int SetParameterRange(int index, double parmin, double parmax);

		/**
		 * Sets the a priori probability for a model.
		 * @param model The model
		 * @param probability The a priori probability */
		void SetModelAPrioriProbability(double probability)
			{ fModelAPriori = probability; };

		/**
		 * Sets the a posteriori probability for a model.
		 * @param model The model
		 * @param probability The a posteriori probability */
		void SetModelAPosterioriProbability(double probability)
			{ fModelAPosteriori = probability; };

		/**
		 * Sets the normalization of the likelihood.
		 * The normalization is the integral of likelihood over all parameters.
		 * @param norm The normalization of the likelihood */
		void SetNormalization(double norm)
			{ fNormalization = norm; };

		/**
		 * Sets the data set.
		 * @param dataset A data set */
		void SetDataSet(BCDataSet* dataset)
			{ fDataSet = dataset; fNormalization = -1.0; };

		/**
		 * Sets a single data point as data set.
		 * @param datapoint A data point */
		void SetSingleDataPoint(BCDataPoint * datapoint);

		void SetSingleDataPoint(BCDataSet * dataset, unsigned int index);

		/**
		 * Sets the minimum number of data points. */
		void SetNDataPointsMinimum(unsigned int minimum)
			{ fNDataPointsMinimum = minimum; };

		/**
		 * Sets the maximum number of data points. */
		void SetNDataPointsMaximum(unsigned int maximum)
			{ fNDataPointsMaximum = maximum; };

		void SetDataBoundaries(unsigned int index, double lowerboundary, double upperboundary, bool fixed=false);

		/**
		 * Sets the error band flag to continuous function */
		void SetErrorBandContinuous(bool flag);

		/**
		 * Set the number of bins for the marginalized distribution of a parameter.
		 * @param parname The name of the parameter in the parameter set
		 * @param nbins   Number of bins (default = 100) */
		void SetNbins(const char * parname, int nbins);

		/**
		 * Set prior for a parameter. 
		 * @param index The parameter index
		 * @param f A pointer to a function describing the prior
		 * @return An error code.
		 */ 
		int SetPrior(int index, TF1* f); 

		/**
		 * Set prior for a parameter. 
		 * @param name The parameter name
		 * @param f A pointer to a function describing the prior
		 * @return An error code.
		 */ 
		int SetPrior(const char* name, TF1* f); 

		/**
		 * Set Gaussian prior for a parameter. 
		 * @param index The parameter index
		 * @param mean The mean of the Gaussian
		 * @param sigma The sigma of the Gaussian
		 * @return An error code.
		 */ 
		int SetPriorGauss(int index, double mean, double sigma); 

		/**
		 * Set Gaussian prior for a parameter. 
		 * @param name The parameter name
		 * @param mean The mean of the Gaussian
		 * @param sigma The sigma of the Gaussian
		 * @return An error code.
		 */ 
		int SetPriorGauss(const char* name, double mean, double sigma); 

		/**
		 * Set Gaussian prior for a parameter with two different widths.
		 * @param index The parameter index
		 * @param mean The mean of the Gaussian
		 * @param sigmadown The sigma (down) of the Gaussian
		 * @param sigmaup The sigma (up)of the Gaussian
		 * @return An error code.
		 */ 
		int SetPriorGauss(int index, double mean, double sigmadown, double sigmaup);

		/**
		 * Set Gaussian prior for a parameter with two different widths.
		 * @param name The parameter name
		 * @param mean The mean of the Gaussian
		 * @param sigmadown The sigma (down) of the Gaussian
		 * @param sigmaup The sigma (up)of the Gaussian
		 * @return An error code.
		 */ 
		int SetPriorGauss(const char* name, double mean, double sigmadown, double sigmaup);

		/**
		 * Set constant prior for this parameter
		 * @param index the index of the parameter
		 * @return An error code
		 */
		int SetPriorConstant(int index);

      /**
       * Set constant prior for this parameter
       * @param name the name of the parameter
       * @return An error code
       */
		int SetPriorConstant(const char* name);

      /**
       * Enable caching the constant value of the prior, so LogAPrioriProbability
       * is called only once. Note that the prior for ALL parameters is
       * assumed to be constant. The value is computed from
       * the parameter ranges, so make sure these are defined before this method is
       * called.
       * @return An error code
       */
      int SetPriorConstantAll();

		/* @} */

		/** \name Member functions (miscellaneous methods) */
		/* @{ */

		/**
		 * Adds a parameter to the parameter set
		 * @param name The name of the parameter
		 * @param lowerlimit The lower limit of the parameter values
		 * @param upperlimit The upper limit of the parameter values
		 * @see AddParameter(BCParameter* parameter); */
		int AddParameter(const char * name, double lowerlimit, double upperlimit);

		/**
		 * Adds a parameter to the model.
		 * @param parameter A model parameter
		 * @see AddParameter(const char * name, double lowerlimit, double upperlimit); */
		int AddParameter(BCParameter* parameter);

		/**
		 * Returns the prior probability.
		 * @param parameters A set of parameter values
		 * @return The prior probability p(parameters)
		 * @see GetPrior(std::vector <double> parameters) */
		double APrioriProbability(std::vector <double> parameters)
			{ return exp( this->LogAPrioriProbability(parameters) ); };

		/**
		 * Returns natural logarithm of the prior probability.
		 * Method needs to be overloaded by the user.
		 * @param parameters A set of parameter values
		 * @return The prior probability p(parameters)
		 * @see GetPrior(std::vector <double> parameters) */
		virtual double LogAPrioriProbability(std::vector <double> parameters); 

		/**
		 * Returns the likelihood
		 * @param parameters A set of parameter values
		 * @return The likelihood */
		double Likelihood(std::vector <double> parameter)
			{ return exp( this->LogLikelihood(parameter) ); };

		/**
		 * Calculates natural logarithm of the likelihood.
		 * Method needs to be overloaded by the user.
		 * @param parameters A set of parameter values
		 * @return Natural logarithm of the likelihood */
		virtual double LogLikelihood(std::vector <double> parameter);

		/**
		 * Returns the likelihood times prior probability given a set of parameter values
		 * @param parameters A set of parameter values
		 * @return The likelihood times prior probability */
		double ProbabilityNN(std::vector <double> parameter)
			{ return exp( this->LogProbabilityNN(parameter) ); };

		/**
		 * Returns the natural logarithm of likelihood times prior probability given
		 * a set of parameter values
		 * @param parameters A set of parameter values
		 * @return The likelihood times prior probability */
		double LogProbabilityNN(std::vector <double> parameter);

		/**
		 * Returns the a posteriori probability given a set of parameter values
		 * @param parameters A set of parameter values
		 * @return The a posteriori probability */
		double Probability(std::vector <double> parameter)
			{ return exp( this->LogProbability(parameter) ); };

		/**
		 * Returns natural logarithm of the  a posteriori probability given a set of parameter values
		 * @param parameters A set of parameter values
		 * @return The a posteriori probability */
		double LogProbability(std::vector <double> parameter);

		/**
		 * Returns a conditional probability.
		 * Method needs to be overloaded by the user.
		 * @param datapoint A data point
		 * @param parameters A set of parameter values
		 * @return The conditional probability p(datapoint|parameters)
		 * @see GetConditionalEntry(BCDataPoint* datapoint, std::vector <double> parameters) */
		double ConditionalProbabilityEntry(BCDataPoint * datapoint, std::vector <double> parameters)
			{ return exp( this->LogConditionalProbabilityEntry(datapoint, parameters) ); };

		/**
		 * Returns a natural logarithm of conditional probability.
		 * Method needs to be overloaded by the user.
		 * @param datapoint A data point
		 * @param parameters A set of parameter values
		 * @return The conditional probability p(datapoint|parameters)
		 * @see GetConditionalEntry(BCDataPoint* datapoint, std::vector <double> parameters) */
		virtual double LogConditionalProbabilityEntry(BCDataPoint * /*datapoint*/, std::vector <double> /*parameters*/)
			{ flag_ConditionalProbabilityEntry = false; return 0.; };

		/**
		 * Sampling function used for importance sampling.
		 * Method needs to be overloaded by the user.
		 * @param parameters A set of parameter values
		 * @return The probability density at the parameter values */
		virtual double SamplingFunction(std::vector <double> parameters);

		/**
		 * Overloaded function to evaluate integral. */
		double Eval(std::vector <double> parameters)
			{ return exp( this->LogEval(parameters) ); };

		/**
		 * Overloaded function to evaluate integral. */
		double LogEval(std::vector <double> parameters);

		/**
		 * Overloaded function to evaluate integral. */
		double EvalSampling(std::vector <double> parameters);

		/**
		 * Integrates over the un-normalized probability and updates fNormalization. */
		double Normalize();

		/**
		 * Checks if a set of parameters values is within the given range.
		 * @param parameters A set of parameter values
		 * @return Error code (0: OK, -1 length of parameters not correct, -2 values not within range)
		 */
		int CheckParameters(std::vector <double> parameters);

		/**
		 * Do the mode finding using a method set via SetOptimizationMethod.
		 * Default is Minuit. The mode can be extracted using the GetBestFitParameters() method.
		 *
		 * A starting point for the mode finding can be specified for Minuit. If not
		 * specified, Minuit default will be used (center of the parameter space).
		 *
		 * If running mode finding after the MCMC it is a good idea to specify the
		 * mode obtained from MCMC as a starting point for the Minuit minimization.
		 * MCMC will find the location of the global mode and Minuit will
		 * converge to the mode precisely. The commands are:
			<pre>
			model -> MarginalizeAll();
			model -> FindMode( model -> GetBestFitParameters() );
			</pre>
		 * @start startinf point of Minuit minimization
		 * */
		void FindMode(std::vector<double> start = std::vector<double>(0));

		/**
		 * Does the mode finding using Minuit. If starting point is not specified,
		 * finding will start from the center of the parameter space.
		 * @param start point in parameter space from which the mode finding is started.
		 * @param printlevel The print level. */
		void FindModeMinuit(std::vector<double> start = std::vector<double>(0), int printlevel = 1);

		/**
		 * Write mode into file */
		void WriteMode(const char * file);

		/**
		 * Read mode from file created by WriteMode() call */
		int ReadMode(const char * file);

		/**
		 * Read */
		int ReadMarginalizedFromFile(const char * file);

		/**
		 * Read */
		int ReadErrorBandFromFile(const char * file);

		/**
		 * Marginalize all probabilities wrt. single parameters and all combinations
		 * of two parameters. The individual distributions can be retrieved using
		 * the GetMarginalized method.
		 * @return Total number of marginalized distributions */
		int MarginalizeAll();

		/**
		 * If MarginalizeAll method was used, the individual marginalized distributions
		 * with respect to one parameter can be retrieved using this method.
		 * @param parameter Model parameter
		 * @return 1D marginalized probability */
		BCH1D * GetMarginalized(BCParameter * parameter);

		BCH1D * GetMarginalized(const char * name)
			{ return this -> GetMarginalized(this -> GetParameter(name)); };

		/**
		 * If MarginalizeAll method was used, the individual marginalized distributions
		 * with respect to otwo parameters can be retrieved using this method.
		 * @param parameter1 First parameter
		 * @param parameter2 Second parameter
		 * @return 2D marginalized probability */
		BCH2D * GetMarginalized(BCParameter * parameter1, BCParameter * parameter2);

		BCH2D * GetMarginalized(const char * name1, const char * name2)
			{ return this -> GetMarginalized(this -> GetParameter(name1), this -> GetParameter(name2)); };

		/**
		 *   */
		int PrintAllMarginalized1D(const char * filebase);
		int PrintAllMarginalized2D(const char * filebase);
		int PrintAllMarginalized(const char * file, unsigned int hdiv=1, unsigned int ndiv=1);

		/**
		 * Constrains a data point
		 * @param x A vector of double */
		virtual void CorrelateDataPointValues(std::vector<double> &x);

		/**
		 * Calculate p-value from Chi2 distribution for Gaussian problems
		 * @param par Parameter set for the calculation of the likelihood
		 * @param sigma_index Index of the sigma/uncertainty for the data points
		 *        (for data in format "x y erry" the index would be 2) */
		double GetPvalueFromChi2(std::vector<double> par, int sigma_index);

		/**
		 * Calculate p-value from asymptotic Chi2 distribution for arbitrary problems
		 * using the definition (3) from
		 * Johnson, V.E. A Bayesian chi2 Test for Goodness-of-Fit. The Annals of Statistics 32, 2361-2384(2004).
		 *
		 * @param par Parameter set for the calculation of the likelihood */
		double GetPvalueFromChi2Johnson(std::vector<double> par);

      /**
       * Calculate p-value from Kolmogorov-Smirnov test statistic
       * for 1D - datasets.
       *
       * @param par Parameter set for the calculation of the likelihood
       * @param index Index of the data point in the BCDataSet
       *        (for data in format "x y erry" the index would be 1) */
      double GetPvalueFromKolmogorov(const std::vector<double>& par, int index);


		/**
		 * Calculate  Chi2  (also called R^{B}) for arbitrary problems with binned data
		 * using the definition (3) from
		 * Johnson, V.E. A Bayesian chi2 Test for Goodness-of-Fit. The Annals of Statistics 32, 2361-2384(2004).
		 *
		 * @param par Parameter set for the calculation of the likelihood
		 * @param nBins how many bins to use for the data, for negative an adapted rule \
		 * of thumb by  Wald(1942) is used, with at least three bins*/
		double GetChi2Johnson(std::vector<double> par, const int nBins=-1);

		/**
		 * Calculate the A-value, a summary statistic. It computes the frequency
		 * that a Chi2 value determined from the data by Johnson's binning prescription is
		 * larger than a value sampled from the reference chi2 distribution. They out
		 * from one chain is used. A=1/2 provides
		 * no evidence against the null hypothesis. Compare
		 * Johnson, V.E. A Bayesian chi2 Test for Goodness-of-Fit. The Annals of Statistics 32, 2361-2384(2004).
		 *
		 * @param par tree contains the samples of posterior of the parameters
		 * @param par histogram filled by function with distribution of p values*/
		double GetAvalueFromChi2Johnson(TTree* tree, TH1D* distribution=0);

		double GetPvalueFromChi2NDoF(std::vector<double> par, int sigma_index);

		BCH1D * CalculatePValue(std::vector<double> par, bool flag_histogram=false);

		/*
		 * @return The p-value */
		double GetPValue()
			{ return fPValue; };

		double GetPValueNDoF()
			{ return fPValueNDoF; };

		double GetChi2NDoF()
			{ return fChi2NDoF; };

		/**
       * For a Gaussian problem, calculate the chi2 of the longest run of consecutive
       * values above/below the expected values
       * @param dataIndex component of datapoint with the observed value
       * @param sigmaIndex component of datapoint with uncertainty */
		std::vector<double> GetChi2Runs(int dataIndex, int sigmaIndex);

		/*
		 * Set maximum number of iterations in the MCMC pre-run of the p-value
		 * evaluation using MCMC */
		void SetGoFNIterationsMax(int n)
			{ fGoFNIterationsMax=n; };

		/*
		 * Set number of iterations in the MCMC normal run of the p-value
		 * evaluation using MCMC */
		void SetGoFNIterationsRun(int n)
			{ fGoFNIterationsRun=n; };

		/*
		 * Set number of chains in the MCMC of the p-value
		 * evaluation using MCMC */
		void SetGoFNChains(int n)
			{ fGoFNChains=n; };

		/**
		 * Calculates the matrix element of the Hessian matrix
		 * @param parameter1 The parameter for the first derivative
		 * @param parameter2 The parameter for the first derivative
		 * @return The matrix element of the Hessian matrix */
		double HessianMatrixElement(BCParameter * parameter1, BCParameter * parameter2, std::vector<double> point);

		/**
		 * Prints a summary on the screen. */
		void PrintSummary();

		/**
		 * Prints a summary of the Markov Chain Monte Carlo to a file. */
		void PrintResults(const char * file);

		/**
		 * Prints a short summary of the fit results on the screen. */
		void PrintShortFitSummary(int chi2flag=0);

		/**
		 * Prints matrix elements of the Hessian matrix
		 * @param parameters The parameter values at which point to evaluate the matrix */
		void PrintHessianMatrix(std::vector<double> parameters);

		void FixDataAxis(unsigned int index, bool fixed);


		/**
		 * 1dim cumulative distribution function of the probability
		 * to get the data f(x_i|param) for a single measurement, assumed to
		 * be of identical functional form for all measurements
		 * @param parameters The parameter values at which point to compute the cdf
		 * @param index The data point index starting at 0,1...N-1
		 * @param lower only needed for discrete distributions!
		 * Return the CDF for the count one less than actually observed, e.g.
		 * in Poisson process, if 3 actually observed, then CDF(2) is returned */
		virtual double CDF(const std::vector<double>& /*parameters*/,  int /*index*/, bool /*lower=false*/)
		{return 0.0;}


		/** 
		 * Reset all results. 
		 * @return An error code. */ 
		int ResetResults();

	/* @} */

	protected:

		/**
		 * Index of the model. */
		int fIndex;

		/**
		 * Name of the model. */
		std::string fName;

		/**
		 * The model prior probability. */
		double fModelAPriori;

		/**
		 * The model a posteriori probability. */
		double fModelAPosteriori;

		/**
		 * A model parameter container. */
		BCParameterSet * fParameterSet;

		/**
		 * A data set */
		BCDataSet * fDataSet;

		/**
		 * Minimum number of data points */
		unsigned int fNDataPointsMinimum;

		/**
		 * Maximum number of data points */
		unsigned int fNDataPointsMaximum;

		/**
		 * A flag for overloading ConditionalProbabilityEntry */
		bool flag_ConditionalProbabilityEntry;

		/**
		 * The p-value */
		double fPValue;

		double fChi2NDoF;
		double fPValueNDoF;

		/**
		* true for a discrete probability, false for continuous pdf  */
		bool flag_discrete;

		/*
		 * Maximum number of iterations in the MCMC pre-run of the p-value
		 * evaluation using MCMC */
		int fGoFNIterationsMax;

		/*
		 * Number of iterations in the MCMC normal run of the p-value
		 * evaluation using MCMC */
		int fGoFNIterationsRun;

		/*
		 * Number of chains in the MCMC of the p-value
		 * evaluation using MCMC */
		int fGoFNChains;

		/*
		 * A vector of prior functions. */ 
		std::vector<TF1*> fPriorContainer;

		/**
		 * Flag to indicate that all parameters have constant prior.
		 */
		bool fPriorConstantAll;

		/**
		 * The value of the product of constant priors of
		 * individual parameters.
		 */
		double fPriorConstantValue;

		/**
		 * List the parameters whose prior is constant
		 */
		std::vector<bool> fPriorContainerConstant;

	private:

		/**
		 * Converts a vector of doubles into a BCDataPoint */
		BCDataPoint * VectorToDataPoint(std::vector<double> data);

		/**
		 * Compares to strings */
		int CompareStrings(const char * string1, const char * string2);

		/**
		 * The Likelihood normalization. */
		double fNormalization;

		/**
		 * rule of thumb for good number of bins (Wald1942, Johnson2004) to group observations
		 * updated so minimum is three bins (for 1-5 observations)!
		 * @param */

		int NumberBins()
			{ return (int)(exp(0.4 * log(this -> GetNDataPoints())) + 2); }

};

// ---------------------------------------------------------

typedef std::vector<BCModel*> BCModelContainer;

// ---------------------------------------------------------

#endif
