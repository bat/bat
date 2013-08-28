#ifndef __BCMODEL__H
#define __BCMODEL__H

/*!
 * \class BCModel
 * \brief The base class for all user-defined models.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a model. It contains a container of prior distributions and the likelihood. The methods that implement the prior and the likelihood
 * have to be overloaded by the user in the user defined model class
 * derived from this class.
 */

/**
 * Copyright (C) 2008-2012, Daniel Kollar, Kevin Kroeninger, and Daniel Greenwald.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCIntegrate.h"

#include <string>

// ROOT classes
class TNamed;
class TH1;
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
      /** @{ */

      /**
       * A constructor.
       * @param name The name of the model */
      BCModel(const char * name="model");

      /**
       * The copy constructor. */
      BCModel(const BCModel & bcmodel);

      /**
       * The default destructor. */
      virtual ~BCModel();

      /** @} */
      /** \name Assignment operators */
      /** @{ */

      /**
       * Defaut assignment operator */
      BCModel & operator = (const BCModel & bcmodel);

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The name of the model. */
      const std::string & GetName() const
         { return fName; }

      /**
       * @return The a priori probability. */
      double GetModelAPrioriProbability() const
         { return fModelAPriori; }

      /**
       * @return The a posteriori probability. */
      double GetModelAPosterioriProbability() const
         { return fModelAPosteriori; }

      /**
       * @return The normalization factor of the probability */
      double GetNormalization() const
         { return fNormalization; }

      /**
       * @return The data set. */
      BCDataSet* GetDataSet() const
         { return fDataSet; }

      /**
       * @return The lower boundaries of possible data values. */
      BCDataPoint * GetDataPointLowerBoundaries() const
         { return fDataPointLowerBoundaries; }

      /**
       * @return The upper boundaries of possible data values. */
      BCDataPoint* GetDataPointUpperBoundaries() const
         { return fDataPointUpperBoundaries; }

      /**
       * @param index The index of the variable.
       * @return The lower boundary of possible data values for a particular variable. */
      double GetDataPointLowerBoundary(unsigned int index) const;

      /**
       * @param index The index of the variable.
       * @return The upper boundary of possible data values for a particular variable. */
      double GetDataPointUpperBoundary(unsigned int index) const;

      /**
       * Checks if the boundaries have been defined
       * @return true, if the boundaries have been set, false otherwise */
      bool GetFlagBoundaries() const;

      /**
       * @return The number of data points in the current data set. */
      unsigned GetNDataPoints() const;

      /**
       * @param index The index of the data point.
       * @return The data point in the current data set at index */
      BCDataPoint * GetDataPoint(unsigned int index) const;

      /**
      const BCParameter * GetParameter(const char * name);
       * @return The 2-d histogram of the error band. */
      TH2D * GetErrorBandXY() const
         { return fErrorBandXY; }

      TH2D * GetErrorBandXY_yellow(double level=.68, int nsmooth=0) const;

      /**
       * Returns a vector of y-values at a certain probability level.
       * @param level The level of probability
       * @return vector of y-values */
      std::vector<double> GetErrorBand(double level) const;

      TGraph * GetErrorBandGraph(double level1, double level2) const;

      TGraph * GetFitFunctionGraph(const std::vector<double> &parameters);

      TGraph * GetFitFunctionGraph()
         { return GetFitFunctionGraph(GetBestFitParameters()); }

      TGraph * GetFitFunctionGraph(const std::vector<double> &parameters, double xmin, double xmax, int n=1000);

      bool GetFixedDataAxis(unsigned int index) const;

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

      /**
       * Sets the name of the model.
       * @param name Name of the model */
      void SetName(const char * name)
         { fName = name; }

      /**
       * Sets the a priori probability for a model.
       * @param model The model
       * @param probability The a priori probability */
      void SetModelAPrioriProbability(double probability)
         { fModelAPriori = probability; }

      /**
       * Sets the a posteriori probability for a model.
       * @param model The model
       * @param probability The a posteriori probability */
      void SetModelAPosterioriProbability(double probability)
         { fModelAPosteriori = probability; }

      /**
       * Sets the normalization of the posterior. The normalization
       * is the integral of likelihood x prior over all parameters.
       * @param norm The normalization */
      void SetNormalization(double norm)
         { fNormalization = norm; }

      /**
       * Due to the combination of overloading and overriding, the method
       * AddParameter(const char * name, double min, double max)
       * from the base class is hidden to the user of BCModel, so explicitly tell the compiler that we want to use both
       */
      using BCEngineMCMC::AddParameter;

      /**
       * Adds a parameter to the model.
       * @param parameter A model parameter
       * @see AddParameter(const char * name, double lowerlimit, double upperlimit); */
      virtual int AddParameter(BCParameter * parameter);

      /**
       * Sets the data set.
       * @param dataset A data set */
      void SetDataSet(BCDataSet* dataset)
         { fDataSet = dataset; fNormalization = -1.0; }

      /**
       * Sets a single data point as data set.
       * @param datapoint A data point */
      void SetSingleDataPoint(BCDataPoint * datapoint);

      void SetSingleDataPoint(BCDataSet * dataset, unsigned int index);

      void SetDataBoundaries(unsigned int index, double lowerboundary, double upperboundary, bool fixed=false);

      /**
       * Sets the data point containing the lower boundaries of possible
       * data values */
      void SetDataPointLowerBoundaries(BCDataPoint * datasetlowerboundaries)
      { fDataPointLowerBoundaries = datasetlowerboundaries; }

      /**
       * Sets the data point containing the upper boundaries of possible
       * data values */
      void SetDataPointUpperBoundaries(BCDataPoint * datasetupperboundaries)
      { fDataPointUpperBoundaries = datasetupperboundaries; }

      /**
       * Sets the lower boundary of possible data values for a particular
       * variable */
      void SetDataPointLowerBoundary(int index, double lowerboundary);

      /**
       * Sets the upper boundary of possible data values for a particular
       * variable */
      void SetDataPointUpperBoundary(int index, double upperboundary);

      /**
       * Sets the error band flag to continuous function */
      void SetErrorBandContinuous(bool flag);

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
       * Set delta-function prior for a parameter. Note: this sets the
       * parameter range to the specified value. The old parameter range
       * is lost.
       * @param index The parameter index
       * @param value The position of the delta function.
       * @return An error code.
       */
      int SetPriorDelta(int index, double value);

      /**
       * Set delta-function prior for a parameter. Note: this sets the
       * parameter range to the specified value. The old parameter range
       * is lost.
       * @param name The parameter name
       * @param value The position of the delta function.
       * @return An error code.
       */
      int SetPriorDelta(const char* name, double value);

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
       * Set prior for a parameter.
       * @param index parameter index
       * @param h pointer to a histogram describing the prior
       * @param flag whether or not to use linear interpolation
       * @return An error code.
       */
      int SetPrior(int index, TH1 * h, bool flag=false);

      /**
       * Set prior for a parameter.
       * @param name parameter name
       * @param h pointer to a histogram describing the prior
       * @param flag whether or not to use linear interpolation
       * @return An error code.
       */
      int SetPrior(const char* name, TH1 * h, bool flag=false);

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
      int SetPriorConstant(const char* name)
      { return SetPriorConstant(fParameters.Index(name)); }

      /**
       * Enable caching the constant value of the prior, so LogAPrioriProbability
       * is called only once. Note that the prior for ALL parameters is
       * assumed to be constant. The value is computed from
       * the parameter ranges, so make sure these are defined before this method is
       * called.
       * @return An error code
       */
      int SetPriorConstantAll();

      /** @} */

      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Copy from object
       * @param bcmodel Object to copy from. */
      void Copy(const BCModel & bcmodel);

      /**
       * Returns the prior probability.
       * @param parameters A set of parameter values
       * @return The prior probability p(parameters)
       * @see GetPrior(std::vector<double> parameters) */
      double APrioriProbability(const std::vector<double> &parameters);

      /**
       * Returns natural logarithm of the prior probability.
       * Method needs to be overloaded by the user.
       * @param parameters A set of parameter values
       * @return The prior probability p(parameters)
       * @see GetPrior(std::vector<double> parameters) */
      virtual double LogAPrioriProbability(const std::vector<double> &parameters);

      /**
       * Returns the likelihood
       * @param params A set of parameter values
       * @return The likelihood */
      virtual double Likelihood(const std::vector<double> &params);

      /**
       * Calculates natural logarithm of the likelihood.
       * Method needs to be overloaded by the user.
       * @param params A set of parameter values
       * @return Natural logarithm of the likelihood */
      virtual double LogLikelihood(const std::vector<double> &params) = 0;

      /**
       * Returns the likelihood times prior probability given a set of parameter values
       * @param params A set of parameter values
       * @return The likelihood times prior probability */
      double ProbabilityNN(const std::vector<double> &params);

      /**
       * Returns the natural logarithm of likelihood times prior probability given
       * a set of parameter values
       * @param parameters A set of parameter values
       * @return The likelihood times prior probability */
      double LogProbabilityNN(const std::vector<double> &parameters)
      { return LogLikelihood(parameters) + LogAPrioriProbability(parameters); }

      /**
       * Returns the a posteriori probability given a set of parameter values
       * @param parameters A set of parameter values
       * @return The a posteriori probability */
      double Probability(const std::vector<double> &parameter);

      /**
       * Returns natural logarithm of the  a posteriori probability given a set of parameter values
       * @param parameters A set of parameter values
       * @return The a posteriori probability */
      double LogProbability(const std::vector<double> &parameter);

      /**
       * Sampling function used for importance sampling.
       * Method needs to be overloaded by the user.
       * @param parameters A set of parameter values
       * @return The probability density at the parameter values */
      virtual double SamplingFunction(const std::vector<double> &parameters);

      /**
       * Overloaded function to evaluate integral. */
      double Eval(const std::vector<double> &parameters);

      /**
       * Overloaded function to evaluate integral. */
      virtual double LogEval(const std::vector<double> &parameters);

      /**
       * Integrates over the un-normalized probability and updates fNormalization. */
      double Normalize();

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
       * Obtain the individual marginalized distributions
       * with respect to one parameter.
       * @note The most efficient method is to access by index.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param parameter Model parameter
       * @return 1D marginalized probability */
      BCH1D * GetMarginalized(const BCParameter * parameter);

      /**
       * Obtain the individual marginalized distributions
       * with respect to one parameter.
       * @note The most efficient method is to access by index.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param name The parameter's name
       * @return 1D marginalized probability */
      BCH1D * GetMarginalized(const char * name)
         { return GetMarginalized(fParameters.Index(name)); }

			/**
       * Obtain the individual marginalized distributions
       * with respect to one parameter.
       * @note The most efficient method is to access by index.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param index The parameter index
       * @return 1D marginalized probability */
      BCH1D * GetMarginalized(unsigned index);

      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters.
       * @note The most efficient method is to access by indices.
       * @note Ownership of the returned heap object is conferred to the caller.

       * @param parameter1 First parameter
       * @param parameter2 Second parameter
       * @return 2D marginalized probability */
      BCH2D * GetMarginalized(const BCParameter * parameter1, const BCParameter * parameter2);

      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters.
       * @note The most efficient method is to access by indices.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param name1 Name of first parameter
       * @param name2 Name of second parameter
       * @return 2D marginalized probability */
      BCH2D * GetMarginalized(const char * name1, const char * name2)
      { return GetMarginalized(fParameters.Index(name1), fParameters.Index(name2)); }

      /**
       * Obtain the individual marginalized distributions
       * with respect to two parameters.
       * @note The most efficient method is to access by indices.
       * @note Ownership of the returned heap object is conferred to the caller.
       * @param index1 Index of first parameter
       * @param index2 Index of second parameter
       * @return 2D marginalized probability */
      BCH2D * GetMarginalized(unsigned index1, unsigned index2);

      /**
       * Returns a one-dimensional slice of the pdf at the point and along a specific direction.
       * @param parameter The model parameter along which the slice is calculated.
       * @param parameters The point at which the other parameters are fixed.
       * @param nbins The number of bins of the 1D-histogram.
       * @return The 1D slice. */
      BCH1D* GetSlice(const BCParameter* parameter, const std::vector<double> parameters = std::vector<double>(0), int bins=0);

      /**
       * Returns a one-dimensional slice of the pdf at the point and along a specific direction.
       * @param name The name of the model parameter along which the slice is calculated.
       * @param parameters The point at which the other parameters are fixed.
       * @param nbins The number of bins of the 1D-histogram.
       * @return The 1D slice. */
      BCH1D* GetSlice(const char * name, const std::vector<double> parameters = std::vector<double>(0), int nbins=0)
      { return GetSlice(GetParameter(name), parameters, nbins); }

      /**
       * Returns a two-dimensional slice of the pdf at the point and along two specified directions.
       * @param parameter1 The first model parameter along which the slice is calculated.
       * @param parameter2 The second model parameter along which the slice is calculated.
       * @param parameters The point at which the other parameters are fixed.
       * @param nbins The number of bins of the 2D-histogram.
       * @return The 2D slice. */
      BCH2D* GetSlice(const BCParameter* parameter1, const BCParameter* parameter2, const std::vector<double> parameters = std::vector<double>(0), int bins=0);

      /**
       * Returns a two-dimensional slice of the pdf at the point and along two specified directions.
       * @param parameter1 The name of the first model parameter along which the slice is calculated.
       * @param parameter2 The name of the second model parameter along which the slice is calculated.
       * @param parameters The point at which the other parameters are fixed.
       * @param nbins The number of bins of the 2D-histogram.
       * @return The 2D slice. */
      BCH2D* GetSlice(const char* name1, const char* name2, const std::vector<double> parameters = std::vector<double>(0), int nbins=0);

      /**
       * Returns a two-dimensional slice of the pdf at the point and along two specified directions.
       * @param parameter1 The name of the first model parameter along which the slice is calculated.
       * @param parameter2 The name of the second model parameter along which the slice is calculated.
       * @param parameters The point at which the other parameters are fixed.
       * @param nbins The number of bins of the 2D-histogram.
       * @return The 2D slice. */
      BCH2D* GetSlice(unsigned index1, unsigned index2, const std::vector<double> parameters = std::vector<double>(0), int nbins=0);

      /**
       *   */
      int PrintAllMarginalized1D(const char * filebase);
      int PrintAllMarginalized2D(const char * filebase);
      int PrintAllMarginalized(const char * file, std::string options1d="BTciB3CS1D0pdf0Lmeanmode", std::string options2d="BTfB3CS1meangmode", unsigned int hdiv=1, unsigned int ndiv=1);

      /**
       * Constrains a data point
       * @param x A vector of double */
      virtual void CorrelateDataPointValues(std::vector<double> &x);

      /**
       * Calculate p-value from Chi2 distribution for Gaussian problems
       * @param par Parameter set for the calculation of the likelihood
       * @param sigma_index Index of the sigma/uncertainty for the data points
       *        (for data in format "x y erry" the index would be 2) */
      double GetPvalueFromChi2(const std::vector<double> &par, int sigma_index);

      /**
       * Calculate p-value from Kolmogorov-Smirnov test statistic
       * for 1D - datasets.
       *
       * @param par Parameter set for the calculation of the likelihood
       * @param index Index of the data point in the BCDataSet
       *        (for data in format "x y erry" the index would be 1) */
      double GetPvalueFromKolmogorov(const std::vector<double>& par, int index);

      double GetPvalueFromChi2NDoF(std::vector<double> par, int sigma_index);

      BCH1D * CalculatePValue(std::vector<double> par, bool flag_histogram=false);

      /**
       * @return The p-value */
      double GetPValue()
         { return fPValue; }

      double GetPValueNDoF()
         { return fPValueNDoF; }

      double GetChi2NDoF()
         { return fChi2NDoF; }

      /**
       * For a Gaussian problem, calculate the chi2 of the longest run of consecutive
       * values above/below the expected values
       * @param dataIndex component of datapoint with the observed value
       * @param sigmaIndex component of datapoint with uncertainty */
      std::vector<double> GetChi2Runs(int dataIndex, int sigmaIndex);

      /**
       * Set maximum number of iterations in the MCMC pre-run of the p-value
       * evaluation using MCMC */
      void SetGoFNIterationsMax(int n)
         { fGoFNIterationsMax=n; }

      /**
       * Set number of iterations in the MCMC normal run of the p-value
       * evaluation using MCMC */
      void SetGoFNIterationsRun(int n)
         { fGoFNIterationsRun=n; }

      /**
       * Set number of chains in the MCMC of the p-value
       * evaluation using MCMC */
      void SetGoFNChains(int n)
         { fGoFNChains=n; }

      /**
       * Calculates the matrix element of the Hessian matrix
       * @param parameter1 The parameter for the first derivative
       * @param parameter2 The parameter for the first derivative
       * @return The matrix element of the Hessian matrix */
      double HessianMatrixElement(const BCParameter * parameter1, const BCParameter * parameter2, std::vector<double> point);

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

   /** @} */

   protected:
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
       * A data set */
      BCDataSet * fDataSet;

      /**
       * data point containing the lower boundaries of possible data values */
      BCDataPoint * fDataPointLowerBoundaries;

      /**
       * data point containing the upper boundaries of possible data values */
      BCDataPoint * fDataPointUpperBoundaries;

      std::vector<bool> fDataFixedValues;

      /**
       * The p-value */
      double fPValue;

      double fChi2NDoF;
      double fPValueNDoF;

      /**
      * true for a discrete probability, false for continuous pdf  */
      bool flag_discrete;

      /**
       * Maximum number of iterations in the MCMC pre-run of the p-value
       * evaluation using MCMC */
      int fGoFNIterationsMax;

      /**
       * Number of iterations in the MCMC normal run of the p-value
       * evaluation using MCMC */
      int fGoFNIterationsRun;

      /**
       * Number of chains in the MCMC of the p-value
       * evaluation using MCMC */
      int fGoFNChains;

      /**
       * A vector of prior functions/histograms/graphs. */
      std::vector<TNamed*> fPriorContainer;

      /**
       * Flag to indicate that all parameters have constant prior. */
      bool fPriorConstantAll;

      /**
       * List the parameters whose prior is constant */
      std::vector<bool> fPriorContainerConstant;

      /**
       * List the parameters for which the histogram prior should be interpolated */
      std::vector<bool> fPriorContainerInterpolate;

   private:

      /**
       * The posterior normalization (evidence). */
      double fNormalization;

      /**
       * Compares to strings */
      int CompareStrings(const char * string1, const char * string2);

      /**
       * Converts a vector of doubles into a BCDataPoint */
      BCDataPoint * VectorToDataPoint(const std::vector<double> &data);
};

// ---------------------------------------------------------

typedef std::vector<BCModel*> BCModelContainer;

// ---------------------------------------------------------

#endif
