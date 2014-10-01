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

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
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
         { fDataSet = dataset; }

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
