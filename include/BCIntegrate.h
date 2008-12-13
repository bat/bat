#ifndef __BCINTEGRATE__H
#define __BCINTEGRATE__H

/*!
 * \class BCIntegrate
 * \brief A class for handling numerical operations for models.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This is a base class for a model class. It contains
 * numerical methods to carry out the integration, marginalization,
 * peak finding etc.
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>

#include <math.h>

#include "BCEngineMCMC.h"
#include "BCParameter.h"
#include "BCDataPoint.h"

// ROOT classes
class TH1D;
class TH2D;
class TRandom3;
class TMinuit;
class TTree;


#define DEBUG 0

// ---------------------------------------------------------

class BCIntegrate : public BCEngineMCMC
{

	public:

		/** \name Enumerators */
		/* @{ */

		/**
		 * An enumerator for the integration algorithm */
		enum BCIntegrationMethod { kIntMonteCarlo, kIntImportance, kIntMetropolis, kIntCuba };

		/**
		 * An enumerator for the marginalization algorithm */
		enum BCMarginalizationMethod { kMargMonteCarlo, kMargMetropolis };

		/**
		 * An enumerator for the mode finding algorithm */
		enum BCOptimizationMethod { kOptSimulatedAnnealing, kOptMetropolis, kOptMinuit };

		/* @} */

		/** \name Constructors and destructors */
		/* @{ */

		/**
		 * The default constructor */
		BCIntegrate();

		/**
		 * A constructor */
		BCIntegrate(BCParameterSet * par);

		/**
		 * The default destructor */
		virtual ~BCIntegrate();

		/* @} */

		/** \name Member functions (get) */
		/* @{ */

		/**
		 * @return The integration method */
		BCIntegrate::BCIntegrationMethod GetIntegrationMethod()
			{ return fIntegrationMethod; };

		/**
		 * @return The marginalization method */
		BCIntegrate::BCMarginalizationMethod GetMarginalizationMethod()
			{ return fMarginalizationMethod; };

		/**
		 * @return The mode finding method */
		BCIntegrate::BCOptimizationMethod GetOptimizationMethod()
			{ return fOptimizationMethod; };

		/**
		 * Fills a vector of random numbers between 0 and 1 into a vector
		 * @param A vector of doubles */
		void GetRandomVector(std::vector <double> &x);

		virtual void GetRandomVectorMetro(std::vector <double> &x);

		/**
		 * Fills a vector of (flat) random numbers in the limits of the parameters and returns
		 * the probability at that point
		 * @param x A vector of doubles
		 * @return The (unnormalized) probability at the random point */
		double GetRandomPoint(std::vector <double> &x);

		/**
		 * Fills a vector of random numbers in the limits of the parameters sampled by the sampling
		 * function and returns the probability at that point
		 * @param x A vector of doubles
		 * @return The (unnormalized) probability at the random point */
		double GetRandomPointImportance(std::vector <double> &x);

		/**
		 * Fills a vector of random numbers in the limits of the parameters sampled by the probality
		 * function and returns the probability at that point (Metropolis)
		 * @param x A vector of doubles */
		void GetRandomPointMetro(std::vector <double> &x);

		/**
		 * Fills a vector of random numbers in the limits of the parameters sampled by the sampling
		 * function and returns the probability at that point (Metropolis)
		 * @param x A vector of doubles */
		void GetRandomPointSamplingMetro(std::vector <double> &x);

		/**
		 * @return The number of iterations per dimension for the Monte Carlo integration */
		int GetNiterationsPerDimension()
			{ return fNiterPerDimension; };

		/**
		 * @return Number of samples per 2D bin per variable in the Metropolis marginalization. */
		int GetNSamplesPer2DBin()
			{ return fNSamplesPer2DBin; };

		/**
		 * @return The number of variables to integrate over */
		int GetNvar()
			{ return fNvar; };

		/**
		 * @return The number of maximum iterations for Monte Carlo integration */
		int GetNIterationsMax()
			{ return fNIterationsMax; };

		/**
		 * @return The number of iterations for the most recent Monte Carlo integration */
		int GetNIterations()
			{ return fNIterations; };

		/**
		 * @return The relative precision for numerical integration */
		double GetRelativePrecision()
			{ return fRelativePrecision; };

		/*
		 * @return The uncertainty in the most recent Monte Carlo integration */
		double GetError()
			{ return fError; };

		/*
		 * @return number of bins per dimension for the marginalized distributions */
		int GetNbins()
			{ return fNbins; };

		/*
		 * @return Minuit used for mode finding */
		TMinuit * GetMinuit();

		int GetMinuitErrorFlag()
			{ return fMinuitErrorFlag; };

		/*
		 * @return The ROOT tree containing the Markov chain */
		TTree * GetMarkovChainTree()
			{ return fMarkovChainTree; };

		/*
		 * Generates a vector x according to the Simulated Annealing algorithm
		 * given the temperature and the stepsize relative to the range
		 * @param x Vector of doubles
		 * @param T temperature used for the stepping probability calculation
		 *  according to exp ( - (p1-p0) / T )
		 * @param step maximum stepsize relative to the range */
		void GetRandomPointSA(std::vector <double> &x, double T, double step);

		/*
		 * Returns the actual point in the markov chain */
		std::vector<double> * GetMarkovChainPoint()
			{ return &fXmetro1; };

		/**
		 * Returns the iteration of the MCMC */
		int * GetMCMCIteration()
			{ return &fMCMCIteration; };

		/**
		 * Returns the value of the loglikelihood at the point fXmetro1 */
		double * GetMarkovChainValue()
			{ return &fMarkovChainValue; };

		/* @} */

		/** \name Member functions (set) */
		/* @{ */

		void SetMinuitArlist(double * arglist)
			{ fMinuitArglist[0] = arglist[0];
			  fMinuitArglist[1] = arglist[1]; };

		/**
		 * @param par The parameter set which gets translated into array
		 * needed for the Monte Carlo integration */
		void SetParameters(BCParameterSet * par);

		/**
		 * @param varlist A list of parameters */
		void SetVarList(int * varlist);

		/**
		 * @param index The index of the variable to be set */
		void SetVar(int index){fVarlist[index]=1;};

		/**
		 * @param method The integration method */
		void SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method)
			{ fIntegrationMethod = method; };

		/**
		 * @param method The marginalization method */
		void SetMarginalizationMethod(BCIntegrate::BCMarginalizationMethod method)
			{ fMarginalizationMethod = method; };

		/**
		 * @param method The mode finding method */
		void SetOptimizationMethod(BCIntegrate::BCOptimizationMethod method)
			{ fOptimizationMethod = method; };

		/**
		 * @param niterations Number of iterations per dimension for Monte Carlo integration. */
		void SetNiterationsPerDimension(int niterations)
			{ fNiterPerDimension = niterations; };

		/**
		 * @param n Number of samples per 2D bin per variable in the Metropolis marginalization.
		 * Default is 100. */
		void SetNSamplesPer2DBin(int n)
			{ fNSamplesPer2DBin = n; };

		/**
		 * @param niterations The maximum number of iterations for Monte Carlo integration */
		void SetNIterationsMax(int niterations)
			{ fNIterationsMax = niterations; };

		/**
		 * @param relprecision The relative precision envisioned for Monte
		 * Carlo integration */
		void SetRelativePrecision(double relprecision)
			{ fRelativePrecision = relprecision; };

		/**
		 * @param n Number of bins per dimension for the marginalized
		 * distributions.  Default is 100. Minimum number allowed is 2. */

		/**
		 * Set the number of bins for the marginalized distribution of a
		 * parameter.
		 * @param nbins Number of bins (default = 100)
		 * @param index Index of the parameter. */
		void SetNbins(int nbins, int index = -1);

//		void SetNbins(int n);
//		void SetNbinsX(int n);
//		void SetNbinsY(int n);

		/*
		 * Turn on or off the filling of the error band during the MCMC run.
		 * @param flag set to true for turning on the filling */
		void SetFillErrorBand(bool flag = true)
			{ fFillErrorBand=flag; };

		/*
		 * Turn off filling of the error band during the MCMC run.
		 * This method is equivalent to SetFillErrorBand(false) */
		void UnsetFillErrorBand()
			{ fFillErrorBand=false; };

		/**
		 * Sets index of the x values in function fits.
		 * @param index Index of the x values */
		void SetFitFunctionIndexX(int index)
			{ fFitFunctionIndexX = index; };

		/**
		 * Sets index of the y values in function fits. @param index Index
		 * of the y values */
		void SetFitFunctionIndexY(int index)
			{ fFitFunctionIndexY = index; };

		void SetFitFunctionIndices(int indexx, int indexy)
			{ this -> SetFitFunctionIndexX(indexx);
				this -> SetFitFunctionIndexY(indexy); };

		/**
		 * Sets the data point containing the lower boundaries of possible
		 * data values */
		void SetDataPointLowerBoundaries(BCDataPoint * datasetlowerboundaries)
			{ fDataPointLowerBoundaries = datasetlowerboundaries; };

		/**
		 * Sets the data point containing the upper boundaries of possible
		 * data values */
		void SetDataPointUpperBoundaries(BCDataPoint * datasetupperboundaries)
			{ fDataPointUpperBoundaries = datasetupperboundaries; };

		/**
		 * Sets the lower boundary of possible data values for a particular
		 * variable */
		void SetDataPointLowerBoundary(int index, double lowerboundary)
			{ fDataPointLowerBoundaries -> SetValue(index, lowerboundary); };

		/**
		 * Sets the upper boundary of possible data values for a particular
		 * variable */
		void SetDataPointUpperBoundary(int index, double upperboundary)
			{ fDataPointUpperBoundaries -> SetValue(index, upperboundary); };

		/**
		 * Flag for writing Markov chain to ROOT file (true) or not (false)
		 */
		void WriteMarkovChain(bool flag)
			{ fFlagWriteMarkovChain = flag; fMCMCFlagWriteChainToFile = flag;};

		/**
		 * Sets the ROOT tree containing the Markov chain */
		void SetMarkovChainTree(TTree * tree)
			{ fMarkovChainTree = tree; };

		/**
		 * Sets the initial position for the Markov chain */
		void SetMarkovChainInitialPosition(std::vector<double> position);

		/**
		 * Sets the step size for Markov chains */
		void SetMarkovChainStepSize(double stepsize)
			{ fMarkovChainStepSize = stepsize; };

		/**
		 * Sets the number of iterations in the markov chain */
		void SetMarkovChainNIterations(int niterations)
			{ fMarkovChainNIterations = niterations;
			  fMarkovChainAutoN = false; }

		/**
		 * Sets a flag for automatically calculating the number of iterations */
		void SetMarkovChainAutoN(bool flag)
			{ fMarkovChainAutoN = flag; };

		/**
		 * Sets mode */
		void SetMode(std::vector <double> mode);

		/**
		 * Sets errorband histogram */
		void SetErrorBandHisto(TH2D * h)
			{ fErrorBandXY = h; };

		/* @} */

		/** \name Member functions (miscellaneous methods) */
		/* @{ */

		/**
		 * Frees the memory for integration variables */
		void DeleteVarList();

		/**
		 * Sets all values of the variable list to a particular value
		 * @v The value */
		void ResetVarlist(int v);

		/**
		 * Set value of a particular integration variable to 0.
		 * @param index The index of the variable */
		void UnsetVar(int index)
			{ fVarlist[index] = 0; };

		/**
		 * Evaluate the unnormalized probability at a point in parameter space.
		 * Method needs to be overloaded by the user.
		 * @param x The point in parameter space
		 * @return The unnormalized probability */
		virtual double Eval(std::vector <double> x);

		/**
		 * Evaluate the natural logarithm of the Eval function. For better numerical
		 * stability, this method should also be overloaded by the user.
		 * @param x The point in parameter space
		 * @return log(Eval(x)) */
		virtual double LogEval(std::vector <double> x);

		/**
		 * Evaluate the sampling function at a point in parameter space.
		 * Method needs to be overloaded by the user.
		 * @param x The point in parameter space
		 * @return The value of the sampling function */
		virtual double EvalSampling(std::vector <double> x);

		/**
		 * Evaluate the natural logarithm of the EvalSampling function.
		 * Method needs to be overloaded by the user.
		 * @param x The point in parameter space
		 * @return log(EvalSampling(x)) */
		double LogEvalSampling(std::vector <double> x);

		/**
		 * Evaluate the un-normalized probability at a point in parameter space
		 * and prints the result to the log.
		 * @param x The point in parameter space
		 * @return The un-normalized probability
		 * @see Eval(std::vector <double> x) */
		double EvalPrint(std::vector <double> x);

		/**
		 * Defines a fit function.
		 * @param parameters A set of parameter values
		 * @param x A vector of x-values
		 * @return The value of the fit function at the x-values given a set of parameters */
		virtual double FitFunction(std::vector <double> x, std::vector <double> parameters)
			{ return 0.; };

		/**
		 * Does the integration over the un-normalized probability.
		 * @return The normalization value */
		double Integrate();

		/**
		 * Perfoms the Monte Carlo integration. For details see documentation.
		 * @param x An initial point in parameter space
		 * @param varlist A list of variables
		 * @return The integral */
		double IntegralMC(std::vector <double> x, int * varlist);

		double IntegralMC(std::vector <double> x);

		/**
		 * Perfoms the Metropolis Monte Carlo integration. For details see documentation.
		 * @param x An initial point in parameter space
		 * @return The integral */
		double IntegralMetro(std::vector <double> x);

		/**
		 * Perfoms the importance sampling Monte Carlo integration. For details see documentation.
		 * @param x An initial point in parameter space
		 * @return The integral */
		double IntegralImportance(std::vector <double> x);

		/**
		 * Calculate integral using the Cuba library. For details see documentation.
		 * @param method A short cut for the method
		 * @param parameters_double A vector of parameters (double)
		 * @param parameters_int A vector of parameters (int)
		 * @return The integral */
		double CubaIntegrate(int method, std::vector<double> parameters_double, std::vector<double> parameters_int);

		double CubaIntegrate();

		/**
		 * Integrand for the Cuba library.
		 * @param ndim The number of dimensions to integrate over
		 * @param xx The point in parameter space to integrate over (scaled to 0 - 1 per dimension)
		 * @param ncomp The number of components of the integrand (usually 1)
		 * @param ff The function value
		 * @return The integral */
		static void CubaIntegrand(const int * ndim, const double xx[], const int * ncomp, double ff[]);

		/**
		 * Performs the marginalization with respect to one parameter.
		 * @param parameter The parameter w.r.t. which the marginalization is performed
		 * @return A histogram which contains the marginalized probability distribution (normalized to 1) */
		TH1D * Marginalize(BCParameter * parameter);

		/**
		 * Performs the marginalization with respect to two parameters.
		 * @param parameter1 The first parameter w.r.t. which the marginalization is performed
		 * @param parameter2 The second parameter w.r.t. which the marginalization is performed
		 * @return A histogram which contains the marginalized probability distribution (normalized to 1) */
		TH2D * Marginalize(BCParameter * parameter1, BCParameter * parameter2);

		/**
		 * Performs the marginalization with respect to one parameter using
		 * the simple Monte Carlo technique.
		 * @param parameter The parameter w.r.t. which the marginalization is performed
		 * @return A histogram which contains the marginalized probability distribution (normalized to 1) */
		TH1D * MarginalizeByIntegrate(BCParameter * parameter);

		/**
		 * Performs the marginalization with respect to two parameters using
		 * the simple Monte Carlo technique.
		 * @param parameter1 The first parameter w.r.t. which the marginalization is performed
		 * @param parameter2 The second parameter w.r.t. which the marginalization is performed
		 * @return A histogram which contains the marginalized probability distribution (normalized to 1) */
		TH2D * MarginalizeByIntegrate(BCParameter * parameter1, BCParameter * parameter2);

		/**
		 * Performs the marginalization with respect to one parameter using
		 * the Metropolis algorithm.
		 * @param parameter The parameter w.r.t. which the marginalization is performed
		 * @return A histogram which contains the marginalized probability distribution (normalized to 1) */
		TH1D * MarginalizeByMetro(BCParameter * parameter);

		/**
		 * Performs the marginalization with respect to two parameters using the Metropolis algorithm.
		 * @param parameter1 The first parameter w.r.t. which the marginalization is performed
		 * @param parameter2 The second parameter w.r.t. which the marginalization is performed
		 * @return A histogram which contains the marginalized probability distribution (normalized to 1) */
		TH2D * MarginalizeByMetro(BCParameter * parameter1, BCParameter * parameter2);

		/**
		 * Performs the marginalization with respect to every single parameter as well as with respect
		 * all combinations to two parameters using the Metropolis algorithm.
		 * @param name Basename for the histograms (e.g. model name)
		 * @return Total number of marginalized distributions */
		int MarginalizeAllByMetro(const char * name);

		/**
		 * @param parIndex1 Index of parameter
		 * @return Pointer to 1D histogram (TH1D) of marginalized distribution wrt. parameter with given index.
		 */
		TH1D * GetH1D(int parIndex);

		/**
		 * @param parIndex1 Index of first parameter
		 * @param parIndex2 Index of second parameter, with parIndex2>parIndex1
		 * @return Index of the distribution in the vector of 2D distributions, which corresponds
		 * to the combination of parameters with given indeces */
		int GetH2DIndex(int parIndex1, int parIndex2);

		/**
		 * @param parIndex1 Index of first parameter
		 * @param parIndex2 Index of second parameter, with parIndex2>parIndex1
		 * @return Pointer to 2D histogram (TH2D) of marginalized distribution wrt. parameters with given indeces.
		 */
		TH2D * GetH2D(int parIndex1, int parIndex2);

		/**
		 * Initializes the Metropolis algorithm (for details see manual) */
		void InitMetro();

		/**
		 * Does the mode finding */
//		void FindMode();

		/**
		 * Does the mode finding using Simulated Annealing (SA) algorithm */
		void FindModeSA();

		/**
		 * Does the mode finding using Minuit. If starting point is not specified,
		 * finding will start from the center of the parameter space.
		 * @param start point in parameter space from which the mode finding is started.
		 * @param printlevel The print level. */
		void FindModeMinuit(std::vector<double> start = std::vector<double>(0), int printlevel = 1);

		/**
		 * Does the mode finding using Markov Chain Monte Carlo */
		void FindModeMCMC();
//		void FindModeMCMC(int flag_run = 0);

		static void FCNLikelihood(int &npar, double * grad, double &fval, double * par, int flag);

		/*
		 * Method executed for every iteration of the MCMC. User's code should be
		 * provided via overloading in the derived class*/
		virtual void MCMCUserIterationInterface()
			{};

		/* @} */

	private:

		/**
		 * Set of parameters for the integration. */
		BCParameterSet * fx;

		/**
		 * Array containing the lower boundaries of the variables to integrate over. */
		double * fMin;

		/**
		 * Array containing the upper boundaries of the variables to integrate over. */
		double * fMax;

		/**
		 * List of variables containing a flag whether to integrate over them or not. */
		int * fVarlist;

		/**
		 * Number of iteration per dimension for Monte Carlo integration. */
		int fNiterPerDimension;

		/**
		 * Current integration method */
		BCIntegrate::BCIntegrationMethod fIntegrationMethod;

		/**
		 * Current marginalization method */
		BCIntegrate::BCMarginalizationMethod fMarginalizationMethod;

		/**
		 * Current mode finding method */
		BCIntegrate::BCOptimizationMethod fOptimizationMethod;

		/**
		 * Maximum number of iterations */
		int fNIterationsMax;

		/**
		 * Number of iterations in the most recent Monte Carlo integation */
		int fNIterations;

		/**
		 * Relative precision aimed at in the Monte Carlo integation */
		double fRelativePrecision;

		/**
		 * The uncertainty in the most recent Monte Carlo integration */
		double fError;

		/**
		 * The number of iterations in the Metropolis integration */
		int fNmetro;
		int fNacceptedMCMC;

		/**
		 * A vector of points in parameter space used for the Metropolis algorithm */
		std::vector <double> fXmetro0;

		/**
		 * A vector of points in parameter space used for the Metropolis algorithm */
		std::vector <double> fXmetro1;

		/**
		 * A double containing the log likelihood value at the point fXmetro1 */
		double fMarkovChainValue;

		/*
		 * Method executed for every iteration of the MCMC, overloaded from BCEngineMCMC. */
		void MCMCIterationInterface();

		/*
		 * Fill error band histogram for curreent iteration. This method is called from MCMCIterationInterface() */
		void MCMCFillErrorBand();

	protected:

		/**
		 * Number of variables to integrate over. */
		int fNvar;

		/**
		 * Number of bins per dimension for the marginalized distributions */
		int fNbins;

		/**
		 * Number of samples per 2D bin per variable in the Metropolis
		 * marginalization. */
		int fNSamplesPer2DBin;

		/**
		 * Step size in the Markov chain relative to min and max */
		double fMarkovChainStepSize;

		int fMarkovChainNIterations;

		bool fMarkovChainAutoN;

		/**
		 * data point containing the lower boundaries of possible data values */
		BCDataPoint * fDataPointLowerBoundaries;

		/**
		 * data point containing the upper boundaries of possible data values */
		BCDataPoint * fDataPointUpperBoundaries;

		std::vector <bool> fDataFixedValues;

		/**
		 * A ROOT random number generator */
		TRandom3 * fRandom;

		/**
		 * A vector of best fit parameters estimated from the global probability */
		std::vector <double> fBestFitParameters;

		/**
		 * A vector of best fit parameters estimated from the marginalized probability */
		std::vector <double> fBestFitParametersMarginalized;

		/**
		 * Vector of TH1D histograms for marginalized probability distributions */
		std::vector <TH1D *> fHProb1D;

		/**
		 * Vector of TH2D histograms for marginalized probability distributions */
		std::vector <TH2D *> fHProb2D;

		/**
		 * Flag whether or not to fill the error band */
		bool fFillErrorBand;

		/**
		 * The indeces for function fits */
		int fFitFunctionIndexX;
		int fFitFunctionIndexY;

		/**
		 * A flag for single point evaluation of the error "band" */
		bool fErrorBandContinuous;
		std::vector <double> fErrorBandX;

		/**
		 * The error band histogram */
		TH2D * fErrorBandXY;

		/**
		 * Nnumber of X bins of the error band histogram */
		int fErrorBandNbinsX;

		/**
		 * Nnumber of Y bins of the error band histogram */
		int fErrorBandNbinsY;

		/**
		 * Minuit */
		TMinuit * fMinuit;

		double fMinuitArglist[2];
		int fMinuitErrorFlag;

		/**
		 * Flag for writing Markov chain to file */
		bool fFlagWriteMarkovChain;

		/**
		 * ROOT tree containing the Markov chain */
		TTree * fMarkovChainTree;

		/**
		 * Iteration of the MCMC */
		int fMCMCIteration;

};

// ---------------------------------------------------------

#endif
