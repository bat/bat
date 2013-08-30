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

/**
 * Copyright (C) 2008-2013, Daniel Kollar, Kevin Kroeninger,
 * Daniel Greenwald, and Frederik Beaujean.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCEngineMCMC.h"

// ROOT classes
class TH1;
class TH1D;
class TH2D;
class TMinuit;
class TTree;

/**
 * Collect all the useful options to tune integration with the CUBA library.
 * Option names match those used in the CUBA manual, where a detailed
 * description is given.
 */
namespace BCCubaOptions
{
struct General
{
   int ncomp, flags, nregions, neval, fail;
   double error, prob;
   General();
protected:
   ~General();
};

struct Vegas : public General
{
   int nstart, nincrease, nbatch, gridno;

   Vegas();
};

struct Suave : public General
{
   int nnew;
   double flatness;

   Suave();
};

struct Divonne : public General
{
   int key1, key2, key3, maxpass;
   double border, maxchisq, mindeviation;

   Divonne();
};

struct Cuhre : public General
{
   int key;

   Cuhre();
};
}
// ---------------------------------------------------------

class BCIntegrate : public BCEngineMCMC
{

public:

   /** \name Enumerators */
   /** @{ */

   /**
    * An enumerator for integration algorithms */
   enum BCIntegrationMethod {
      kIntMonteCarlo,
      kIntImportance,
      kIntCuba,
      NIntMethods };

   /**
    * An enumerator for marginalization algorithms */
   enum BCMarginalizationMethod {
      kMargMCMC,
      kMargMonteCarlo,
      kMargSlice,
      kMargDefault,
      NMargMethods };

   /**
    * An enumerator for the Simulated Annealing schedule */
   enum BCSASchedule { kSACauchy, kSABoltzmann, kSACustom, NSAMethods };

   /**
    * An enumerator for Cuba integration methods */
   enum BCCubaMethod { kCubaVegas, kCubaSuave, kCubaDivonne, kCubaCuhre, NCubaMethods};

   /** @} */

   /** \name Function pointer types */
   /** @{ */

   /**
    * A pointer for a function that chooses a next random point */
   typedef void (BCIntegrate::*tRandomizer)(std::vector<double> &) const;

   /**
    * A pointer for a function that evaluates at a point */
   typedef double (BCIntegrate::*tEvaluator)(std::vector<double> &, const std::vector<double> &, bool &);

   /**
    * A pointer for a function that updates the integral and absolute precision */
   typedef void (*tIntegralUpdater)(const std::vector<double> &, const int &, double &, double &);

   /** @} */


   /** \name Constructors and destructors */
   /** @{ */

   /**
    * A constructor */
   BCIntegrate();

   /**
    * The copy constructor */
   BCIntegrate(const BCIntegrate & bcintegrate);

   /**
    * The default destructor */
   virtual ~BCIntegrate();

   /** @} */
   /** \name Assignment operators */
   /** @{ */

   /**
    * Defaut assignment operator */
   BCIntegrate & operator = (const BCIntegrate & bcintegrate);

   /** @} */
   /** \name Member functions (get) */
   /** @{ */

      /**
       * Read */
      int ReadMarginalizedFromFile(const char * file);

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
       *   */
      int PrintAllMarginalized1D(const char * filebase);
      int PrintAllMarginalized2D(const char * filebase);
      int PrintAllMarginalized(const char * file, std::string options1d="BTciB3CS1D0pdf0Lmeanmode", std::string options2d="BTfB3CS1meangmode", unsigned int hdiv=1, unsigned int ndiv=1);


      /**
       * @return The integral. */
      double GetIntegral() const
         { return fIntegral; }

   /**
    * @return The integration method */
   BCIntegrate::BCIntegrationMethod GetIntegrationMethod() const
   { return fIntegrationMethod; }

   /**
    * @return The marginalization method */
   BCIntegrate::BCMarginalizationMethod GetMarginalizationMethod() const
   { return fMarginalizationMethod; }

   /**
    * @return The Simulated Annealing schedule */
   BCIntegrate::BCSASchedule GetSASchedule() const
   { return fSASchedule; }

   /**
    * Fills a vector of random numbers between 0 and 1 into a vector
    * @param A vector of doubles */
   void GetRandomVectorUnitHypercube(std::vector<double> &x) const;

   /**
    * Fills a vector of random numbers x[i] between fMin[i] and fMax[i] into a vector
    * @param A vector of doubles */
   void GetRandomVectorInParameterSpace(std::vector<double> &x) const;

   /**
    * Fills a vector of (flat) random numbers in the limits of the parameters and returns
    * the probability at that point
    * @param x A vector of doubles
    * @return The (unnormalized) probability at the random point */
   double GetRandomPoint(std::vector<double> &x);

   /**
    * @return The number of minimum iterations for integration */
   int GetNIterationsMin() const
   { return fNIterationsMin; }

   /**
    * @return The number of maximum iterations for integration */
   int GetNIterationsMax() const
   { return fNIterationsMax; }

   /**
    * @return The interval for checking precision in integration */
   int GetNIterationsPrecisionCheck() const
   { return fNIterationsPrecisionCheck; }

   /**
    * @return The interval for outputting during integration */
   int GetNIterationsOutput() const
   { return fNIterationsOutput; }

   /**
    * @return The number of iterations for the most recent Monte Carlo integration */
   int GetNIterations() const
   { return fNIterations; }

   /**
    * @return Requested relative precision of the numerical integation */
   double GetRelativePrecision() const
   { return fRelativePrecision; }

   /**
    * @return Requested absolute precision of the numerical integration */
   double GetAbsolutePrecision() const
   { return fAbsolutePrecision; }

   /**
    * @return Cuba Integration method */
   BCCubaMethod GetCubaIntegrationMethod() const
   { return fCubaIntegrationMethod; }

   /**
    * @return Options used for integration with CUBA
    */
   const BCCubaOptions::Vegas & GetCubaVegasOptions() const
   { return fCubaVegasOptions; }

   const BCCubaOptions::Suave & GetCubaSuaveOptions() const
   { return fCubaSuaveOptions; }

   const BCCubaOptions::Divonne & GetCubaDivonneOptions() const
   { return fCubaDivonneOptions; }

   const BCCubaOptions::Cuhre & GetCubaCuhreOptions() const
   { return fCubaCuhreOptions; }

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
    * @return The uncertainty in the most recent Monte Carlo integration */
   double GetError() const
   { return fError; }

   /**
    * @return Minuit used for mode finding */
   TMinuit * GetMinuit();

   /**
    * @return Error flag from Minuit run */
   int GetMinuitErrorFlag() const
   { return fMinuitErrorFlag; }

   /**
    * Returns the Simulated Annealing starting temperature. */
   double GetSAT0() const
   { return fSAT0; }

   /**
    * Returns the Simulated Annealing threshhold temperature. */
   double GetSATmin() const
   { return fSATmin; }

   /**
    * Returns the value of a parameter (defined by index) at
    * the global mode of the posterior pdf.
    * @param index index of the parameter.
    * @return best fit value of the parameter or -1e+111 on error or center of the range if mode finding not yer run */
   double GetBestFitParameter(unsigned index) const;

   /**
    * Returns the error on the value of a parameter (defined by index) at
    * the global mode of the posterior pdf.
    * @param index index of the parameter.
    * @return error on the best fit value of the parameter or -1 if undefined */
   double GetBestFitParameterError(unsigned index) const;

   /**
    * Returns the set of values of the parameters at the global mode of
    * the posterior pdf.
    * @return The best fit parameters */
   const std::vector<double> & GetBestFitParameters() const
   { return fBestFitParameters; }

   /**
    * Returns the set of errors on the values of the parameters at the global mode */
   const std::vector<double> & GetBestFitParameterErrors() const
	{ return fBestFitParameterErrors; }

   /** @} */

   /** \name Member functions (set) */
   /** @{ */

   /**
    * @arglist pointer to list of doubles to be passed as arguments to Minuit */
   void SetMinuitArlist(double * arglist)
   { fMinuitArglist[0] = arglist[0];
   fMinuitArglist[1] = arglist[1]; }

   /**
    * @flag Flag whether or not to ignore result of previous mode finding */
   void SetFlagIgnorePrevOptimization(bool flag)
   { fFlagIgnorePrevOptimization = flag; }

   /**
    * @param method The integration method */
   void SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method);

   /**
    * @param method The marginalization method */
   void SetMarginalizationMethod(BCIntegrate::BCMarginalizationMethod method)
   { fMarginalizationMethod = method; }

   /**
    * @param method The Simulated Annealing schedule */
   void SetSASchedule(BCIntegrate::BCSASchedule schedule)
   { fSASchedule = schedule; }

   /**
    * @param niterations The maximum number of iterations for integration */
   void SetNIterationsMin(int niterations)
   { fNIterationsMin = niterations; }

   /**
    * @param niterations The maximum number of iterations for integration */
   void SetNIterationsMax(int niterations)
   { fNIterationsMax = niterations; }

   /**
    * @param niterations interval for checking precision in integration routines */
   void SetNIterationsPrecisionCheck(int niterations)
   { fNIterationsPrecisionCheck = niterations; }

   /**
    * @param niterations interval for outputting during integration. If negative, frequency is autogenerated. */
   void SetNIterationsOutput(int niterations)
   { fNIterationsOutput = niterations; }

   /**
    * @param relprecision The relative precision envisioned for Monte
    * Carlo integration */
   void SetRelativePrecision(double relprecision)
   { fRelativePrecision = relprecision; }

   /**
    * Set absolute precision of the numerical integation */
   void SetAbsolutePrecision(double absprecision)
   { fAbsolutePrecision = absprecision; }

   /**
    * Set Cuba integration method */
   void SetCubaIntegrationMethod(BCCubaMethod type);

   /**
    * Set options for individual cuba methods
    */
   void SetCubaOptions(const BCCubaOptions::Vegas & options)
   {  fCubaVegasOptions = options; }

   void SetCubaOptions(const BCCubaOptions::Suave & options)
   {  fCubaSuaveOptions = options; }

   void SetCubaOptions(const BCCubaOptions::Divonne & options)
   {  fCubaDivonneOptions = options; }

   void SetCubaOptions(const BCCubaOptions::Cuhre & options)
   {  fCubaCuhreOptions = options; }

   /**
    * @param T0 new value for Simulated Annealing starting temperature. */
   void SetSAT0(double T0)
   { fSAT0 = T0; }

   /**
    * @param Tmin new value for Simulated Annealing threshold temperature. */
   void SetSATmin(double Tmin)
   { fSATmin = Tmin; }

   void SetFlagWriteSAToFile(bool flag)
   { fFlagWriteSAToFile = flag; }

   /**
    * Sets the tree containing the Simulated Annealing  chain. */
   void SetSATree(TTree * tree)
   { fTreeSA = tree; }

   /**
    * Getter for the tree containing the  Simulated Annealing  chain. */
   TTree * GetSATree()
   { return fTreeSA; }

   /**
    * Initialization of the tree for the Simulated Annealing */
   void InitializeSATree();

   /** @} */

   /** \name Member functions (miscellaneous methods) */
   /** @{ */

   /**
    * Evaluate the unnormalized probability at a point in parameter space.
    * Method needs to be overloaded by the user.
    * @param x The point in parameter space
    * @return The unnormalized probability */
   virtual double Eval(const std::vector<double> &x);

   /**
    * Evaluate the natural logarithm of the Eval function. For better numerical
    * stability, this method should also be overloaded by the user.
    * @param x The point in parameter space
    * @return log(Eval(x)) */
   virtual double LogEval(const std::vector<double> &x);

   /**
    * Evaluate the sampling function at a point in parameter space.
    * Method needs to be overloaded by the user.
    * @param x The point in parameter space
    * @return The value of the sampling function */
   virtual double EvalSampling(const std::vector<double> &x);

   /**
    * Performs integration. */
   double Normalize()
   { return Integrate(); };

   /**
    * Does the integration over the un-normalized probability.
    * @return The normalization value */
   double Integrate()
   { return Integrate(fIntegrationMethod);}

   /**
    * Does the integration over the un-normalized probability.
    * @param type The integration method to used
    * @return The normalization value */
   double Integrate(BCIntegrationMethod type);

   /**
    * Does the integration over the un-normalized probability.
    * @param type The integration method to used (for printing out status updates by name)
    * @param randomizer Pointer to function to choose next random point
    * @param evaluator Pointer to function to evaluate point
    * @param updater Pointer to function to update integral and precision
    * @param sums Vector of doubles holding values used in integral calculation
    * @param x An initial point to start integration routine at
    * @return The integral value */
   double Integrate(BCIntegrationMethod type, tRandomizer randomizer, tEvaluator evaluator, tIntegralUpdater updater, std::vector<double> &sums);

   // todo document
   double EvaluatorMC(std::vector<double> &sums, const std::vector<double> &point, bool &accepted);
   static void IntegralUpdaterMC(const std::vector<double> &sums, const int &nIterations, double &integral, double &absprecision);
   double EvaluatorImportance(std::vector<double> &sums, const std::vector<double> &point, bool &accepted);
   static void IntegralUpdaterImportance(const std::vector<double> &sums, const int &nIterations, double &integral, double &absprecision);

   /**
    * Calculate integral using the Cuba library. For details see documentation.
    * @return The integral */
   double IntegrateCuba()
   { return IntegrateCuba(fCubaIntegrationMethod); }

   /**
    * Calculate integral using the Cuba library. For details see documentation.
    * @param Cuba integration method to use
    * @return The integral */
   double IntegrateCuba(BCCubaMethod cubatype);

   /**
    * Integrand for the Cuba library.
    * @param ndim The number of dimensions to integrate over
    * @param xx The point in parameter space to integrate over (scaled to 0 - 1 per dimension)
    * @param ncomp The number of components of the integrand (usually 1)
    * @param ff The function value
    * @return An error code */
   static int CubaIntegrand(const int * ndim, const double xx[], const int * ncomp, double ff[], void *userdata);

   TH1D * Marginalize(BCIntegrationMethod type, unsigned index);

   TH2D * Marginalize(BCIntegrationMethod type, unsigned index1, unsigned index2);

   bool Marginalize(TH1* hist, BCIntegrationMethod type, const std::vector<unsigned> &index);

   /**
    * Marginalize all probabilities wrt. single parameters and all combinations
    * of two parameters. The individual distributions can be retrieved using
    * the GetMarginalized method.
    * @return Total number of marginalized distributions */
   int MarginalizeAll();

   /**
    * Method executed for before marginalization. User's code should
    * be provided via overloading in the derived class */
   virtual void MarginalizePreprocess()
   {};

   /**
    * Method executed after marginalization. User's code should be
    * provided via overloading in the derived class*/
   virtual void MarginalizePostprocess()
   {};

   /**
    * Initializes the Simulated Annealing algorithm (for details see manual) */
   void SAInitialize();

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
   * @return The mode found.
    * @note The result may not coincide with the result of @code GetBestFitParameters()
    * if a previous optimization found a better value. */
    std::vector<double> FindMode(std::vector<double> start = std::vector<double>(0));

   /**
    * Does the mode finding using Minuit. If starting point is not specified,
    * finding will start from the center of the parameter space.
    * @param start point in parameter space from which the mode finding is started.
    * @param printlevel The print level.
    * @return The mode found.
    * @note The result may not coincide with the result of @code GetBestFitParameters()
    * if a previous optimization found a better value. */
   std::vector<double> FindModeMinuit(std::vector<double> start = std::vector<double>(0), int printlevel = 1);

   /**
    * Does the mode finding using Markov Chain Monte Carlo (prerun only!)
    * @return The mode.
    * @note The result may not coincide with the result of @code GetBestFitParameters()
    * if a previous optimization found a better value. */
   std::vector<double> FindModeMCMC();

   /**
    * Does the mode finding using Simulated Annealing. If starting point
    * is not specified, finding will start from the center of the
    * parameter space.
    * @param start point in parameter space from thich the mode finding is started.
    * @return The mode.
    * @note The result may not coincide with the result of @code GetBestFitParameters()
    * if a previous optimization found a better value.*/
   std::vector<double> FindModeSA(std::vector<double> start = std::vector<double>(0));

   /**
    * Temperature annealing schedule for use with Simulated Annealing.
    * Delegates to the appropriate method according to
    * fSASchedule.
    * @param t iterator for lowering the temperature over time. */
   double SATemperature(double t);

   /**
    * Temperature annealing schedule for use with Simulated Annealing.
    * This method is used for Boltzmann annealing schedule.
    * @param t iterator for lowering the temperature over time. */
   double SATemperatureBoltzmann(double t);

   /**
    * Temperature annealing schedule for use with Simulated Annealing.
    * This method is used for Cauchy annealing schedule.
    * @param t iterator for lowering the temperature over time. */
   double SATemperatureCauchy(double t);

   /**
    * Temperature annealing schedule for use with Simulated Annealing.
    * This is a virtual method to be overridden by a user-defined
    * custom temperature schedule.
    * @param t iterator for lowering the temperature over time. */
   virtual double SATemperatureCustom(double t);

   /**
    * Generates a new state in a neighbourhood around x that is to be
    * accepted or rejected by the Simulated Annealing algorithm.
    * Delegates the generation to the appropriate method according
    * to fSASchedule.
    * @param x last state.
    * @param t time iterator to determine current temperature. */
   std::vector<double> GetProposalPointSA(const std::vector<double> &x, int t);

   /**
    * Generates a new state in a neighbourhood around x that is to be
    * accepted or rejected by the Simulated Annealing algorithm.
    * This method is used for Boltzmann annealing schedule.
    * @param x last state.
    * @param t time iterator to determine current temperature. */
   std::vector<double> GetProposalPointSABoltzmann(const std::vector<double> &x, int t);

   /**
    * Generates a new state in a neighbourhood around x that is to be
    * accepted or rejected by the Simulated Annealing algorithm.
    * This method is used for Cauchy annealing schedule.
    * @param x last state.
    * @param t time iterator to determine current temperature. */
   std::vector<double> GetProposalPointSACauchy(const std::vector<double> &x, int t);

   /**
    * Generates a new state in a neighbourhood around x that is to be
    * accepted or rejected by the Simulated Annealing algorithm.
    * This is a virtual method to be overridden by a user-defined
    * custom proposal function.
    * @param x last state.
    * @param t time iterator to determine current temperature. */
   virtual std::vector<double> GetProposalPointSACustom(const std::vector<double> &x, int t);

   /**
    * Generates a uniform distributed random point on the surface of
    * a fNvar-dimensional Hypersphere.
    * Used as a helper to generate proposal points for Cauchy annealing. */
   std::vector<double> SAHelperGetRandomPointOnHypersphere();

   /**
    * Generates the radial part of a n-dimensional Cauchy distribution.
    * Helper function for Cauchy annealing. */
   double SAHelperGetRadialCauchy();

   /**
    * Returns the Integral of sin^dim from 0 to theta.
    * Helper function needed for generating Cauchy distributions. */
   double SAHelperSinusToNIntegral(int dim, double theta);


   static void FCNLikelihood(int &npar, double * grad, double &fval, double * par, int flag);

   /**
    * Reset all information on the best fit parameters. */
   virtual void ResetResults();

   /**
    * Return string with the name for a given integration type.
    * @param type code for the integration type
    * @return string containing the name of the integration type */
   std::string DumpIntegrationMethod(BCIntegrationMethod type);

   /**
    * Return string with the name for the currently set integration type.
    * @return string containing the name of the integration type */
   std::string DumpIntegrationMethod()
   { return DumpIntegrationMethod(fIntegrationMethod); }

   /**
    * Return string with the name for a given marginalization type.
    * @param type code for the marginalization type
    * @return string containing the name of the marginalization type */
   std::string DumpMarginalizationMethod(BCMarginalizationMethod type);

   /**
    * Return string with the name for the currently set marginalization type.
    * @return string containing the name of the marginalization type */
   std::string DumpMarginalizationMethod()
   { return DumpMarginalizationMethod(fMarginalizationMethod); }

   /**
    * Return string with the name for a given optimization type.
    * @param type code for the optimization type
    * @return string containing the name of the optimization type */
   std::string DumpOptimizationMethod(BCOptimizationMethod type);

   /**
    * Return string with the name for the currently set optimization type.
    * @return string containing the name of the optimization type */
   std::string DumpOptimizationMethod()
   { return DumpOptimizationMethod(fOptimizationMethod); }

   /**
    * Return string with the name for the optimization type used to find the current mode.
    * @return string containing the name of the optimization type */
   std::string DumpUsedOptimizationMethod()
   { return DumpOptimizationMethod(fOptimizationMethodMode); }

   /**
    * Return string with the name for a given Cuba integration type.
    * @param type code for the Cuba integration type
    * @return string containing the name of the Cuba integration type */
   std::string DumpCubaIntegrationMethod(BCCubaMethod type);

   /**
    * Return string with the name for the currently set Cuba integration type.
    * @return string containing the name of the Cuba integration type */
   std::string DumpCubaIntegrationMethod()
   { return DumpCubaIntegrationMethod(fCubaIntegrationMethod); }

   /**
    * Set best fit parameters values*/
   void SetBestFitParameters(const std::vector<double> &x)
   { fBestFitParameters = x; }

   /**
    * Set best fit parameters if best fit
    * @param new_value is the value of the function at x
    * @param old_value is the old best fit value, updated to new_value, if it is larger */
   void SetBestFitParameters(const std::vector<double> &x, const double &new_value, double &old_value);

   /**
    * Get number of variables that are varied in the integration
    * @return fNvar minus the number of fixed variables */
   unsigned GetNIntegrationVariables();

   /**
    * Calculate the integration volume
    * @return integration volume */
   double CalculateIntegrationVolume();

   /**
    * Check availability of integration routine for marginalization */
   bool CheckMarginalizationAvailability(BCMarginalizationMethod type);

   /**
    * Check that indices of parameters to marginalize w/r/t are correct */
   bool CheckMarginalizationIndices(TH1* hist, const std::vector<unsigned> &index);

   /** @} */

protected:

   /**
    * An identification number in case several models exist .*/
   int fID;

   /**
    * A vector of best fit parameters estimated from the global
    * probability and the estimate on their uncertainties */
   std::vector<double> fBestFitParameterErrors;

   /**
    * Minuit */
   TMinuit * fMinuit;

   double fMinuitArglist[2];
   int fMinuitErrorFlag;

   /**
    * Flag for ignoring older results of minimization */
   bool fFlagIgnorePrevOptimization;

   /**
    * Starting temperature for Simulated Annealing */
   double fSAT0;

   /**
    * Minimal/Threshold temperature for Simulated Annealing */
   double fSATmin;

   /**
    * Tree for the Simulated Annealing */
   TTree * fTreeSA;

   /**
    * Flag deciding whether to write SA to file or not. */
   bool fFlagWriteSAToFile;

   int fSANIterations;
   double fSATemperature;
   double fSALogProb;
   std::vector<double> fSAx;

   /**
    * Set of marginalized distributions. */
   std::vector<BCH1D*> fMarginalized1D;
   
   /**
    * Set of marginalized distributions. */
   std::vector<BCH2D*> fMarginalized2D;

 protected:
   /**
    * Determine frequency of output during integration */
   unsigned IntegrationOutputFrequency() const;

   /**
    * Helper methods to unify output for integration methods
    * @param type
    * @param cubatype
    */
   void LogOutputAtStartOfIntegration(BCIntegrationMethod type, BCCubaMethod cubatype);
   void LogOutputAtIntegrationStatusUpdate(BCIntegrationMethod type, double integral, double absprecision, int nIterations);
   void LogOutputAtEndOfIntegration(double integral, double absprecision, double relprecision, int nIterations);

   /**
    * Copy into object
    * @param bcintegrate BCIntegrate object to copy values from */
   void Copy(const BCIntegrate & bcintegrate);

private:
   /**
    * Current integration method */
   BCIntegrate::BCIntegrationMethod fIntegrationMethod;

   /**
    * Current marginalization method */
   BCIntegrate::BCMarginalizationMethod fMarginalizationMethod;

   /**
    * Current Simulated Annealing schedule */
   BCIntegrate::BCSASchedule fSASchedule;

   /**
    * Maximum number of iterations */
   unsigned fNIterationsMin;

   /**
    * Maximum number of iterations */
   unsigned fNIterationsMax;

   /**
    * Maximum number of iterations */
   unsigned fNIterationsPrecisionCheck;

   /**
    * Output frequency during integration */
   unsigned fNIterationsOutput;

   /**
    * Number of iterations in the most recent Monte Carlo integration */
   int fNIterations;

   /**
    * The integral. */
   double fIntegral;
   
   /** Requested relative precision of the integration */
   double fRelativePrecision;

   /** Requested relative precision of the integration */
   double fAbsolutePrecision;

   /**
    * The uncertainty in the most recent Monte Carlo integration */
   double fError;

   /** Cuba integration method */
   BCCubaMethod fCubaIntegrationMethod;
   BCCubaOptions::Vegas fCubaVegasOptions;
   BCCubaOptions::Suave fCubaSuaveOptions;
   BCCubaOptions::Divonne fCubaDivonneOptions;
   BCCubaOptions::Cuhre fCubaCuhreOptions;
};

// ---------------------------------------------------------

#endif
