#ifndef __BCINTEGRATE__H
#define __BCINTEGRATE__H

/**
 * @class BCIntegrate
 * @brief A class for handling numerical operations for models.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This is a base class for a model class. It contains
 * numerical methods to carry out the integration, marginalization,
 * peak finding etc.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCEngineMCMC.h"

#include <TMinuitMinimizer.h>

#include <string>

// forward declarations
class BCIntegrate;

class TFile;
class TH1;
class TH2;
class TTree;

/**
 * Collect all the useful options to tune integration with the CUBA library.
 * Option names match those used in the CUBA manual, where a detailed
 * description is given. Default values are taken from the demo that ships with CUBA.
 */
namespace BCCubaOptions
{

struct General {
    int ncomp, flags, nregions, neval, fail;
    double error, prob;
    General();
protected:
    ~General();
};

struct Vegas : public General {
    int nstart, nincrease, nbatch, gridno;

    Vegas();
};

struct Suave : public General {
    int nmin, nnew;
    double flatness;

    Suave();
};

struct Divonne : public General {
    int key1, key2, key3, maxpass;
    double border, maxchisq, mindeviation;

    Divonne();
};

struct Cuhre : public General {
    int key;

    Cuhre();
};
}

namespace BCMinimizer
{

class Adapter : public ROOT::Math::IMultiGenFunction
{
public:
    Adapter(BCIntegrate& m);
    virtual unsigned int NDim() const;
    virtual ROOT::Math::IMultiGenFunction* Clone() const;

    BCIntegrate* m;
    mutable std::vector<double> par;
private:
    virtual double DoEval (const double* x) const;
};

/**
 * Wrapper to approximate RAII for TMinuitMinimizer which is not
 * copyable unfortunately.
 */
class Wrapper
{
public:
    Wrapper(BCIntegrate& m);
    void Init();
    void Init(const std::vector<double>& start, int printlevel);
    void Reset(BCIntegrate& m);

    TMinuitMinimizer min;
    Adapter adapt;
};
}

// ---------------------------------------------------------

class BCIntegrate : public BCEngineMCMC
{

public:

    /** \name Enumerators */
    /** @{ */

    /**
     * An enumerator for the mode finding algorithm */
    enum BCOptimizationMethod {
        kOptEmpty,                                ///< No optimization method set.
        kOptSimAnn,                               ///< Simulated annealing
        kOptMetropolis,                           ///< Metropolis Hastings
        kOptMinuit,                               ///< ROOT's Minuit
        kOptDefault,                              ///< Default
        NOptMethods                               ///< number of available optimization methods
    };

    /**
     * An enumerator for integration algorithms */
    enum BCIntegrationMethod {
        kIntEmpty,                                ///< No integration method set
        kIntMonteCarlo,                           ///< Sample mean method
        kIntCuba,                                 ///< Use CUBA interface
        kIntGrid,                                 ///< Integration by gridding of parameter space
        kIntLaplace,                              ///< Laplace approximation
        kIntDefault,                              ///< Default
        NIntMethods                               ///< number of available integration methods
    };

    /**
     * An enumerator for marginalization algorithms */
    enum BCMarginalizationMethod {
        kMargEmpty,                               ///< No marginalization method set
        kMargMetropolis,                          ///< Metropolis Hastings
        kMargMonteCarlo,                          ///< Sample mean Monte Carlo
        kMargGrid,                                ///< Marginalization by gridding of parameter space
        kMargDefault,                             ///< Default
        NMargMethods                              ///< number of available marginalization methods
    };

    /**
     * An enumerator for the Simulated Annealing schedule */
    enum BCSASchedule {
        kSACauchy,                                ///< Cauchy scheduler
        kSABoltzmann,                             ///< Boltzman scheduler
        kSACustom,                                ///< Custom scheduler
        NSAMethods                                ///< number of available schedulers
    };

    /**
     * An enumerator for Cuba integration methods */
    enum BCCubaMethod {
        kCubaVegas,                               ///< Vegas
        kCubaSuave,                               ///< Suave
        kCubaDivonne,                             ///< Divonne
        kCubaCuhre,                               ///< Cuhre
        kCubaDefault,                             ///< Default
        NCubaMethods                              ///< number of available CUBA methods
    };

    /** @} */

    /** \name Function pointer types */
    /** @{ */

    /**
     * A pointer for a function that chooses a next random point */
    typedef void (BCIntegrate::*tRandomizer)(std::vector<double>&) const;

    /**
     * A pointer for a function that evaluates at a point */
    typedef double (BCIntegrate::*tEvaluator)(std::vector<double>&, const std::vector<double>&, bool&);

    /**
     * A pointer for a function that updates the integral and absolute precision */
    typedef void (*tIntegralUpdater)(const std::vector<double>&, const int&, double&, double&);

    /** @} */


    /** \name Constructors and destructors */
    /** @{ */

    /**
     * Default constructor */
    BCIntegrate(const std::string& name = "model");

    /**
     * Read in MCMC constructor.
     * @param filename Path of file holding model.
     * @param name Name of model (file should contain TTree's [name]_mcmc and [name]_parameters.\n
     * if empty string is given, properly matching TTrees are searched for in the file.
     * @param loadObservables Flag for whether to load observables for file (true; default) or to let user respecify observables.*/
    BCIntegrate(const std::string& filename, const std::string& name, bool loadObservables = true);

    /**
     * Copy constructor */
    BCIntegrate(const BCIntegrate& other);

    // No assignment operator for abstract class

    /**
     * Destructor */
    virtual ~BCIntegrate() {};

    /** @} */

    /** \name swap*/
    /** @{ */

    friend void swap(BCIntegrate& A, BCIntegrate& B);

    /** @} */

    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return The integral. */
    double GetIntegral() const
    { return fIntegral; }

    /**
     * @return The current optimization method */
    BCIntegrate::BCOptimizationMethod GetOptimizationMethod() const
    { return fOptimizationMethodCurrent; }

    /**
     * @return The current integration method */
    BCIntegrate::BCIntegrationMethod GetIntegrationMethod() const
    { return fIntegrationMethodCurrent; }

    /**
     * @return The current marginalization method */
    BCIntegrate::BCMarginalizationMethod GetMarginalizationMethod() const
    { return fMarginalizationMethodCurrent; }

    /**
     * @return The Simulated Annealing schedule */
    BCIntegrate::BCSASchedule GetSASchedule() const
    { return fSASchedule; }

    /**
     * Fills a vector of random numbers x[i] between fMin[i] and fMax[i] into a vector
     * @param x A vector of doubles to fill*/
    void GetRandomVectorInParameterSpace(std::vector<double>& x) const;

    /**
     * Fills a vector of (flat) random numbers in the limits of the parameters and returns
     * the probability at that point
     * @param x A vector of doubles to fill
     * @return The (unnormalized) probability at the random point */
    double GetRandomPoint(std::vector<double>& x);

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
     * @return Options used for integration with CUBA's Vegas */
    const BCCubaOptions::Vegas& GetCubaVegasOptions() const
    { return fCubaVegasOptions; }

    /**
     * @return Options used for integration with CUBA's Suave */
    const BCCubaOptions::Suave& GetCubaSuaveOptions() const
    { return fCubaSuaveOptions; }

    /**
     * @return Options used for integration with CUBA's Divonne */
    const BCCubaOptions::Divonne& GetCubaDivonneOptions() const
    { return fCubaDivonneOptions; }

    /**
     * @return Options used for integration with CUBA's Cuhre */
    const BCCubaOptions::Cuhre& GetCubaCuhreOptions() const
    { return fCubaCuhreOptions; }

    /**
     * Returns a one-dimensional slice of the pdf at the point and along a specific direction.
     * @param name The name of the model parameter along which the slice is calculated.
     * @param nIterations Add the number of posterior evaluations performed.
     * @param log_max_val Stores the log of the maximum value before normalizing
     * @param parameters The point at which the other parameters are fixed.
     * @param nbins The number of bins of the 1D-histogram.
     * @param normalize Flag for turning on normalization of histogram.
     * @return The slice histogram. */
    TH1* GetSlice(std::vector<unsigned> indices, unsigned& nIterations, double& log_max_val, const std::vector<double> parameters = std::vector<double>(0), int nbins = 0, bool normalize = true);

    /**
     * Returns a one-dimensional slice of the pdf at the point and along a specific direction.
     * @param name The name of the model parameter along which the slice is calculated.
     * @param nIterations Add the number of posterior evaluations performed.
     * @param log_max_val Stores the log of the maximum value before normalizing
     * @param parameters The point at which the other parameters are fixed.
     * @param nbins The number of bins of the 1D-histogram.
     * @param normalize Flag for turning on normalization of histogram.
     * @return The 1D slice. */
    TH1* GetSlice(const std::string& name, unsigned& nIterations, double& log_max_val, const std::vector<double> parameters = std::vector<double>(0), int nbins = 0, bool normalize = true)
    { return GetSlice(fParameters.Index(name), nIterations, log_max_val, parameters, nbins, normalize); }

    /**
     * Returns a one-dimensional slice of the pdf at the point and along a specific direction.
     * @param name The name of the model parameter along which the slice is calculated.
     * @param nIterations Add the number of posterior evaluations performed.
     * @param log_max_val Stores the log of the maximum value before normalizing
     * @param parameters The point at which the other parameters are fixed.
     * @param nbins The number of bins of the 1D-histogram.
     * @param normalize Flag for turning on normalization of histogram.
     * @return The 1D slice. */
    TH1* GetSlice(unsigned index, unsigned& nIterations, double& log_max_val, const std::vector<double> parameters = std::vector<double>(0), int nbins = 0, bool normalize = true)
    { return GetSlice(std::vector<unsigned>(1, index), nIterations, log_max_val, parameters, nbins, normalize); }

    /**
     * Returns a two-dimensional slice of the pdf at the point and along two specified directions.
     * @param name1 The name of the first model parameter along which the slice is calculated.
     * @param name2 The name of the second model parameter along which the slice is calculated.
     * @param nIterations Add the number of posterior evaluations performed.
     * @param log_max_val Stores the log of the maximum value before normalizing
     * @param parameters The point at which the other parameters are fixed.
     * @param nbins The number of bins on each axis of the 2D-histogram.
     * @param normalize Flag for turning on normalization of histogram.
     * @return The 2D slice. */
    TH2* GetSlice(const std::string& name1, const std::string& name2, unsigned& nIterations, double& log_max_val, const std::vector<double> parameters = std::vector<double>(0), int nbins = 0, bool normalize = true)
    { return GetSlice(fParameters.Index(name1), fParameters.Index(name2), nIterations, log_max_val, parameters, nbins, normalize); }

    /**
     * Returns a two-dimensional slice of the pdf at the point and along two specified directions.
     * @param name1 The name of the first model parameter along which the slice is calculated.
     * @param name2 The name of the second model parameter along which the slice is calculated.
     * @param nIterations Add the number of posterior evaluations performed.
     * @param log_max_val Stores the log of the maximum value before normalizing
     * @param parameters The point at which the other parameters are fixed.
     * @param nbins The number of bins on each axis of the 2D-histogram.
     * @param normalize Flag for turning on normalization of histogram.
     * @return The 2D slice. */
    TH2* GetSlice(unsigned index1, unsigned index2, unsigned& nIterations, double& log_max_val, const std::vector<double> parameters = std::vector<double>(0), int nbins = 0, bool normalize = true);

    /**
     * @return The uncertainty in the most recent Monte Carlo integration */
    double GetError() const
    { return fError; }

    /**
     * @return Minuit minimizer used for mode finding. The object is in a
     * valid but unusable state if no mode finding has been
     * performed. */
    TMinuitMinimizer& GetMinuit()
    { return fMinimizer.min; }

    /**
     * Returns the Simulated Annealing starting temperature. */
    double GetSAT0() const
    { return fSAT0; }

    /**
     * Returns the Simulated Annealing threshhold temperature. */
    double GetSATmin() const
    { return fSATmin; }

    /**
     * @return vector of parameter and observable values at global mode. */
    virtual const std::vector<double>& GetBestFitParameters() const;

    /**
     * Returns the set of errors on the values of the parameters at the global mode */
    const std::vector<double>& GetBestFitParameterErrors() const
    { return fBestFitParameterErrors; }

    /**
     * Returns the posterior at the mode.
     * @return the posterior. */
    double GetLogMaximum() const
    { return fLogMaximum; };

    /** @} */

    /** \name Member functions (set) */
    /** @{ */
    /**
     * @param flag Flag whether or not to ignore result of previous mode finding */
    void SetFlagIgnorePrevOptimization(bool flag)
    { fFlagIgnorePrevOptimization = flag; }

    /**
     * @param method The current optimization method */
    void SetOptimizationMethod(BCIntegrate::BCOptimizationMethod method)
    { fOptimizationMethodCurrent = method; }

    /**
     * @param method The current integration method */
    void SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method);

    /**
     * @param method The current marginalization method */
    void SetMarginalizationMethod(BCIntegrate::BCMarginalizationMethod method)
    { fMarginalizationMethodCurrent = method; }

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
     * Set options for CUBA's Vegas.
     * @param options Options for CUBA*/
    void SetCubaOptions(const BCCubaOptions::Vegas& options)
    { fCubaVegasOptions = options; }

    /**
     * Set options for CUBA's Suave.
     * @param options Options for CUBA*/
    void SetCubaOptions(const BCCubaOptions::Suave& options)
    { fCubaSuaveOptions = options; }

    /**
     * Set options for CUBA's Divonne.
     * @param options Options for CUBA*/
    void SetCubaOptions(const BCCubaOptions::Divonne& options)
    { fCubaDivonneOptions = options; }

    /**
     * Set options for CUBA's Cuhre.
     * @param options Options for CUBA*/
    void SetCubaOptions(const BCCubaOptions::Cuhre& options)
    { fCubaCuhreOptions = options; }

    /**
     * Set starting temperature for Simulated Annealing
     * @param T0 starting temperature. */
    void SetSAT0(double T0)
    { fSAT0 = T0; }

    /**
     * Set threshold temperature for Simulated Annealing
     * @param Tmin threshold temperature. */
    void SetSATmin(double Tmin)
    { fSATmin = Tmin; }

    /**
     * Turn on/off writing of simulated annealing to root file.
     * If setting true, use function with filename arguments.
     * @param flag Flag for writing simulated annealing to ROOT file (true) or not (false). */
    void WriteSAToFile(bool flag);

    /**
     * Turn on writing of simulated annealing to root file.
     * @param filename Name of file to.
     * @param file-open options (TFile), must be "NEW", "CREATE", "RECREATE", or "UPDATE" (i.e. writeable).
     * @param autoclose Toggle autoclosing of file after simulated annealing. */
    void WriteSAToFile(const std::string& filename, const std::string& option, bool autoclose = true);

    /**
     * Close SA output file. */
    void CloseSAOutputFile();

    /**
     * Getter for the tree containing the  Simulated Annealing  chain. */
    TTree* GetSATree()
    { return fSATree; }

    /**
     * Initialization of the tree for the Simulated Annealing
     * @param replacetree Whether to delete and recreate tree object if already existing.
     * @param replacefile Whether to delete and recreate file object if already existing. */
    void InitializeSATree(bool replacetree = false, bool replacefile = false);

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Evaluate the unnormalized probability at a point in parameter space.
     * Method needs to be overloaded by the user.
     * @param x The point in parameter space
     * @return The unnormalized probability */
    virtual double Eval(const std::vector<double>& x) = 0;

    /**
     * Evaluate the natural logarithm of the Eval function. For better numerical
     * stability, this method should also be overloaded by the user.
     * @param x The point in parameter space
     * @return log(Eval(x)) */
    virtual double LogEval(const std::vector<double>& x)
    { return log(Eval(x)); }

    /**
     * Performs integration. */
    double Normalize()
    { return Integrate(); };

    /**
     * Does the integration over the un-normalized probability.
     * @param intmethod The integration method to used
     * @return The normalization value */
    double Integrate(BCIntegrationMethod intmethod);

    /**
     * Perform the integration
     * @return the integral */
    double Integrate();

    /**
     * Does the integration over the un-normalized probability.
     * @param type The integration method to used (for printing out status updates by name)
     * @param randomizer Pointer to function to choose next random point
     * @param evaluator Pointer to function to evaluate point
     * @param updater Pointer to function to update integral and precision
     * @param sums Vector of doubles holding values used in integral calculation
     * @param x An initial point to start integration routine at
     * @return The integral value */
    double Integrate(BCIntegrationMethod type, tRandomizer randomizer, tEvaluator evaluator, tIntegralUpdater updater, std::vector<double>& sums);

    /**
     * Evaluates integrator */
    double EvaluatorMC(std::vector<double>& sums, const std::vector<double>& point, bool& accepted);

    /**
     * Updates info about integrator */
    static void IntegralUpdaterMC(const std::vector<double>& sums, const int& nIterations, double& integral, double& absprecision);

    /**
     * Marginalize all probabilities wrt. single parameters and all combinations
     * of two parameters. The individual distributions can be retrieved using
     * the GetMarginalized method.
     * @return Total number of marginalized distributions */
    int MarginalizeAll();

    /**
     * Marginalize all probabilities wrt. single parameters and all combinations
     * of two parameters. The individual distributions can be retrieved using
     * the GetMarginalized method.
     * @param margmethod the marginalization method.
     * @return Total number of marginalized distributions */
    int MarginalizeAll(BCMarginalizationMethod margmethod);

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
     * Do the mode finding using a method set via SetOptimizationMethod.
     * Default is Minuit. The mode can be extracted using the GetBestFitParameters() method.
     *
     * A starting point for the mode finding can be specified for Minuit. If not
     * specified, the previously found maximum (typically from marginalization)
     * is used as an initial point. If that is not available,
     * then the Minuit default will be used (center of the parameter space).
     * @return The mode found.
     * @note The result may not coincide with the result of GetBestFitParameters()
     * if a previous optimization found a better value. */
    std::vector<double> FindMode(std::vector<double> start = std::vector<double>());

    /**
     * Find mode using a specific method. The original method will be reset.
     * @param optmethod the optimization method
     * @param start the starting point for the optimization algorithm
     * @return the mode
     * @see std::vector<double> FindMode(std::vector<double> start = std::vector<double>(0)); */
    std::vector<double> FindMode(BCIntegrate::BCOptimizationMethod optmethod, std::vector<double> start = std::vector<double>());

    /**
     * Temperature annealing schedule for use with Simulated Annealing.
     * Delegates to the appropriate method according to
     * fSASchedule.
     * @param t iterator for lowering the temperature over time. */
    double SATemperature(double t) const;

    /**
     * Temperature annealing schedule for use with Simulated Annealing.
     * This method is used for Boltzmann annealing schedule.
     * @param t iterator for lowering the temperature over time. */
    double SATemperatureBoltzmann(double t) const;

    /**
     * Temperature annealing schedule for use with Simulated Annealing.
     * This method is used for Cauchy annealing schedule.
     * @param t iterator for lowering the temperature over time. */
    double SATemperatureCauchy(double t) const;

    /**
     * Temperature annealing schedule for use with Simulated Annealing.
     * This is a virtual method to be overridden by a user-defined
     * custom temperature schedule.
     * @param t iterator for lowering the temperature over time. */
    virtual double SATemperatureCustom(double t) const;

    /**
     * Generates a new state in a neighbourhood around x that is to be
     * accepted or rejected by the Simulated Annealing algorithm.
     * Delegates the generation to the appropriate method according
     * to fSASchedule.
     * @param x last state.
     * @param t time iterator to determine current temperature. */
    std::vector<double> GetProposalPointSA(const std::vector<double>& x, int t) const;

    /**
     * Generates a new state in a neighbourhood around x that is to be
     * accepted or rejected by the Simulated Annealing algorithm.
     * This method is used for Boltzmann annealing schedule.
     * @param x last state.
     * @param t time iterator to determine current temperature. */
    std::vector<double> GetProposalPointSABoltzmann(const std::vector<double>& x, int t) const;

    /**
     * Generates a new state in a neighbourhood around x that is to be
     * accepted or rejected by the Simulated Annealing algorithm.
     * This method is used for Cauchy annealing schedule.
     * @param x last state.
     * @param t time iterator to determine current temperature. */
    std::vector<double> GetProposalPointSACauchy(const std::vector<double>& x, int t) const;

    /**
     * Generates a new state in a neighbourhood around x that is to be
     * accepted or rejected by the Simulated Annealing algorithm.
     * This is a virtual method to be overridden by a user-defined
     * custom proposal function.
     * @param x last state.
     * @param t time iterator to determine current temperature. */
    virtual std::vector<double> GetProposalPointSACustom(const std::vector<double>& x, int t) const;

    /**
     * Generates a uniform distributed random point on the surface of
     * a fNvar-dimensional Hypersphere.
     * Used as a helper to generate proposal points for Cauchy annealing. */
    std::vector<double> SAHelperGetRandomPointOnHypersphere() const;

    /**
     * Generates the radial part of a n-dimensional Cauchy distribution.
     * Helper function for Cauchy annealing. */
    double SAHelperGetRadialCauchy() const;

    /**
     * Returns the Integral of sin^dim from 0 to theta.
     * Helper function needed for generating Cauchy distributions. */
    double SAHelperSinusToNIntegral(int dim, double theta) const;

    /**
     * Reset all information on the best fit parameters. */
    virtual void ResetResults();

    /**
     * Return string with the name for a given integration type.
     * @param type code for the integration type
     * @return string containing the name of the integration type */
    std::string DumpIntegrationMethod(BCIntegrationMethod type) const;

    /**
     * Return string with the name for the currently set integration type.
     * @return string containing the name of the integration type */
    std::string DumpCurrentIntegrationMethod() const
    { return DumpIntegrationMethod(fIntegrationMethodCurrent); }

    /**
     * Return string with the name for the previously used integration type.
     * @return string containing the name of the integration type */
    std::string DumpUsedIntegrationMethod() const
    { return DumpIntegrationMethod(fIntegrationMethodUsed); }

    /**
     * Return string with the name for a given marginalization type.
     * @param type code for the marginalization type
     * @return string containing the name of the marginalization type */
    std::string DumpMarginalizationMethod(BCMarginalizationMethod type) const;

    /**
     * Return string with the name for the currently set marginalization type.
     * @return string containing the name of the marginalization type */
    std::string DumpCurrentMarginalizationMethod() const
    { return DumpMarginalizationMethod(fMarginalizationMethodCurrent); }

    /**
     * Return string with the name for the marginalization type used.
     * @return string containing the name of the marginalization type */
    std::string DumpUsedMarginalizationMethod() const
    { return DumpMarginalizationMethod(fMarginalizationMethodUsed); }

    /**
     * Return string with the name for a given optimization type.
     * @param type code for the optimization type
     * @return string containing the name of the optimization type */
    std::string DumpOptimizationMethod(BCOptimizationMethod type) const;

    /**
     * Return string with the name for the currently set optimization type.
     * @return string containing the name of the optimization type */
    std::string DumpCurrentOptimizationMethod() const
    { return DumpOptimizationMethod(fOptimizationMethodCurrent); }

    /**
     * Return string with the name for the optimization type used to find the current mode.
     * @return string containing the name of the optimization type */
    std::string DumpUsedOptimizationMethod() const
    { return DumpOptimizationMethod(fOptimizationMethodUsed); }

    /**
     * Return string with the name for a given Cuba integration type.
     * @param type code for the Cuba integration type
     * @return string containing the name of the Cuba integration type */
    std::string DumpCubaIntegrationMethod(BCCubaMethod type) const;

    /**
     * Return string with the name for the currently set Cuba integration type.
     * @return string containing the name of the Cuba integration type */
    std::string DumpCubaIntegrationMethod() const
    { return DumpCubaIntegrationMethod(fCubaIntegrationMethod); }

    /**
     * Set best fit parameters values
     * @param x Parameter set to designate as best fit parameters. */
    void SetBestFitParameters(const std::vector<double>& x);

    /**
     * Set best fit parameters if best fit
     * @param new_value is the value of the function at x
     * @param old_value is the old best fit value, updated to new_value, if it is larger */
    void SetBestFitParameters(const std::vector<double>& x, const double& new_value, double& old_value);

    /**
     * Check availability of integration routine for marginalization
     * @return availability of marginalization method */
    bool CheckMarginalizationAvailability(BCMarginalizationMethod type);

    /**
     * Check that indices of parameters to marginalize w/r/t are correct */
    bool CheckMarginalizationIndices(TH1* hist, const std::vector<unsigned>& index);

    /** @} */

protected:
    /**
     * Determine frequency of output during integration */
    unsigned IntegrationOutputFrequency() const;

    /**
     * Print best fit to log */
    virtual void PrintBestFitSummary() const;

    /**
     * Get string summarizing best fit for single variable.
     * @param i Index of variable to summarize. */
    virtual std::string GetBestFitSummary(unsigned i) const;

    /**
     * Print marginalization to log. */
    virtual void PrintMarginalizationSummary() const;

    /**
     * Flag for ignoring older results of optimization */
    bool fFlagIgnorePrevOptimization;

    /**
     * Starting temperature for Simulated Annealing */
    double fSAT0;

    /**
     * Minimal/Threshold temperature for Simulated Annealing */
    double fSATmin;

    /**
     * Tree for the Simulated Annealing */
    TTree* fSATree;

    /**
     * Flag deciding whether to write simulated annealing to file or not. */
    bool fFlagWriteSAToFile;

    /**
     * Number of iterations for simualted annealing. */
    int fSANIterations;

    /**
     * Current temperature of simulated annealing algorithm. */
    double fSATemperature;

    /**
     * Log probability of current simulated annealing iteration. */
    double fSALogProb;

    /**
     * Current simulated annealing parameter point. */
    std::vector<double> fSAx;

    /**
     * Helper method to output at beginning of integration. */
    void LogOutputAtStartOfIntegration(BCIntegrationMethod type, BCCubaMethod cubatype);

    /**
     * Helper method to output integration status. */
    void LogOutputAtIntegrationStatusUpdate(BCIntegrationMethod type, double integral, double absprecision, int nIterations);

    /**
     * Helper method to output at end of integration. */
    void LogOutputAtEndOfIntegration(double integral, double absprecision, double relprecision, int nIterations);

    /**
     * flag indicating if the model was marginalized */
    bool fFlagMarginalized;

    /**
     * Output file for writing SA Tree. */
    TFile* fSAOutputFile;

    /**
     * Output filename for writing SA Tree. */
    std::string fSAOutputFilename;

    /**
     * Output file open option for for writing SA Tree. */
    std::string fSAOutputFileOption;

    /**
     * flag for autoclosing SA output file. */
    bool fSAOutputFileAutoclose;

private:

    ///> Wrapper to run minuit
    BCMinimizer::Wrapper fMinimizer;

    /**
     * Does the mode finding using Minuit. If starting point is not specified,
     * finding will start from the center of the parameter space.
     * @param start point in parameter space from which the mode finding is started.
     * @param printlevel The print level.
     * @param mode a reference to a vector holding the mode
     * @param errors a reference to a vector holding the errors
     * @return The mode found.
     * @note The result may not coincide with the result of GetBestFitParameters()
     * if a previous optimization found a better value. */
    std::vector<double> FindModeMinuit(std::vector<double>& mode, std::vector<double>& errors, std::vector<double> start = std::vector<double>(0), int printlevel = -1);

    /**
     * Does the mode finding using Markov Chain Monte Carlo (prerun only!)
     * @param mode a reference to a vector holding the mode
     * @param errors a reference to a vector holding the errors
     * @return The mode.
     * @note The result may not coincide with the result of GetBestFitParameters()
     * if a previous optimization found a better value. */
    std::vector<double> FindModeMCMC(std::vector<double>& mode, std::vector<double>& errors);

    /**
     * Does the mode finding using Simulated Annealing. If starting point
     * is not specified, finding will start from the center of the
     * parameter space.
     * @param mode a reference to a vector holding the mode
     * @param errors a reference to a vector holding the errors
     * @param start point in parameter space from thich the mode finding is started.
     * @return The mode.
     * @note The result may not coincide with the result of GetBestFitParameters()
     * if a previous optimization found a better value.*/
    std::vector<double> FindModeSA(std::vector<double>& mode, std::vector<double>& errors, std::vector<double> start = std::vector<double>(0));

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
    static int CubaIntegrand(const int* ndim, const double xx[], const int* ncomp, double ff[], void* userdata);

    /**
     * Integrate using the slice method
     * @return the integral; */
    double IntegrateSlice();

    /**
     * Integrate using the Laplace approximation.
     *
     * Take the result of a previous successful minuit run to estimate
     * the covariance matrix. Else it runs minuit (again).
     *
     * @return the integral; */
    double IntegrateLaplace();

    /**
     * Current mode finding method */
    BCIntegrate::BCOptimizationMethod fOptimizationMethodCurrent;

    /**
     * Method with which the global mode was found (can differ from
     * fOptimization method in case more than one algorithm is used). */
    BCIntegrate::BCOptimizationMethod fOptimizationMethodUsed;

    /**
     * Current integration method */
    BCIntegrate::BCIntegrationMethod fIntegrationMethodCurrent;

    /**
     * Integration method used for the current results */
    BCIntegrate::BCIntegrationMethod fIntegrationMethodUsed;

    /**
     * Current marginalization method */
    BCIntegrate::BCMarginalizationMethod fMarginalizationMethodCurrent;

    /**
     * Marginalization method used for the current results */
    BCIntegrate::BCMarginalizationMethod fMarginalizationMethodUsed;

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
     * Number of iterations in the most recent Monte Carlo integration */
    int fNIterations;

    /**
     * A vector of best fit parameters found by MCMC */
    std::vector<double> fBestFitParameters;

    /**
     * A vector of estimates on the uncertainties */
    std::vector<double> fBestFitParameterErrors;

    /**
     * The function value at the mode on the @em log scale */
    double fLogMaximum;

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

    BCCubaOptions::Vegas fCubaVegasOptions;			///< CUBA Vegas options
    BCCubaOptions::Suave fCubaSuaveOptions;			///< CUBA Suave options
    BCCubaOptions::Divonne fCubaDivonneOptions;		///< CUBA Divonne options
    BCCubaOptions::Cuhre fCubaCuhreOptions;			///< CUBA Cuhre options
};

#endif
