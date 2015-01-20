#ifndef __BCMODEL__H
#define __BCMODEL__H

/*!
 * \class BCModel
 * \brief The base class for all user-defined models.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Daniel Greenwald
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
#include "BCDataSet.h"

#include <string>

// ROOT classes
class TNamed;
class TH1;
class TF1;

//BAT classes
class BCDataPoint;
class BCParameter;
class BCH1D;
class BCH2D;
class BCPriorModel;


const int MAXNDATAPOINTVALUES = 20;

// ---------------------------------------------------------

class BCModel : public BCIntegrate
{

   public:

	/** \name Enumerators  */
	/** @{ */

	/** An enumerator for the knowledge update drawing style presets. */
	enum BCKnowledgeUpdateDrawingStyle {
		kKnowledgeUpdateDefaultStyle      = 0,
		kKnowledgeUpdateDetailedPosterior = 1,
		kKnowledgeUpdateDetailedPrior     = 2
	};

	/** @} */
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
			 * Read in MCMC constructor.
			 * @param filename Path of file holding model.
			 * @param name Name of model (file should contain TTree's [name]_mcmc and [name]_parameters.\n
			 * if empty string is given, properly matching TTrees are searched for in the file.
			 * @param reuseObservables Flag for whether to load observables for file (true; default) or to let user respecify observables.*/
	    BCModel(std::string filename, std::string name, bool reuseObservables=true);

      /**
       * The default destructor. */
      virtual ~BCModel();

      /** @} */
      /** \name Assignment operators */
      /** @{ */

      /**
       * Defaut assignment operator */
      BCModel & operator = (const BCModel & bcmodel)
	      { Copy(bcmodel); return *this; }

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

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
       * @return The number of data points in the current data set. */
	    unsigned GetNDataPoints() const
	      { return (fDataSet) ? fDataSet->GetNDataPoints() : 0; }

			/**
			 * @return The number of degrees of freedom:
			 * (number of data points) - (number of free parameters) */
			int GetNDoF() const
			  { return GetNDataPoints() - fParameters.GetNFreeParameters(); }

	    /**
	     * @return BCPriorModel. */
	    virtual BCPriorModel * GetPriorModel(bool prepare=true, bool call_likelihood=false);

	    /**
	     * @return BCH1D object for controlling drawing options of priors in knowledge update plots. */
	    BCH1D * GetBCH1DPriorDrawingOptions()
	    { return fBCH1DPriorDrawingOptions; }

	    /**
	     * @return BCH2D object for controlling drawing options of priors in knowledge update plots. */
	    BCH2D * GetBCH2DPriorDrawingOptions()
	    { return fBCH2DPriorDrawingOptions; }

	    /**
	     * @return BCH1D object for controlling drawing options of posteriors in knowledge update plots. */
	    BCH1D * GetBCH1DPosteriorDrawingOptions()
	    { return fBCH1DPosteriorDrawingOptions; }

	    /**
	     * @return BCH2D object for controlling drawing options of posteriors in knowledge update plots. */
	    BCH2D * GetBCH2DPosteriorDrawingOptions()
	    { return fBCH2DPosteriorDrawingOptions; }

	/**
	 * @return flag for order of drawing prior first, posterior second (true); or the other way around (false). */
	bool GetDrawPriorPosteriorNormalOrder()
	{ return fPriorPosteriorNormalOrder; }

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

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
       * Sets the data set. Model does not own data set!
       * @param dataset A data set */
      void SetDataSet(BCDataSet* dataset)
         { fDataSet = dataset; }

      /**
       * Set prior for a parameter.
       * @param index The parameter index
       * @param f A pointer to a function describing the prior
       * @return success of action. */
      bool SetPrior(unsigned index, TF1* f)
	       { return GetParameter(index) ? GetParameter(index)->SetPrior(f) : false; }

      /**
       * Set prior for a parameter.
       * @param name The parameter name
       * @param f A pointer to a function describing the prior
       * @return success of action. */
      bool SetPrior(const char* name, TF1* f)
	      { return SetPrior(fParameters.Index(name), f); }

      /**
       * Set delta-function prior for a parameter. Note: this sets the
       * parameter range to the specified value. The old parameter range
       * is lost.
       * @param index The parameter index
       * @param value The position of the delta function.
       * @return success of action. */
	    bool SetPriorDelta(unsigned index, double value)
	      { return GetParameter(index) ? GetParameter(index)->Fix(value) : false; }

      /**
       * Set delta-function prior for a parameter. Note: this sets the
       * parameter range to the specified value. The old parameter range
       * is lost.
       * @param name The parameter name
       * @param value The position of the delta function.
       * @return success of action. */
      bool SetPriorDelta(const char* name, double value)
	       { return SetPriorDelta(fParameters.Index(name),value); }

      /**
       * Set Gaussian prior for a parameter.
       * @param index The parameter index
       * @param mean The mean of the Gaussian
       * @param sigma The sigma of the Gaussian
       * @return success of action. */
      bool SetPriorGauss(unsigned index, double mean, double sigma)
	      { return GetParameter(index) ? GetParameter(index)->SetPriorGauss(mean,sigma) : false; }

      /**
       * Set Gaussian prior for a parameter.
       * @param name The parameter name
       * @param mean The mean of the Gaussian
       * @param sigma The sigma of the Gaussian
       * @return success of action. */
      bool SetPriorGauss(const char* name, double mean, double sigma)
	       { return SetPriorGauss(fParameters.Index(name), mean, sigma); }

      /**
       * Set Gaussian prior for a parameter with two different widths.
       * @param index The parameter index
       * @param mean The mean of the Gaussian
       * @param sigma_below Standard deviation below mean.
       * @param sigma_above Standard deviation above mean.
       * @return success of action. */
     bool SetPriorGauss(unsigned index, double mean, double sigma_below, double sigma_above)
	      { return GetParameter(index) ? GetParameter(index)->SetPriorGauss(mean,sigma_below,sigma_above) : false; }

      /**
       * Set Gaussian prior for a parameter with two different widths.
       * @param name The parameter name
       * @param mean The mean of the Gaussian
       * @param sigmadown The sigma (down) of the Gaussian
       * @param sigmaup The sigma (up)of the Gaussian
       * @return success of action. */
      int SetPriorGauss(const char* name, double mean, double sigmadown, double sigmaup)
	       {	return SetPriorGauss(fParameters.Index(name), mean, sigmadown, sigmaup); }

      /**
       * Set prior for a parameter.
       * @param index parameter index
       * @param h pointer to a histogram describing the prior
       * @param interpolate whether or not to use linear interpolation
       * @return success of action. */
      bool SetPrior(unsigned index, TH1 * h, bool interpolate=false)
	       { return GetParameter(index) -> SetPrior(h,interpolate); }

      /**
       * Set prior for a parameter.
       * @param name parameter name
       * @param h pointer to a histogram describing the prior
       * @param interpolate whether or not to use linear interpolation
       * @return success of action. */
	    bool SetPrior(const char* name, TH1 * h, bool interpolate=false)
	       { return SetPrior(fParameters.Index(name),h,interpolate); }

      /**
       * Set constant prior for this parameter
       * @param index the index of the parameter
       * @return success of action. */
      bool SetPriorConstant(unsigned index)
	       { return GetParameter(index) -> SetPriorConstant(); }

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
       * @return success of action. */
	    bool SetPriorConstantAll()
	       { return fParameters.SetPriorConstantAll(); }

	    /**
	     * Set default drawing options for knowledge update plots. */
	    void SetKnowledgeUpdateDrawingStyle(BCModel::BCKnowledgeUpdateDrawingStyle style=BCModel::kKnowledgeUpdateDefaultStyle);

	    /**
	     * Set drawing of prior first, posterior second (true), or reverse (false) for knowledge update plots. */
	    void SetDrawPriorPosteriorNormalOrder(bool b=true)
	      { fPriorPosteriorNormalOrder = b; }

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
      double APrioriProbability(const std::vector<double> &parameters)
	       { return exp(this->LogAPrioriProbability(parameters)); }
	
      /**
       * Returns natural logarithm of the prior probability.
			 * Assumes prior factorizes into independent parts set by user for each parameter.
       * Method can be overloaded by the user (especially for nonfactorizable priors)
       * @param parameters A set of parameter values
       * @return The prior probability p(parameters)
       * @see GetPrior(std::vector<double> parameters) */
      virtual double LogAPrioriProbability(const std::vector<double> &parameters)
	      { return fParameters.GetLogPrior(parameters); }

      /**
       * Returns the likelihood
       * @param params A set of parameter values
       * @return The likelihood */
      virtual double Likelihood(const std::vector<double> &params)
	      { return exp(LogLikelihood(params)) ; }
	
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
      double ProbabilityNN(const std::vector<double> &params)
	      { return exp(LogProbabilityNN(params)); }

      /**
       * Returns the natural logarithm of likelihood times prior probability given
       * a set of parameter values
       * @param parameters A set of parameter values
       * @return The likelihood times prior probability */
	    double LogProbabilityNN(const std::vector<double> &parameters);


      /**
       * Returns the a posteriori probability given a set of parameter values
       * @param parameters A set of parameter values
       * @return The a posteriori probability */
      double Probability(const std::vector<double> &parameter)
	      { return exp(LogProbability(parameter)); }

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
      virtual double Eval(const std::vector<double> &parameters)
			{ return exp(LogEval(parameters)); }

      /**
       * Overloaded function to evaluate integral. */
      virtual double LogEval(const std::vector<double> &parameters)
	      { return LogProbabilityNN(parameters); }
	
	    /**
       * Initialize the trees containing the Markov chains and parameter info. */
	    virtual void InitializeMarkovChainTree(bool replacetree=false, bool replacefile=false);

      /**
       * Calculates the matrix element of the Hessian matrix
       * @param parameter1 The parameter for the first derivative
       * @param parameter2 The parameter for the first derivative
       * @return The matrix element of the Hessian matrix */
      double HessianMatrixElement(const BCParameter * parameter1, const BCParameter * parameter2, std::vector<double> point);

      /**
       * Prints a short summary of the fit results on the screen. */
      void PrintShortFitSummary(int chi2flag=0);

      /**
       * Prints matrix elements of the Hessian matrix
       * @param parameters The parameter values at which point to evaluate the matrix */
      void PrintHessianMatrix(std::vector<double> parameters);

      /**
       * Draw a comparison of the prior knowledge to the posterior
       * knowledge for each parameter.
       * @return An error flag. */
	    virtual int DrawKnowledgeUpdatePlot1D(unsigned index, bool flag_slice_post=false, bool flag_slice_prior=false);
	
      /**
       * Draw a comparison of the prior knowledge to the posterior.
       * @return An error flag. */
	    virtual int DrawKnowledgeUpdatePlot2D(unsigned index1, unsigned index2, bool flag_slice=false);

      /**
       * Print a comparison of the prior knowledge to the posterior
       * knowledge for each parameter.
       * @return An error flag. */
    	virtual int PrintKnowledgeUpdatePlots(const char * filename = "update.pdf", unsigned hdiv=1, unsigned vdiv=1, bool flag_slice=false, bool call_likelihood=false);


   /** @} */

   protected:
      /**
       * The model prior probability. */
      double fModelAPriori;

      /**
       * The model a posteriori probability. */
      double fModelAPosteriori;

      /**
       * A data set */
      BCDataSet * fDataSet;

      std::vector<bool> fDataFixedValues;

	    /**
	     * BCPriorModel object for drawing of knowledge update, and saving of samples according to prior.*/
	    BCPriorModel * fPriorModel;

			/**
			 * knowledge update plot 1D prior options. */
			BCH1D * fBCH1DPriorDrawingOptions;
			
			/**
			 * knowledge update plot 2D prior options. */
			BCH2D * fBCH2DPriorDrawingOptions;
			
			/**
			 * knowledge update plot 1D posterior options. */
			BCH1D * fBCH1DPosteriorDrawingOptions;

			/**
			 * knowledge update plot 2D posterior options. */
			BCH2D * fBCH2DPosteriorDrawingOptions;

			/**
			 * flag for ordering of drawing of prior and posterior in knowledge update plots. */
			bool fPriorPosteriorNormalOrder;

};

// ---------------------------------------------------------

#endif
