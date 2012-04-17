#ifndef __BCTEMPLATEFITTER__H
#define __BCTEMPLATEFITTER__H

/*!
 * \class BCTemplateFitter
 * This class is DEPRECETED and it will be removed in the future version
 * of BAT. Use BCMTF instead.
 *
 * This class can be used for fitting several template
 * histograms to a data histogram. The templates are assumed to have
 * no statistical uncertainty whereas the data are assumed to have
 * Poissonian fluctuations in each bin. Several methods to judge the
 * validity of the model are available.
 * \brief DEPRECETED class for fitting several templates to a data set.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \date 10.04.2010
 */

/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>

#include <TH1D.h>
#include <TF1.h>

#include "../../BAT/BCModel.h"

class TH2D;

// ---------------------------------------------------------

class BCTemplateFitter : public BCModel
{
   public:

      /**
       * The default constructor.
       */
      BCTemplateFitter();

      /**
       * A constructor.
       */
      BCTemplateFitter(const char * name);

      /**
       * The default destructor.
       */
      ~BCTemplateFitter();

      /**
       * Return the number of templates.
       */
      int GetNTemplates()
         { return int(fTemplateParIndexContainer.size()); }

      /**
       * Return the number of templates of a certain type
       * @param templatetype the type of the template
       */
      int GetNTemplatesType(int templatetype);

      /**
       * Return the number of sources of systematic uncertainties.
       */
      int GetNSystErrors()
         { return int(fSystErrorParIndexContainer.size()); }

      /**
       * Return the number of degrees of freedom.
       */
      int GetNDF()
         { return fHistData.GetNbinsX() - GetNParameters(); }

      /**
       * Return the number of ratios to calculate.
       */
      int GetNRatios()
         { return int(fHistRatios1D.size()); }

      /**
       * Return the vector of histgrams containing the ratios.
       */
      std::vector<TH1D> GetHistRatios1D()
         { return fHistRatios1D; }

      /**
       * Return the container of template histograms.
       */
      std::vector<TH1D> GetTemplateHistogramContainer()
         { return fTemplateHistogramContainer; };

      /**
       * Return the container of template functions.
       */
      std::vector<TF1> GetTemplateFunctionContainer()
         { return fTemplateFunctionContainer; };

      /**
       * Return a ratio histogram
       */
      TH1D GetHistRatio1D(int index)
         { return fHistRatios1D.at(index); }

      TH2D* GetUncertaintyHistogramExp()
      { return fUncertaintyHistogramExp; };

      TH2D* GetUncertaintyHistogramObsPosterior()
      { return fUncertaintyHistogramObsPosterior; };

      /**
       * Return the index of a template.
       * @param name The template name.
       */
      int GetIndexTemplate(const char * name);

      /**
       * Return the index of a source of systematic uncertainty.
       * @param name The name of the source.
       */
      int GetIndexSystError(const char * name);

      /**
       * Return the parameter index corresponding to a template.
       * @param name The template name.
       */
      int GetParIndexTemplate(const char * name);

      /**
       * Return the parameter index corresponding to a template.
       * @param index The template index.
       */
      int GetParIndexTemplate(int index);

      /**
       * Return the parameter index corresponding to an efficiency.
       * @param name The name of the template associated with the efficiency.
       */
      int GetParIndexEff(const char * name);

      /**
       * Return the parameter index corresponding to an efficiency.
       * @param index The index of the template associated with the efficiency.
       */
      int GetParIndexSystError(const char * name);

      /**
       * Return a template histogram.
       * @param name The template name.
       */
      //      TH1D GetTemplate(const char * name);

      /**
       * Return the histogram containing the data.
       */
      TH1D GetData()
         { return fHistData; }

      /**
       * Get histogram index of template.
       * @templatenumber The number of the template
       */
      int GetHistogramIndex(int templatenumber)
         { return fHistogramIndexContainer.at(templatenumber); };

      /**
       * Get function index of template.
       * @templatenumber The number of the template
       */
      int GetFunctionIndex(int templatenumber)
         { return fFunctionIndexContainer.at(templatenumber); };

      /**
       * Return container of parameter indeces for templates
       */
      std::vector<int> GetTemplateParIndexContainer()
         { return fTemplateParIndexContainer; };

      /**
       * Return the expectation value in a certain bin.
       * @param binnumber The number of the bin.
       * @param parameters The parameter values.
       */
      double Expectation(int binnumber, const std::vector<double>& parameters);

      /**
       * Return the efficiency of a template in a certain bin.
       * @param templatenumber The number of the template.
       * @param binnumber The bin number.
       * @param parameters The parameter values.
       */
      double TemplateEfficiency(int templatenumber, int binnumber, const std::vector<double>& parameters);

      /**
       * Return the probability of a template in a certain bin.
       * @param templateindex The number of the template.
       * @param binnumber The bin number.
       * @param parameters The parameter values.
       */
      double TemplateProbability(int templatenumber, int binnumber, const std::vector<double>& parameters);


      /**
       * Set a flag for having physical limits (expectation values greater
       * or equal to 0).
       * @param flag The flag.
       */
      void SetFlagPhysicalLimits(bool flag)
         { fFlagPhysicalLimits = flag; }

      /**
       * Set the flag for using a fixed normalization (true) or floating
       * normalization (false).
       */
      void SetFlagFixNorm(bool flag)
         { fFlagFixNorm = flag; }

      /**
       * Set normalization constant.
       */
      void SetNorm(double norm)
         { fNorm = norm; }

      /**
       * Set the histogram containing the data.
       * @param hist The data histogram.
       * @return An error code.
       */
      int SetData(const TH1D& hist);

      int FixTemplateFunctions(const std::vector<double>& parameters);

      /**
       * Initialize the fitting procedure.
       * @return An error code.
       */
      int Initialize();

      /**
       * Add a template histogram. The templates do not have to be
       * normalized. The histogram has to have the same number of bins
       * and cover the same region as the data histogram.
       * @param hist The template histogram
       * @param name The name of the template
       * @param Nmin The lower limit of the normalization.
       * @param Nmax The upper limit of the normalization.
       */
      int AddTemplate(TH1D hist, const char * name, double Nmin=0, double Nmax=0);

      int AddTemplate(TF1 func, const char * name, double Nmin=0, double Nmax=0, const char* options="");

      /**
       * Describe the efficiency and the uncertainty for all bins.
       * @param name The template name.
       * @param effmean The mean value of the efficiency.
       * @param errsigma The uncertainty on the efficiency.
       */
      int SetTemplateEfficiency(const char * name, double effmean = 1., double effsigma = 0.);

      /**
       * Describe the efficiency and the uncertainty for each bin.
       * @param name The template name.
       * @param eff A histogram describing the efficieny.
       * @param efferr A histogram describing the uncertainty on the efficiency.
       * @return An error code.
       */
      int SetTemplateEfficiency(const char * name, TH1D eff, TH1D efferr);

      /**
       * Add a source of systematic uncertainty.
       * @param errorname The name of the source.
       * @param errortype The shape of the uncertainty in each bin.
       * @return An error code.
       */
      int AddSystError(const char * errorname, const char * errtype = "gauss");

      /**
       * Set the systematic uncertainty on a template.
       * @param errorname The name of the source.
       * @param templatename The template name.
       * @param parerror A histogram describing the uncertainty.
       * @return An error code.
       */
      int SetTemplateSystError(const char * errorname, const char * templatename, TH1D parerror);

      /**
       * Add a correlation among two sources of systematic uncertainty.
       * @param errorname1 The name of the first source.
       * @param errorname2 The name of the second source.
       * @param corr The correlation coefficiency.
       * @return An error code.
       */
      //      int AddSystErrorCorrelation(const char * errorname1, const char * errnorame2, double corr);

      /**
       * Constrains a sum of contributions. Assume a Gaussian prior.
       * @param indices The vector of indicies for the contributions.
       * @param mean The mean value of the prior.
       * @param rms The standard deviation of the prior.
       * @return An error code.
       */
      int ConstrainSum(std::vector<int> indices, double mean, double rms);

      /**
       * Add the calculation of a certain ratio.
       * @param index Index in the numerator.
       * @param indices Vector of indices in the denominator.
       * @param rmin The minimum ratio
       * @param rmax The maximum ratio
       * @return An error code.
       */
      int CalculateRatio(int index, std::vector<int> indices, double rmin = -1.0, double rmax = 1.0);

      /**
       * Checks if two histograms have the same properties.
       * @return 0 if not, 1 if they have the same properties.
       */
      int CompareHistogramProperties(TH1D hist1, TH1D hist2);

      /**
       * Calculates and returns the chi2 value. The chi2 is calculated
       * using the expectation value for the uncertainty.
       */
      double CalculateChi2(const std::vector<double> & parameters);

      /**
       * Calculates and returns the chi2-probability.
       */
      double CalculateChi2Prob(const std::vector<double> & parameters);

      /**
       * Calculates and returns the Likelihood at the global mode.
       */
      double CalculateMaxLike();

      /**
       * Calculates and returns the p-value for the global mode.
       */
      double CalculatePValue();

      /**
       * Calculates and returns the Kolmogorov-Smirnov-probability.
       */
      double CalculateKSProb();

      /**
       * Overloaded from BCIntegrate.
       */
      void MCMCUserIterationInterface();

      /**
       * Prints the stack of templates scaled with the global mode. The
       * data is plotted on top. The following options are available:\n\n
       * "L"  : adds a legend\n\n
       * "E0" : symmetric errorbars on the data points with a length of sqrt(obs).\n\n
       * "E1" : symmetric errorbars on the sum of templates with a length of sqrt(exp).\n\n
       * "E2" : asymmetric errorbars on the sum of templates from error
       * propagation of the parameters.\n\n
       * "E3" : asymmetric errorbars on the data points. The errorbars
       * mark the 16% and 85% quantiles of the Poisson distribution
       * rounded to the lower or upper integer, respectively.
       * @param filename The name of the file the output is printed to.
       * @param options Plotting options
       */
      void PrintStack(const char * filename = "stack.ps", const char * options="LE2E0D");

      /**
       * Print the ratios and the norm.
       * @param filename The filename.
       * @param option Plot options
       * @param pvalue Value which goes with plot options (see BAT manual).
       */
      void PrintRatios(const char * filename = "ratio.ps", int option = 0, double ovalue = 0.);

      /**
       * Print a template to a file.
       * @param name The template name.
       * @param filename The filename.
       * @return An error code.
       */
      int PrintTemplate(const char * name, const char * filename);

      /**
       * Temporary entry. Used for debugging.
       */
      void PrintTemp();

      /**
       * Create a histogram with specified uncertainties.
       * @param hist The histogram to be copied.
       * @param histerr The uncertainties of the new histogram.
       * @return A histogram with uncertainties.
       */
      TH1D CreateErrorHist(TH1D hist, TH1D histerr);

      /**
       * Calculates and returns the log of the Likelihood at a given point
       * in parameter space.
       */
      double LogLikelihood(const std::vector<double> & parameters);

      /**
       * Perform the template fit.
       * @return An error code.
       */
      int PerformFit();

      /**
       * Cobine all sources of systematic uncertainties (including
       * efficiency) by summing the squares.
       * @param name The name of the template
       * @return A histogram with the total uncertainty
       */
      TH1D CombineUncertainties(const char * name);

    protected:

      /**
       * The data histogram.
       */
      TH1D fHistData;

      // histogram container

      /**
       * A container of template histograms. */
      std::vector<TH1D> fTemplateHistogramContainer;

      /**
       * A container of template functions. */
      std::vector<TF1> fTemplateFunctionContainer;

      /**
       * A container of efficiency histograms. */
      std::vector<TH1D> fEffHistogramContainer;

      /**
       * A container of efficiency uncertainty histograms. */
      std::vector<TH1D> fEffErrHistogramContainer;

      /**
       * A matrix of histograms describing systematic uncertainties. */
      std::vector< std::vector<TH1D> > fSystErrorHistogramContainer;

      // name container

      /**
       * A container of template names. */
      std::vector<std::string> fTemplateNameContainer;

      /**
       * A container of names of sources of systematic uncertainties. */
      std::vector<std::string> fSystErrorNameContainer;

      // index container

      /**
       * A container of parameter indeces for templates. */
      std::vector<int> fTemplateParIndexContainer;

      /**
       * A container of histogram indeces for templates. */
      std::vector<int> fHistogramIndexContainer;

      /**
       * A container of function indeces for templates. */
      std::vector<int> fFunctionIndexContainer;

      /**
       * A container of parameter indeces for efficiencies. */
      std::vector<int> fEffParIndexContainer;

      /**
       * A container of parameter indeces for sources of systematic
       * uncertainty. */
      std::vector<int> fSystErrorParIndexContainer;

      // error type container

      /**
       * A container of error types. */
      std::vector<std::string> fSystErrorTypeContainer;


      /**
       * A container for template types. \n
       * 0: histogram
       * 1: function */
      std::vector<int> fTemplateTypeContainer;

      /**
       * A container for constrained sums: indices.
       */
      std::vector< std::vector<int> > fConstraintSumIndices;

      /**
       * A container for constrained sums: mean values.
       */
      std::vector< double > fConstraintSumMean;

      /**
       * A container for constrained sums: mean values.
       */
      std::vector< double > fConstraintSumRMS;

      /**
       * Histogram containing the overall number of expected events.
       */
      TH1D fHistNorm;

      /**
       * Vector of indices for the calculation of ratios.
       */
      std::vector< std::vector<int> > fIndicesRatios1D;

      /**
       * 1-D histograms containing ratios.
       */
      std::vector<TH1D> fHistRatios1D;

      /**
       * Flag for fixing the normalization or not
       */
      bool fFlagFixNorm;

      /**
       * Flag for having physical limits (expectation values greater or
       * equal to 0).
       */
      bool fFlagPhysicalLimits;

      /**
       * A 2-d histogram for calculating the error bars
       */
      TH2D * fUncertaintyHistogramExp;

      /**
       * A 2-d histogram for calculating the error bars
       */
      TH2D * fUncertaintyHistogramObsPosterior;

      /**
       * Normalization constant
       */
      double fNorm;

      /**
       * The number of bins in the data. */
      int fNBins;

      /**
       * The minimum value of the data range. */
      double fXmin;

      /**
       * The maximum value of the data range. */
      double fXmax;

};

// ---------------------------------------------------------

#endif

