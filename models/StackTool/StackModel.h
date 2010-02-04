#ifndef __STACKMODEL__H
#define __STACKMODEL__H

/*!
 * \class StackModel
 * This class can be used for fitting several template
 * histograms to a data histogram. The templates are assumed to have
 * no statistical uncertainty whereas the data are assumed to have
 * Poissonian fluctuations in each bin. Several methods to judge the
 * validity of the model are available.
 * \brief A class for fitting several templates to a data set.
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 10.2009
 */

/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <BAT/BCModel.h>

#include <TH1D.h>
#include <TH2D.h>

// ---------------------------------------------------------

class StackModel : public BCModel
{
 public:

	/**
	 * The default constructor.
	 */
	StackModel();

	/**
	 * A constructor.
	 */
	StackModel(const char * name);

	/**
	 * The default destructor.
	 */
	~StackModel();

	/**
	 * Calculates and returns the log of the prior probability at a
	 * given point in parameter space.
	 */
	double LogAPrioriProbability(std::vector <double> parameters);

	/**
	 * Calculates and returns the log of the Likelihood at a given point
	 * in parameter space.
	 */
	double LogLikelihood(std::vector <double> parameters);

	/**
	 * Returns the number of degrees of freedom of the model. These
	 * are the number of bins - the number of templates.
	 */
	int GetNDF()
		{ return fDataHistogram.GetNbinsX() - GetNParameters(); };

	TH1D* GetHistNorm()
		{ return fHistNorm; }; 

	/**
	 * Set the histogram containing the data.
	 */
	int SetDataHistogram(TH1D hist);

	/**
	 * Return the histogram containing the data. 
	 */ 
	TH1D GetDataHistogram()
	{ return fDataHistogram; }; 

	/**
	 * Set the flag for using a fixed normalization (true) or floating
	 * normalization (false).
	 */ 
	void SetFlagFixNorm(bool flag)
	{ fFlagFixNorm = flag; }; 

	/**
	 * Set normalization constant.
	 */ 
	void SetNorm(double norm)
	{ fNorm = norm; }; 

	/**
	 * Add a template histogram. The templates do not have to be
	 * normalized. The histogram has to have the same number of bins
	 * and cover the same region as the data histogram.
	 * @param hist The template histogram
	 * @param name The name of the template
	 * @param Nmin The lower limit of the normalization.
	 * @param Nmax The upper limit of the normalization.
	 */
	int AddTemplateHistogram(TH1D hist, const char * name, double Nmin=0, double Nmax=0);

	/**
	 * Add a template histogram. The templates do not have to be
	 * normalized. The histogram has to have the same number of bins
	 * and cover the same region as the data histogram.
	 * @param hist The template histogram
	 * @param name The name of the template
	 * @param prior The prior probability
	 */
	int AddTemplateHistogram(TH1D hist, const char * name, TH1D prior);

	/**
	 * Constrains a sum of contributions. Assume a Gaussian prior. 
	 * @param indices The vector of indicies for the contributions. 
	 * @param mean The mean value of the prior. 
	 * @param rms The standard deviation of the prior. 
	 */ 
	int ConstrainSum(std::vector <int> indices, double mean, double rms); 

	/**
	 * Set a Gaussian prior on the expectation value. 
	 * @param index The index of the contribution. 
	 * @param mean The mean value of the prior. 
	 * @param rms The standard deviation of the prior. 
	 * @param adjust If true, the 5-sigma region around the mean value will be 
	 * set as the parameter range for the efficiency. 
	 */ 
	int SetTemplatePrior(int index, double mean, double rms, bool adjust=false); 

	/**
	 * Set a Gaussian prior on the efficiencies. If the err is negative,
	 * no uncertainty is assumed. 
	 * @param index The index of the contribution.
	 * @param eff The mean value of the efficiency prior. 
	 * @param err The standard deviation of the efficiency prior.
	 * @param adjust If true, the 5-sigma region around the mean value will be 
	 * set as the parameter range for the efficiency. 
	 */ 
	int SetTemplateEfficiency(int index, double eff, double err, bool adjust=false); 

	/**
	 * Add the calculation of a certain ratio. 
	 * @param index Index in the numerator.
	 * @param indices Vector of indices in the denominator.
	 */ 
	int CalculateRatio(int index, std::vector<int> indices); 

	/**
	 * Set a flag for having physical limits (expectation values greater
	 * or equal to 0).
	 * @param flag The flag.
	 */ 
	void SetFlagPhysicalLimits(bool flag)
	{ fFlagPhysicalLimits = flag; }; 

	/**
	 * Checks if two histograms have the same properties.
	 * @return 0 if not, 1 if they have the same properties.
	 */
	int CompareHistogramProperties(TH1D hist1, TH1D hist2);

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
	 */ 
	void PrintRatios(const char * filename = "ratio.ps");
	
	/**
	 * Calculates and returns the chi2 value. The chi2 is calculated
	 * using the expectation value for the uncertainty.
	 */
	double CalculateChi2();

	/**
	 * Calculates and returns the chi2-probability.
	 */
	double CalculateChi2Prob();

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
	 * Temporary entry. Used for debugging.
	 */
	void PrintTemp();

 protected:

	/**
	 * The data histogram.
	 */
	TH1D fDataHistogram;

	/**
	 * A container for the template histograms.
	 */
	std::vector <TH1D> fTemplateHistogramContainer;

	/**
	 * A container for the template names.
	 */
	std::vector <std::string> fTemplateNameContainer;

	/**
	 * A constainer of the efficiency. 
	 */ 
	std::vector<double> fTemplateEff; 

	/**
	 * A constainer of the efficiency uncertainty.
	 */ 
	std::vector<double> fTemplateEffErr; 

	/**
	 * A constainer of the efficiency. 
	 */ 
	std::vector<double> fTemplatePriorMean; 

	/**
	 * A constainer of the efficiency uncertainty.
	 */ 
	std::vector<double> fTemplatePriorSigma; 

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
	TH1D * fHistNorm; 

	/**
	 * 1-D histograms containing the prior for each template.
	 */ 
	std::vector <TH1D> fHistPrior; 

	/**
	 * Vector of indices for the calculation of ratios.
	 */ 
	std::vector< std::vector<int> > fIndicesRatios1D; 

	/**
	 * 1-D histograms containing ratios. 
	 */ 
	std::vector <TH1D> fHistRatios1D; 

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
	 * Normalization constant
	 */ 
	double fNorm; 

	/**
	 * A 2-d histogram for calculating the error bars
	 */
	TH2D * fUncertaintyHistogramExp;

	/**
	 * A 2-d histogram for calculating the error bars
	 */
	TH2D * fUncertaintyHistogramObsPosterior;

};

// ---------------------------------------------------------

#endif

