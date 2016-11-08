#ifndef __BCMTF__H
#define __BCMTF__H

/**
 * @class BCMTF
 * @brief A class for fitting several templates to a data set.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.1
 * @date 06.2012
 * @details This class can be used for fitting several template
 * histograms to a data histogram. The templates are assumed to have
 * no statistical uncertainty whereas the data are assumed to have
 * Poissonian fluctuations in each bin. Several methods to judge the
 * validity of the model are available.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "../../BAT/BCModel.h"

#include <TH1D.h>

class BCMTFChannel;
class BCMTFProcess;
class BCMTFSystematic;
class TF1;

// ---------------------------------------------------------
class BCMTF : public BCModel
{

public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * A constructor.
     * @param name The name of the model */
    BCMTF(const std::string& name = "multi_template_fitter");

    /**
     * The default destructor. */
    ~BCMTF();

    /** @} */
    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return The number of channels. */
    int GetNChannels() const
    { return fNChannels; };

    /**
     * @return The number of processes. */
    int GetNProcesses() const
    { return fNProcesses; };

    /**
     * @return The number of systematics. */
    int GetNSystematics() const
    { return fNSystematics; };

    /**
     * @param name The name of the channel.
     * @return The channel index. */
    int GetChannelIndex(const std::string& name) const;

    /**
     * @param name The name of the process.
     * @return The process index. */
    int GetProcessIndex(const std::string& name) const;

    /**
     * @param name The name of the systematic.
     * @return The systematic uncertainty index. */
    int GetSystematicIndex(const std::string& name) const;

    /**
     * @param index The parameter index (mtf counting) .
     * @return The parameter number corresponding to the parameter index (BAT counting). */
    int GetParIndexProcess(int index) const
    { return fProcessParIndexContainer.at(index); };

    /**
     * @param name The systematic uncertainty index (mtf counting).
     * @return The parameter number corresponding to the systematic uncertainty (BAT counting). */
    int GetParIndexSystematic(int index) const
    { return fSystematicParIndexContainer.at(index); };

    /**
     * @param name The channel index.
     * @return The channel object. */
    BCMTFChannel* GetChannel(int index)
    { return fChannelContainer.at(index); };

    /**
     * @param name The process index.
     * @return The process object. */
    BCMTFProcess* GetProcess(int index)
    { return fProcessContainer.at(index); };

    /**
     * @param name The systematic ucnertainty index.
     * @return The systematic uncertainty object. */
    BCMTFSystematic* GetSystematic(int index)
    { return fSystematicContainer.at(index); };

    /** @} */

    /** \name Member functions (set) */
    /** @{ */

    /**
     * Set the data histogram in a particular channel.
     * @param channelname The name of the channel.
     * @param hist The TH1D histogram.
     * @param minimum The minimum number of expected events (used for calculation of uncertainty bands).
     * @param maximum The maximum number of expected events (used for calculation of uncertainty bands). */
    void SetData(const std::string& channelname, TH1D hist, double minimum = -1, double maximum = -1);

    /**
     * Set the template for a specific process in a particular channel.
     * @param channelname The name of the channel.
     * @param processname The name of the process.
     * @param hist The TH1D histogram.
     * @param efficiency The efficiency of this process in this channel. */
    void SetTemplate(const std::string& channelname, const std::string& processname, TH1D hist, double efficiency = 1., double norm = 1.);

    /**
     * Set the template for a specific process in a particular
     * channel. This is an alternative way to describe the number of
     * expected events. It is used in the rare case that processes
     * cannot be summed directly, but interference effects have to
     * be taken into account. The expected number of events is then
     * parametrized as a function for each bin.
     * @param channelname The name of the channel.
     * @param processname The name of the process.
     * @param funccont A vector of pointers of TF1 functions.
     * @param nbins The number of bins used for the histogram.
     * @param efficiency The efficiency of this process in this channel.
     *
     * @see SetTemplate(const std::string& channelname, const std::string& processname, TH1D hist, double efficiency = 1.) */
    void SetTemplate(const std::string& channelname, const std::string& processname, std::vector<TF1*>* funccont, int nbins, double efficiency = 1.);

    /**
     * Set an expectation function.
     * @param parindex The index of the parameter
     * @param func The pointer to a TF1 function.
     * @see SetTemplate(const std::string& channelname, const std::string& processname, std::vector<TF1 *> * funccont, int nbins, double efficiency = 1.) */
    void SetExpectationFunction(int parindex, TF1* func)
    { fExpectationFunctionContainer[parindex] = func; };

    /**
     * Set the impact of a source of systematic uncertainty for a
     * particular source of systematic uncertainty and process in a
     * given channel. The impact is the relative deviation between
     * the varied and the nominal template, i.e., if the variation
     * is set to 0.05 then the efficiency is multiplied by
     * (1+0.05*eps), where eps is the corresponding nuisance
     * parameter.
     * @param channelname The name of the channel.
     * @param processname The name of the process.
     * @param systematicname The name of the source of systematic uncertainty.
     * @param variation_up The relative shift between the up-variation and the nominal template: (up-nom)/nom.
     * @param variation_down The relative shift between the down-variation and the nominal template: (nom-down)/nom. */
    void SetSystematicVariation(const std::string& channelname, const std::string& processname,  const std::string& systematicname, double variation_up, double variation_down);

    /**
     * Set the impact of a source of systematic uncertainty for a
     * particular source of systematic uncertainty and process in a
     * given channel. The variation depends on x and is given as a
     * histogram.
     * @param channelname The name of the channel.
     * @param processname The name of the process.
     * @param systematicname The name of the source of systematic uncertainty.
     * @param hist_up The TH1D histogram defining the relative shift between the up-variation and the nominal template: (up-nom)/nom.
     * @param hist_down The TH1D histogram defining the relative shift between the down-variation and the nominal template: (nom-down)/nom.
     *
     * @see SetSystematicVariation(const std::string& channelname, const std::string& processname,  const std::string& systematicname, double variation_up, double variation_down) */
    void SetSystematicVariation(const std::string& channelname, const std::string& processname,  const std::string& systematicname, TH1D hist_up, TH1D hist_down);

    /**
     * Set the impact of a source of systematic uncertainty for a
     * particular source of systematic uncertainty and process in a
     * given channel. The variation depends on x. The histograms are
     * the raw histograms after the shift and the variation wrt the
     * nominal template will be calculated automatically.
     * @param channelname The name of the channel.
     * @param processname The name of the process.
     * @param systematicname The name of the source of systematic uncertainty.
     * @param hist The histogram with the nominal template
     * @param hist_up The TH1D histogram after up-scaling of the systematic uncertainty.
     * @param hist_down The TH1D histogram after down-scaling of the systematic uncertainty.
     *
     * @see SetSystematicVariation(const std::string& channelname, const std::string& processname,  const std::string& systematicname, double variation_up, double variation_down) */
    void SetSystematicVariation(const std::string& channelname, const std::string& processname,  const std::string& systematicname, TH1D hist, TH1D hist_up, TH1D hist_down);

    /**
     * Set a flag for the efficiency: if true then the total
     * efficiency including all systematic uncertainties has to be
     * between 0 and 1 at all times. Larger (smaller) values are set
     * to 1 (0) during the calculation.
     * @param flag The flag */
    void SetFlagEfficiencyConstraint(bool flag)
    { fFlagEfficiencyConstraint = flag; };

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Add a channel
     * @param name The channel name. */
    void AddChannel(const std::string& name);

    /**
     * Add a process and the associated BAT parameter.
     * @param name The process name
     * @param nmin The minimum number of expected events (lower limit on the BAT parameter values).
     * @param nmax The maximum number of expected events (upper limit on the BAT parameter values).
     * @param color The histogram color
     * @param fillstyle The histogram fill style
     * @param linestyle The histogram line style */
    void AddProcess(const std::string& name, double nmin = 0., double nmax = 1., int color = -1, int fillstyle = -1, int linestyle = -1);

    /**
     * Add a source of systematic uncertainty and the associated BAT (nuisance) parameter.
     * @param name The systematic uncertainty name.
     * @param min The lower limit on the BAT parameter values, typically -5 sigma if Gaussian constraint is used.
     * @param max The upper limit on the BAT parameter values, typically +5 sigma if Gaussian constraint is used. */
    void AddSystematic(const std::string& name, double min = -5., double max = 5.);

    /**
     * Return the expected number of events for a channel and bin.
     * @param channelindex The channel index.
     * @param binindex The bin index.
     * @param parameters A reference to the parameters used to calculate the expectation.
     * @return The expectation value */
    double Expectation(int channelindex, int binindex, const std::vector<double>& parameters);

    /**
     * Return the function value of the expectation function for a
     * parameter, channel and process.
     * @param parindex The parameter index.
     * @param channelindex The channel index.
     * @param processindex The process index.
     * @param A reference to the parameters used to calculate the expectation.
     * @return The expectation function value. */
    double ExpectationFunction(int parindex, int channelindex, int processindex, const std::vector<double>& parameters);

    /**
     * Return the efficiency for a process in a channel and for a particular bin.
     * @param channelindex The channel index.
     * @param processindex The process index.
     * @param binindex The bin index.
     * @param parameters A reference to the parameters used to calculate the efficiency.
     * @return The efficiency. */
    double Efficiency(int channelindex, int processindex, int binindex, const std::vector<double>& parameters);

    // avoid accidental hiding
    using BCModel::Probability;

    /**
     * Return the probability for a process in a channel and for a
     * particular bin. This corresponds to the (normalized) bin
     * content of the template.
     * @param channelindex The channel index.
     * @param processindex The process index.
     * @param binindex The bin index.
     * @param parameters A reference to the parameters used to calculate the probability.
     * @return The probability. */
    double Probability(int channelindex, int processindex, int binindex, const std::vector<double>& parameters);

    /**
     * Calculate a chi2 for a single channel given a set of parameters.
     * @param channelindex The channel index.
     * @param parameters A reference to the parameters used to calculate the chi2.
     * @return A chi2 value.
     * @see CalculateChi2(const std::vector<double> & parameters) */
    double CalculateChi2(int channelindex, const std::vector<double>& parameters);

    /**
     * Calculate a chi2 for all channels together given a set of parameters.
     * @param parameters A reference to the parameters used to calculate the chi2.
     * @return A chi2 value.
     * @see CalculateChi2(int channelindex, const std::vector<double> & parameters) */
    double CalculateChi2(const std::vector<double>& parameters);

    /**
     * Calculate the Cash statistic for a single channel
     * @param channelindex The channel index.
     * @param parameters A reference to the parameters used to calculate the Cash statistic.
     * @return The Cash statistic.
     * @see CalculateCash(const std::vector<double> & parameters) */
    double CalculateCash(int channelindex, const std::vector<double>& parameters);

    /**
     * Calculate the Cash statistic for all channels
     * @param parameters A reference to the parameters used to calculate the Cash statistic.
     * @return The Cash statistic.
     * @see CalculateCash(int channelindex, const std::vector<double> & parameters) */
    double CalculateCash(const std::vector<double>& parameters);

    /**
     * Calculates and returns the fast p-value for the total likelihood as test statistic.
     * @see BCMath::CorrectPValue for correcting the fitting bias
     * @param channelindex The channel index.
     * @param parameters A reference to the parameters for which the model expectations are computed.
     * @return the uncorrected p-value
     */
    double CalculatePValue(int channelindex, const std::vector<double>& parameters);

    /**
     * Calculates and returns the fast p-value for the total likelihood as test statistic.
     *
     * @note Obtain the results with @see BCModel::GetPValue and the value corrected for degrees
     * of freedom with @see BCModel::GetPValueNDoF
     *
     * @param parameters A reference to the parameters for which the model expectations are computed.
     */
    double CalculatePValue(const std::vector<double>& parameters);

    /** @} */

    /** \name Member functions (output methods) */
    /** @{ */

    /**
     * Print a summary of the fit to the log.*/
    void PrintFitSummary();

    /**
     * Print the stack of templates together with the data in a
     * particular channel. Several plot options are available:\n
     * "logx" : plot the x-axis on a log scale \n
     * "logy" : plot the y-axis on a log scale \n
     * "bw"   : plot in black and white \n
     * "sum"  : draw a line corresponding to the sum of all templates \n
     * "stack" : draw the templates as a stack \n
     * "e0" : do not draw error bars \n
     * "e1" : draw error bars corresponding to sqrt(n) \n
     * "b0" : draw an error band on the expectation corresponding to the central 68% probability \n
     * "b1" : draw bands showing the probability to observe a certain number of events given the expectation. The green (yellow, red) bands correspond to the central 68% (95%, 99.8%) probability \n
     * @param channelindex The channel index.
     * @param parameters A reference to the parameters used to scale the templates.
     * @param options The plotting options. */
    void PrintStack(int channelindex, const std::vector<double>& parameters, const std::string& filename = "stack.pdf", const std::string& options = "e1b0stack");

    /**
     * Print the stack of templates together with the data in a
     * particular channel.
     * @param channelname The name of the channel.
     * @param parameters A reference to the parameters used to scale the templates.
     * @param options The plotting options.
     *
     * @see PrintStack(int channelindex, const std::vector<double> & parameters, const std::string& filename = "stack.pdf", const std::string& options = "e1b0stack") */
    void PrintStack(const std::string& channelname, const std::vector<double>& parameters, const std::string& filename = "stack.pdf", const std::string& options = "e1b0stack")
    {
        return PrintStack(GetChannelIndex(channelname), parameters, filename, options);
    }

    /** @} */

    /** \name Member functions (overloaded from BCModel) */
    /** @{ */

    /**
     * Calculates natural logarithm of the likelihood.
     * Method needs to be overloaded by the user.
     * @param params A set of parameter values
     * @return Natural logarithm of the likelihood */
    double LogLikelihood(const std::vector<double>& parameters);

    /**
     * Method executed for every iteration of the MCMC. User's code should be
     * provided via overloading in the derived class*/
    void MCMCUserIterationInterface();

    /** @} */

private:

    /**
     * A container of channels. */
    std::vector<BCMTFChannel*> fChannelContainer;

    /**
     * A container of processes. */
    std::vector<BCMTFProcess*> fProcessContainer;

    /**
     * A container of sources of systematic uncertainty. */
    std::vector<BCMTFSystematic*> fSystematicContainer;

    /**
     * The number of channels. */
    int fNChannels;

    /**
     * The number of processes. */
    int fNProcesses;

    /**
     * The number of systematics uncertainties. */
    int fNSystematics;

    /**
     * A container of parameter indeces for the process
     * normalization. */
    std::vector<int> fProcessParIndexContainer;

    /**
     * A container of parameter indeces for the systematics. */
    std::vector<int> fSystematicParIndexContainer;

    /**
     * A flag for the efficiency constraint: efficiency within 0 and
     * 1 (true) or not (false). */
    bool fFlagEfficiencyConstraint;

    /**
     * A container of functions for the expectation. */
    std::vector<TF1*> fExpectationFunctionContainer;

    /**
     * P value. */
    double fPValue;

    /**
     * P value accounting for degrees of freedom. */
    double fPValueNDoF;

};
// ---------------------------------------------------------

#endif
