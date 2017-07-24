Predefined Models {#cha-predefined-models}
=================

[TOC]

@todo Only describe how to use BAT, not to extend it

@section sec-native-data Using BAT's Native Data Structures

@section sec-eff-fitter Efficiency Fitter

@section sec-graph-fitter Graph Fitter
@section sec-hist-fitter Histogram Fitter

@section sec-mtf Multi-template Fitter

@todo let Kevin check if info up to date, this is just copied from introduction.tex Since this section is long, it deserves its own chapter

The Multi-Template Fitter (MTF) is a tool which allows to fit several
template histograms to a data histogram. The content of the bins in the
templates are assumed to fluctuate independently according to Poisson
distributions. Several channels can be fitted simultaneously.

@subsection sec-mtf-math  Mathematical Formulation

The multi-template fitter is formulated in terms of Bayesian reasoning.
The posterior probability is proportional to the product of the
Likelihood and the prior probability. The latter can be freely chosen by
the user whereas the Likelihood is predefined. It is a binned Likelihood
which assumes that the fluctuations in each bin are of Poisson nature
and independent of each other. All channels, processes and sources of
systematic uncertainties are assumed to be uncorrelated.

The parameters of the model are thus the expectation values of the
different processes, \f$\lambda_{k}\f$, and the nuisance parameters,
\f$\delta_{l}\f$.

@subsubsection sec-mtf-exclsyst Excluding Systematic Uncertainies

In case no sources of systematic uncertainty are taken into account the
Likelihood is defined as
\f{eqnarray}{
L = \prod_{i=1}^{N_{\mathrm{ch}}} \prod_{j=1}^{N_{\mathrm{bin}}} \frac{\lambda_{ij}^{n_{ij}}}{n_{ij}!} e^{-\lambda_{ij}} \, ,
\label{eqn:likelihood}
\f}
 where \f$N_{\mathrm{ch}}\f$ and
\f$N_{\mathrm{bin}}\f$ are the number of channels and bins, respectively.
\f$n_{ij}\f$ and \f$\lambda_{ij}\f$ are the observed and expected number of
events in the \f$j\f$th bin of the \f$i\f$th channel. The expected number of
events are calculated via
\f{eqnarray}{
\lambda_{ij} & = & \sum_{k=1}^{N_{\mathrm{p}}} \lambda_{ijk} \\
            & = & \sum_{k=1}^{N_{\mathrm{p}}} \lambda_{k} \cdot f_{ijk} \cdot \epsilon_{ik} \, ,
\label{eqn:expectation}
\f}
 where \f$f_{ij}\f$ is the bin content
of the \f$j\f$th bin in the normalized template of the \f$k\f$th process in the
\f$i\f$th channel. \f$\epsilon_{ik}\f$ is the efficiency of the \f$k\f$th process in
the \f$i\f$th channel specified when setting the template. \f$\lambda_{k}\f$ is
the contribution of the \f$k\f$th process and is a free parameter of the
fit.

@subsubsection sec-mtf-inclsyst Including Systematic Uncertainies

In case sources of systematic uncertainties are taken into account, the
efficiency \f$\epsilon_{ik}\f$ is modified according to a nuisance
parameter:
\f{eqnarray}{
\epsilon_{ik} \rightarrow \epsilon_{ik} \cdot (1 + \sum_{l=1}^{N_{\mathrm{syst}}} \delta_{l}\cdot \Delta\epsilon_{ijkl}) \, ,
\label{eqn:systematic}
\f}
 where \f$\delta_{l}\f$ is the nuisance
parameter associated with the source of systematic uncertainty and
\f$\Delta\epsilon_{ijkl}\f$ is the change in efficiency due to the \f$l\f$th
source of systematic uncertainty in the \f$i\f$th channel and \f$j\f$th bin for
the \f$k\f$th process.

@subsection sec-mtf-basics  Creating the Fitter

The main MTF class `BCMTF` is derived from the `BCModel` class. A new
instance can be created via

    @code{.cpp}
    BCMTF::BCMTF()
    BCMTF::BCMTF(const char * name)
    @endcode

where the name of the MTF can be specified via argument `name`.

@subsection sec-mtf-addchannel Adding a Channel

The MTF fits several channels simultaneously. These channels can be
physics channels, e.g., \f$Z^{0}\rightarrow e^{+}e^{-}\f$ and
\f$Z^{0}\rightarrow \mu^{+}\mu^{-}\f$, samples with disjunct jet
multiplicity or entirely different classes altogether. A new channel can
be added using the following method:

    @code{.cpp}
    int BCMTF::AddChannel(const char * name)
    @endcode

where `name` is the name of the process, and the return value is an
error code. Note that at least one channel has to be added.

@subsection sec-mtf-adddata Adding a Data Set

Each channel added to the MTF has a unique data set which comes in form
of a (`TH1D`) histogram. It can be defined using the following method:

    @code{.cpp}
    int BCMTF::SetData(const char * channelname, TH1D hist)
    @endcode

where `channelname` is the name of the channel and `hist` is the
histogram representing the data. The return value is an error code.

@subsection sec-mtf-addprocess Adding a Process

Each template that is fit to the data set corresponds to a process,
where one process can occur in several channels. The fit then defines
the contribution of the process and thus each process comes with one
model parameter. A process can be added using the following method:

    @code{.cpp}
    int BCMTF::AddProcess(const char * name,
                          double nmin = 0.,
                          double nmax = 1.) ,
    @endcode

where `name` is the name of the process and `nmin` and `nmax` are the
lower and upper bound of the parameter associated with the contribution
of the process. The parameter is denoted
\f$\lambda_{k} \, (k=1\dots N_{\mathrm{p}})\f$ in the @ref sec-mtf-math.
Note that at least one process has to be added. A prior needs to be
defined for each process, using the default `BCModel` methods.

It is likely that a single process will have different shapes in
different channels. Thus, templates for a process need to be defined for
each channel separately using the following method:

    @code{.cpp}
    int BCMTF::SetTemplate(const char * channelname,
                           const char * processname,
                           TH1D hist,
                           double efficiency = 1.) ,
    @endcode

where `channelname` and `processname` are the names of the channel and
the process, respectively. The parameter `hist` is the (`TH1D`)
histogram (or template) which represents the process. The histogram will
be normalized to unity and the entries in the normalized histogram are
the probabilities to find an event of a process \f$k\f$ and channel \f$i\f$ in
bin \f$j\f$. This probability is denoted
\f$f_{ijk}\,(i=1\dots N_{\mathrm{ch}}, \, j=1\dots N_{\mathrm{b}}, \,
k=1\dots N_{\mathrm{p}})\f$ in the @ref sec-mtf-math. The last
parameter, `efficiency`, is the efficiency of the process in that
channel and is used to scale to template during the fit. This is needed
if a process contributes with different amounts in two separate
channels. The efficiency is denoted \f$\epsilon_{ik} \, (i=1\dots
N_{\mathrm{ch}}, \, k=1\dots N_{\mathrm{p}})\f$ in
the @ref sec-mtf-math. The return value is an error code. Note that
templates do not have to be set if the process does not contribute to a
particular channel.

@subsection sec-mtf-addunc Adding Systematic Uncertainties

Systematic uncertainties can alter the shape of a template. Sources of
systematic uncertainty can be included in the fit using nuisance
parameters. This nuisance parameter is assumed to alter the original
template linearly, where values of -1, 0, and 1 correspond to the
“downwards” shifted, nominal and “upwards shifted” template,
respectively. The nuisance parameters are denoted \f$\delta_{l} \,
(l=1\dots N_{\mathrm{syst}})\f$ in the @ref sec-mtf-math. Shifted
refers to a change of one standard deviation. An example for a nuisance
parameter could be the jet energy scale (JES). With a nominal JES of 1
and and uncertainty of 5%, the scaled templates correspond to a JES of
0.95 and 1.05, respectively. A prior needs to be defined for each
nuisance parameter which is usually chosen to be a standard normal
distribution. A source of systematic uncertainty can be added using the
following method:

    @code{.cpp}
    int BCMTF::AddSystematic(const char * name,
                             double min = -5.,
                             double max = 5.) ,
    @endcode

where `name` is the name of the source of systematic uncertainty and
`min` and `max` are the lower and upper bound of the nuisance parameter,
respectively. The return value is an error code.

Since the different sources of systematic uncertainty have an individual
impact on each process and in each channel, these need to be specified.
Two method can be used to define the impact:

    @code{.cpp}
    int BCMTF::SetSystematicVariation(const char * channelname,
                                      const char * processname,
                                      const char * systematicname,
                                      TH1D hist_up,
                                      TH1D hist_down) ,
    @endcode

where `channelname`, `processname` and `systematicname` are the names of
the channel, the process and the source of systematic uncertainty. The
(`TH1D`) histograms `hist_up` and `hist_down` are the histograms
corresponding to an “up”- and “down”-scaling of the systematic
uncertainty of one standard deviation, i.e., for each bin entry \f$y\f$ they
are calculated as
\f{eqnarray}{
\label{eqn:deltay}
\Delta_{\mathrm{up}}   & = & (y_{\mathrm{up}} - y_{\mathrm{nominal}})/y_{\mathrm{nominal}} \, , \\
\Delta_{\mathrm{down}} & = & (y_{\mathrm{nominal}} - y_{\mathrm{down}})/y_{\mathrm{nominal}} \, .
\f}

Note the sign of the down-ward fluctuation. These histograms define the
change of the bins in each template in the efficiency which is denoted
\f$\Delta\epsilon_{ijkl} \, (i=1\dots N_{\mathrm{ch}}, \,
j=1\dots N_{\mathrm{b}}, \, k=1\dots N_{\mathrm{p}}, \, l=1\dots
N_{\mathrm{syst}})\f$. For example, if the value for a particular bin of
`hist_up` is 0.05, i.e., if the systematic uncertaintiy is 5% in that
bin, then the efficiency of the process in that channel will be
multiplied by \f$(1+0.05)\f$. The return value is an error code.

The second variant does not take the difference in efficiency, but
calculates it internally from the absolute values:

    @code{.cpp}
    int BCMTF::SetSystematicVariation(const char * channelname,
                                      const char * processname,
                                      const char * systematicname,
                                      TH1D hist,
                                      TH1D hist_up,
                                      TH1D hist_down) ,
    @endcode

where `channelname`, `processname` and `systematicname` are the names of
the channel, the process and the source of systematic uncertainty. The
(`TH1D`) histograms `hist`, `hist_up` and `hist_down` are the nominal
histogram and the histograms corresponding to an “up”- and
“down”-scaling of the systematic uncertainty of one standard deviation.
In this case, the histograms are not the relative differences but the
absolute values. The return value is an error code.

@subsection sec-mtf-runfit Running the Fit

The fit can be started using one of the standard `BCModel` fitting
methods, e.g.

    @code{.cpp}
    BCMTF::MarginalizeAll() ,
    BCMTF::FindMode() .
    @endcode

@subsection sec-mtf-mtfoutput  Output
The MTF produces several outputs:

1. `PrintAllMarginalized(const char* name)` prints the marginalized
distributions in 1D and 2D for all parameters, i.e., the processes
and nuisance parameters into a PostScript file `name`.

2. `PrintResults(const char* name)` writes a summary of the fit into a
text file `name`.

3. To print a stacked histogram of the templates and the data histogram in
the file `name` using a set of parameters `parameters`, do
@code{.cpp}
PrintStack(int channelindex,
           const std::vector<double> & parameters,
           const char * filename = "stack.pdf",
           const char * options = "")
PrintStack(const char * channelname,
           const std::vector<double> & parameters,
           const char * filename = "stack.pdf",
           const char * options = "")
@endcode

For example, these could be the best fit results. Several options can
be specified:

-   `logx`: uses a log-scale for the x-axis.

-   `logy`: uses a log-scale for the y-axis.

-   `logx`: plot the x-axis on a log scale

-   `logy`: plot the y-axis on a log scale

-   `bw`: plot in black and white

-   `sum`: draw a line corresponding to the sum of all templates

-   `stack`: draw the templates as a stack

-   `e0`: do not draw error bars

-   `e1`: draw error bars corresponding to sqrt(n)

-   `b0`: draw an error band on the expectation corresponding to the
    central 68% probability

-   `b1`: draw bands showing the probability to observe a certain
    number of events given the expectation. The green (yellow, red)
    bands correspond to the central 68% (95%, 99.8%) probability

@subsection sec-mtf-mtfsettings  Settings
Several settings can be changed which impact the fit.

-   `SetFlagEfficiencyConstraint` sets a flag if the overall efficiency
    (calculated from the value given when setting a template and the
    corresponding systematic uncertainties) is constrained to be
    between 0 and 1 or not. The default value is `true`.

@subsection sec-mtf-facility  Analysis Facility
The analysis facility allows to perform a variety of analyses and
ensemble tests for a given MTF. It can be created using the constructor:

    @code{.cpp}
    BCMTFAnalysisFacility::BCMTFAnalysisFacility(BCMTF * mtf)
    @endcode

where `mtf` is the corresponding MTF object.

@subsection sec-mtf-perf_ensemble Performing Ensemble Tests

Ensemble testing is done in two steps: first, ensembles are generated
according to the processes defined in the MTF. The ensembles are stored
in root files. In a second step, the ensembles are analyzed using the
MTF specified.

@subsection sec-mtf-create_ensem Creating Ensembles

Ensembles can be generated using several methods. A single ensemble can
be generated using the following method:

    @code{.cpp}
    std::vector<TH1D> BCMTFAnalysisFacility::BuildEnsemble(const std::vector<double> & parameters)
    @endcode

where `parameters` is a set of parameters which corresponds to those in
the template fitter, i.e., the process contributions and nuisance
parameters. For most applications, the best fit parameters of the data
set at hand is used. The return value is a set of histograms
corresponding to a pseudo data set for the different channels.

A similar method is used to generate multiple ensembles:

    @code{.cpp}
    std::vector<TH1D> BCMTFAnalysisFacility::BuildEnsembles(const std::vector<double> & parameters, int nensembles)
    @endcode

where `nensembles` is the number of ensembles to be generated. The
return value is a pointer to a `TTree` object in which the ensembles are
stored. The entries in the tree are the parameters and the number of
entries in each bin of the data histograms.

The third method is based on a tree where the tree contains a set of
parameters for each ensemble. This option is preferred if, e.g., the
ensembles should be varied accoring to the prior probabilities. The
method used to generate ensembles is

    @code{.cpp}
    std::vector<TH1D> BCMTFAnalysisFacility::BuildEnsembles(TTree * tree, int nensembles)
    @endcode

where `tree` is the input tree. Note that the ensembles are randomized,
i.e., the first event in the tree does not correspond to the first
ensemble. This is done to avoid biases if the tree itself is the output
of a Markov Chain.

@subsection sec-mtf-analyze_ensem Analyzing Ensembles

Ensemble tests can be performed usign the ensembles defined earlier or
using a set of parameters. In the former case, the method is:

    @code{.cpp}
    TTree * BCMTFAnalysisFacility::PerformEnsembleTest(TTree * tree, int nensembles)
    @endcode

where `tree` is the tree of ensembles and `nensembles` is the number of
ensembles to be analyzed. The return value is a tree containing the
information about the analyzed ensemble. The list of variables is

-   `parameter_i`: the \f$i\f$th parameter value used at the generation of
    the ensemble.

-   `mode_global_i`: the \f$i\f$th global mode.

-   `std_global_i`: the \f$i\f$th standard deviation evaluated with the
    global mode.

-   `chi2_generated_i`: the \f$\chi^{2}\f$ calculated using the parameters
    at generation of the ensemble for channel \f$i\f$.

-   `chi2_mode_i`: the \f$\chi^{2}\f$ calculated using the global mode
    parameters for channel \f$i\f$.

-   `cash_generated_i`: the Cash statistic (Likelihood ratio) calculated
    using the parameters at generation of the ensemble for channel \f$i\f$.

-   `cash_mode_i`: the Cash statistic (Likelihood ratio) calculated
    using the global mode parameters for channel \f$i\f$.

-   `n_events_i`: the number of events in the ensemble in channel \f$i\f$.

-   `chi2_generated_total`: the total \f$\chi^{2}\f$ calculated using the
    parameters at generation of the ensemble.

-   `chi2_mode_total`: the total \f$\chi^{2}\f$ calculated using the global
    mode parameters.

-   `cash_generated_total`: the total Cash statistic calculated using
    the parameters at generation of the ensemble.

-   `cash_mode_total`: the total Cash statistic calculated using the
    global mode parameters.

-   `n_events_total`: the total number of events in the ensemble.

Ensemble tests can also be performed using the following methd:

    @code{.cpp}
    TTree * BCMTFAnalysisFacility::PerformEnsembleTest(const std::vector<double> & parameters, int nensembles)
    @endcode

in which case the ensembles are generated internally using the
parameters and are then analyzed.

By default the log messages for both the screen and the log-file are
suppressed while performing the ensemble test. This can be changed using

    @code{.cpp}
    void BCMTFAnalysisFacility::SetLogLevel(BCLog::LogLevel level)
    @endcode

@subsection sec-mtf-autom_analyses Performing Automated Analyses

The analysis facility also allows to perform an automated analysis over
individual channels and or systematic uncertainties.

@subsubsection sec-mtf-perf_single_channel Performing Single-Channel Analyses

The current data set can be analyzed automatically for each channel
separately using the analysis facility method

    @code{.cpp}
    int BCMTFAnalysisFacility::PerformSingleChannelAnalyses(const char * dirname, const char * options = "")
    @endcode

where `dirname` is the name of a directory which will be created and
into which all plots will be copied. If `mcmc` is specified in the
`options` then the MCMC will be run for each channel. The method creates
all marginalized distributions and results as well as an overview plot.
If the option `nosyst` is specified, the systematic uncertainties are
all switched off.

@subsubsection sec-mtf-perf_single_channel_syst Performing Single Systematic Analyses

Similarly, the method

    @code{.cpp}
    int BCMTFAnalysisFacility::PerformSingleSystematicAnalyses(const char * dirname, const char * options = "")
    @endcode

can be used to perform a set of analyses for each systematic uncertainty
separately.

@subsubsection sec-mtf-perf_calib Performing Calibration Analyses

Ensemble tests for different sets of parameters can be automized by
using the method

    @code{.cpp}
    int BCMTFAnalysisFacility::PerformCalibrationAnalysis(const char * dirname,
                                   const std::vector<double> & default_parameters,
                                   int index,
                                   const std::vector<double> & parametervalues,
                                   int nensembles = 1000)
    @endcode

which can be used to easily generate calibration curves. The ensembles
are generated for a set of parameters, `default_parameters` where one of
the parameters, `index`, can vary. The paramter values are defined by
`parametervalues`. `nensembles` defines the number of pseudo data sets
used for each ensemble.
