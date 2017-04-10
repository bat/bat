Predefined Models {#cha-predefined-models}
=================

@todo Only describe how to use BAT, not to extend it

@section sec-native-data Using BAT's Native Data Structures

@section sec-eff-fitter Efficiency Fitter

@section sec-graph-fitter Graph Fitter
@section sec-hist-fitter Histogram Fitter

<!-- ------------------- -->
@section sec-mtf Multitemplate Fitter
<!-- ------------------- -->

@todo let Kevin check if info up to date, this is just copied from introduction.tex

The Multi Template Fitter (MTF) is a tool that allows to fit several
template histograms to a data histogram. The content of the bins in
the templates are assumed to fluctuate independently according to
Poisson distributions. Several channels can be fitted simultaneously.

<!-- % -------------------------------------------------------- -->
<!-- % Mathematical formulation -->
<!-- % -------------------------------------------------------- -->
@subsection sec-mtf-math Mathematical formulation

The multi-template fitter is formulated in terms of Bayesian reasoning. The
posterior probability is proportional to the product of the Likelihood
and the prior probability. The latter can be freely chosen by the user
whereas the Likelihood is predefined. It is a binned Likelihood which
assumes that the fluctuations in each bin are of Poisson nature and
independent of each other. All channels, processes and sources of
systematic uncertainties are assumed to be uncorrelated.

The parameters of the model are thus the expectation values of the
different processes, \f$\lambda_{k}\f$, and the nuisance parameters,
\f$\delta_{l}\f$.

@subsubsection Excluding systematic uncertainies

In case no sources of systematic uncertainty are taken into account
the likelihood is defined as
\f{eqnarray}
L = \prod_{i=1}^{N_{\mathrm{ch}}} \prod_{j=1}^{N_{\mathrm{bin}}} \frac{\lambda_{ij}^{n_{ij}}}{n_{ij}!} e^{-\lambda_{ij}} \, ,
\label{eqn:likelihood}
\f}
%
where \f$N_{\mathrm{ch}}\f$ and \f$N_{\mathrm{bin}}\f$ are the number of
channels and bins, respectively. \f$n_{ij}\f$ and \f$\lambda_{ij}\f$ are the
observed and expected number of events in the \f$j\f$th bin of the \f$i\f$th
channel. The expected number of events are calculated via
%
\f{eqnarray}
\lambda_{ij} & = & \sum_{k=1}^{N_{\mathrm{p}}} \lambda_{ijk} \\
            & = & \sum_{k=1}^{N_{\mathrm{p}}} \lambda_{k} \cdot f_{ijk} \cdot \epsilon_{ik} \, ,
\label{eqn:expectation}
\f}
%
where \f$f_{ij}\f$ is the bin content of the \f$j\f$th bin in the normalized
template of the \f$k\f$th process in the \f$i\f$th channel. \f$\epsilon_{ik}\f$ is
the efficiency of the \f$k\f$th process in the \f$i\f$th channel specified
when setting the template. \f$\lambda_{k}\f$ is the contribution of the
\f$k\f$th process and is a free parameter of the fit.

\paragraph{Including systematic uncertainties.}
In case sources of systematic uncertainties are taken into account,
the efficiency \f$\epsilon_{ik}\f$ is modified according to a nuisance
parameter:
%
\f{eqnarray}
\epsilon_{ik} \rightarrow \epsilon_{ik} \cdot (1 + \sum_{l=1}^{N_{\mathrm{syst}}} \delta_{l}\cdot \Delta\epsilon_{ijkl}) \, ,
\label{eqn:systematic}
\f}
%
where \f$\delta_{l}\f$ is the nuisance parameter associated with the
source of systematic uncertainty and \f$\Delta\epsilon_{ijkl}\f$ is the
change in efficiency due to the \f$l\f$th source of systematic uncertainty
in the \f$i\f$th channel and \f$j\f$th bin for the \f$k\f$th process.

% --------------------------------------------------------
% Basics
% --------------------------------------------------------
\subsubsection{Creating the fitter}
\label{section:basics}

The main MTF class \verb|mtf| is derived from the \verb|BCModel| class.
A new instance can be created via
%
\begin{verbatim}
BCMTF::BCMTF()
BCMTF::BCMTF(const char * name) ,
\end{verbatim}
%
where the name of the MTF can be specified via argument \verb|name|.

\subsubsection{Adding a channel}

The MTF fits several channels simultaneously. These channels can be
physics channels, e.g., \f$Z^{0}\rightarrow e^{+}e^{-}\f$ and %\f$
\f$Z^{0}\rightarrow \mu^{+}\mu^{-}\f$, samples with disjunct jet
multiplicity or entirely different classes altogether. A new channel
can be added using the following method:
%
\begin{verbatim}
int BCMTF::AddChannel(const char * name) ,
\end{verbatim}
%
where \verb|name| is the name of the process, and the return value is
an error code. Note that at least one channel has to be added.

\subsubsection{Adding a data set}

Each channel added to the MTF has a unique data set which comes in
form of a (\verb|TH1D|) histogram. It can be defined using the
following method:
%
\begin{verbatim}
int BCMTF::SetData(const char * channelname,
                   TH1D hist) ,
\end{verbatim}
%
where \verb|channelname| is the name of the channel and \verb|hist| is
the histogram representing the data. The return value is an error
code.

\subsubsection{Adding a process}

Each template that is fit to the data set corresponds to a process,
where one process can occur in several channels. The fit then defines
the contribution of the process and thus each process comes with one
model parameter. A process can be added using the following method:
%
\begin{verbatim}
int BCMTF::AddProcess(const char * name,
                      double nmin = 0.,
                      double nmax = 1.) ,
\end{verbatim}
%
where \verb|name| is the name of the process and \verb|nmin| and
\verb|nmax| are the lower and upper bound of the parameter associated
with the contribution of the process. The parameter is denoted
\f$\lambda_{k} \, (k=1\dots N_{\mathrm{p}})\f$ in
section~\ref{section:math}. Note that at least one process has to be
added. A prior needs to be defined for each process, using the default
\verb|BCModel| methods.

It is likely that a single process will have different shapes in
different channels. Thus, templates for a process need to be defined
for each channel separately using the following method:
%
\begin{verbatim}
int BCMTF::SetTemplate(const char * channelname,
                       const char * processname,
                       TH1D hist,
                       double efficiency = 1.) ,
\end{verbatim}
%
where \verb|channelname| and \verb|processname| are the names of the
channel and the process, respectively. The parameter \verb|hist| is
the (\verb|TH1D|) histogram (or template) which represents the
process. The histogram will be normalized to unity and the entries in
the normalized histogram are the probabilities to find an event of a
process~\f$k\f$ and channel~\f$i\f$ in bin~\f$j\f$. This probability is denoted
\f$f_{ijk}\,(i=1\dots N_{\mathrm{ch}}, \, j=1\dots N_{\mathrm{b}}, \,
k=1\dots N_{\mathrm{p}})\f$ in section~\ref{section:math}. The last
parameter, \verb|efficiency|, is the efficiency of the process in that
channel and is used to scale to template during the fit. This is
needed if a process contributes with different amounts in two separate
channels. The efficiency is denoted \f$\epsilon_{ik} \, (i=1\dots
N_{\mathrm{ch}}, \, k=1\dots N_{\mathrm{p}})\f$ in
section~\ref{section:math}. The return value is an error code. Note
that templates do not have to be set if the process does not
contribute to a particular channel.

\subsubsection{Adding systematic uncertainties}

Systematic uncertainties can alter the shape of a template. Sources of
systematic uncertainty can be included in the fit using nuisance
parameters. This nuisance parameter is assumed to alter the original
template linearly, where values of~-1,~0, and~1 correspond to the
``downwards'' shifted, nominal and ``upwards shifted'' template,
respectively. The nuisance parameters are denoted \f$\delta_{l} \,
(l=1\dots N_{\mathrm{syst}})\f$ in section~\ref{section:math}. Shifted
refers to a change of one standard deviation. An example for a
nuisance parameter could be the jet energy scale (JES). With a nominal
JES of~1 and and uncertainty of 5\%, the scaled templates correspond
to a JES of 0.95 and 1.05, respectively. A prior needs to be defined
for each nuisance parameter which is usually chosen to be a standard
normal distribution. A source of systematic uncertainty can be added
using the following method:
%
\begin{verbatim}
int BCMTF::AddSystematic(const char * name,
                         double min = -5.,
                         double max = 5.) ,
\end{verbatim}
%
where \verb|name| is the name of the source of systematic uncertainty
and \verb|min| and \verb|max| are the lower and upper bound of the
nuisance parameter, respectively. The return value is an error code.

Since the different sources of systematic uncertainty have an
individual impact on each process and in each channel, these need to
be specified. Two method can be used to define the impact:
%
\begin{verbatim}
int BCMTF::SetSystematicVariation(const char * channelname,
                                  const char * processname,
                                  const char * systematicname,
                                  TH1D hist_up,
                                  TH1D hist_down) ,
\end{verbatim}
%
where \verb|channelname|, \verb|processname| and \verb|systematicname|
are the names of the channel, the process and the source of systematic
uncertainty. The (\verb|TH1D|) histograms \verb|hist_up| and
\verb|hist_down| are the histograms corresponding to an ``up''- and
``down''-scaling of the systematic uncertainty of one standard
deviation, i.e., for each bin entry \f$y\f$ they are calculated as
%
\f{eqnarray}
\label{eqn:deltay}
\Delta_{\mathrm{up}}   & = & (y_{\mathrm{up}} - y_{\mathrm{nominal}})/y_{\mathrm{nominal}} \, , \\
\Delta_{\mathrm{down}} & = & (y_{\mathrm{nominal}} - y_{\mathrm{down}})/y_{\mathrm{nominal}} \, .
\f}
%
Note the sign of the down-ward fluctuation. These histograms define
the change of the bins in each template in the efficiency which is
denoted \f$\Delta\epsilon_{ijkl} \, (i=1\dots N_{\mathrm{ch}}, \,
j=1\dots N_{\mathrm{b}}, \, k=1\dots N_{\mathrm{p}}, \, l=1\dots
N_{\mathrm{syst}})\f$. For example, if the value for a particular bin of
\verb|hist_up| is 0.05, i.e., if the systematic uncertaintiy is 5\% in
that bin, then the efficiency of the process in that channel will be
multiplied by \f$(1+0.05)\f$. The return value is an error code.

The second variant does not take the difference in efficiency, but
calculates it internally from the absolute values:
%
\begin{verbatim}
int BCMTF::SetSystematicVariation(const char * channelname,
                                  const char * processname,
                                  const char * systematicname,
                                  TH1D hist,
                                  TH1D hist_up,
                                  TH1D hist_down) ,
\end{verbatim}
%
where \verb|channelname|, \verb|processname| and \verb|systematicname|
are the names of the channel, the process and the source of systematic
uncertainty. The (\verb|TH1D|) histograms \verb|hist|, \verb|hist_up|
and \verb|hist_down| are the nominal histogram and the histograms
corresponding to an ``up''- and ``down''-scaling of the systematic
uncertainty of one standard deviation. In this case, the histograms
are not the relative differences but the absolute values. The return
value is an error code.

\subsubsection{Running the fit}

The fit can be started using one of the standard \verb|BCModel| fitting methods,
e.g.
%
\begin{verbatim}
BCMTF::MarginalizeAll() ,
BCMTF::FindMode() .
\end{verbatim}
