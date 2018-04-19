Getting started {#cha-basics}
=================
[TOC]

To demonstrate the basic usage of BAT, we will build an example analysis step by
step. Step one, naturally, is to install BAT---please refer to the installation
chapter in this [manual](INSTALL.md) or <a
href="https://github.com/bat/bat/blob/master/INSTALL.md">online</a>. To check if
you have BAT installed and accessible, run the following command

    bat-config

in your terminal. It should output a usage statement. This program
outputs information about your BAT installation:
<table>
<tr><th> `bat-config` Flag   <th> Returned Information
<tr><td> `--prefix`          <td> Path to BAT installation
<tr><td> `--version`         <td> BAT version identifier
<tr><td> `--libs`            <td> Linker flags for compiling code using BAT
<tr><td> `--cflags`          <td> C++ compiler flags for compiling code using BAT
<tr><td> `--bindir`          <td> Path to BAT binaries directory
<tr><td> `--incdir`          <td> Path to BAT include directory
<tr><td> `--libdir`          <td> Path to BAT libraries directory
</table>

BAT has a second executable, which creates for you the files necessary
to start a basic analysis. We will use this executable to
initialize our tutorial project:

    bat-project MyTut MyMod

This will create a directory called `MyTut` that contains
<table>
<tr><td> `Makefile`      <td> a makefile to compile our tutorial project
<tr><td> `runMyTut.cxx`  <td> the C++ source to an executable to run our tutorial project
<tr><td> `MyMod.h`       <td> the C++ header for our tutorial model `MyMod`
<tr><td> `MyMod.cxx`     <td> the C++ source for our tutorial model `MyMod`
</table>

If BAT is installed correctly, you can compile and run this project already:

    cd MyTut
    make
    ./runMyTut

BAT will issue a series of errors telling you your model has no
parameters. Because of course we haven't actually put anything into
our model yet. Let's do that.

@section sec-basics-define-a-model Defining a model

To define a valid BAT model we must make three additions to the empty
model that `bat-project` has created:
1. we must add parameters to our
model;
2. we must implement a log-likelihood function;
3. and we must implement or state our priors.

We will start with a simple model that fits a normal distribution to
data, and so has three parameters: the mode (\f$\mu\f$) and standard
deviation (\f$\sigma\f$) of our distribution, and a scaling factor
("height"). We will start with flat priors for all.

How you store and access your data is entirely up to you. For this
example, we are going to fit to a binned data set that we store as a
ROOT histogram in a private member of our model class. In the header,
add to the class

@code{.cpp}

#include <TH1D.h>

...

class MyMod : public BCModel
{

    ...

private:
    TH1D fDataHistogram

};
@endcode

And in the source file, we initialize our `fDataHistogram` in the
constructor; let's also fill it with some random data, which ROOT can
do for us using a TF1:

@code{.cpp}
#include "MyMod.h"
#include <TF1.h>

...

// ---------------------------------------------------------
MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", ";mass [GeV];count", 100, 5.0, 5.6)
{
    // create function to fill data according to
    TF1 data_func("data_func", "exp(-0.5*((x - 5.27926) / 0.04)^2)", 5.0, 5.6);

    // fill data histogram randomly from data_func 1,000 times
    fDataHistogram.FillRandom("data_func", 1e3);
}
@endcode

@subsection subsec-basics-add-parameters Adding parameters

To add a parameter to a model, call its member function `BCModel::AddParameter`:
@code{.cpp}
bool AddParameter(const std::string& name, double min, double max, const std::string& latexname = "", const std::string& unitstring = "")
@endcode

You must indicate a parameter's name and its allowed range, via `min`
and `max`. Each parameter must be added with a unique name. BAT will
also create a "safe version" of the name which removes all non
alpha-numeric characters except for the underscore, which is needed
for naming of internal storage objects. Each parameter name should
also convert to a unique safe name; BAT will complain if this is not
the case (and `AddParameter` will return `false`).

Optionally, you may add a display name for the parameter and a unit
string, both of which are used when BAT creates plots.

The most logical place to add parameters to a model is in its
constructor. Let us now edit `MyMod.cxx` to add our mode, standard
deviation, and height parameters inside the constructor:

@code{.cpp}
MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", ";mass [GeV];count", 100, 5.0, 5.6)
{
    ...

    AddParameter("mu",    5.27, 5.29, "#mu", "[GeV]");
    AddParameter("sigma", 25e-3, 45e-3, "#sigma", "[GeV]");
    AddParameter("height", 0, 10, "", "[events]");
}
@endcode

I have chosen the ranges because I have prior information about the
data: the mode of my distribution will be between 5.27 and 5.29; the
standard deviation will be between 0.025 and 0.045; and the height
will be between 0 and 10. When we provide a prior for a parameter,
\f$p_0(\lambda)\f$, the prior BAT uses is

\f{equation}{
    P_0(\lambda) = \begin{cases} p_0(\lambda),& \text{if }
    \lambda\in[\lambda_{\text{min}}, \lambda_{\text{max}}],\\ 0,&
    \text{otherwise.}
    \end{cases}
\f}

So keep mind that the ranges your provide for parameters become part
of the prior: BAT will not explore parameter space outside of the
range limits you provide.

If you have written the code correctly for adding parameters to your
model, your code should compile. BAT will again, though, issue an
error if you try to run `runMyTut`, since we are still missing priors.

@subsection subsec-basics-setting-priors Setting prior distributions

There are two ways we may set the prior, \f$P_0(\vec\lambda)\f$ for a parameter
point: We can override the function that returns the log _a priori_ probability
for a model (`BCModel::LogAPrioriProbability`) and code anything we can dream of
in C++:

@code{.cpp}
double MyMod::LogAPrioriProbability(const std::vector<double>& pars)
{
   ...
}
@endcode

In this case, you have to make sure the prior is properly normalized if you wish to compare different models as BAT will not do this for you.

Or we can set individual priors for each parameter, and BAT will multiply them
together for us including the proper normalization. If each parameter has a
prior that factorizes from all the others, then this is the much better option.
In this case, we **must not** override `LogAPrioriProbability`. Instead we tell
BAT what the factorized prior is for each parameter, by adding a `BCPrior`
object to each parameter after we create it:

@code{.cpp}
#include <BAT/BCGaussianPrior.h>

...

MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", ";mass [GeV];count", 100, 5.0, 5.6)
{
    ...

    // add parameters for Gaussian distribution
    AddParameter("mu",    5.27, 5.29, "#mu", "[GeV]");
    GetParameters().Back().SetPrior(new BCGaussianPrior(5.28, 2e-3));

    AddParameter("sigma", 25e-3, 45e-3, "#sigma", "[GeV]");
    GetParameters().Back().SetPrior(new BCGaussianPrior(35e-3, 3e-3));

    AddParameter("height", 0, 10, "", "[events]");
    GetParameters().Back().SetPriorConstant();
}
@endcode

We could create a `BCConstantPrior` object for `height` just as we
created a `BCGaussianPrior` object for `mu` and `sigma`, but BAT has a
convenient function that creates one for us.

We can now compile our code and run it. The results will be
meaningless, though, since we have yet to define our likelihood.

@subsection subsec-basics-loglikelihood Defining a likelihood

The heart of our model is the likelihood function---more specifically,
since BAT works with the natural logarithm of functions, our
log-likelihood function.

Given a histogrammed data set containing numbers of events in each
bin, our statistical model is the product of Poisson probabilities for
each bin---the probability of the events observed in the bin given the
expectation from our model function:

\f{equation}{
    \mathcal{L}(\vec\lambda) \equiv \prod_{i} \text{Poisson}(n_i|\nu_i),
\f}

where \f$n_i\f$ is the number of events in bin \f$i\f$ and \f$\nu_i\f$
is the expected number of events given our model, which we will take
as the value of our model function at the center of the bin,
\f$m_i\f$,

\f{equation}{
    \nu_i \equiv f(m_i) = \frac{1}{\sqrt{2\pi}\sigma}\exp{\left(-\frac{(m_i - \mu)^2}{2\sigma^2}\right)}.
\f}

Working with the log-likelihood, this transforms into a sum:

\f{equation}{
    \log\mathcal{L}(\vec\lambda) \equiv \sum_{i} \log\text{Poisson}(n_i|\nu_i).
\f}

BAT conveniently has a function to calculate the logarithm of the
Poisson distribution (with observed \f$x\f$ and expected
\f$\lambda\f$) for you:

@code{.cpp}
double BCMath::LogPoisson(double x, double lambda)
@endcode

Let us code this into our log-likelihood function:

@code{.cpp}

#include <BAT/BCMath.h>

#include <TMath.h>

...

// ---------------------------------------------------------
double MyMod::LogLikelihood(const std::vector<double>& pars)
{
    // store our log-likelihood as we loop through bins
    double LL = 0.;

    // loop over bins of our data
    for (int i = 1; i <= fDataHistogram.GetNbinsX(); ++i) {

        // retrieve observed number of events
        double x = fDataHistogram.GetBinContent(i);

        // retrieve bin center
        double m = fDataHistogram.GetBinCenter(i);

        // calculate expected number of events, using ROOT Gaus function
        double nu = TMath::Gaus(m, pars[0], pars[1], true);

        // add to log-likelihood sum
        LL += BCMath::LogPoisson(x, nu);

    }

    // return log-likelihood
    return LL;
}
@endcode

@section sec-basics-output Looking at the output

Our model class is now ready to go. Let's just make one edit to the
`runMyTut.cxx` file to change our model name from the default
"name_me" to something sensible:

@code{.cpp}
// create new MyMod object
MyMod m("gaus_mod");
@endcode

We can now compile and run our project:
@code
make
./runMyTut
@endcode

This will sample from the posterior probability distribution and
marginalize the results, saving plots to `gaus_mod_plots.pdf`. In that
file you should see the 1D and 2D marginalizations of our three
parameters:

@image html gaus_mod_plots-1.png "1D Posteriors of Gaussian distribution model."\
@image latex gaus_mod_plots-1.pdf "1D Posteriors of Gaussian distribution model." width=\textwidth

@image html gaus_mod_plots-2.png "2D Posteriors of Gaussian distribution model."\
@image latex gaus_mod_plots-2.pdf "2D Posteriors of Gaussian distribution model." width=\textwidth

In each plot, we see the global mode and the marginalized mean and standard
deviation of the posterior distribution; and three credibility
intervals.

BAT has also printed a summary of the results to the log file and
command line. This includes the global mode

    Summary :  Global mode:
    Summary :  0)  Parameter "mu"          : 5.279741 +- 0.001072076
    Summary :  1)  Parameter "sigma"       : 0.040267 +- 0.00089553
    Summary :  2)  Parameter "height"      : 6 +- 0.1897

and marginalized posteriors:

    Summary :   (0) Parameter "mu" :
    Summary :       Mean +- sqrt(Variance):         5.279748 +- 0.001065279
    Summary :       Median +- central 68% interval: 5.279749 +  0.00105881 - -0.00106418
    Summary :       (Marginalized) mode:            5.2797
    Summary :        5% quantile:                   5.278005
    Summary :       10% quantile:                   5.278381
    Summary :       16% quantile:                   5.278685
    Summary :       84% quantile:                   5.280808
    Summary :       90% quantile:                   5.281112
    Summary :       95% quantile:                   5.281498
    Summary :       Smallest interval containing 69.8% and local mode:
    Summary :       (5.2786, 5.2808) (local mode at 5.2797 with rel. height 1; rel. area 1)

@section sec-basics-observables Adding an observable

We can store the posterior distribution of any function of the
parameters of our model using BAT's `BCObservable`'s. Let us suppose
we want to know the posterior distribution for the total number of
events our model predicts---the signal yield; and we want to know the
standard deviation's relation to the mean: \f$\sigma/\mu\f$. In our
constructor, we add the observables in the same way we added our
parameters:

@code{.cpp}
// ---------------------------------------------------------
MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", "mass [GeV];count", 100, 5.0, 5.6)
{
    ...

    AddObservable("SignalYield", 900, 1100, "Y_{S}", "[events]");
    AddObservable("Resolution",
                  100. * GetParameter("sigma").GetLowerLimit() / GetParameter("mu").GetUpperLimit(),
                  100. * GetParameter("sigma").GetUpperLimit() / GetParameter("mu").GetLowerLimit(),
                  "#sigma / #mu", "[%]");
}
@endcode

Note that an observable does not need a prior since it is not a
parameter; but it does need a range for setting the marginalized
histogram's limits. Since we know already that the answer will be 1000
events with an uncertainty of \f$\sqrt{1000}\f$, we set the range for
the yield to \f$[900, 1100]\f$.

And we need to calculate this observable in our
`CalculateObservables(...)` member function. Uncomment it in the
header file `MyMod.h`:

@code{.cpp}
    void CalculateObservables(const std::vector<double> & pars);
@endcode

and implement it in the source:
@code{.cpp}
// ---------------------------------------------------------
void MyMod::CalculateObservables(const std::vector<double>& pars)
{
    // store total of number events expected
    double nu = 0;

    // loop over bins of our data
    for (int i = 1; i <= fDataHistogram.GetNbinsX(); ++i)
        // calculate expected number of events in that bin
        // and add to total expectation
        nu += pars[2] * TMath::Gaus(fDataHistogram.GetBinCenter(i), pars[0], pars[1], true);

    // store in the observable
    GetObservable(0) = nu;

    // Store sigma as percentage of mu:
    GetObservable(1) = 100. * pars[1] / pars[0];
}
@endcode

Compile and run `runMyTut.cxx` and you will see new marginalized
distributions and text output for the observables:

    Summary :   (3) Observable "SignalYield" :
    Summary :       Mean +- sqrt(Variance):         1000.9 +- 31.53
    Summary :       Median +- central 68% interval: 1000.6 +  31.937 - -31.194
    Summary :       (Marginalized) mode:            1003
    Summary :        5% quantile:                   949.52
    Summary :       10% quantile:                   960.52
    Summary :       16% quantile:                   969.44
    Summary :       84% quantile:                   1032.6
    Summary :       90% quantile:                   1041.9
    Summary :       95% quantile:                   1053.5
    Summary :       Smallest intervals containing 68.7% and local modes:
    Summary :       (970, 1032) (local mode at 1003 with rel. height 1; rel. area 0.97769)
    Summary :       (1032, 1034) (local mode at 1033 with rel. height 0.57407; rel. area 0.022311)

As we expected, the mean of the total yield posterior is just the
number of events in our data set and the standard deviation is its
square root.

@section sec-basics-other-output Further Output

BAT has a few more output options than the two mentioned above. The
code for turning them on is already included in the `runMyTut.cxx`
generated by `bat-project`. You can write the Markov-chain samples to
a ROOT `TTree` for further postprocessing by uncommenting the following line:

@code{.cpp}
m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");
@endcode

You can also modify the plotting output to include more plots per
page. Without specifying, we used the default of one plot per
page. Let's instead plot 4 plots (in a 2-by-2 grid):

@code{.cpp}
m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf", 2, 2);
@endcode

There are four more graphical outputs from BAT. Turn them on by uncommenting the following lines
@code{.cpp}
m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf", 3, 2);
@endcode
(We have edited the last line to specify 3-by-2 printing.)

The parameter plot graphically summarizes the output for all
parameters (and in a separate page, all observables) in a single
image:

@image html gaus_mod_parameters-1.png "Summary of parameter marginalizations."
@image latex gaus_mod_parameters-1.pdf "Summary of parameter marginalizations." width=\textwidth

The correlation plot and the correlation matrix summarize graphically
the correlations among parameters and observables:

@image html gaus_mod_correlation.png "Parameter and observable correlation plots."
@image latex gaus_mod_correlation.pdf "Parameter and observable correlation plots." width=\textwidth

@image html gaus_mod_correlationMatrix.png "Parameter and observable correlation matrix."
@image latex gaus_mod_correlationMatrix.pdf "Parameter and observable correlation matrix." width=\textwidth

The knowledge update plots show the marginalized priors and
marginalized posteriors together in one plot for each variable (and
also for 2D marginalizations):

@image html gaus_mod_update-1.png "Knowledge update plots."
@image latex gaus_mod_update-1.pdf "Knowledge update plots." width=\textwidth

Note that this really shows the _marginalized_ prior: Here, having
used factorized priors, we see exactly the factorized priors. But had
we used a multivariate prior, then we'd see what this looks like for
each parameter. We also see what the prior of our observables look
like given the priors of the variables they are functions of. These
are very useful things to know.
