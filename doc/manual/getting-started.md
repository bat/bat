Getting started {#cha-basics}
=================
[TOC]

To demonstrate the basic usage of BAT, we will build an example
analysis step by step. Step one, naturally, is to install BAT---please
refer to chapter in this manual concerning installation. To check if
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
to start a basic analysis with. We will use this executable to
initialize our tutorial project:

    bat-project MyTut MyMod

This will create a directory called `MyTut` that contains
<table>
<tr><td> `Makefile`     <td> a makefile to compile our tutorial project
<tr><td> `runMyTut.cxx` <td> the C++ source to an executable to run our tutorial project
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
model that `bat-project` has created: we must add parameters to our
model; we must implement a log-likelihood function; and we must
implement or state our priors.

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
contructor; let's also fill it with some random data, which ROOT can
do for us using a TF1:
@code{.cpp}
#include <TF1.h>

// ---------------------------------------------------------
MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", ";mass [GeV];count", 100, 5.0, 5.6)
{
    // create function to fill data according to
    TF1 data_func("data_func", "exp(-0.5*((x - 5.28) / 0.04)^2)", 5.0, 5.6);

    // fill data histogram randomly from data_func 1,000 times
    fDataHistogram.FillRandom("data_func", 1e3);

    ...

}
@endcode

@subsection subsec-basics-add-parameters Adding parameters

To add a parameter to a model, call its member function `#AddParameter`:
@code{.cpp}
bool AddParameter(const std::string& name, double min, double max, const std::string& latexname = "", const std::string& unitstring = "")
@endcode

You must indicate a parameter's name and its allowed range, via `min`
and `max`. Each parameter must be added with a unique name. BAT will
also create a "safe version" of the name which removes all non
alpha-numeric characters except for the underscore, which is needed
for naming of internal storage objects. Each parameters name should
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
    AddParameter("sigma", 35e-3, 45e-3, "#sigma", "[GeV]");
    AddParameter("height", 0, 10, "", "[events]");
}
@endcode

I have chosen the ranges because I have prior information about the
data: the mode of my distribution will be between 5.27 and 5.29; the
standard deviation will be between 0.035 and 0.045; and the height
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

There are two ways we may set the prior, \f$P_0(\vec\lambda)\f$ for a
parameter point: We can override the function that returns the log
_a priori_ probability for a model and code anything we can dream of in C++:

@code{.cpp}
double MyMod::LogAPrioriProbability(const std::vector<double>& pars)
{
   ...
}
@endcode

Or we can set individual priors for each parameter, and BAT will
multiply them together for us. If each parameter has a prior that
factorizes from all the others, then this is the much better
option. In this case, we **must not** override
`LogAPrioriProbability`. Instead we tell BAT what the factorized prior
is for each parameter. Let us start with each parameter having a flat
prior, which we can tell BAT about with one line added to our
constructor:

@code{.cpp}
MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", ";mass [GeV];count", 100, 5.0, 5.6)
{
    ...
    
    // add parameters for Gaussian distribution
    AddParameter("mu",    5.2, 5.4, "#mu", "[GeV]");
    AddParameter("sigma", 10e-3, 100e-3, "#sigma", "[GeV]");
    AddParameter("yield", 0, 1e5, "", "[events]");

    // Set all priors flat
    SetPriorConstantAll();
}
@endcode

Naturally, this must be done after we have added the parameters.

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

BAT convienently has a function to calculate the logarithm of the
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

@section sec-basics-output Looking at the ouput

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

@image html MyTut/gaus_mod_plots_1D.png "1D Posteriors of Gaussian distribution model."\
@image latex MyTut/gaus_mod_plots_1D.pdf "1D Posteriors of Gaussian distribution model." width=\textwidth

@image html MyTut/gaus_mod_plots_2D.png "2D Posteriors of Gaussian distribution model."\
@image latex MyTut/gaus_mod_plots_2D.pdf "2D Posteriors of Gaussian distribution model." width=\textwidth

In each plot, we see the global mode and the marginalized mean and standard
deviation of the posterior distribution; and three credibility
intervals.

BAT has also printed a summary of the results to the log file and
command line. This includes the global mode

    Summary :  Global mode:
    Summary :  0)  Parameter "mu"          : 5.280326 +- 0.00126801
    Summary :  1)  Parameter "sigma"       : 0.040206 +- 0.0008942
    Summary :  2)  Parameter "height"      : 6 +- 0.1897

and marginalized posteriors:

	Summary :   (0) Parameter "mu" :
	Summary :       Mean +- sqrt(Variance):         5.280331 +- 0.001277385
	Summary :       Median +- central 68% interval: 5.280331 +  0.001264596 - -0.001275553
	Summary :       (Marginalized) mode:            5.2803
	Summary :        5% quantile:                   5.278232
	Summary :       10% quantile:                   5.278694
	Summary :       16% quantile:                   5.279056
	Summary :       84% quantile:                   5.281596
	Summary :       90% quantile:                   5.28197
	Summary :       95% quantile:                   5.282441
	Summary :       Smallest interval containing 69.2% and local mode:
	Summary :       (5.279, 5.2816) (local mode at 5.2803 with rel. height 1; rel. area 1)

@section sec-basics-observables Adding an observable

We can store the posterior distribution of any function of the
parameters of our model using BAT's `BCObversable`'s. Let us suppose
we want to know the posterior distribution for the total number of
events our model predicts---the signal yield. In our constructor, we
add the observable in the same way we added our parameters:

@code{.cpp}
// ---------------------------------------------------------
MyMod::MyMod(const std::string& name)
    : BCModel(name),
      fDataHistogram("data", "mass [GeV];count", 100, 5.0, 5.6)
{
    ...

    AddParameter("mu",    5.27, 5.29, "#mu", "[GeV]");
    ...
    
    AddObservable("SignalYield", 900, 1100, "Y_{S}", "[events]");
}
@endcode

Note that an observable does not need a prior since it is not a parameter.

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
}
@endcode

Compile and run `runMyTut.cxx` and you will see new marginalized
distributions for the yield and a summary of the yield observable in
the text output:

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

As we should expect, the total mean of the total yield posterior is
just the number of events in our data set and the standard deviation
is the square root of the number of events.

@todo This code should be one of the examples and tested to be
compilable. Does doxygen support including some lines of a source
file? Follow the DRY principle

@todo Dan, copy over from manual.tex

<!-- %%%%%%%%%%%%%%%%%%%%%%%%% -->
@section sec-bat-project bat-project
