// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// a graph using a function defined by the user. TGraphErrors
// with defined y uncertainties has to be defined. In the fit the
// uncertainties are considered to be Gaussian.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x graphFitterAdvancedExample.C
//
// or
//
//    root[1] .L graphFitterAdvancedExample.C
//    root[2] graphFitterAdvancedExample()
//
// or from the command line
//
//    $ root graphFitterAdvancedExample.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x graphFitterAdvancedExample.C+
//
// See ROOT documentation for details.
//
//
// Below are the includes needed for compilation of the macro.
// The #if ... #endif directives around the includes allow to
// run the macro in both normal and compiled mode.
#if defined(__MAKECINT__) || defined(__ROOTCLING__) || (!defined(__CINT__) && !defined(__CLING__))

#include <BAT/BCLog.h>
#include <BAT/BCGraphFitter.h>

#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include <string>
#include <vector>

#endif

// The macro performs fits with four different functions (models):
//   2nd order polynomial
//   gaussian peak + constant
//   gaussian peak + straight line
//   gaussian peak + 2nd order polynomial

// ----------------------------------------------------------------------
void graphFitterAdvancedExample()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // -------------------------
    // Prepare fitting functions
    // (all parameter limits must be set!)

    // 2nd order polynomial
    TF1 f1("f1", "[0] + [1]*x + [2]*x^2", 0., 100.);
    f1.SetParNames("offset", "slope", "quad");
    f1.SetParLimits(0,   0.,    5.);  // offset
    f1.SetParLimits(1,   0.,    1.2); // slope
    f1.SetParLimits(2,  -0.1,   0.1); // quad

    // gaussian + constant
    TF1 f2("f2", "[0]/sqrt(2*pi)/[2] * exp(-0.5*((x-[1])/[2])^2) + [3]", 0., 100.);
    f2.SetParNames("A_gauss", "mean", "sigma", "offset");
    f2.SetParLimits(0,   0.,  200.); // A_gauss
    f2.SetParLimits(1,   2.,   18.); // mean
    f2.SetParLimits(2,   .2,    4.); // sigma
    f2.SetParLimits(3,   0.,   10.); // offset

    // gaussian + line
    TF1 f3("f3", "[0]/sqrt(2*pi)/[2] * exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x", 0., 100.);
    f3.SetParNames("A_gauss", "mean", "sigma", "offset", "slope");
    f3.SetParLimits(0,   0.,  200.); // A_gauss
    f3.SetParLimits(1,   2.,   18.); // mean
    f3.SetParLimits(2,   .2,    4.); // sigma
    f3.SetParLimits(3,   0.,   10.); // offset
    f3.SetParLimits(4,   0.,    2.); // slope

    // gaussian + 2nd order polynomial
    TF1 f4("f4", "[0]/sqrt(2*pi)/[2] * exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2", 0., 100.);
    f4.SetParNames("A_gauss", "mean", "sigma", "offset", "slope", "quad");
    f4.SetParLimits(0,   0.,  200.); // A_gauss
    f4.SetParLimits(1,   2.,   18.); // mean
    f4.SetParLimits(2,   .2,    4.); // sigma
    f4.SetParLimits(3,   0.,   10.); // offset
    f4.SetParLimits(4,   0.,    2.); // slope
    f4.SetParLimits(5,   0.,    .5); // quad

    // -------------------------

    //-------------------------
    // Create data
    TGraphErrors gr;

    if (true) { // set false to generate data anew (below)

        // read data from the file supplied with BAT distribution
        // is data/datax.txt and contains the data which were
        // used in the BAT paper
        gr = TGraphErrors("data.txt", "%lg %lg %lg");

    } else {

        // generate data anew from f4 defined above, with parameters:
        f4.SetParameters(150, 5, 0.5, 0, 50, 2);

        // and smear from value by Gaussian with sigma:
        double sigmay = 40.;

        // initialize random number generator
        TRandom3 random(1234);

        // define range to generate data in
        double xmin = 0.1;
        double xmax = 19.9;

        // define number of points to generate
        unsigned n = 100;

        double dx = (xmax - xmin) / (n - 1);

        for (unsigned i = 0; i < n; ++i) {

            // get x value
            double x = xmin + i * dx;

            // smear from f4 value
            gr.SetPoint(i, x, random.Gaus(f4.Eval(x), sigmay));
            // set uncertainty to tenth of smearing sigma
            gr.SetPointError(i, 0, sigmay / 10);
        }
    }
    // -------------------------


    // setup all graph fitters
    // old-fashioned approach for compatibility with CINT from root 5
    const unsigned nmodels = 4;
    BCGraphFitter* models[] = {
        new BCGraphFitter(gr, f1, f1.GetName()),
        new BCGraphFitter(gr, f2, f2.GetName()),
        new BCGraphFitter(gr, f3, f3.GetName()),
        new BCGraphFitter(gr, f4, f4.GetName())
    };

    // perform the analysis on all models
    for (unsigned i = 0; i < nmodels; ++i) {
        // Set precision to medium or high because there are multiple
        // modes for f3 and f4 and the chains need time to explore both.
        // Only the multivariate proposal does a good job here.
        models[i]->SetPrecision(BCEngineMCMC::kMedium);

        // run the fitting (MCMC + Minuit)
        models[i]->Fit();

        // print all marginalized distributions into a PostScript file
        models[i]->PrintAllMarginalized("plots-" + models[i]->GetSafeName() + ".pdf");

        // draw summary plots and tables
        models[i]->PrintParameterPlot("summary_pars-" + models[i]->GetSafeName() + ".pdf");
        // models[i]->PrintCorrelationPlot("summary_corr-" + models[i]->GetSafeName() + ".pdf");
        // models[i]->PrintCorrelationMatrix("summary_corr_matrix-" + models[i]->GetSafeName() + ".pdf");
        // models[i]->PrintKnowledgeUpdatePlots("summary_update-" + models[i]->GetSafeName() + ".pdf");
        // models[i]->PrintParameterLatex("summary_pars-" + models[i]->GetSafeName() + ".tex");
    }

    TCanvas c;

    // draw fit including the error band
    c.Divide(2, 2);
    for (unsigned i = 0; i < nmodels; ++i) {
        c.cd(i + 1);
        models[i]->DrawFit();
    }
    c.Print("data-all-band.pdf");

    // draw data and all fits in the same plot (w/o error bands)
    c.Clear();
    c.cd(1);

    gr.SetMarkerStyle(20);
    gr.SetMarkerSize(.5);
    gr.Draw("ap");

    for (unsigned i = 0; i < nmodels; ++i) {
        models[i]->GetFitFunction().SetLineColor(i + 1);
        models[i]->GetFitFunction().SetLineWidth(2);
        models[i]->GetFitFunction().Draw("l same");
    }
    c.Print("data-all.pdf");

    // clean up
    for (unsigned i = 0; i < nmodels; ++i)
        delete models[i];
}
