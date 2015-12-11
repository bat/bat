/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "BCGraphFitter.h"

#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCLog.h>
#include <BAT/BCMath.h>

#include "BCGraphFitter.h"

#include <BAT/BCDataPoint.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCLog.h>
#include <BAT/BCMath.h>

#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <Math/ProbFuncMathCore.h>
#include <TPad.h>
#include <TString.h>

// ---------------------------------------------------------
BCGraphFitter::BCGraphFitter(const TGraphErrors& graph, const TF1& func, const std::string& name)
    : BCFitter(func, name),
      fGraph(graph)
{
    // why not just one point?
    if (fGraph.GetN() <= 1) {
        throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": TGraphErrors needs at least two points.");
    }

    double* x  = fGraph.GetX();
    double* y  = fGraph.GetY();
    double* ex = fGraph.GetEX();
    double* ey = fGraph.GetEY();

    if (!ey) {
        throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) + ": TGraphErrors has NO errors set on Y. Not able to fit.");
    }

    // fill the dataset
    for (int i = 0; i < fGraph.GetN(); ++i) {
        // create the data point
        GetDataSet()->AddDataPoint(BCDataPoint(4));
        GetDataSet()->Back()[0] = x[i];
        GetDataSet()->Back()[1] = y[i];
        GetDataSet()->Back()[2] = ex ? ex[i] : 0;
        GetDataSet()->Back()[3] = ey[i];
    }
    // adjust bounds for 1 sigma of the uncertainties
    GetDataSet()->AdjustBoundForUncertainties(0, 1, 2); // x +- 1*errx
    GetDataSet()->AdjustBoundForUncertainties(1, 1, 3); // y +- 1*erry

    SetFitFunctionIndices(0, 1);

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
BCGraphFitter::~BCGraphFitter() {}

// ---------------------------------------------------------
double BCGraphFitter::LogLikelihood(const std::vector<double>& params)
{
    TF1& f = GetFitFunction();

    // set the parameters of the function
    // passing the pointer to first element of the vector is
    // not completely safe as there might be an implementation where
    // the vector elements are not stored consecutively in memory.
    // however it is much faster than copying the contents, especially
    // for large number of parameters
    f.SetParameters(&params[0]);

    // initialize probability
    double logl = 0.;

    // loop over all data points
    for (unsigned i = 0; i < GetNDataPoints(); i++)
        // calculate log of probability assuming
        // a Gaussian distribution for each point
        logl += BCMath::LogGaus(GetDataSet()->GetDataPoint(i)[1], // y value of point
                                f.Eval(GetDataSet()->GetDataPoint(i)[0]), // f(x value of point)
                                GetDataSet()->GetDataPoint(i)[3], // uncertainty on y value of point
                                true); // include normalization factor

    return logl;
}

// ---------------------------------------------------------
void BCGraphFitter::Fit()
{
    // check setup
    BCLog::OutDetail(Form("Fitting %d data points with function of %d parameters", GetNDataPoints(), GetNParameters()));
    if (GetNDataPoints() <= GetNParameters()) {
        BCLog::OutWarning(Form("Number of parameters (%d) lower than or equal to number of points (%d).", GetNParameters(), GetNDataPoints()));
        BCLog::OutWarning("Fit doesn't have much meaning.");
    }

    // perform marginalization
    MarginalizeAll();

    // maximize posterior probability, using the best-fit values close
    // to the global maximum from the MCMC
    BCIntegrate::BCOptimizationMethod method_temp = GetOptimizationMethod();
    SetOptimizationMethod(BCIntegrate::kOptMinuit);
    FindMode(GetBestFitParameters());
    SetOptimizationMethod(method_temp);

    // p value with (approximate) correction for degrees of freedom
    CalculatePValue(GetBestFitParameters(), true);

    // print summary to screen
    PrintShortFitSummary();
}

// ---------------------------------------------------------
void BCGraphFitter::DrawFit(const std::string& options, bool flaglegend)
{
    // check wheather options contain "same"
    TString opt = options;
    opt.ToLower();

    // if not same, draw the histogram first to get the axes
    if (!opt.Contains("same"))
        fGraph.Draw("ap");

    // draw the error band as central 68% probability interval
    TGraph* errorBand = GetErrorBandGraph(0.16, 0.84);
    fObjectTrash.Put(errorBand);
    errorBand->Draw("f same");

    // draw the fit function on top
    TGraph* graphFitFunction = GetFitFunctionGraph();
    fObjectTrash.Put(graphFitFunction);
    graphFitFunction->SetLineColor(kRed);
    graphFitFunction->SetLineWidth(2);
    graphFitFunction->Draw("l same");

    // now draw the histogram again since it was covered by the band and the best fit
    fGraph.Draw("p same");

    // draw legend
    if (flaglegend) {
        TLegend* legend = new TLegend(0.25, 0.75, 0.55, 0.9);
        fObjectTrash.Put(legend);
        legend->SetBorderSize(0);
        legend->SetFillColor(kWhite);
        legend->AddEntry(&fGraph, "Data", "PE");
        legend->AddEntry(graphFitFunction, "Best fit", "L");
        legend->AddEntry(errorBand, "Error band", "F");
        legend->Draw();
    }

    gPad->RedrawAxis();
}

// ---------------------------------------------------------
double BCGraphFitter::CalculateChi2(const std::vector<double>& pars)
{
    TF1& f = GetFitFunction();

    // set pars into fit function
    f.SetParameters(&pars[0]);

    double chi, chi2 = 0;
    for (unsigned i = 0; i < GetDataSet()->GetNDataPoints(); ++i) {
        chi = (GetDataSet()->GetDataPoint(i)[1] - f.Eval(GetDataSet()->GetDataPoint(i)[0])) / GetDataSet()->GetDataPoint(i)[3];
        chi2 += chi * chi;
    }

    return chi2;
}

// ---------------------------------------------------------
double BCGraphFitter::CalculatePValue(const std::vector<double>& pars, bool ndf)
{
    const double chi2 = CalculateChi2(pars);

    if (chi2 < 0) {
        BCLOG_ERROR("chi2 is negative.");
        fPValue = -1;
    }

    else if (ndf) {
        if (GetNDoF() <= 0) {
            BCLOG_ERROR("number of degrees of freedom is not positive.");
            fPValue = -1;
        }
        fPValue = TMath::Prob(chi2, GetNDoF());

    } else if (GetNDataPoints() == 0) {
        BCLog::OutError("BCGraphFitter::CalculatePValue : number of data points is zero.");
        fPValue = -1;

    } else
        fPValue = TMath::Prob(chi2, GetNDataPoints());

    return fPValue;
}
