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
BCGraphFitter::BCGraphFitter(std::string name)
    : BCFitter(name),
      fGraph(0)
{
    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
BCGraphFitter::BCGraphFitter(TGraphErrors* graph, TF1* func, std::string name)
    : BCFitter(func, name),
      fGraph(0)
{
    SetGraph(graph);

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
int BCGraphFitter::SetGraph(TGraphErrors* graph)
{
    if (!graph) {
        BCLog::Out(BCLog::error, BCLog::error, "BCGraphFitter::SetGraph() : TGraphErrors not created.");
        return 0;
    }

    if (graph->GetN() <= 0) {
        BCLog::Out(BCLog::error, BCLog::error, "BCGraphFitter::SetGraph() : TGraphErrors is empty.");
        return 0;
    }
    if (graph->GetN() == 1) {
        BCLog::Out(BCLog::error, BCLog::error, "BCGraphFitter::SetGraph() : TGraphErrors has only one point. Not able to fit.");
        return 0;
    }

    fGraph = graph;

    double* x  = fGraph->GetX();
    double* y  = fGraph->GetY();
    double* ex = fGraph->GetEX();
    double* ey = fGraph->GetEY();

    if (!ey) {
        BCLog::Out(BCLog::error, BCLog::error, "BCGraphFitter::SetGraph() : TGraphErrors has NO errors set on Y. Not able to fit.");
        return 0;
    }

    SetDataSet(new BCDataSet());

    // fill the dataset
    for (int i = 0; i < fGraph->GetN(); ++i) {
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

    return GetNDataPoints();
}

// ---------------------------------------------------------
BCGraphFitter::~BCGraphFitter()
{
    // todo memory leak
}

// ---------------------------------------------------------
double BCGraphFitter::LogLikelihood(const std::vector<double>& params)
{
    unsigned c = MCMCGetCurrentChain();

    if (!fFitFunction.at(c) or !GetDataSet() or GetDataSet()->GetNDataPoints() == 0)
        return std::numeric_limits<double>::quiet_NaN();

    // set the parameters of the function
    // passing the pointer to first element of the vector is
    // not completely safe as there might be an implementation where
    // the vector elements are not stored consecutively in memory.
    // however it is much faster than copying the contents, especially
    // for large number of parameters
    fFitFunction[c]->SetParameters(&params[0]);

    // initialize probability
    double logl = 0.;

    // loop over all data points
    for (unsigned i = 0; i < GetNDataPoints(); i++)
        // calculate log of probability assuming
        // a Gaussian distribution for each point
        logl += BCMath::LogGaus(GetDataSet()->GetDataPoint(i)[1], // y value of point
                                fFitFunction[c]->Eval(GetDataSet()->GetDataPoint(i)[0]), // f(x value of point)
                                GetDataSet()->GetDataPoint(i)[3], // uncertainty on y value of point
                                true); // include normalization factor

    return logl;
}

// ---------------------------------------------------------
bool BCGraphFitter::Fit(TGraphErrors* graph, TF1* func)
{
    // set graph
    if (!SetGraph(graph)) {
        BCLog::OutError("BCEfficiencyFitter::Fit : Graph not defined.");
        return 0;
    }

    // set function
    if (!SetFitFunction(func)) {
        BCLog::OutError("BCEfficiencyFitter::Fit : Fit function not defined.");
        return 0;
    }

    return Fit();
}

// ---------------------------------------------------------
bool BCGraphFitter::Fit()
{
    // set graph
    if (!fGraph) {
        BCLog::OutError("BCEfficiencyFitter::Fit : Graph not defined.");
        return 0;
    }

    // set function
    if (fFitFunction.empty()) {
        BCLog::OutError("BCEfficiencyFitter::Fit : Fit function not defined.");
        return 0;
    }

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
    FindMode(GetGlobalMode());
    SetOptimizationMethod(method_temp);

    // p value with (approximate) correction for degrees of freedom
    CalculatePValue(GetGlobalMode(), true);

    // print summary to screen
    PrintShortFitSummary();

    return 1;
}

// ---------------------------------------------------------
void BCGraphFitter::DrawFit(const char* options, bool flaglegend)
{
    if (!fGraph) {
        BCLog::OutError("BCGraphFitter::DrawFit() : TGraphErrors not defined.");
        return;
    }

    if (fFitFunction.empty()) {
        BCLog::OutError("BCGraphFitter::DrawFit() : Fit function not defined.");
        return;
    }

    // check wheather options contain "same"
    TString opt = options;
    opt.ToLower();

    // if not same, draw the histogram first to get the axes
    if (!opt.Contains("same"))
        fGraph->Draw("ap");

    // draw the error band as central 68% probability interval
    fErrorBand = GetErrorBandGraph(0.16, 0.84);
    fErrorBand->Draw("f same");

    // draw the fit function on top
    fGraphFitFunction = GetFitFunctionGraph();
    fGraphFitFunction->SetLineColor(kRed);
    fGraphFitFunction->SetLineWidth(2);
    fGraphFitFunction->Draw("l same");

    // now draw the histogram again since it was covered by the band and the best fit
    fGraph->Draw("p same");

    // draw legend
    if (flaglegend) {
        TLegend* legend = new TLegend(0.25, 0.75, 0.55, 0.9);
        legend->SetBorderSize(0);
        legend->SetFillColor(kWhite);
        legend->AddEntry(fGraph, "Data", "PE");
        legend->AddEntry(fGraphFitFunction, "Best fit", "L");
        legend->AddEntry(fErrorBand, "Error band", "F");
        legend->Draw();
    }

    gPad->RedrawAxis();
}

// ---------------------------------------------------------
double BCGraphFitter::CDF(const std::vector<double>& parameters,  int index, bool /*lower*/)
{

    //format: x y error_x error_y
    std::vector<double> values = fDataSet->GetDataPoint(index).GetValues();

    if (values[2])
        BCLog::OutWarning("BCGraphFitter::CDF: Non-zero errors in x-direction are ignored!");

    // get the observed value
    double yObs = values[1];

    // expectation value
    double yExp = FitFunction(values, parameters);

    return ROOT::Math::normal_cdf(yObs, values[3], yExp);
}

// ---------------------------------------------------------
double BCGraphFitter::CalculateChi2(const std::vector<double>& pars)
{
    if (fFitFunction.empty() or !GetDataSet())
        return std::numeric_limits<double>::quiet_NaN();

    const unsigned c = MCMCGetCurrentChain();

    // set pars into fit function
    fFitFunction[c]->SetParameters(&pars[0]);

    double chi, chi2 = 0;
    for (unsigned i = 0; i < GetDataSet()->GetNDataPoints(); ++i) {
        chi = (GetDataSet()->GetDataPoint(i)[1] - fFitFunction[c]->Eval(GetDataSet()->GetDataPoint(i)[0])) / GetDataSet()->GetDataPoint(i)[3];
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
