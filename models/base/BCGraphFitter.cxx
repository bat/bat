/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

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
    : BCFitter(name)
    , fGraph(0)
    , fFitFunction(0)
    , fErrorBand(0)
    , fGraphFitFunction(0)
{
    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
BCGraphFitter::BCGraphFitter(TGraphErrors* graph, TF1* func, std::string name)
    : BCFitter(name)
    , fGraph(0)
    , fFitFunction(0)
    , fErrorBand(0)
    , fGraphFitFunction(0)
{
    SetGraph(graph);
    SetFitFunction(func);

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
int BCGraphFitter::SetFitFunction(TF1* func)
{
    if (!func) {
        BCLog::Out(BCLog::error, BCLog::error, "BCGraphFitter::SetFitFunction() : TF1 not created.");
        return 0;
    }

    // get the new number of parameters
    if (func->GetNpar() <= 0) {
        BCLog::Out(BCLog::error, BCLog::error, "BCGraphFitter::SetFitFunction() : TF1 has zero parameters. Not able to fit.");
        return 0;
    }

    // set the function
    fFitFunction = func;

    // update the model name to contain the function name
    if (fName == "model")
        SetName(std::string("GraphFitter with ") + fFitFunction->GetName());

    // reset parameters
    fParameters = BCParameterSet();

    // add parameters
    for (int i = 0; i < fFitFunction->GetNpar(); ++i) {
        double xmin;
        double xmax;
        fFitFunction->GetParLimits(i, xmin, xmax);
        AddParameter(fFitFunction->GetParName(i), xmin, xmax);
    }

    // set flat prior for all parameters by default
    SetPriorConstantAll();

    return GetNParameters();
}

// ---------------------------------------------------------
BCGraphFitter::~BCGraphFitter()
{
    // todo memory leak
}

// ---------------------------------------------------------
double BCGraphFitter::LogLikelihood(const std::vector<double>& params)
{
    if (!fFitFunction or !GetDataSet() or GetDataSet()->GetNDataPoints() == 0)
        return -std::numeric_limits<double>::infinity();

    // set the parameters of the function
    // passing the pointer to first element of the vector is
    // not completely safe as there might be an implementation where
    // the vector elements are not stored consecutively in memory.
    // however it is much faster than copying the contents, especially
    // for large number of parameters
    fFitFunction->SetParameters(&params[0]);

    // initialize probability
    double logl = 0.;

    // loop over all data points
    for (unsigned i = 0; i < GetNDataPoints(); i++)
        // calculate log of probability assuming
        // a Gaussian distribution for each point
        logl += BCMath::LogGaus(GetDataSet()->GetDataPoint(i)[1], // y value of point
                                fFitFunction->Eval(GetDataSet()->GetDataPoint(i)[0]), // f(x value of point)
                                GetDataSet()->GetDataPoint(i)[3], // uncertainty on y value of point
                                true); // include normalization factor

    return logl;
}

// ---------------------------------------------------------

double BCGraphFitter::FitFunction(const std::vector<double>& x, const std::vector<double>& params)
{
    // set the parameters of the function
    // passing the pointer to first element of the vector is
    // not completely safe as there might be an implementation where
    // the vector elements are not stored consecutively in memory.
    // however it is much faster than copying the contents, especially
    // for large number of parameters
    fFitFunction->SetParameters(&params[0]);
    return fFitFunction->Eval(x[0]);
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
    if (!fFitFunction) {
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
    FindMode( GetGlobalMode());
    SetOptimizationMethod(method_temp);

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

    if (!fFitFunction) {
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
double BCGraphFitter::GetChi2(const std::vector<double>& pars)
{
    if (!fFitFunction)
        return -std::numeric_limits<double>::infinity();
    if (!GetDataSet())
        return 0;

    // set pars into fit function
    fFitFunction->SetParameters(&pars[0]);

    double chi2 = 0;
    for (unsigned i = 0; i < GetDataSet()->GetNDataPoints(); ++i)
        chi2 += BCMath::LogGaus(GetDataSet()->GetDataPoint(i)[1], // y value of point
                                fFitFunction->Eval(GetDataSet()->GetDataPoint(i)[0]), // f(x value of point)
                                GetDataSet()->GetDataPoint(i)[3], // uncertainty on y value of point
                                false); // forego normalization factor
    return chi2;
}


// ---------------------------------------------------------
double BCGraphFitter::GetPValue(const std::vector<double>& pars, bool ndf)
{
    double chi2 = GetChi2(pars);

    if (chi2 < 0)
        fPValue = -1;

    else if (ndf) {
        if (GetNDoF() <= 0) {
            BCLog::OutError("BCGraphFitter::GetPValue : number of degrees of freedom is less than or equal to zero.");
            fPValue = -1;
        }
        fPValue = TMath::Prob(chi2, GetNDoF());

    } else if (GetNDataPoints() == 0) {
        BCLog::OutError("BCGraphFitter::GetPValue : number of data points is zero.");
        fPValue = -1;

    } else
        fPValue = TMath::Prob(chi2, GetNDataPoints());

    return fPValue;
}
