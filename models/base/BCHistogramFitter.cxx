/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <config.h>

#include "BCHistogramFitter.h"

#include <BAT/BCDataPoint.h>
#include <BAT/BCLog.h>
#include <BAT/BCMath.h>

#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <Math/ProbFuncMathCore.h> //for ROOT::Math namespace
#include <TPad.h>
#include <TRandom3.h>
#include <TString.h>

// ---------------------------------------------------------
BCHistogramFitter::BCHistogramFitter(const TH1& hist, const TF1& func, const std::string& name) :
    BCFitter(func, name)
{
    if (hist.GetDimension() != 1)
        throw std::invalid_argument("Only 1D histograms supported");

    hist.Copy(fHistogram);

    // create data points and add them to the data set.
    // the x value is the lower edge of the bin, and
    // the y value is the bin count
    int nbins = hist.GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        GetDataSet()->AddDataPoint(BCDataPoint(2));
        GetDataSet()->Back()[0] = hist.GetBinLowEdge(i);
        GetDataSet()->Back()[1] = hist.GetBinContent(i);
    }

    // set the data boundaries for x and y values
    GetDataSet()->SetBounds(0, hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
    GetDataSet()->SetBounds(1, std::max<double>(hist.GetMinimum() - sqrt(hist.GetMinimum()) / 2, 0),
                            hist.GetMaximum() + sqrt(hist.GetMaximum()) / 2);

    // set the indices for fitting
    SetFitFunctionIndices(0, 1);

    SetNIterationsRun(2000);
    SetFillErrorBand(true);
    fFlagIntegration = true;

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
double BCHistogramFitter::LogLikelihood(const std::vector<double>& params)
{
    // initialize probability
    double loglikelihood = 0;

    // loop over all bins
    for (int ibin = 1; ibin <= fHistogram.GetNbinsX(); ++ibin) {
        // get bin edges and integrate for expectation
        const double xmin = fHistogram.GetXaxis()->GetBinLowEdge(ibin);
        const double xmax = fHistogram.GetXaxis()->GetBinUpEdge(ibin);
        double yexp = Integral(params, xmin, xmax);

        // get the number of observed events
        double y = fHistogram.GetBinContent(ibin);

        // get the value of the Poisson distribution
        loglikelihood += BCMath::LogPoisson(y, yexp);
    }

    return loglikelihood;
}

// ---------------------------------------------------------
void BCHistogramFitter::Fit()
{
    // perform marginalization
    MarginalizeAll();

    // maximize posterior probability, using the best-fit values close
    // to the global maximum from the MCMC
    BCIntegrate::BCOptimizationMethod method_temp = GetOptimizationMethod();
    SetOptimizationMethod(BCIntegrate::kOptMinuit);
    FindMode(GetBestFitParameters());
    SetOptimizationMethod(method_temp);

    // calculate the p-value using the fast MCMC algorithm
    if (!CalculatePValueFast(GetBestFitParameters()))
        BCLOG_WARNING("Could not use the fast p-value evaluation.");

    // print summary to screen
    PrintShortFitSummary();
}

// ---------------------------------------------------------
void BCHistogramFitter::DrawFit(const std::string& options, bool flaglegend)
{
    if (GetBestFitParameters().empty()) {
        BCLOG_ERROR("Fit not performed yet.");
        return;
    }

    // check whether options contain "same"
    TString opt = options;
    opt.ToLower();

    // if not same, draw the histogram first to get the axes
    if (!opt.Contains("same"))
        fHistogram.Draw(Form("hist%s", opt.Data()));

    // draw the error band as central 68% probability interval
    TGraph* errorBand = GetErrorBandGraph(0.16, 0.84);
    fObjectTrash.Put(errorBand);
    errorBand->Draw("f same");

    // now draw the histogram again since it was covered by the band
    fHistogram.Draw(Form("hist%ssame", opt.Data()));

    // draw the fit function on top
    TGraph* graphFitFunction = GetFitFunctionGraph();
    fObjectTrash.Put(graphFitFunction);
    graphFitFunction->SetLineColor(kRed);
    graphFitFunction->SetLineWidth(2);
    graphFitFunction->Draw("l same");

    // draw legend
    if (flaglegend) {
        TLegend* legend = new TLegend(0.25, 0.75, 0.55, 0.9);
        fObjectTrash.Put(legend);
        legend->SetBorderSize(0);
        legend->SetFillColor(kWhite);
        legend->AddEntry(&fHistogram, "Data", "L");
        legend->AddEntry(graphFitFunction, "Best fit", "L");
        legend->AddEntry(errorBand, "Error band", "F");
        legend->Draw();
    }

    gPad->RedrawAxis();
}

// ---------------------------------------------------------
double BCHistogramFitter::CalculatePValueFast(const std::vector<double>& pars, unsigned nIterations)
{
    // check size of parameter vector
    if (pars.size() != GetNParameters()) {
        BCLOG_ERROR("Number of parameters is inconsistent.");
        return false;
    }

    // copy observed and expected values
    std::vector<unsigned> observed(fHistogram.GetNbinsX());
    std::vector<double> expected(observed.size());

    for (int ibin = 0 ; ibin < fHistogram.GetNbinsX(); ++ibin) {
        observed[ibin] = (unsigned int) fHistogram.GetBinContent(ibin + 1);
        expected[ibin] = Integral(pars, fHistogram.GetBinLowEdge(ibin + 1), fHistogram.GetBinLowEdge(ibin + 2));
    }

    fPValue = BCMath::FastPValue(observed, expected, nIterations, fRandom.GetSeed());
    return fPValue;
}

// ---------------------------------------------------------
double BCHistogramFitter::CalculatePValueLikelihood(const std::vector<double>& pars)
{
    // initialize test statistic -2*lambda
    double logLambda = 0.0;

    for (int ibin = 1; ibin <= fHistogram.GetNbinsX(); ++ibin) {

        // get the number of observed events
        const double y = fHistogram.GetBinContent(ibin);

        // get the number of expected events
        const double yexp = Integral(pars, fHistogram.GetBinLowEdge(ibin), fHistogram.GetBinLowEdge(ibin + 1));

        // get the contribution from this datapoint
        if (y == 0)
            logLambda += yexp;
        else
            logLambda += yexp - y + y * log(y / yexp);
    }

    // rescale
    logLambda *= 2.0;

    //p value from chi^2 distribution, returns zero if logLambda < 0
    fPValue = TMath::Prob(logLambda, GetNDataPoints());
    return fPValue;
}

// ---------------------------------------------------------
double BCHistogramFitter::CalculatePValueLeastSquares(const std::vector<double>& pars, bool weightExpect)
{
    // initialize test statistic chi^2
    double chi2 = 0.0;

    for (int ibin = 1; ibin <= fHistogram.GetNbinsX(); ++ibin) {

        // get the number of observed events
        double y = fHistogram.GetBinContent(ibin);

        // get the number of expected events
        const double yexp = Integral(pars, fHistogram.GetBinLowEdge(ibin), fHistogram.GetBinLowEdge(ibin + 1));

        //convert 1/0.0 into 1 for weighted sum
        double weight;
        if (weightExpect)
            weight = (yexp > 0) ? yexp : 1.0;
        else
            weight = (y > 0) ? y : 1.0;

        // get the contribution from this datapoint
        chi2 += (y - yexp) * (y - yexp) / weight;
    }

    // p value from chi^2 distribution
    fPValue = TMath::Prob(chi2, GetNDataPoints());
    return fPValue;
}

// ---------------------------------------------------------
double BCHistogramFitter::GraphCorrection(unsigned ibin) const
{
    return fHistogram.GetXaxis()->GetBinUpEdge(ibin) - fHistogram.GetXaxis()->GetBinLowEdge(ibin);
}
