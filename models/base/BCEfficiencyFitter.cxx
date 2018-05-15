/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include <config.h>

#include "BCEfficiencyFitter.h"

#include <BAT/BCDataPoint.h>
#include <BAT/BCH1D.h>
#include <BAT/BCLog.h>
#include <BAT/BCMath.h>

#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TString.h>

// ---------------------------------------------------------
BCEfficiencyFitter::BCEfficiencyFitter(const TH1& trials, const TH1& successes, const TF1& func, const std::string& name)
    : BCFitter(func, name),
      fHistogramBinomial("hist_binomial", "", 1000, 0., 1.),
      fDataPointType(kDataPointSmallestInterval)
{
    // check compatibility of both histograms : number of bins
    if (trials.GetNbinsX() != successes.GetNbinsX())
        throw std::invalid_argument("BCEfficiencyFitter: Histograms do not have the same number of bins.");

    // check compatibility of both histograms : bin edges incl. right-most edge = left edge of overflow bin
    for (int i = 1; i <= trials.GetNbinsX() + 1; ++i) {
        if (trials.GetBinLowEdge(i) != successes.GetBinLowEdge(i)) {
            BCLOG_WARNING("Histograms of trials and successes don't have same bin boundaries.");
            break;
        }
    }

    // check that there aren't more successes than trials
    for (int i = 1; i <= trials.GetNbinsX(); ++i) {
        if (trials.GetBinContent(i) < successes.GetBinContent(i)) {
            throw std::invalid_argument("Trials histogram has fewer entries than successes histogram.");
        }
    }

    trials.Copy(fTrials);
    successes.Copy(fSuccesses);

    // create data points and add them to the data set.
    for (int i = 0; i < fTrials.GetNbinsX(); ++i)
        fFitterDataSet.AddDataPoint(BCDataPoint(2));

    // set the data boundaries for x and y values.
    GetDataSet()->SetBounds(0, fTrials.GetXaxis()->GetXmin(), fTrials.GetXaxis()->GetXmax());
    GetDataSet()->SetBounds(1, 0.0, 1.0);

    // set the indices for fitting.
    SetFitFunctionIndices(0, 1);

    fFlagIntegration = false;

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
BCEfficiencyFitter::~BCEfficiencyFitter() {}

// ---------------------------------------------------------
double BCEfficiencyFitter::LogLikelihood(const std::vector<double>& params)
{
    TF1& f = GetFitFunction();

    // initialize probability
    double loglikelihood = 0;

    // set the parameters of the function
    // passing the pointer to first element of the vector is
    // not completely safe as there might be an implementation where
    // the vector elements are not stored consecutively in memory.
    // however it is much faster than copying the contents, especially
    // for large number of parameters
    f.SetParameters(&params[0]);

    // get the number of bins
    int nbins = fTrials.GetNbinsX();

    // loop over all bins
    for (int ibin = 1; ibin <= nbins; ++ibin) {
        // get n
        int n = int(fTrials.GetBinContent(ibin));

        // get k
        int k = int(fSuccesses.GetBinContent(ibin));

        // get bin edges
        const double xmin = fTrials.GetXaxis()->GetBinLowEdge(ibin);
        const double xmax = fTrials.GetXaxis()->GetBinUpEdge(ibin);

        const double eff = Integral(params, xmin, xmax);

        // get the value of the Binomial distribution
        loglikelihood += BCMath::LogApproxBinomial(n, k, eff);
    }

    return loglikelihood;
}

// ---------------------------------------------------------
void BCEfficiencyFitter::Fit()
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
    CalculatePValueFast(GetBestFitParameters());

    // print summary to screen
    PrintShortFitSummary();
}

// ---------------------------------------------------------
void BCEfficiencyFitter::DrawData(bool flaglegend)
{
    BCAux::RootSideEffectGuard g;

    // create new histogram
    TH2D* hist_axes = new TH2D("hist_axes",
                               (std::string(fTrials.GetTitle()) + ";" + fTrials.GetXaxis()->GetTitle() + ";" + fTrials.GetYaxis()->GetTitle()).data(),
                               100, fTrials.GetXaxis()->GetXmin(), fTrials.GetXaxis()->GetXmax(),
                               100, std::min(fTrials.GetMinimum(), fSuccesses.GetMinimum()), 1.1 * std::max(fTrials.GetMaximum(), fSuccesses.GetMaximum()));
    hist_axes->SetStats(false);

    fObjectTrash.Put(hist_axes);
    hist_axes->Draw();

    fTrials.SetLineWidth(1);
    fTrials.SetLineColor(kGray + 3);
    fTrials.SetFillColor(fTrials.GetLineColor());
    fTrials.SetFillStyle(1001);

    fSuccesses.SetMarkerColor(kRed);
    fSuccesses.SetMarkerStyle(8);
    fSuccesses.SetMarkerSize(0.7);

    fTrials.Draw("same");
    fSuccesses.Draw("Psame");

    // draw values when trial is nonzero but success is zero
    TMarker M;
    M.SetMarkerColor(fSuccesses.GetMarkerColor());
    M.SetMarkerStyle(fSuccesses.GetMarkerStyle());
    M.SetMarkerSize(fSuccesses.GetMarkerSize());
    for (int i = 1; i <= fTrials.GetNbinsX(); ++i)
        if (fTrials.GetBinContent(i) > 0 && fSuccesses.GetBinContent(i) == 0)
            M.DrawMarker(fTrials.GetBinCenter(i), 0);

    // draw legend
    if (flaglegend) {
        TLegend* legend = new TLegend(0.135, 0.955, 0.85, 0.98);
        fObjectTrash.Put(legend);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetFillStyle(0);
        legend->SetTextAlign(12);
        legend->SetTextSize(0.03);
        legend->SetNColumns(2);
        legend->SetColumnSeparation(5e-2);

        legend->AddEntry(&fTrials, fTrials.GetTitle(), "F");
        legend->AddEntry(&fSuccesses, fSuccesses.GetTitle(), "P");

        legend->Draw();
    }
}

// ---------------------------------------------------------
void BCEfficiencyFitter::DrawFit(const std::string& options, bool flaglegend)
{
    // create efficiency graph
    TGraphAsymmErrors* graphRatio = new TGraphAsymmErrors();
    fObjectTrash.Put(graphRatio);
    graphRatio->SetLineColor(kGray + 3);
    graphRatio->SetMarkerColor(graphRatio->GetLineColor());
    graphRatio->SetMarkerStyle(8);
    graphRatio->SetMarkerSize(0.7);

    // set points
    for (int i = 1; i <= fTrials.GetNbinsX(); ++i) {
        // calculate central value and uncertainties
        double xexp, xmin, xmax;
        if (fTrials.GetBinContent(i) > 0) {
            GetUncertainties( static_cast<int>(fTrials.GetBinContent(i)), static_cast<int>(fSuccesses.GetBinContent(i)), 0.68, xexp, xmin, xmax);
            graphRatio->SetPoint(graphRatio->GetN(), fTrials.GetBinCenter(i), xexp);
            graphRatio->SetPointError(graphRatio->GetN() - 1, 0, 0, xexp - xmin, xmax - xexp);
        }
    }

    // check wheather options contain "same"
    TString opt = options;
    opt.ToLower();

    // if not same, draw the histogram first to get the axes
    if (!opt.Contains("same")) {
        BCAux::RootSideEffectGuard g;

        // create new histogram
        TH2D* hist_axes = new TH2D("hist_axes",
                                   Form(";%s;ratio", fTrials.GetXaxis()->GetTitle()),
                                   100, fTrials.GetXaxis()->GetXmin(), fTrials.GetXaxis()->GetXmax(),
                                   100, 0., 1.05);

        fObjectTrash.Put(hist_axes);
        hist_axes->SetStats(kFALSE);
        hist_axes->Draw();

        graphRatio->Draw(Form("%sp", opt.Data()));
    }

    // draw the error band as central 68% probability interval
    TGraph* errorBand = GetErrorBandGraph(0.16, 0.84);
    fObjectTrash.Put(errorBand);
    errorBand->Draw("f same");

    // now draw the histogram again since it was covered by the band
    graphRatio->Draw(Form("%ssamepX", opt.Data()));

    // draw the fit function on top
    TGraph* graphFitFunction = GetFitFunctionGraph();
    fObjectTrash.Put(graphFitFunction);
    graphFitFunction->SetLineColor(kRed);
    graphFitFunction->SetLineWidth(2);
    graphFitFunction->Draw("l same");

    // draw legend
    if (flaglegend) {
        TLegend* legend = new TLegend(0.135, 0.955, 0.85, 0.98);
        fObjectTrash.Put(legend);
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetFillStyle(0);
        legend->SetTextAlign(12);
        legend->SetTextSize(0.03);
        legend->SetNColumns(3);
        legend->SetColumnSeparation(5e-2);

        legend->AddEntry(graphRatio, "Data", "PE");
        legend->AddEntry(graphFitFunction, "Best fit", "L");
        legend->AddEntry(errorBand, "Error band", "F");

        legend->Draw();

    }
    gPad->RedrawAxis();
}

// ---------------------------------------------------------
double BCEfficiencyFitter::CalculatePValueFast(const std::vector<double>& par, unsigned nIterations)
{
    //use NULL pointer for no callback
    return CalculatePValueFast(par, NULL, nIterations);
}

// ---------------------------------------------------------
double BCEfficiencyFitter::CalculatePValueFast(const std::vector<double>& pars, BCEfficiencyFitter::ToyDataInterface* callback, unsigned nIterations)
{
    // define temporary variables
    int nbins = fTrials.GetNbinsX();

    std::vector<unsigned> histogram(nbins, 0);
    std::vector<double> expectation(nbins, 0);

    double logp = 0;
    double logp_start = 0;
    int counter_pvalue = 0;

    // define starting distribution
    for (int ibin = 0; ibin < nbins; ++ibin) {
        // get bin boundaries
        const double xmin = fTrials.GetXaxis()->GetBinLowEdge(ibin + 1);
        const double xmax = fTrials.GetXaxis()->GetBinUpEdge(ibin + 1);

        // get the number of expected events
        double yexp = Integral(pars, xmin, xmax);

        // get n
        int n = int(fTrials.GetBinContent(ibin));

        // get k
        int k = int(fSuccesses.GetBinContent(ibin));

        histogram[ibin]   = k;
        expectation[ibin] = n * yexp;

        // calculate p;
        logp += BCMath::LogApproxBinomial(n, k, yexp);
    }
    logp_start = logp;

    // loop over iterations
    for (unsigned iiter = 0; iiter < nIterations; ++iiter) {
        // loop over bins
        for (int ibin = 0; ibin < nbins; ++ibin) {
            // get n
            int n = int(fTrials.GetBinContent(ibin));

            // get k
            int k = histogram.at(ibin);

            // get expectation
            double yexp = 0;
            if (n > 0)
                yexp = expectation.at(ibin) / n;

            // random step up or down in statistics for this bin
            double ptest = fRandom.Rndm() - 0.5;

            // continue, if efficiency is at limit
            if (!(yexp > 0. || yexp < 1.0))
                continue;

            // increase statistics by 1
            if (ptest > 0 && (k < n)) {
                // calculate factor of probability
                double r = (double(n) - double(k)) / (double(k) + 1.) * yexp / (1. - yexp);

                // walk, or don't (this is the Metropolis part)
                if (fRandom.Rndm() < r) {
                    histogram[ibin] = k + 1;
                    logp += log(r);
                }
            }

            // decrease statistics by 1
            else if (ptest <= 0 && (k > 0)) {
                // calculate factor of probability
                double r = double(k) / (double(n) - (double(k) - 1)) * (1 - yexp) / yexp;

                // walk, or don't (this is the Metropolis part)
                if (fRandom.Rndm() < r) {
                    histogram[ibin] = k - 1;
                    logp += log(r);
                }
            }

        } // end of looping over bins

        //one new toy data set is created
        //call user interface to calculate arbitrary statistic's distribution
        if (callback)
            (*callback)(expectation, histogram);

        // increase counter
        if (logp < logp_start)
            counter_pvalue++;

    }

    // pvalue estimated as simple fraction
    fPValue = double(counter_pvalue) / double(nIterations);
    return fPValue;
}

// ---------------------------------------------------------
bool BCEfficiencyFitter::GetUncertainties(int n, int k, double p, double& xexp, double& xmin, double& xmax)
{
    if (n == 0) {
        xexp = 0.;
        xmin = 0.;
        xmax = 0.;
        return false;
    }

    BCLog::OutDebug(Form("Calculating efficiency data-point of type %d for (n,k) = (%d,%d)", fDataPointType, n, k));

    // create a clean histogram with binomial distribution
    fHistogramBinomial.Reset();

    // loop over bins and fill histogram
    for (int i = 1; i <= fHistogramBinomial.GetNbinsX(); ++i)
        fHistogramBinomial.SetBinContent(i, BCMath::ApproxBinomial(n, k, fHistogramBinomial.GetBinCenter(i)));
    // normalize
    fHistogramBinomial.Scale(1. / fHistogramBinomial.Integral());

    switch (fDataPointType) {

        case kDataPointRMS: {
            xexp = fHistogramBinomial.GetMean();
            double rms = fHistogramBinomial.GetRMS();
            xmin = xexp - rms;
            xmax = xexp + rms;
            BCLog::OutDebug(Form(" - mean = %f , rms = %f)", xexp, rms));
            break;
        }

        case kDataPointSmallestInterval: {
            xexp = (double)k / (double)n;
            BCH1D fbh(&fHistogramBinomial);
            BCH1D::BCH1DSmallestInterval SI = fbh.GetSmallestIntervals(p);
            if (SI.intervals.empty()) {
                xmin = xmax = xexp = 0.;
            } else {
                xmin = SI.intervals.front().xmin;
                xmax = SI.intervals.front().xmax;
            }

            break;
        }

        case kDataPointCentralInterval: {
            // calculate quantiles
            int nprobSum = 3;
            double q[3];
            double probSum[3];
            probSum[0] = (1. - p) / 2.;
            probSum[1] = 0.5;
            probSum[2] = (1. + p) / 2.;

            fHistogramBinomial.GetQuantiles(nprobSum, q, probSum);

            xmin = q[0];
            xexp = q[1];
            xmax = q[2];
            break;
        }

        default: {
            BCLog::OutError("BCEfficiencyFitter::GetUncertainties - invalid DataPointType specified.");
            xmin = xmax = xexp = 0.;
        }
    }
    BCLog::OutDebug(Form(" - efficiency = %f , range (%f - %f)", xexp, xmin, xmax));

    return true;
}
