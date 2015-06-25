/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCHistogramFitter.h"
#include "BAT/BCLog.h"
#include "BAT/BCDataSet.h"
#include "BAT/BCDataPoint.h"
#include "BAT/BCMath.h"

#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TMath.h>
#include <Math/ProbFuncMathCore.h> //for ROOT::Math namespace
#include <TString.h>

// ---------------------------------------------------------
BCHistogramFitter::BCHistogramFitter()
    : BCFitter()
    , fHistogram(0)
    , fFitFunction(0)
    , fHistogramExpected(0)
{
    // set default options and values
    MCMCSetNIterationsRun(2000);
    SetFillErrorBand();
    fFlagIntegration = true;

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
BCHistogramFitter::BCHistogramFitter(std::string name)
    : BCFitter(name)
    , fHistogram(0)
    , fFitFunction(0)
    , fHistogramExpected(0)
{
    // set default options and values
    MCMCSetNIterationsRun(2000);
    SetFillErrorBand();
    fFlagIntegration = true;

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
BCHistogramFitter::BCHistogramFitter(TH1D* hist, TF1* func)
    : BCFitter()
    , fHistogram(0)
    , fFitFunction(0)
    , fHistogramExpected(0)
{
    SetHistogram(hist);
    SetFitFunction(func);

    MCMCSetNIterationsRun(2000);
    SetFillErrorBand();
    fFlagIntegration = true;

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
BCHistogramFitter::BCHistogramFitter(std::string name, TH1D* hist, TF1* func)
    : BCFitter(name)
    , fHistogram(0)
    , fFitFunction(0)
    , fHistogramExpected(0)
{
    SetHistogram(hist);
    SetFitFunction(func);

    MCMCSetNIterationsRun(2000);
    SetFillErrorBand();
    fFlagIntegration = true;

    // set MCMC for marginalization
    SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------
bool BCHistogramFitter::SetHistogram(TH1D* hist)
{
    // check if histogram exists
    if (!hist) {
        BCLog::OutError("BCHistogramFitter::SetHistogram : TH1D not created.");
        return false;
    }

    // set pointer to histogram
    fHistogram = hist;

    // create a data set. this is necessary in order to calculate the
    // error band. the data set contains as many data points as there
    // are bins.
    SetDataSet(new BCDataSet());

    // create data points and add them to the data set.
    // the x value is the lower edge of the bin, and
    // the y value is the bin count
    int nbins = fHistogram->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        GetDataSet()->AddDataPoint(BCDataPoint(2));
        GetDataSet()->Back()[0] = fHistogram->GetBinLowEdge(i);
        GetDataSet()->Back()[1] = fHistogram->GetBinContent(i);
    }

    // set the data boundaries for x and y values.
    GetDataSet()->SetBounds(0, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    GetDataSet()->SetBounds(1, std::max<double>(hist->GetMinimum() - sqrt(hist->GetMinimum()) / 2, 0),
                            hist->GetMaximum() + sqrt(hist->GetMaximum()) / 2);

    // set the indeces for fitting.
    SetFitFunctionIndices(0, 1);

    // no error
    return false;
}

// ---------------------------------------------------------
bool BCHistogramFitter::SetHistogramExpected(const std::vector<double>& parameters)
{
    // clear up previous memory
    if (fHistogramExpected)
        delete fHistogramExpected;
    // copy all properties from the data histogram
    fHistogramExpected = new TH1D(*fHistogram);

    // get the number of bins
    int nBins = fHistogramExpected->GetNbinsX();

    //set the parameters of fit function
    fFitFunction->SetParameters(&parameters[0]);

    // loop over all bins, skip underflow
    for (int ibin = 1; ibin <= nBins; ++ibin) {
        // get bin edges
        double xmin = fHistogramExpected->GetXaxis()->GetBinLowEdge(ibin);
        double xmax = fHistogramExpected->GetXaxis()->GetBinUpEdge(ibin);

        //expected count
        double yexp = 0.;

        if (fFlagIntegration)				// use ROOT's TH1D::Integral method
            yexp = fFitFunction->Integral(xmin, xmax);
        else												// use linear interpolation
            yexp = (fFitFunction->Eval(xmin) + fFitFunction->Eval(xmax)) * (xmax - xmin) / 2;

        //write the expectation for the bin
        fHistogramExpected->SetBinContent(ibin, yexp);

        //avoid automatic error as sqrt(yexp), used e.g. in Kolmogorov correction factor
        fHistogramExpected->SetBinError(ibin, 0.0);

        // but the data under this model have that sqrt(yexp) uncertainty
        fHistogram->SetBinError(ibin, sqrt(yexp));

    }
    return true;
}

// ---------------------------------------------------------
bool BCHistogramFitter::SetFitFunction(TF1* func)
{
    // check if function exists
    if (!func) {
        BCLog::OutError("BCHistogramFitter::SetFitFunction : TF1 not created.");
        return false;
    }

    // set the function
    fFitFunction = func;

    // update the model name to contain the function name
    if (fName == "model")
        SetName(std::string("HistogramFitter with ") + fFitFunction->GetName());

    // reset parameters
    fParameters = BCParameterSet();

    // get the new number of parameters
    int n = func->GetNpar();

    // add parameters
    for (int i = 0; i < n; ++i) {
        double xmin;
        double xmax;

        func->GetParLimits(i, xmin, xmax);

        AddParameter(func->GetParName(i), xmin, xmax);
    }

    // set flat prior for all parameters by default
    SetPriorConstantAll();

    return true;
}

// ---------------------------------------------------------
BCHistogramFitter::~BCHistogramFitter()
{
    // todo memory leak, many members not removed
    delete fHistogramExpected;
}

// ---------------------------------------------------------
double BCHistogramFitter::LogLikelihood(const std::vector<double>& params)
{
    // initialize probability
    double loglikelihood = 0;

    // set the parameters of the function
    fFitFunction->SetParameters(&params[0]);

    // get the number of bins
    int nbins = fHistogram->GetNbinsX();

    // loop over all bins
    for (int ibin = 1; ibin <= nbins; ++ibin) {
        // get bin edges
        double xmin = fHistogram->GetXaxis()->GetBinLowEdge(ibin);
        double xmax = fHistogram->GetXaxis()->GetBinUpEdge(ibin);

        // get the number of observed events
        double y = fHistogram->GetBinContent(ibin);

        double yexp = 0.;

        if (fFlagIntegration)				// use ROOT's TH1D::Integral method
            yexp = fFitFunction->Integral(xmin, xmax);
        else												// use linear interpolation
            yexp = (fFitFunction->Eval(xmin) + fFitFunction->Eval(xmax)) * (xmax - xmin) / 2.;

        // get the value of the Poisson distribution
        loglikelihood += BCMath::LogPoisson(y, yexp);
    }

    return loglikelihood;
}

// ---------------------------------------------------------
double BCHistogramFitter::FitFunction(const std::vector<double>& x, const std::vector<double>& params)
{
    // set the parameters of the function
    fFitFunction->SetParameters(&params[0]);

    return fFitFunction->Eval(x[0]) * fHistogram->GetBinWidth(fHistogram->FindBin(x[0]));
}

// ---------------------------------------------------------
bool BCHistogramFitter::Fit(TH1D* hist, TF1* func)
{
    // set histogram
    if (hist)
        SetHistogram(hist);
    else {
        BCLog::OutError("BCHistogramFitter::Fit : Histogram not defined.");
        return false;
    }

    // set function
    if (func)
        SetFitFunction(func);
    else {
        BCLog::OutError("BCHistogramFitter::Fit : Fit function not defined.");
        return false;
    }

    return Fit();
}

// ---------------------------------------------------------
bool BCHistogramFitter::Fit()
{
    if (!fHistogram) {
        BCLog::OutError("BCHistogramFitter::Fit : Histogram not defined.");
        return false;
    }

    if (!fFitFunction) {
        BCLog::OutError("BCHistogramFitter::Fit : Fit function not defined.");
        return false;
    }

    // perform marginalization
    MarginalizeAll();

    // maximize posterior probability, using the best-fit values close
    // to the global maximum from the MCMC
    BCIntegrate::BCOptimizationMethod method_temp = GetOptimizationMethod();
    SetOptimizationMethod(BCIntegrate::kOptMinuit);
    FindMode(GetGlobalMode());
    SetOptimizationMethod(method_temp);

    // calculate the p-value using the fast MCMC algorithm
    if ( !CalculatePValueFast(GetGlobalMode()))
        BCLog::OutWarning("BCHistogramFitter::Fit : Could not use the fast p-value evaluation.");

    // print summary to screen
    PrintShortFitSummary();

    // no error
    return true;
}

// ---------------------------------------------------------
void BCHistogramFitter::DrawFit(const char* options, bool flaglegend)
{
    if (!fHistogram) {
        BCLog::OutError("BCHistogramFitter::DrawFit : Histogram not defined.");
        return;
    }

    if (!fFitFunction) {
        BCLog::OutError("BCHistogramFitter::DrawFit : Fit function not defined.");
        return;
    }

    if (!fErrorBandXY or GetGlobalMode().empty()) {
        BCLog::OutError("BCHistogramFitter::DrawFit : Fit not performed yet.");
        return;
    }

    // check wheather options contain "same"
    TString opt = options;
    opt.ToLower();

    // if not same, draw the histogram first to get the axes
    if (!opt.Contains("same"))
        fHistogram->Draw(opt.Data());

    // draw the error band as central 68% probability interval
    fErrorBand = GetErrorBandGraph(0.16, 0.84);
    fErrorBand->Draw("f same");

    // now draw the histogram again since it was covered by the band
    fHistogram->Draw(Form("%ssame", opt.Data()));

    // draw the fit function on top
    fGraphFitFunction = GetFitFunctionGraph();
    fGraphFitFunction->SetLineColor(kRed);
    fGraphFitFunction->SetLineWidth(2);
    fGraphFitFunction->Draw("l same");

    // draw legend
    if (flaglegend) {
        TLegend* legend = new TLegend(0.25, 0.75, 0.55, 0.95);
        legend->SetBorderSize(0);
        legend->SetFillColor(kWhite);
        legend->AddEntry(fHistogram, "Data", "L");
        legend->AddEntry(fGraphFitFunction, "Best fit", "L");
        legend->AddEntry(fErrorBand, "Error band", "F");
        legend->Draw();
    }

    gPad->RedrawAxis();
}

// ---------------------------------------------------------
bool BCHistogramFitter::CalculatePValueFast(const std::vector<double>& par, unsigned nIterations)
{
    // check size of parameter vector
    if (par.size() < GetNParameters()) {
        BCLog::OutError("BCHistogramFitter::CalculatePValueFast : Number of parameters is inconsistent.");
        return false;
    }

    // check if histogram exists
    if (!fHistogram) {
        BCLog::OutError("BCHistogramFitter::CalculatePValueFast : Histogram not defined.");
        return false;
    }

    //update the expected number of events for all bins
    SetHistogramExpected(par);

    // copy observed and expected values
    std::vector<unsigned> observed(fHistogram->GetNbinsX());
    std::vector<double> expected(fHistogramExpected->GetNbinsX());

    for (int ibin = 0 ; ibin < fHistogram->GetNbinsX(); ++ibin) {
        observed[ibin] = (unsigned int) fHistogram->GetBinContent(ibin + 1);
        expected[ibin] = fHistogramExpected->GetBinContent(ibin + 1);
    }

    // create pseudo experiments
    fPValue = BCMath::FastPValue(observed, expected, nIterations, fRandom.GetSeed());

    // correct for fitting bias
    fPValueNDoF = BCMath::CorrectPValue(fPValue, GetNParameters(), fHistogram->GetNbinsX());

    // no error
    return true;
}

// ---------------------------------------------------------
bool BCHistogramFitter::CalculatePValueLikelihood(const std::vector<double>& par)
{
    // initialize test statistic -2*lambda
    double logLambda = 0.0;

    //Calculate expected counts to compare with
    SetHistogramExpected(par);

    for (int ibin = 1; ibin <= fHistogram->GetNbinsX(); ++ibin) {

        // get the number of observed events
        double y = fHistogram->GetBinContent(ibin);

        // get the number of expected events
        double yexp = fHistogramExpected->GetBinContent(ibin);

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
    fPValueNDoF = TMath::Prob(logLambda, GetNDoF());

    // no error
    return true;
}

// ---------------------------------------------------------
bool BCHistogramFitter::CalculatePValueLeastSquares(const std::vector<double>& par, bool weightExpect)
{
    // initialize test statistic chi^2
    double chi2 = 0.0;

    //Calculate expected counts to compare with
    SetHistogramExpected(par);

    for (int ibin = 1; ibin <= fHistogram->GetNbinsX(); ++ibin) {

        // get the number of observed events
        double y = fHistogram->GetBinContent(ibin);

        // get the number of expected events
        double yexp = fHistogramExpected->GetBinContent(ibin);

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
    fPValueNDoF = TMath::Prob(chi2, GetNDataPoints() - GetNParameters());

    // no error
    return true;
}

// ---------------------------------------------------------
bool BCHistogramFitter::CalculatePValueKolmogorov(const std::vector<double>& par)
{
    if (!fHistogramExpected or !fHistogram) {
        BCLog::OutError("BCHistogramFitter::CalculatePValueKolmogorov: "
                        "Please define the reference distribution by calling \n"
                        "BCHistogramFitter::SetHistogramExpected() first!");
        return false;
    }

    //update expected counts with current parameters
    SetHistogramExpected(par);

    //perform the test
    fPValue = fHistogramExpected->KolmogorovTest(fHistogram);
    fPValue = BCMath::CorrectPValue(fPValue, GetNParameters(), GetNDataPoints());

    // no error
    return true;
}

// ---------------------------------------------------------
double BCHistogramFitter::CDF(const std::vector<double>& parameters, int index, bool lower)
{

    if (parameters.empty())
        BCLog::OutWarning("BCHistogramFitter::CDF: parameter vector empty!");
    //histogram stores underflow in bin 0, so datapoint 0 is in bin 1
    index++;

    // get the number of observed events (should be integer)
    double yObs = fHistogram->GetBinContent(index);

    // get function value of lower bin edge
    double xmin = fHistogram->GetXaxis()->GetBinLowEdge(index);
    double xmax = fHistogram->GetXaxis()->GetBinUpEdge(index);

    // expectation value of this bin
    double yExp = 0.0;

    // use ROOT's TH1D::Integral method
    if (fFlagIntegration) {
        TF1 tmpF(*fFitFunction);
        tmpF.SetParameters(&parameters[0]);
        yExp = tmpF.Integral(xmin, xmax);
    } else													// use linear interpolation
        yExp = (fFitFunction->Eval(xmin) + fFitFunction->Eval(xmax)) * (xmax - xmin) / 2;

    // usually Poisson bins do not agree with fixed probability bins
    // so choose where it should belong

    if (lower) {
        if ((int) yObs >= 1)
            return ROOT::Math::poisson_cdf((unsigned int) (yObs - 1), yExp);
        else
            return 0.0;
    }
    // what if yObs as double doesn't represent a whole number? exception?
    if ((double) (unsigned int) yObs != yObs) {
        BCLog::OutWarning(Form(
                              "BCHistogramFitter::CDF: Histogram values should be integer!\n"
                              " Bin %d = %f", index, yObs));


        //convert randomly to integer
        // ex. yObs = 9.785 =>
        //      yObs --> 10 with P = 78.5%
        //      yObs --> 9  with P = 21.5%
        double U = fRandom.Rndm();
        double yObsLower = (unsigned int) yObs;
        if (U > (yObs - yObsLower))
            yObs = yObsLower;
        else
            yObs = yObsLower + 1;
    }

    return ROOT::Math::poisson_cdf((unsigned int) yObs, yExp);
}

// ---------------------------------------------------------
