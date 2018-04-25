/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TCanvas.h>
#include <TF2.h>
#include <TH2D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <BAT/BCH2D.h>

#include <stdexcept>
#include <iostream>

#include "include/PerfTest2DFunction.h"

//______________________________________________________________________________
PerfTest2DFunction::PerfTest2DFunction(const std::string& name, TF2* func)
    : PerfTestMCMC(name)
    , fFunction(func)
{
    // set test type
    fTestType = PerfTest::kFunction2D;

    // manipulate function
    fFunction->SetNDF(100000);

    // get limits
    double xmin = fFunction->GetXmin();
    double xmax = fFunction->GetXmax();
    double ymin = fFunction->GetYmin();
    double ymax = fFunction->GetYmax();

    // add parameters
    AddParameter("x", xmin, xmax);
    AddParameter("y", ymin, ymax);

    GetParameter(0).SetNbins(25);
    GetParameter(1).SetNbins(25);

    DefineSubtests();
}

//______________________________________________________________________________
PerfTest2DFunction::~PerfTest2DFunction()
{
    if (fFunction)
        delete fFunction;
}

//______________________________________________________________________________
int PerfTest2DFunction::PostTest()
{
    PerfTestMCMC::PostTest();

    // get histogram
    TH2D* hist_marg = (TH2D*) GetMarginalized(0u, 1u).GetHistogram()->Clone();
    hist_marg->SetContour(20);
    hist_marg->Scale(hist_marg->GetEntries() / hist_marg->Integral());
    TH2D* hist_diff = (TH2D*) GetMarginalized(0u, 1u).GetHistogram()->Clone();
    hist_diff->SetContour(20);
    TH2D* hist_func = (TH2D*) GetMarginalized(0u, 1u).GetHistogram()->Clone();
    hist_func->SetContour(20);

    TH1D* hist_pull = new TH1D("", ";#Deltaf/sqrt(f);N", 50, -5.0, 5.0);

    // calculate ndf
    int nbinsx = hist_marg->GetNbinsX();
    int nbinsy = hist_marg->GetNbinsY();
    int ndf = nbinsx * nbinsy;

    // calculate norms
    double norm_hist = hist_marg->GetEntries();

    // get bin widths
    double binwidthx = hist_marg->GetXaxis()->GetBinWidth(1);
    double binwidthy = hist_marg->GetYaxis()->GetBinWidth(1);

    // fill histograms
    for (int i = 1; i <= nbinsx; ++i) {
        for (int j = 1; j <= nbinsy; ++j) {
            double e = fFunction->Integral(hist_marg->GetXaxis()->GetBinCenter(i) - 0.5 * binwidthx,
                                           hist_marg->GetXaxis()->GetBinCenter(i) + 0.5 * binwidthx,
                                           hist_marg->GetYaxis()->GetBinCenter(j) - 0.5 * binwidthy,
                                           hist_marg->GetYaxis()->GetBinCenter(j) + 0.5 * binwidthy
                                          );
            hist_func->SetBinContent(i, j, e);
        }
    }
    hist_func->Scale(norm_hist / hist_func->Integral());

    // calculate chi2
    double chi2 = 0;
    for (int i = 1; i <= nbinsx; ++i) {
        for (int j = 1; j <= nbinsy; ++j) {
            // use if histograms are properly normalized as 2D PDF
            //            double n = norm_hist * binwidthx * binwidthy * hist_marg->GetBinContent(i, j);
            double n = hist_marg->GetBinContent(i, j);
            double e = hist_func->GetBinContent(i, j);

            if (e >= 10)
                chi2 += (n - e) * (n - e) / e;
            else
                ndf--;

            // fill histograms
            hist_func->SetBinContent(i, j, e);
            hist_diff->SetBinContent(i, j, (n - e) / sqrt(e));
            hist_pull->Fill((n - e) / sqrt(e));
        }
    }
    if (ndf == 0)
        throw std::runtime_error(PerfTest::GetName() + ": Zero degrees of freedom");

    // define test results
    GetSubtest("chi2")->SetTargetValue(ndf);
    GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kGood,   3.0 * sqrt(2.0 * ndf));
    GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kAcceptable, 5.0 * sqrt(2.0 * ndf));
    GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kBad,    7.0 * sqrt(2.0 * ndf));
    GetSubtest("chi2")->SetTestValue(chi2);
    GetSubtest("chi2")->SetTestUncertainty(sqrt(2.0 * ndf));

    // add canvases
    TCanvas* c_func1 = new TCanvas();
    fFunction->Draw();
    AddCanvas(c_func1);
    AddCanvasDescription("The analytic function drawn with contours.");

    TCanvas* c_func2 = new TCanvas();
    hist_func->Draw("COLZ");
    AddCanvas(c_func2);
    AddCanvasDescription("The histogrammed analytic function. Each bin contains the integral of the analytic function over the bin.");

    TCanvas* c_marg = new TCanvas();
    hist_marg->Draw("COLZ");
    AddCanvas(c_marg);
    AddCanvasDescription("The distribution from MCMC.");

    TCanvas* c_func3 = new TCanvas();
    c_func3->SetLogz(kTRUE);
    hist_func->Draw("COLZ");
    AddCanvas(c_func3);
    AddCanvasDescription("The histogrammed analytic function in log-scale. Each bin contains the integral of the analytic function over the bin.");

    TCanvas* c_marg_log = new TCanvas();
    c_marg_log->SetLogz(kTRUE);
    hist_marg->Draw("COLZ");
    AddCanvas(c_marg_log);
    AddCanvasDescription("The distribution from MCMC in log-scale.");

    TCanvas* c_diff = new TCanvas();
    c_diff->cd();
    hist_diff->Draw("COLZ");
    AddCanvas(c_diff);
    AddCanvasDescription("The difference between the distribution from MCMC and the analytic function divided by the square root of the analytic function value in the corresponding bin.");

    TCanvas* c_pull = new TCanvas();
    c_pull->cd();
    hist_pull->Draw();
    TF1* g = new TF1("g", "[0]/sqrt(2.0*TMath::Pi())*exp(-x*x/2.)", -5.0, 5.0);
    g->SetParameter(0, hist_pull->Integral("")*hist_pull->GetBinWidth(1));
    g->Draw("SAMEP");
    AddCanvas(c_pull);
    AddCanvasDescription("The pull between the distribution from MCMC and the analytic function.");

    // debugKK
    // also add projections?

    // no error
    return 1;
}

//______________________________________________________________________________
void PerfTest2DFunction::DefineSubtests()
{
    PerfTestMCMC::DefineSubtests();

    PerfSubTest* subtest = new PerfSubTest("chi2");
    subtest->SetDescription("Calculate &chi;<sup>2</sup> and compare with prediction for dof=number of bins with an expectation >= 10. <br> Tolerance good: |&chi;<sup>2</sup>-E[&chi;<sup>2</sup>]| < 3 &middot; (2 dof)<sup>1/2</sup>, <br> Tolerance acceptable: |&chi;<sup>2</sup>-E[&chi;<sup>2</sup>]| < 5 &middot; (2 dof)<sup>1/2</sup>, <br> Tolerance bad: |&chi;<sup>2</sup>-E[&chi;<sup>2</sup>]| < 7 &middot; (2 dof)<sup>1/2</sup>.");
    AddSubtest(subtest);
}

//______________________________________________________________________________
