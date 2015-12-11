/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "include/PerfTestMCMC.h"

#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TGraph.h>

#include <iostream>

//______________________________________________________________________________
PerfTestMCMC::PerfTestMCMC(const std::string& name)
    : PerfTest(name)
    , BCModel(name.c_str())
    , fCorrelation(std::vector<TGraph * >(0))
    , fHistCorr(std::vector<TH2D * >(0))
    , fXOld(std::vector<std::vector<double> >(0))
{
    // define subtests
    DefineSubtests();
}

//______________________________________________________________________________
PerfTestMCMC::~PerfTestMCMC()
{
    // delete all graphs
    while (!fCorrelation.empty()) {
        TGraph* can = fCorrelation.front();
        fCorrelation.erase(fCorrelation.begin());
        delete can;
    }
    fCorrelation.clear();

    // delete all histograms
    while (!fHistCorr.empty()) {
        TH2D* hist = fHistCorr.front();
        fHistCorr.erase(fHistCorr.begin());
        delete hist;
    }
    fHistCorr.clear();

}

//______________________________________________________________________________
int PerfTestMCMC::SetVarPar(double value, const std::string& name)
{
    if (name == "lag") {
        int n = GetNIterationsRun();
        int lag = GetNLag();
        int iter = int( n * value / double(lag) );
        SetNLag(int(value));
        SetNIterationsRun(iter);
        return 1;
    } else if (name == "iterations") {
        SetNIterationsRun(int(value));
        return 1;
    } else
        return 0;
}

//______________________________________________________________________________
int PerfTestMCMC::PreTest()
{
    //add histograms
    int npar = GetNParameters();

    // clear histograms
    fHistCorr.clear();
    fCorrelation.clear();

    // loop over parameters
    for (int i = 0; i < npar; ++i) {
        double xmin = GetParameter(i).GetLowerLimit();
        double xmax = GetParameter(i).GetLowerLimit();
        TH2D* hist = new TH2D("", "", 100, xmin, xmax, 100, xmin, xmax);
        fHistCorr.push_back(hist);
        fCorrelation.push_back(new TGraph(0));
    }

    return 1;
}

//______________________________________________________________________________
int PerfTestMCMC::PostTest()
{

    // loop over parameters
    int npar = GetNParameters();
    for (int i = 0; i < npar; ++i) {
        TCanvas* c_corr = new TCanvas();
        TH2D* hist = new TH2D("", ";Iteration; correlation coefficient", 1, 0., GetNIterationsRun() / GetNLag() + 1, 1, -1., 1.0);
        hist->SetStats(kFALSE);
        hist->Draw();
        fCorrelation.at(i)->Draw("SAMEP");
        AddCanvas(c_corr);
        AddCanvasDescription("Correlation coefficient of iterations with lag in between.");

        double x;
        double y;
        fCorrelation.at(i)->GetPoint(fCorrelation.at(i)->GetN() - 1, x, y);

        // define test results
        GetSubtest(Form("correlation par %i", i))->SetTargetValue(0.0);
        GetSubtest(Form("correlation par %i", i))->SetStatusRegion(PerfSubTest::kGood,   0.3);
        GetSubtest(Form("correlation par %i", i))->SetStatusRegion(PerfSubTest::kAcceptable, 0.5);
        GetSubtest(Form("correlation par %i", i))->SetStatusRegion(PerfSubTest::kBad,    0.7);
        GetSubtest(Form("correlation par %i", i))->SetTestValue(y);
        GetSubtest(Form("correlation par %i", i))->SetTestUncertainty(fCorrelation.at(i)->GetRMS(2));
        GetSubtest(Form("correlation par %i", i))->SetStatusOff(true);
    }

    // no error
    return 1;
}

//______________________________________________________________________________
int PerfTestMCMC::RunTest()
{
    // define error code
    int err = 1;

    // perform mcmc
    err *= MarginalizeAll(BCIntegrate::kMargMetropolis);

    // return error code
    return err;
}

//______________________________________________________________________________
void PerfTestMCMC::DefineSubtests()
{

    // loop over parameters
    int npar = GetNParameters();
    for (int i = 0; i < npar; ++i) {
        PerfSubTest* subtest = new PerfSubTest(Form("correlation par %i", i));
        subtest->SetDescription("Calculate the auto-correlation among the points.");
        AddSubtest(subtest);
    }
}

//______________________________________________________________________________
int PerfTestMCMC::WriteResults()
{
    PerfTest::WriteResults();

    PrintSummary();

    return 1;
}

//______________________________________________________________________________
void PerfTestMCMC::PrecisionSettings(PerfTest::Precision precision)
{
    unsigned lag = 1;
    if (precision == PerfTest::kCoarse) {
        BCEngineMCMC::SetPrecision(BCEngineMCMC::kLow);
    } else if (precision == PerfTest::kMedium) {
        BCEngineMCMC::SetPrecision(BCEngineMCMC::kMedium);
        SetNIterationsRun(5000);
        lag = 50;
    } else if (precision == PerfTest::kDetail) {
        BCEngineMCMC::SetPrecision(BCEngineMCMC::kMedium);
        lag = 100;
    } else {
        SetNIterationsRun(10000);
    }

    // interpret NIterationsRun as desired #independent samples
    unsigned temp(GetNIterationsRun() * lag);
    unsigned iter = temp / GetNLag();
    SetNLag(lag);
    SetNIterationsRun(iter);
}

//______________________________________________________________________________
void PerfTestMCMC::MCMCUserIterationInterface()
{
    // copy over old point on first call
    if (fXOld.size() < fMCMCx.size()) {
        fXOld = fMCMCx;
        return;
    }

    unsigned iteration = GetCurrentIteration();
    unsigned nlag = GetNLag();

    if ((iteration % nlag) == 0) {

        // loop over parameters
        unsigned npar = GetNParameters();
        unsigned nchains = GetNChains();

        for (unsigned i = 0; i < npar; ++i) {
            TH2D* hist = fHistCorr.at(i);

            for (unsigned j = 0; j < nchains; ++j) {
                hist->Fill(fXOld.at(j).at(i), fMCMCx.at(j).at(i));
            }

            if (iteration / nlag % (GetNIterationsRun() / 100 / nlag) == 0) {
                (fCorrelation[i])->SetPoint( (iteration / nlag) / (GetNIterationsRun() / 100 / nlag),
                                             iteration / nlag,
                                             hist->GetCorrelationFactor());
            }
        }
        // copy old point
        fXOld = fMCMCx;
    }
}

//______________________________________________________________________________
