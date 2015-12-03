/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "include/PerfTestVarPar.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPostScript.h>
#include <TGraphErrors.h>

#include <cmath>
#include <iostream>
#include <fstream>

//______________________________________________________________________________
PerfTestVarPar::PerfTestVarPar(const std::string& name, PerfTestMCMC* test)
    : PerfTest(name)
    , fTest(test)
    , fTargetContainer(std::vector<TGraphErrors * >(0))
    , fTestContainer(std::vector<TGraphErrors * >(0))
{
    // set test type
    fTestType = PerfTest::kVarPar;
}

//______________________________________________________________________________
PerfTestVarPar::~PerfTestVarPar()
{
    // delete all target value graphs
    while (!fTargetContainer.empty()) {
        TGraphErrors* graph = fTargetContainer.front();
        fTargetContainer.erase(fTargetContainer.begin());
        delete graph;
    }
    fTargetContainer.clear();

    // delete all test value graphs
    while (!fTestContainer.empty()) {
        TGraphErrors* graph = fTestContainer.front();
        fTestContainer.erase(fTestContainer.begin());
        delete graph;
    }
    fTestContainer.clear();

    delete fTest;
}

//______________________________________________________________________________
int PerfTestVarPar::AddVarPar(std::vector<double> values, const std::string& name)
{
    // get number of values
    int nval = int(values.size());

    // check values
    if (nval  == 0)
        return 0;

    // copy parameter values
    fVarParValues = values;

    // set parameter name
    fVarParName = name;

    // no error
    return 1;
}

//______________________________________________________________________________
int PerfTestVarPar::PreTest()
{
    // get number of subtests
    int ntests = fTest->GetNSubtests();

    // get number of variation parameter values
    int npar = GetNVarPar();

    for (int i = 0; i < ntests; ++i) {
        TGraphErrors* g1 = new TGraphErrors(npar);
        g1->SetMarkerColor(kRed);
        g1->SetMarkerStyle(20);
        fTargetContainer.push_back(g1);
        TGraphErrors* g2 = new TGraphErrors(npar);
        g2->SetMarkerColor(kBlack);
        g2->SetMarkerStyle(21);
        fTestContainer.push_back(g2);
    }

    // no errors
    return 1;
}

//______________________________________________________________________________
int PerfTestVarPar::PostTest()
{

    // get number of subtests
    int ntests = fTest->GetNSubtests();

    // get number of variation parameter values
    int npar = GetNVarPar();

    // loop over subtests
    for (int i = 0; i < ntests; ++i) {
        TCanvas* c = new TCanvas();
        c->cd();

        double x;
        double ymin;
        double ymax;
        fTargetContainer.at(i)->GetPoint(0, x, ymin);
        ymax = ymin;

        for (int j = 0; j < npar; ++j) {
            double valx;
            double valy;
            fTargetContainer.at(i)->GetPoint(j, valx, valy);
            if (valy < ymin)
                ymin = valy;
            if (valy > ymax)
                ymax = valy;
            fTestContainer.at(i)->GetPoint(j, valx, valy);
            if (valy - fTestContainer.at(i)->GetErrorY(j) < ymin)
                ymin = valy - fTestContainer.at(i)->GetErrorY(j);
            if (valy + fTestContainer.at(i)->GetErrorY(j) > ymax)
                ymax = valy + fTestContainer.at(i)->GetErrorY(j);
        }

        TH2D* hist = new TH2D("", Form(";%s;%s", fVarParName.c_str(), fTest->GetSubtest(i)->GetName().c_str()),
                              1, fVarParValues.at(0) - fabs(0.1 * fVarParValues.at(npar - 1)), fVarParValues.at(npar - 1) + fabs(0.1 * fVarParValues.at(npar - 1)),
                              1, ymin - fabs(0.1 * ymax), ymax + fabs(0.1 * ymax));
        hist->SetStats(kFALSE);
        hist->Draw();
        fTargetContainer.at(i)->Draw("PL");
        fTestContainer.at(i)->Draw("PL");
        AddCanvas(c);
        AddCanvasDescription("Subtest values as a function of the variation parameter. Red: target values, blue: test values and uncertainties.");
    }

    // no error
    return 1;
}

//______________________________________________________________________________
int PerfTestVarPar::RunTest()
{
    // get number of subtests
    int ntests = fTest->GetNSubtests();

    // get number of variation parameter values
    int npar = GetNVarPar();

    // run all tests
    for (int i = 0; i < npar; ++i) {
        // set the corresponding value
        fTest->SetVarPar(fVarParValues.at(i), fVarParName);

        // run the test
        fTest->Run();

        // copy the output
        // fill histograms with target and test value vs. var par value

        for (int j = 0; j < ntests; ++j) {
            // get target and test values
            double target = fTest->GetSubtest(j)->GetTargetValue();
            double test = fTest->GetSubtest(j)->GetTestValue();
            double sigma = fTest->GetSubtest(j)->GetTestUncertainty();

            // fill graphs
            fTargetContainer.at(j)->SetPoint(i, fVarParValues.at(i), target);
            fTestContainer.at(j)->SetPoint(i, fVarParValues.at(i), test);
            fTestContainer.at(j)->SetPointError(i, 0., sigma);
        }

    }

    // no error
    return 1;
}

//______________________________________________________________________________
void PerfTestVarPar::PrecisionSettings(PerfTest::Precision precision)
{
    if (fTest)
        fTest->PrecisionSettings(precision);
}

//______________________________________________________________________________

