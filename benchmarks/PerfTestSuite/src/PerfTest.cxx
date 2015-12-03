/*
 * Copyright (C) 2008-2010, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "include/PerfTest.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TPDF.h>

#include <iostream>
#include <fstream>

//______________________________________________________________________________
PerfTest::PerfTest(const std::string& name)
    : fTestType(PerfTest::kUnknown)
    , fPrecision(PerfTest::kCoarse)
    , fSubtestContainer(std::vector<PerfSubTest * >(0))
    , fCanvasContainer(std::vector<TCanvas * >(0))
    , fCanvasDescriptionContainer(std::vector<std::string>(0))
    , fName(name)
    , fRealTime(0.)
    , fCpuTime(0.)
{
    // define subtests
    //   DefineSubtests();
}

//______________________________________________________________________________
PerfTest::~PerfTest()
{
    // delete all subtests
    while (!fSubtestContainer.empty()) {
        PerfSubTest* test = fSubtestContainer.front();
        fSubtestContainer.erase(fSubtestContainer.begin());
        delete test;
    }
    fSubtestContainer.clear();

    // delete all canvases
    while (!fCanvasContainer.empty()) {
        TCanvas* can = fCanvasContainer.front();
        fCanvasContainer.erase(fCanvasContainer.begin());
        delete can;
    }
    fCanvasContainer.clear();

}

//______________________________________________________________________________
std::string TypeToString(PerfTest::TestType type)
{
    switch (type) {
        case PerfTest::kFunction1D :
            return std::string("Function1D");

        case PerfTest::kFunction2D :
            return std::string("Function2D");

        default :
            return std::string("-");
    }
}

//______________________________________________________________________________
std::string PerfTest::ToString(PerfSubTest::Status status)
{
    PerfSubTest st;
    return st.ToString(status);
}

//______________________________________________________________________________
std::string PerfTest::ToStringHTML(PerfSubTest::Status status)
{
    PerfSubTest st;
    return st.ToStringHTML(status);
}

//______________________________________________________________________________
void PerfTest::SetPrecision(PerfTest::Precision precision)
{
    fPrecision = precision;

    PrecisionSettings(precision);
}

//______________________________________________________________________________
int PerfTest::GetNSubtests(PerfSubTest::Status status)
{
    // get number of sub tests
    int n = GetNSubtests();

    // initialize counter
    int counter = 0;

    // loop over all subtests and compare status
    for (int i = 0; i < n; ++i) {
        if (fSubtestContainer.at(i)->GetStatus() == status)
            counter++;
    }

    // return counter
    return counter;
}

//______________________________________________________________________________
PerfSubTest::Status PerfTest::GetStatus()
{
    // get number of active sub tests
    int n = GetNSubtests() - GetNSubtests(PerfSubTest::kOff);

    // get number of successful sub tests
    int ngood = GetNSubtests(PerfSubTest::kGood);

    // get number of acceptable sub tests
    int nacceptable = GetNSubtests(PerfSubTest::kAcceptable);

    // get number of failed sub tests
    int nbad = GetNSubtests(PerfSubTest::kBad);

    // get number of failed sub tests
    int nfatal = GetNSubtests(PerfSubTest::kFatal);

    // get number of unkown
    int nunknown = GetNSubtests(PerfSubTest::kUnknown);

    // calculate overall status
    if (n == ngood)
        return PerfSubTest::kGood;

    else if (nunknown > 0)
        return PerfSubTest::kUnknown;

    else if (nfatal > 0)
        return PerfSubTest::kFatal;

    else if (nbad > 0)
        return PerfSubTest::kBad;

    else if (nacceptable > 0)
        return PerfSubTest::kAcceptable;

    return PerfSubTest::kUnknown;
}

//______________________________________________________________________________
PerfSubTest* PerfTest::GetSubtest(const std::string& name)
{
    // get number of sub tests
    int n = GetNSubtests();

    // loop over all subtests and compare status
    for (int i = 0; i < n; ++i) {
        if (!name.compare(GetSubtest(i)->GetName()))
            return GetSubtest(i);
    }

    return 0;
}

//______________________________________________________________________________
TCanvas* PerfTest::GetCanvas(int index)
{
    int ncanvases = GetNCanvases();

    // check index
    if (index < 0 || index >= ncanvases)
        return 0;

    // return canvas pointer
    return fCanvasContainer.at(index);
}

//______________________________________________________________________________
std::string PerfTest::GetCanvasDescription(int index)
{
    int ncanvases = int(fCanvasDescriptionContainer.size());

    // check index
    if (index < 0 || index >= ncanvases)
        return std::string("-");

    // return canvas pointer
    return fCanvasDescriptionContainer.at(index);
}

//______________________________________________________________________________
int PerfTest::ReadResults()
{
    /*
    // open file
    std::fstream file;
    file.open((fName+std::string(".tst")).c_str(), std::fstream::in);

    // check if file is open
    if (!file.is_open())
    {
    std::cout << "Could not open file." << std::endl;
    return 0;
    }

    // read data from file
    std::string dummy_string;
    double dummy_double;
    bool dummy_bool;
    int dummy_int;

    file >> dummy_string;
    file >> dummy_int;

    // check name
    if (fName.compare(dummy_string)) {
    std::cout << "Test name and name in file do not agree." << std::endl;
    file.close();
    return 0;
    }

    // check number of subtests
    if (GetNSubtests() != dummy_int) {
    std::cout << "Number of subtests no consistent." << std::endl;
    file.close();
    return 0;
    }

    // loop over all subtests and read results
    for (int i = 0; i < dummy_int; ++i) {

    // get subtest
    PerfSubTest * subtest = GetSubtest(i);

    file >> dummy_string;

    // check name
    if (!(subtest->GetName().compare(dummy_string)==0)) {
    std::cout << "Subtest name and name in file do not agree." << std::endl;
    file.close();
    return 0;
    }

    file >> dummy_double;
    subtest->SetTestValue(dummy_double);
    file >> dummy_double;
    subtest->SetStatusRegion(PerfSubTest::kGood, dummy_double);
    file >> dummy_double;
    subtest->SetStatusRegion(PerfSubTest::kAcceptable, dummy_double);
    file >> dummy_double;
    subtest->SetStatusRegion(PerfSubTest::kBad, dummy_double);
    file >> dummy_bool;
    subtest->SetStatusUnknown(dummy_bool);
    file >> dummy_bool;
    subtest->SetStatusOff(dummy_bool);
    }

    // close file
    file.close();
    */

    // no error
    return 1;
}

//______________________________________________________________________________
int PerfTest::WriteResults()
{
    // open file
    std::fstream file;
    file.open((fName.data() + std::string(".tst")).c_str(), std::fstream::out);

    // check if file is open
    if (!file.is_open()) {
        std::cout << "Could not open file." << std::endl;
        return 0;
    }

    // write to file
    file << fName.data() << std::endl;
    file << GetNSubtests() << std::endl;

    // get number of active sub tests
    int n = GetNSubtests() - GetNSubtests(PerfSubTest::kOff);

    // loop over all subtests and write results
    for (int i = 0; i < n; ++i) {
        file << fSubtestContainer.at(i)->GetName().c_str() << std::endl;
        file << fSubtestContainer.at(i)->GetTestValue() << std::endl;
        file << fSubtestContainer.at(i)->GetTargetValue() << std::endl;
        file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kGood) << std::endl;
        file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kAcceptable) << std::endl;
        file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kBad) << std::endl;
        file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kFatal) << std::endl;
        file << fSubtestContainer.at(i)->GetStatusUnknown() << std::endl;
        file << fSubtestContainer.at(i)->GetStatusOff() << std::endl;
    }

    // close file
    file.close();

    // create postscript
    TPDF* pdf = new TPDF((fName.data() + std::string(".pdf")).c_str());

    // get number of canvases
    int nhist = GetNCanvases();

    pdf->NewPage();

    // loop over histograms
    for (int i = 0; i < nhist; ++i) {
        // get canvas
        TCanvas* c = GetCanvas(i);

        // update post script
        c->Update();
        if (i != nhist - 1)
            pdf->NewPage();
        c->cd();
    }

    // close ps
    pdf->Close();

    // print thumbnials
    for (int i = 0; i < nhist; ++i) {
        // get canvas
        TCanvas* c = GetCanvas(i);
        c->Print( Form("%s_%i.png", GetName().c_str(), i) );
        c->Print( Form("%s_%i.pdf", GetName().c_str(), i) );
    }

    // delete postscript
    delete pdf;

    // no error
    return 1;
}

//______________________________________________________________________________
int PerfTest::Run()
{
    // define error code
    int err = 1;

    // call pre-test
    err *= PreTest();

    // perform mcmc
    err *= RunTest();

    // call post-test
    err *= PostTest();

    // return error code
    return err;
}

//______________________________________________________________________________

