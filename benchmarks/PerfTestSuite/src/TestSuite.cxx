/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "include/TestSuite.h"
#include "include/PerfTest1DFunction.h"

#include <TStopwatch.h>

#include <fstream>
#include <iomanip>
#include <iostream>

//______________________________________________________________________________
TestSuite::TestSuite(bool multivariate, double dof):
    fMultivariate(multivariate),
    fDof(dof),
    fTestContainer(std::vector<PerfTest * >(0)),
    fNPlotColumns(3),
    fThumbSize(300),
    fIncludeHtmlHeader(true),
    fIncludeHtmlFooter(true),
    fLinkPrefix(""),
    fFileLinkPrefix(""),
    fHtmlFileExtension(".html")
{
    DefineTests();
}

//______________________________________________________________________________
TestSuite::~TestSuite()
{
    // delete all tests
    while (!fTestContainer.empty()) {
        PerfTest* test = fTestContainer.front();
        fTestContainer.erase(fTestContainer.begin());
        delete test;
    }
    fTestContainer.clear();

}

//______________________________________________________________________________
void TestSuite::SetPrecision(PerfTest::Precision precision)
{
    // loop over all tests and set precision
    for (unsigned i = 0; i < GetNTests(); ++i) {
        fTestContainer.at(i)->SetPrecision(precision);
    }

}

//______________________________________________________________________________
int TestSuite::GetNTests(PerfSubTest::Status status)
{
    // get number of sub tests
    unsigned n = GetNTests();

    // initialize counter
    int counter = 0;

    // loop over all tests and compare status
    for (unsigned i = 0; i < n; ++i) {
        if (fTestContainer.at(i)->GetStatus() == status)
            counter++;
    }

    // return counter
    return counter;
}

//______________________________________________________________________________
PerfTest* TestSuite::GetTest(const std::string& name)
{
    // get number of sub tests
    int n = GetNTests();

    // loop over all subtests and compare status
    for (int i = 0; i < n; ++i) {
        if (!name.compare(GetTest(i)->GetName()))
            return GetTest(i);
    }

    return 0;
}

//______________________________________________________________________________
int TestSuite::AddTest(PerfTest* test)
{
    // add test to container
    fTestContainer.push_back(test);

    // no error
    return 1;
}

//______________________________________________________________________________
int TestSuite::RunTests()
{
    // get number of sub tests
    int n = GetNTests();

    // initialize error code
    int err = 1;

    // loop over all tests, run and write the output to file
    for (int i = 0; i < n; ++i) {
        TStopwatch* sw = new TStopwatch();
        sw->Start();
        err *= fTestContainer.at(i)->Run();
        sw->Stop();
        fTestContainer.at(i)->SetCpuTime(sw->CpuTime());
        fTestContainer.at(i)->SetRealTime(sw->RealTime());
        err *= fTestContainer.at(i)->WriteResults();
        delete sw;
    }

    // return error code
    return err;
}

//______________________________________________________________________________
void TestSuite::PrintResultsScreen()
{
    std::cout << std::endl;
    std::cout << " Test results: " << std::endl << std::endl;

    std::cout << " Overview: " << std::endl;
    std::cout << " Number of tests               : " << GetNTests() << std::endl;
    std::cout << " Number of successful tests    : " << GetNTests(PerfSubTest::kGood) << std::endl;
    std::cout << " Number of acceptable tests        : " << GetNTests(PerfSubTest::kAcceptable) << std::endl;
    std::cout << " Number of bad tests           : " << GetNTests(PerfSubTest::kBad) << std::endl;
    std::cout << " Number of fatal tests         : " << GetNTests(PerfSubTest::kFatal) << std::endl;
    std::cout << " Number of tests unkown status : " << GetNTests(PerfSubTest::kUnknown) << std::endl;
    std::cout << std::endl;

    // loop over tests
    std::cout << " Details: " << std::endl;
    int n = GetNTests();
    for (int i = 0; i < n; ++i) {
        std::cout << " Test \"" << (GetTest(i)->GetName()).data() << "\" : " << GetTest(i)->GetStatusString().data() << std::endl;

        // loop over subtests
        int nsub = GetTest(i)->GetNSubtests();
        for (int j = 0; j < nsub; ++j)
            std::cout << "  Subtest \"" << (GetTest(i)->GetSubtest(j))->GetName() << "\" : " << (GetTest(i)->GetSubtest(j))->GetStatusString().data() << std::endl;
        std::cout << std::endl;
    }

}

//______________________________________________________________________________
void TestSuite::PrintResultsHTML(const std::string& filename)
{
    // open file
    std::ofstream file_main;
    file_main.open(filename.c_str());

    if (fIncludeHtmlHeader) {
        file_main << "<html>" << std::endl;
        file_main << "<head>" << std::endl;
        file_main << "<title>Test suite results</title>" << std::endl;
        file_main << "</head>" << std::endl << std::endl;
        file_main << "<body>" << std::endl << std::endl;
    }

    file_main << "<h3>Overview</h3>" << std::endl << std::endl;
    file_main << "<table border=\"0\"  width=\"600\">" << std::endl;
    file_main << " <tr> <td align=\"left\"> Number of tests </td> <td>" << GetNTests() << " </td> </tr>" << std::endl;
    file_main << " <tr> <td align=\"left\"> Number of successful tests </td> <td>" << GetNTests(PerfSubTest::kGood) << " </td> </tr>" << std::endl;
    file_main << " <tr> <td align=\"left\"> Number of acceptable tests </td> <td>" << GetNTests(PerfSubTest::kAcceptable) << " </td> </tr>" << std::endl;
    file_main << " <tr> <td align=\"left\"> Number of bad tests </td> <td>" << GetNTests(PerfSubTest::kBad) << " </td> </tr>" << std::endl;
    file_main << " <tr> <td align=\"left\"> Number of fatal tests </td> <td>" << GetNTests(PerfSubTest::kFatal) << " </td> </tr>" << std::endl;
    file_main << " <tr> <td align=\"left\"> Number of tests unkown status </td> <td>" << GetNTests(PerfSubTest::kUnknown) << " </td> </tr>" << std::endl;
    file_main << " <tr> <td align=\"left\">" << (fMultivariate ? "Multivariate" : "Factorized") <<  " proposal</td> <td>" << (fMultivariate ? Form("%g degrees of freedom", fDof) : "") << "</td> </tr>" << std::endl;
    file_main << "</table>" << std::endl;
    file_main << std::endl;

    file_main << "<h3>Function1D</h3>" << std::endl << std::endl;
    file_main << "<table border=\"0\"  width=\"70\">" << std::endl;
    file_main << "<tr>" << std::endl;
    file_main << "  <th align=\"left\" width=\"15%\"> Test </th>" << std::endl;
    file_main << "  <th align=\"left\"> Status </th>" << std::endl;
    file_main << "  <th align=\"left\"> Subtests </th>" << std::endl;
    file_main << "  <th align=\"left\"> Good </th>" << std::endl;
    file_main << "  <th align=\"left\"> Acceptable </th>" << std::endl;
    file_main << "  <th align=\"left\"> Bad </th>" << std::endl;
    file_main << "  <th align=\"left\"> Fatal </th>" << std::endl;
    file_main << "  <th align=\"left\"> Unknown </th>" << std::endl;
    file_main << "</tr>" << std::endl;
    for (unsigned i = 0; i < GetNTests(); ++i) {
        if (GetTest(i)->GetTestType() == PerfTest::kFunction1D) {
            file_main << " <tr>";
            file_main << " <td align=\"left\"><a href=\"" << fLinkPrefix << (GetTest(i)->GetName()).c_str() << fHtmlFileExtension.c_str() << "\">" << (GetTest(i)->GetName()).c_str() << "</a> </td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetStatusStringHTML().c_str() << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests() - GetTest(i)->GetNSubtests(PerfSubTest::kOff) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kGood) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kAcceptable) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kBad) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kFatal) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kUnknown) << "</td></tr>" << std::endl;
        }
    }
    file_main << "</table>" << std::endl;

    file_main << "<h3>Function2D</h3>" << std::endl << std::endl;
    file_main << "<table border=\"0\"  width=\"400\">" << std::endl;
    file_main << "<tr>" << std::endl;
    file_main << "  <th align=\"left\" width=\"70\"> Test </th>" << std::endl;
    file_main << "  <th align=\"left\"> Status </th>" << std::endl;
    file_main << "  <th align=\"left\"> Subtests </th>" << std::endl;
    file_main << "  <th align=\"left\"> Good </th>" << std::endl;
    file_main << "  <th align=\"left\"> Acceptable </th>" << std::endl;
    file_main << "  <th align=\"left\"> Bad </th>" << std::endl;
    file_main << "  <th align=\"left\"> Fatal </th>" << std::endl;
    file_main << "  <th align=\"left\"> Unknown </th>" << std::endl;
    file_main << "</tr>" << std::endl;
    for (unsigned i = 0; i < GetNTests(); ++i) {
        if (GetTest(i)->GetTestType() == PerfTest::kFunction2D) {
            file_main << " <tr>";
            file_main << " <td align=\"left\"><a href=\"" << fLinkPrefix << (GetTest(i)->GetName()).c_str() << fHtmlFileExtension.c_str() << "\">" << (GetTest(i)->GetName()).c_str() << "</a> </td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetStatusStringHTML().c_str() << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests() - GetTest(i)->GetNSubtests(PerfSubTest::kOff) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kGood) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kAcceptable) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kBad) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kFatal) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kUnknown) << "</td> </tr>" << std::endl;
        }
    }
    file_main << "</table>" << std::endl;

    file_main << "<h3>Varying parameters</h3>" << std::endl << std::endl;
    file_main << "<table border=\"0\"  width=\"70\">" << std::endl;
    file_main << "<tr>" << std::endl;
    file_main << "  <th align=\"left\" width=\"15%\"> Test </th>" << std::endl;
    file_main << "  <th align=\"left\"> Status </th>" << std::endl;
    file_main << "  <th align=\"left\"> Subtests </th>" << std::endl;
    file_main << "  <th align=\"left\"> Good </th>" << std::endl;
    file_main << "  <th align=\"left\"> Acceptable </th>" << std::endl;
    file_main << "  <th align=\"left\"> Bad </th>" << std::endl;
    file_main << "  <th align=\"left\"> Fatal </th>" << std::endl;
    file_main << "  <th align=\"left\"> Unknown </th>" << std::endl;
    file_main << "</tr>" << std::endl;
    for (unsigned i = 0; i < GetNTests(); ++i) {
        if (GetTest(i)->GetTestType() == PerfTest::kVarPar) {
            file_main << " <tr>";
            file_main << " <td align=\"left\"><a href=\"" << fLinkPrefix << (GetTest(i)->GetName()).c_str() << fHtmlFileExtension.c_str() << "\">" << (GetTest(i)->GetName()).c_str() << "</a> </td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetStatusStringHTML().c_str() << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests() - GetTest(i)->GetNSubtests(PerfSubTest::kOff) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kGood) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kAcceptable) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kBad) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kFatal) << "</td>" << std::endl;
            file_main << " <td align=\"left\">" << GetTest(i)->GetNSubtests(PerfSubTest::kUnknown) << "</td></tr>" << std::endl;
        }
    }
    file_main << "</table>" << std::endl;

    // loop over tests
    for (unsigned i = 0; i < GetNTests(); ++i) {

        // open file
        std::ofstream file;
        file.open((GetTest(i)->GetName() + fHtmlFileExtension).c_str());

        if (fIncludeHtmlHeader) {
            file << "<html>" << std::endl;
            file << "<head>" << std::endl;
            file << "<title>Test \"" << (GetTest(i)->GetName()).data() << "\" </title>" << std::endl;
            file << "</head>" << std::endl << std::endl;
            file << "<body>" << std::endl << std::endl;
        }

        file << " <h3>Test \"" << (GetTest(i)->GetName()).data() << "\" </h3>" << std::endl;

        file << "<table border=\"0\" width=\"30%\">" << std::endl;
        file << "<tr>";

        file << "  <th align=\"left\">Results</th>" << std::endl;
        file << "</tr>" << std::endl;

        file << " <tr> <td>Status</td> <td>" << GetTest(i)->GetStatusStringHTML().data() << " </td> </tr>" << std::endl;
        file << " <tr> <td>CPU time</td> <td>" << std::setprecision(4) << GetTest(i)->GetCpuTime() << " s </td></tr>" << std::endl;
        file << " <tr> <td>Real time</td> <td>" << std::setprecision(4) << GetTest(i)->GetRealTime() << " s </td></tr>" << std::endl;
        file << " <tr> <td>Plots</td> <td>" << "<a href=\"" << fFileLinkPrefix << GetTest(i)->GetName().data() << ".pdf" << "\">" << GetTest(i)->GetName().data() << ".pdf</a>" << " </td> </tr>" << std::endl;
        file << " <tr> <td>Log</td> <td>" << "<a href=\"" << fFileLinkPrefix << GetTest(i)->GetName().data() << ".log" << "\">" << GetTest(i)->GetName().data() << ".log</a>" << " </td> </tr>" << std::endl;
        file << "</table>" << std::endl;
        file << "</br>" << std::endl;

        if (GetTest(i)->GetTestType() ==  PerfTest::kFunction1D || GetTest(i)->GetTestType() ==  PerfTest::kFunction2D) {
            const PerfTest1DFunction* m = (PerfTest1DFunction*) GetTest(i);
            const bool converged = (m->GetNIterationsConvergenceGlobal() != -1);

            file << "<table border=\"0\" width=\"30%\">" << std::endl;
            file << "<tr>";
            file << "  <th align=\"left\">Settings</th>" << std::endl;
            file << "</tr>" << std::endl;
            file << " <tr> <td>multivariate</td> <td>" << std::boolalpha << m->GetProposeMultivariate() << " </td></tr>"  << std::endl;
            file << " <tr> <td>N chains</td> <td>" << m->GetNChains() << " </td></tr>"  << std::endl;
            file << " <tr> <td>N lag</td> <td>" << m->GetNLag() << " </td></tr>"  << std::endl;
            file << " <tr> <td>Convergence</td> <td>" << (converged ? "<font color=\"#4cc417\">true</font>" : "<font color=\"#FF8000\">false</font>") << " </td></tr>"  << std::endl;
            file << " <tr> <td>prerun iterations </td> <td>" << (converged ? m->GetNIterationsConvergenceGlobal() : m->GetNIterationsPreRunMax()) << " </td></tr>"  << std::endl;
            file << " <tr> <td>iterations (with lag)</td> <td>" << m->GetNIterationsRun() << " ("
                 << m->GetNIterationsRun() / m->GetNLag() << ") </td></tr>"  << std::endl;
            file << "</table>" << std::endl;
            file << "</br>" << std::endl;
        }

        // plots
        file << "<b>Plots</b>" << std::endl;
        int nplots = GetTest(i)->GetNCanvases();
        int rows = nplots / fNPlotColumns;

        file << "<table border=\"0\" width=\"80%\">" << std::endl;
        for (int j = 0; j <= rows; ++j) {
            file << "<tr>";
            for (int k = 0; k < fNPlotColumns; ++k) {
                if (j * fNPlotColumns + k < nplots)
                    file << "<td width=\"20%\"><a href=\"" << fFileLinkPrefix << GetTest(i)->GetName().c_str() << "_" << j* fNPlotColumns + k << ".pdf\"><img src=\"" << fFileLinkPrefix << GetTest(i)->GetName().c_str() << "_" << j* fNPlotColumns + k << ".png\" width=\"" << fThumbSize << "\" height=\"" << fThumbSize << "\"></a></td>";
            }
            file << "</tr>" << std::endl;
            file << "<tr>" << std::endl;
            for (int k = 0; k < fNPlotColumns; ++k) {
                if (j * fNPlotColumns + k < nplots)
                    file << "<td width=\"20%\" valign=\"top\">" << GetTest(i)->GetCanvasDescription(j * fNPlotColumns + k) << "</td>";
            }
            file << "</tr>" << std::endl;
        }
        file << "</table>" << std::endl;
        file << "</br>" << std::endl;

        // loop over subtests
        int nsub = GetTest(i)->GetNSubtests();
        if (nsub > 0) {
            file << "<table border=\"0\" width=\"80%\">" << std::endl;
            file << "<tr>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Subtest</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Status</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Target</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Test</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Uncertainty</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Deviation [%]</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Deviation [sigma]</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Tol. (Good)</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Tol. (Acceptable)</th>" << std::endl;
            file << "  <th align=\"left\" width = \"150\">Tol. (Bad)</th>" << std::endl;
            file << "</tr>" << std::endl;

            for (int j = 0; j < nsub; ++j) {
                double target = (GetTest(i)->GetSubtest(j))->GetTargetValue();
                double test = (GetTest(i)->GetSubtest(j))->GetTestValue();
                double err = (GetTest(i)->GetSubtest(j))->GetTestUncertainty();
                file << "<tr>" << std::endl;
                file << "  <td align=\"left\"> " << (GetTest(i)->GetSubtest(j))->GetName() << " </td> " << std::endl;
                file << "  <td align=\"left\"> " << (GetTest(i)->GetSubtest(j))->GetStatusStringHTML().data() << "</td>" << std::endl;
                file << "  <td align=\"left\"> " << std::setprecision(4) << target  << "</td>" << std::endl;
                file << "  <td align=\"left\"> " << std::setprecision(4) << test << "</td>" << std::endl;
                file << "  <td align=\"left\"> " << std::setprecision(4) << err << "</td>" << std::endl;
                if (target != 0)
                    file << "  <td align=\"left\"> " << std::setprecision(4) << -(target - test) / target * 100. << "</td>" << std::endl;
                else
                    file << "  <td align=\"left\"> " << std::setprecision(4) << " - " << "</td>" << std::endl;
                file << "  <td align=\"left\"> " << std::setprecision(4) << (target - test) / err << "</td>" << std::endl;
                file << "  <td align=\"left\"> " << std::setprecision(4) << (GetTest(i)->GetSubtest(j))->GetStatusRegion(PerfSubTest::kGood) << "</td>" << std::endl;
                file << "  <td align=\"left\"> " << std::setprecision(4) << (GetTest(i)->GetSubtest(j))->GetStatusRegion(PerfSubTest::kAcceptable) << "</td>" << std::endl;
                file << "  <td align=\"left\"> " << std::setprecision(4) << (GetTest(i)->GetSubtest(j))->GetStatusRegion(PerfSubTest::kBad) << "</td>" << std::endl;
                file << "</tr>" << std::endl;
            }
            file << "</table>" << std::endl;
            file << "<br>" << std::endl;

            // loop over subtests
            file << "<table border=\"0\" width=\"100%\">" << std::endl;
            file << "<tr>" << std::endl;
            file << "  <th align=\"left\" valign=\"top\">Subtest</th>" << std::endl;
            file << "  <th align=\"left\">Description</th>" << std::endl;
            file << "</tr>" << std::endl;
            for (int j = 0; j < nsub; ++j) {
                file << "<tr>" << std::endl;
                file << "  <td align=\"left\" width = \"150\"  valign=\"top\"> " << (GetTest(i)->GetSubtest(j))->GetName() << " </td> " << std::endl;
                file << "  <td align=\"left\"> " << (GetTest(i)->GetSubtest(j))->GetDescription().data() << "</td>" << std::endl;
                file << "</tr>" << std::endl;
            }
            file << "</table>" << std::endl;
        }

        if (fIncludeHtmlFooter) {
            file << std::endl;
            file << "<br>" << std::endl;
            file << "<br>" << std::endl;
            file << std::endl;
            file << "</body>" << std::endl << std::endl;
            file << "</html>" << std::endl;
        }

        // close file
        file.close();
    } // end loop tests

    if (fIncludeHtmlFooter) {
        file_main << std::endl;
        file_main << "<br>" << std::endl;
        file_main << "<br>" << std::endl;
        file_main << std::endl;
        file_main << "</body>" << std::endl << std::endl;
        file_main << "</html>" << std::endl;
    }

    // close file
    file_main.close();
}

//______________________________________________________________________________
void TestSuite::DefineTests()
{
    // ...
}

//______________________________________________________________________________
