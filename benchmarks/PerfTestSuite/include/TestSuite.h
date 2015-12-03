/*!
   \class BAT::TestSuite
   \brief A test suite class for BAT

   Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
   All rights reserved.

   For the licensing terms see doc/COPYING.
*/

#ifndef BAT_TESTSUITE
#define BAT_TESTSUITE

#include <string>
#include <vector>

#include <include/PerfTest.h>

class TestSuite
{

public:

    /** \name Enumerators  */
    /* @{ */

    /* @} */
    /** \name Constructors and destructors */
    /* @{ */

    /** The default constructor */
    TestSuite(bool multivariate, double dof);

    /** The default destructor */
    ~TestSuite();

    /* @} */
    /** \name Member functions (Set) */
    /* @{ */

    /** Set the precision of the tests
     * @param the precision. */
    void SetPrecision(PerfTest::Precision precision);

    /** Set number of plots per line in the model html page */
    void SetNPlotColumns(int n)
    { fNPlotColumns = n; };

    /** Set size of the thumbnail in the html page of the model results */
    void SetThumbSize(int n)
    { fThumbSize = n; };

    /** Flag whether or not to include html header in the output html files */
    void IncludeHtmlHeader(bool flag = true)
    { fIncludeHtmlHeader = flag; };

    /** Flag whether or not to include html footer in the output html files */
    void IncludeHtmlFooter(bool flag = true)
    { fIncludeHtmlFooter = flag; };

    /** Prefix to be used for internal html links */
    void SetLinkPrefix(const std::string& prefix)
    { fLinkPrefix = prefix; };

    /** Prefix to be used for html links to external files */
    void SetFileLinkPrefix(const std::string& prefix)
    { fFileLinkPrefix = prefix; };

    /** Set file extension for the html files (.html by default) */
    void SetHtmlFileExtension(const std::string& ext)
    { fHtmlFileExtension = ext; };

    /* @} */
    /** \name Member functions (Get) */
    /* @{ */

    /** Get the number of tests which belong to this test.
     * @return the number of tests. */
    size_t GetNTests()
    { return fTestContainer.size(); };

    /** Get the number of tests which belong to this test with
     * specified status.
     * @param status the status code.
     * @return the number of tests. */
    int GetNTests(PerfSubTest::Status status);

    /** Find a test by index
     * @param index the index of the test.
     * @return the subtest. */
    PerfTest* GetTest(unsigned index)
    { return fTestContainer.at(index); };

    /** Find a test by name
     * @param name the name of the test.
     * @return the test. */
    PerfTest* GetTest(const std::string& name);

    /* @} */
    /** \name Member functions (misc) */
    /* @{ */

    /** Define all tests here. */
    void DefineTests();

    /** Add a test.
     * @param test the test to be added.
     * @return An error code. */
    int AddTest(PerfTest* test);

    /** Create the HTML output of the tests.
     * @param filename the name of the HTML file. */
    void PrintResultsHTML(const std::string& filename = "results.html");

    /** Create an ASCII output of the tests.
     * @param filename the name of the ASCII file. */
    void PrintResultsASCII(const std::string& filename = "results.txt");

    /** Create an ASCII output of the tests.
     * @param filename the name of the ASCII file. */
    void PrintResultsScreen();

    /** Run all tests.
     * @return an error code. */
    int RunTests();

    /* @} */

protected:
    /** Multivariate MCMC proposal */
    bool fMultivariate;

    /** Student's T degree of freedom */
    double fDof;

private:

    /** A container of tests which belong to the test suite. */
    std::vector<PerfTest*> fTestContainer;

    /** Number of plots per line in the html page of the model results */
    int fNPlotColumns;

    /** Size of the thumbnail in the html page of the model results */
    int fThumbSize;

    /** Flag whether or not to include html header in the output html files */
    bool fIncludeHtmlHeader;

    /** Flag whether or not to include html footer in the output html files */
    bool fIncludeHtmlFooter;

    /** Prefix to be used for internal html links */
    std::string fLinkPrefix;

    /** Prefix to be used for html links to external files */
    std::string fFileLinkPrefix;

    /** Set file extension for the html files (.html by default) */
    std::string fHtmlFileExtension;
};

#endif
