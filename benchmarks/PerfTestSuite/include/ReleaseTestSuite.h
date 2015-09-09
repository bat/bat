/*!
 * \class BAT::ReleaseTestSuite
 * \brief A generic test suite for BAT checking basic functionality.
 * It incorporates other test suite. The whole suite hould run in under 10s.
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_RELEASETESTSUITE
#define BAT_RELEASETESTSUITE

#include <include/TestSuite.h>

class ReleaseTestSuite : public TestSuite
{

public:

    /** \name Enumerators  */
    /* @{ */

    /* @} */
    /** \name Constructors and destructors  */
    /* @{ */

    /** The default constructor */
    ReleaseTestSuite(bool multivariate, double dof);

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */

    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /**
     * Prepare all tests.
     * @return An error code. */
    int PrepareTests();

    /** Setup html as needed for BAT webpage */
    void WebpageSetup();

    /* @} */
};

#endif
