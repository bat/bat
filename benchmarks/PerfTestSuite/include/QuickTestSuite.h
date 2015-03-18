/*!
 * \class BAT::QuickTestSuite
 * \brief A generic test suite for BAT checking basic functionality.
 * It incorporates other test suite. The whole suite hould run in under 10s.
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_QUICKTESTSUITE
#define BAT_QUICKTESTSUITE

#include <cpptest.h>



namespace BAT
{

class QuickTestSuite: public Test::Suite
{

public:

    /** \name Enumerators  */
    /* @{ */

    /* @} */
    /** \name Constructors and destructors  */
    /* @{ */

    /** The default constructor */
    QuickTestSuite();

    /** The default destructor */
    ~QuickTestSuite();

    /* @} */
    /** \name Member functions (Set)  */
    /* @{ */

    /* @} */
    /** \name Member functions (Get)  */
    /* @{ */


    /* @} */
    /** \name Member functions (misc)  */
    /* @{ */

    /* @} */

private:

};

} // namespace BAT

#endif
