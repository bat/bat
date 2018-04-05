/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>

#include <BAT/BCParameterSet.h>

using namespace test;

class BCParameterTest :
    public TestCase
{
public:
    BCParameterTest() :
        TestCase("BCParameter test")
    {
    }

    virtual void run() const
    {
        BCParameter a("a", 0.2, 1.7);
        TEST_CHECK_EQUAL(a.GetLowerLimit(), 0.2);
        TEST_CHECK_EQUAL(a.GetUpperLimit(), 1.7);

        // copy ctor
        BCParameter b(a);
        TEST_CHECK_EQUAL(b.GetLowerLimit(), 0.2);
        TEST_CHECK_EQUAL(b.GetUpperLimit(), 1.7);

        // copy assignment
        BCParameter c("c", 12, 145.4);
        c = a;
        TEST_CHECK_EQUAL(c.GetLowerLimit(), 0.2);
        TEST_CHECK_EQUAL(c.GetUpperLimit(), 1.7);
    }

} bcparameter_Test;
