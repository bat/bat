/*
 * Copyright (C) 2012, Danny van Dyk and Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <test.h>

#include <cmath>

using namespace test;

class NoThrowTest :
    public TestCase
{
    public:
        NoThrowTest() :
            TestCase("no_throw_test")
        {
        }

        virtual void run() const
        {
            try
            {
                TEST_CHECK_NO_THROW(throw std::string("failed"));
            }
            catch (TestCaseFailedException & e)
            {
                // as should be
            }
            catch (std::string & e)
            {
                TEST_CHECK_FAILED("failed");
            }
        }
} no_throw_test;

class EqualTest :
    public TestCase
{
    public:
        EqualTest() :
            TestCase("equal_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_NO_THROW(TEST_CHECK_EQUAL(0, 0));
            TEST_CHECK_NO_THROW(TEST_CHECK_EQUAL(std::string("foo"), std::string("foo")));
            TEST_CHECK_NO_THROW(TEST_CHECK_EQUAL(0.0, 0.0));

            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_EQUAL(0, 1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_EQUAL(std::string("foo"), std::string("bar")));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_EQUAL(17.0, 23.0));
        }
} equal_test;

class RelativeErrorTest :
    public TestCase
{
    public:
        RelativeErrorTest() :
            TestCase("relative_error_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK_NO_THROW(TEST_CHECK_RELATIVE_ERROR(1.0, 1.09, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(+1.0, +2.0, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(+1.0, -2.0, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(-1.0, +2.0, 0.1));
            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(-1.0, -2.0, 0.1));

            TEST_CHECK_THROWS(TestCaseFailedException, TEST_CHECK_RELATIVE_ERROR(-0.1, 0.0, 0.2));
        }
} relative_error_test;
