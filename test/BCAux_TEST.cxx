/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>
#include <BAT/BCAux.h>

using namespace test;
using namespace BCAux;

class DefaultToPdfTest :
    public TestCase
{
public:
    DefaultToPdfTest() :
        TestCase("DefaultToPdf test")
    {
    }

    virtual void run() const
    {
        std::string file = "foo.pdf";
        DefaultToPDF(file);
        TEST_CHECK_EQUAL(file, "foo.pdf");

        file = "foo.ps";
        DefaultToPDF(file);
        TEST_CHECK_EQUAL(file, "foo.ps");

        file = "foo";
        DefaultToPDF(file);
        TEST_CHECK_EQUAL(file, "foo.pdf");

        file = "foo.";
        DefaultToPDF(file);
        TEST_CHECK_EQUAL(file, "foo.pdf");

        file = "foo.PDF";
        DefaultToPDF(file);
        TEST_CHECK_EQUAL(file, "foo.PDF");

        file = "foo.PS";
        DefaultToPDF(file);
        TEST_CHECK_EQUAL(file, "foo.PS");

        file = "foo.bar";
        DefaultToPDF(file);
        TEST_CHECK_EQUAL(file, "foo.bar.pdf");
    }
} defaultToPdfTest;
