/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>
#include <BAT/BCAux.h>

#include <limits>

using namespace test;
using namespace BCAux;

class BCAuxTest :
    public TestCase
{
public:
    BCAuxTest() :
        TestCase("BCAux test")
    {
    }

    void defaultToPDFTest() const
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

    void rangeTypeTest() const
    {
        double inf_val = std::numeric_limits<double>::infinity();
        double fin_val = 6.2;
        TEST_CHECK_EQUAL(RangeType(-fin_val, fin_val),  kFiniteRange);
        TEST_CHECK_EQUAL(RangeType(-inf_val, fin_val),  kNegativeInfiniteRange);
        TEST_CHECK_EQUAL(RangeType(fin_val, inf_val),   kPositiveInfiniteRange);
        TEST_CHECK_EQUAL(RangeType(-inf_val, inf_val),  kInfiniteRange);
        TEST_CHECK_EQUAL(RangeType(fin_val, fin_val),   kEmptyRange);
        TEST_CHECK_EQUAL(RangeType(inf_val, inf_val),   kEmptyRange);
        TEST_CHECK_EQUAL(RangeType(-inf_val, -inf_val), kEmptyRange);
        TEST_CHECK_EQUAL(RangeType(+inf_val, -inf_val), kReverseRange);
        TEST_CHECK_EQUAL(RangeType(+inf_val, fin_val),  kReverseRange);
        TEST_CHECK_EQUAL(RangeType(fin_val, -inf_val),  kReverseRange);
        TEST_CHECK_EQUAL(RangeType(fin_val, -fin_val),  kReverseRange);
    }

    virtual void run() const
    {
        defaultToPDFTest();
        rangeTypeTest();
    }

} bcaux_Test;
