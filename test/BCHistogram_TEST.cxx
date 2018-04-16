/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>

#include <BAT/BCH1D.h>

#include <TH1D.h>
#include <TMath.h>

using namespace test;

class BCHistogramTest :
    public TestCase
{
public:
    BCHistogramTest() :
        TestCase("BC histogram test")
    {
    }

    virtual void run() const
    {
        double eps = 1e-6;

        // create a histogram and fill it with Gaus(bin center|0, 1)
        int N = 10;
        int x = 5;
        TH1D h("h", "", (N + 1) * x, -1. * x, 1. * x);
        for (int b = 1; b <= h.GetNbinsX(); ++b)
            h.SetBinContent(b, TMath::Gaus(h.GetBinCenter(b)));

        // create BCH1D
        BCH1D bch(&h);

        double alpha = 0.68;

        /////////////////////////
        // check for one interval with proper values
        BCH1D::BCH1DSmallestInterval smallest_interval = bch.GetSmallestIntervals(alpha);
        TEST_CHECK_EQUAL(smallest_interval.intervals.size(), 1);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].xmin, -1.);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].xmax, +1.);
        TEST_CHECK_NEARLY_EQUAL(smallest_interval.intervals[0].mode, 0., eps);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].relative_height, 1.);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].relative_mass, 1.);

        /////////////////////////
        // set peak bin to be exactly equal to threshold for inclusion,
        // to test that still only one contiguous interval is found
        int peak_bin = bch.GetHistogram()->GetNbinsX() / 2 + 1;
        double peak_val = bch.GetHistogram()->GetBinContent(peak_bin);
        double lower_bound = bch.GetSmallestIntervalBounds(std::vector<double>(1, alpha))[0].second;
        bch.GetHistogram()->SetBinContent(peak_bin, lower_bound);
        for (int b = 1; b <= bch.GetHistogram()->GetNbinsX(); ++b)
            bch.GetHistogram()->SetBinContent(b, bch.GetHistogram()->GetBinContent(b) + (peak_val - lower_bound) / bch.GetHistogram()->GetNbinsX());

        smallest_interval = bch.GetSmallestIntervals(alpha);
        TEST_CHECK_EQUAL(smallest_interval.intervals.size(), 1);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].xmin, -1.);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].xmax, +1.);
        TEST_CHECK_NEARLY_EQUAL(smallest_interval.intervals[0].mode, 0., eps);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].relative_height, 1.);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].relative_mass, 1.);

        /////////////////////////
        // set peak to be below threshold
        // check for two intervals
        bch.GetHistogram()->SetBinContent(peak_bin, 0);
        for (int b = 1; b <= bch.GetHistogram()->GetNbinsX(); ++b)
            bch.GetHistogram()->SetBinContent(b, bch.GetHistogram()->GetBinContent(b) + lower_bound / bch.GetHistogram()->GetNbinsX());

        smallest_interval = bch.GetSmallestIntervals(alpha);
        // check there are two intervals
        TEST_CHECK_EQUAL(smallest_interval.intervals.size(), 2);
        // check below-peak interval
        TEST_CHECK_NEARLY_EQUAL(smallest_interval.intervals[0].xmin, -1. - bch.GetHistogram()->GetBinWidth(1), eps);
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].xmax, bch.GetHistogram()->GetXaxis()->GetBinUpEdge(peak_bin - 1));
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].mode, bch.GetHistogram()->GetBinCenter(peak_bin - 1));
        TEST_CHECK_EQUAL(smallest_interval.intervals[0].relative_height, 1.);
        TEST_CHECK_NEARLY_EQUAL(smallest_interval.intervals[0].relative_mass, 0.5, eps);
        // check above-peak interval
        TEST_CHECK_EQUAL(smallest_interval.intervals[1].xmin, bch.GetHistogram()->GetBinLowEdge(peak_bin + 1));
        TEST_CHECK_NEARLY_EQUAL(smallest_interval.intervals[1].xmax, +1. + bch.GetHistogram()->GetBinWidth(1), eps);
        TEST_CHECK_EQUAL(smallest_interval.intervals[1].mode, bch.GetHistogram()->GetBinCenter(peak_bin + 1));
        TEST_CHECK_EQUAL(smallest_interval.intervals[1].relative_height, 1.);
        TEST_CHECK_NEARLY_EQUAL(smallest_interval.intervals[1].relative_mass, 0.5, eps);

    }

} bchistogram_Test;
