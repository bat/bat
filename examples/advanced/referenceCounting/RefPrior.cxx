#include "RefPrior.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
RefPrior::RefPrior(double shape, double rate)
    : BCPrior(),
      fShape(shape),
      fRate(rate),
      fMaxK(1000),
      fPrecisionLimit(1e-4),
      fLogSqrtIZero(std::numeric_limits<double>::quiet_NaN())
{
}

// ---------------------------------------------------------
RefPrior::RefPrior(const RefPrior& other)
    : BCPrior(other),
      fShape(other.fShape),
      fRate(other.fRate),
      fMaxK(other.fMaxK),
      fPrecisionLimit(other.fPrecisionLimit),
      fLogSqrtIZero(other.fLogSqrtIZero)
{
}

// ---------------------------------------------------------
double RefPrior::LogSqrtI(double s) const
{
    double sum = 0;
    double inc = 1.;

    double t2 = F(s, 0);
    for (unsigned n = 0; n <= fMaxK; ++n) {
        double t1 = t2;
        t2 = F(s, n + 1);
        double dsum = t1 * t1 / t2;
        sum += dsum;
        inc = fabs( dsum / sum );

        // check if sum does not change by more than predefined precision
        if (inc < fPrecisionLimit)
            break;
    }

    return 0.5 * log(fabs(pow(fRate / (1 + fRate), fShape) * exp(-s) * sum - 1));
}

// ---------------------------------------------------------
double RefPrior::F(double s, unsigned k) const
{
    if (s <= 0) {
        if (k == 0)
            return 1;
        else
            return (fShape + k - 1) / (k * (1 + fRate)) * F(0, k - 1);
    }
    double logs = log(s);
    double logb = log(1 + fRate);

    double sum = 0;

    for (unsigned n = 0; n <= k; ++n)
        sum += exp(BCMath::LogBinomFactor(fShape + n - 1, n) + (k - n) * logs - n * logb - BCMath::ApproxLogFact(k - n));

    return sum;
}
