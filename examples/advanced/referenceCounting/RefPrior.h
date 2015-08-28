#ifndef __BAT__REFPRIOR__H
#define __BAT__REFPRIOR__H

#include <BAT/BCPrior.h>

#include <limits>

class RefPrior : public BCPrior
{
public:

    // constructor
    RefPrior(double shape = 1, double rate = 1);

    // copy constructor
    RefPrior(const RefPrior& other);

    // destructor
    ~RefPrior()
    { /* empty destructor */ }

    ////////////////////////////////////////
    // methods that must be overloaded from BCPrior:
    double GetLogPrior(double x)
    { return LogSqrtI(x) - LogSqrtIZero(); }

    BCPrior* Clone() const
    { return new RefPrior(*this); }

    bool IsValid() const
    { return (fShape > 0) && (fRate > 0); }
    ////////////////////////////////////////

    double LogSqrtI(double s) const;
    double F(double s, unsigned k) const;

    double LogSqrtIZero()
    { return (std::isfinite(fLogSqrtIZero)) ? fLogSqrtIZero : fLogSqrtIZero = LogSqrtI(0); }

    ////////////////////////////////////////
    // setters
    void SetShape(double shape)
    { fShape = shape; fLogSqrtIZero = std::numeric_limits<double>::quiet_NaN();}

    void SetRate(double rate)
    { fRate = rate; fLogSqrtIZero = std::numeric_limits<double>::quiet_NaN();}

    void SetMaxK(unsigned k)
    { fMaxK = k; fLogSqrtIZero = std::numeric_limits<double>::quiet_NaN();}

    void SetPrecisionLimit(double l)
    { fPrecisionLimit = l; fLogSqrtIZero = std::numeric_limits<double>::quiet_NaN();}
    ////////////////////////////////////////

    ////////////////////////////////////////
    // getters
    double GetShape() const
    { return fShape; }

    double GetRate() const
    { return fRate; }

    unsigned GetMaxK() const
    { return fMaxK; }

    double GetPrecisionLimit() const
    { return fPrecisionLimit; }
    ////////////////////////////////////////

protected:

    double fShape;              // "alpha"
    double fRate;               // "beta"

    unsigned fMaxK;
    double fPrecisionLimit;

    double fLogSqrtIZero;

};

#endif
