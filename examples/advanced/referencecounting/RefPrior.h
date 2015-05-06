#ifndef __BAT__REFPRIOR__H
#define __BAT__REFPRIOR__H

#include <BAT/BCPrior.h>

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
    double GetLogPrior(double x) const
    { return LogSqrtI(x) - LogSqrtI(0); }

    BCPrior* Clone() const
    { return new RefPrior(*this); }

    bool IsValid() const
    { return (fShape > 0) && (fRate > 0); }
    ////////////////////////////////////////

    double LogSqrtI(double s) const;
    double F(double s, unsigned k) const;

    ////////////////////////////////////////
    // setters
    void SetShape(double shape)
    { fShape = shape; }

    void SetRate(double rate)
    { fRate = rate; }

    void SetMaxK(unsigned k)
    { fMaxK = k; }

    void SetPrecisionLimit(double l)
    { fPrecisionLimit = l; }
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

};

#endif
