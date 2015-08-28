#ifndef __BAT__GammaDistPrior__H
#define __BAT__GammaDistPrior__H

#include <BAT/BCPrior.h>

class GammaDistPrior : public BCPrior
{
public:

    // constructor
    GammaDistPrior(double shape = 1, double rate = 1);

    // copy constructor
    GammaDistPrior(const GammaDistPrior& other);

    // destructor
    ~GammaDistPrior()
    { /* empty destructor */ }

    ////////////////////////////////////////
    // methods that must be overloaded from BCPrior:
    double GetLogPrior(double x);

    BCPrior* Clone() const
    { return new GammaDistPrior(*this); }

    bool IsValid() const
    { return (fShape > 0) && (fRate > 0); }
    ////////////////////////////////////////

    void SetShape(double shape)
    { fShape = shape; }

    void SetRate(double rate)
    { fRate = rate; }

    double GetShape() const
    { return fShape; }

    double GetRate() const
    { return fRate; }

protected:

    double fShape;
    double fRate;

};

#endif
