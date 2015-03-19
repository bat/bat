#ifndef __POL2ASYMM__H
#define __POL2ASYMM__H

/*
 * This class derives from Pol1Asymm. It describes a linear
 * correlation relation between measured points. Two parameters
 * are defined within the model, an offset and a slope.
 * The data are points (x,y) with an asymmetric uncertainty on y.
 * The uncertainty is assumed to be two half gaussians with different
 * widths.
 */

// ---------------------------------------------------------

#include <vector>

#include "Pol1Asymm.h"

// ---------------------------------------------------------

class Pol2Asymm : public Pol1Asymm
{
public:
    // constructor
    Pol2Asymm(const char* name);

    // destructor
    virtual ~Pol2Asymm();

    // define parameters of the model
    virtual void DefineParameters();

};

// ---------------------------------------------------------

#endif

