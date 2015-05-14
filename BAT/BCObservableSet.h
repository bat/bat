#ifndef __BCOBSERVABLESET__H
#define __BCOBSERVABLESET__H

/**
 * @class BCObservableSet Wrapper to allow access by name into list of BCObservable.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Observables are owned, and will be deleted by BCObservableSet.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCObservable.h"
#include "BCVariableSet.h"

// ---------------------------------------------------------

class BCObservableSet : public BCVariableSet<BCObservable>
{
public:

    /**
     * Constructor */
    BCObservableSet();

    /*
     * Destructor */
    virtual ~BCObservableSet();

};
#endif
