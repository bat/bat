#ifndef __BCOBSERVABLESET__H
#define __BCOBSERVABLESET__H

/**
 * @class BCObservableSet
 * @brief Wrapper to allow access by name into list of BCObservable.
 * @author Frederik Beaujean
 * @author Daniel Greenwald
 * @note Observables are owned, and will be deleted by BCObservableSet.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
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

    /**
     * Copy constructor
     * @param share_pointers whether to share pointers in new observables*/
    BCObservableSet(const BCObservableSet& other, bool share_pointers = false);

    /*
     * Destructor */
    virtual ~BCObservableSet();

};
#endif
