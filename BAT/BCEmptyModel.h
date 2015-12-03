#ifndef __BCEMPTYMODEL__H
#define __BCEMPTYMODEL__H

/*!
 * \class BCEmptyModel
 * \brief An empty model, used for reading in a chain
 * \author Daniel Greenwald
 * \version 1.0
 * \date 11.2014
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModel.h"

#include <string>
#include <vector>

// ---------------------------------------------------------

class BCEmptyModel : public BCModel
{
public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * A constructor.
     * @param name The name of the model */
    BCEmptyModel(const std::string& name = "model");

    /**
     * Read in MCMC constructor.
     * @param filename Path of file holding model.
     * @param name Name of model (file should contain TTree's [name]_mcmc and [name]_parameters.\n
     * if empty string is given, properly matching TTrees are searched for in the file.
     * @param loadObservables Flag for whether to load observables for file (true; default) or to let user respecify observables.*/
    BCEmptyModel(const std::string& filename, const std::string& name, bool loadObservables = true);

    /**
     * The default destructor. */
    virtual ~BCEmptyModel() {};

    /** @} */

    /**
     * Calculates natural logarithm of the likelihood.
     * @param params A set of parameter values
     * @return Natural logarithm of the likelihood */
    virtual double LogLikelihood(const std::vector<double>& /*params*/)
    { return 0; }

};
// ---------------------------------------------------------

#endif
