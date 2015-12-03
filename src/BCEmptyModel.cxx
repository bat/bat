/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCEmptyModel.h"

// ---------------------------------------------------------
BCEmptyModel::BCEmptyModel(const std::string& name)
    : BCModel(name)
{
}

// ---------------------------------------------------------
BCEmptyModel::BCEmptyModel(const std::string& filename, const std::string& name, bool reuseObservables)
    : BCModel(filename, name, reuseObservables)
{
}

