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
BCEmptyModel::BCEmptyModel(const char* name)
    : BCModel(name)
{}

// ---------------------------------------------------------
BCEmptyModel::BCEmptyModel(const BCEmptyModel& model)
    : BCModel(model)
{}

// ---------------------------------------------------------
BCEmptyModel::BCEmptyModel(std::string filename, std::string name, bool reuseObservables)
    : BCModel(filename, name, reuseObservables)
{}

// ---------------------------------------------------------
BCEmptyModel::~BCEmptyModel()
{}

// ---------------------------------------------------------
double BCEmptyModel::LogLikelihood(const std::vector<double>& /*params*/)
{
    return 0;
}
