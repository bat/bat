/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMVCPhysicsModel.h"
#include "BCMVCObservable.h"

#include "../../BAT/BCModel.h"
#include "../../BAT/BCMath.h"
#include "../../BAT/BCLog.h"

// ---------------------------------------------------------
BCMVCPhysicsModel::BCMVCPhysicsModel() : BCMVCombination()
{
}

// ---------------------------------------------------------
BCMVCPhysicsModel::~BCMVCPhysicsModel()
{
}

// ---------------------------------------------------------
double BCMVCPhysicsModel::LogLikelihood(const std::vector<double> &parameters)
{
  std::vector<double> observables;

  for (int i = 0; i < fNObservables; ++i)
    observables.push_back(CalculateObservable(i, parameters));

  return BCMVCombination::LogLikelihood(observables);
}

// ---------------------------------------------------------
void BCMVCPhysicsModel::AddObservable(std::string name, double min, double max)
{
  // check if observable exists already
  int index = GetIndexObservable(name);

  if (index >= 0)
    return;

  BCMVCObservable* observable = new BCMVCObservable();
  observable->SetName(name);
  observable->SetMinMax(min, max);
  fObservables.push_back(observable);

  fNObservables++;
}

// ---------------------------------------------------------
