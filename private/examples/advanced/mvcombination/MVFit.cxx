#include "MVFit.h"
#include "Observable.h"

#include <BAT/BCModel.h>
#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MVFit::MVFit() : MVCombination()
{
}

// ---------------------------------------------------------
MVFit::~MVFit()
{
}

// ---------------------------------------------------------
double MVFit::LogLikelihood(const std::vector<double> &parameters)
{
	std::vector<double> observables;
	
	for (int i = 0; i < fNObservables; ++i) 
		observables.push_back(CalculateObservable(i, parameters));
	
	return MVCombination::LogLikelihood(observables);
}

// ---------------------------------------------------------
void MVFit::AddObservable(std::string name, double min, double max)
{
	// check if observable exists already
  int index = GetIndexObservable(name);

	if (index >= 0)
		return;

	Observable* observable = new Observable();
 	observable->SetName(name);
 	observable->SetMinMax(min, max);
 	fObservables.push_back(observable);
	
  fNObservables++;
}

// ---------------------------------------------------------
