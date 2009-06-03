#ifndef __BCBENCHMARKMCMC2D__H
#define __BCBENCHMARKMCMC2D__H

#include "BAT/BCModel.h"

#include <TF2.h>

class BCBenchmarkMCMC2D : public BCModel
{
public:
	BCBenchmarkMCMC2D(const char* name);
	~BCBenchmarkMCMC2D(){};

	double LogLikelihood(std::vector <double> parameters);

	double PerformTest(
			std::vector<double> parameters,
			int index,
			BCH2D * hist,
			bool flag_print = true,
			const char * filename = "test.eps");

	void SetTestFunction(TF2 * testfunction)
	{fTestFunction = testfunction;}

private:
	TF2 * fTestFunction;

};

#endif

