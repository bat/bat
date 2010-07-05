/*!
 * \class BAT::PerfTest2DFunction
 * \brief A performance test class for BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_PERFTEST2DFunction
#define BAT_PERFTEST2DFunction

#include <string>
#include <vector>

#include <TF2.h>

#include <include/PerfTest.h>
#include <include/PerfTestMCMC.h>

class PerfTest2DFunction : public PerfTestMCMC
{
	
 public:
	
	/** \name Constructors and destructors  */
	/* @{ */
	
	/** The default constructor */
	PerfTest2DFunction(std::string name = "unknown", TF2* func = 0);

	/** The default destructor */
	~PerfTest2DFunction();

	/* @} */

	/** Defines the subtests. */ 
	void DefineSubtests(); 

	/** Run after the test
	 * @return an error code. */ 
	int PostTest(); 

	/* @} */

	// inherited methods
	double LogAPrioriProbability(std::vector <double> parameters)
	{return 0;}

	double LogLikelihood(std::vector <double> parameters)
	{ return log(fFunction->Eval(parameters[0], parameters[1])); }

 private:
		
	/** The test function. */
	TF2* fFunction;
};

#endif

