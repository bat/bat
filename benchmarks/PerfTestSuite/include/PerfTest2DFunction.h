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

#include <BAT/BCModel.h>

#include <include/PerfTest.h>

class PerfTest2DFunction : public PerfTest, public BCModel
{
	
 public:
	
	/** \name Constructors and destructors  */
	/* @{ */
	
	/** The default constructor */
	PerfTest2DFunction(std::string name = "unknown", TF2* func = 0);

	/** The default destructor */
	~PerfTest2DFunction();

	/* @} */

	/** Run the test. 
	 * @return an error code. */ 
	int Run(); 

	/** Defines the subtests. */ 
	void DefineSubtests(); 

	/* @} */

	// inherited methods
	double LogAPrioriProbability(std::vector <double> parameters)
	{return 0;}

	double LogLikelihood(std::vector <double> parameters)
	{ return log(fFunction->Eval(parameters[0], parameters[1])); }

	/** Writes the test to file. 
	 * @return an error code. */ 
	int WriteResults(); 

	/** Define precision settings. */ 
	void PrecisionSettings(PerfTest::Precision);

 private:
		
	/** The test function. */
	TF2* fFunction;
};

#endif

