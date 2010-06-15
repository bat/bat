/*!
 * \class BAT::PerfTest1DFunction
 * \brief A performance test class for BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_PERFTEST1DFUNCTION
#define BAT_PERFTEST1DFUNCTION

#include <string>
#include <vector>

#include <TF1.h>

#include <BAT/BCModel.h>

#include <include/PerfTest.h>

class PerfTest1DFunction : public PerfTest, public BCModel
{
	
 public:
	
	/** \name Constructors and destructors  */
	/* @{ */
	
	/** The default constructor */
	PerfTest1DFunction(std::string name = "unknown", TF1* func = 0);

	/** The default destructor */
	~PerfTest1DFunction();

	/* @} */

	/** Run the test. 
	 * @return an error code. */ 
	int Run(); 

	/** Defines the subtests. */ 
	void DefineSubtests(); 

	/** Writes the test to file. 
	 * @return an error code. */ 
	int WriteResults(); 

	/** Define precision settings. */ 
	void PrecisionSettings(PerfTest::Precision);

	/* @} */

	// inherited methods
	double LogAPrioriProbability(std::vector <double> parameters)
	{return 0;}

	double LogLikelihood(std::vector <double> parameters)
	{ return log(fFunction->Eval(parameters[0])); }

 private:
		
	/** The test function. */
	TF1* fFunction;
};

#endif

