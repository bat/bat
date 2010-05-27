/*!
 * \class BAT::PerfTestMCMC
 * \brief A performance test class for BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_PERFTESTMCMC
#define BAT_PERFTESTMCMC

#include <string>
#include <vector>

#include <include/PerfTest.h>

namespace BAT
{

	class PerfTestMCMC : public PerfTest
	{

	public:

		/** \name Constructors and destructors  */
		/* @{ */
		 
		/** The default constructor */
		PerfTestMCMC(std::string name = "unknown");

		/** The default destructor */
		~PerfTestMCMC();

		/* @} */

		/** Run the test. 
		 * @return an error code. */ 
		int Run(); 

		/** Defines the subtests. */ 
		void DefineSubtests(); 

		/* @} */
	};

} // namespace BAT

#endif

