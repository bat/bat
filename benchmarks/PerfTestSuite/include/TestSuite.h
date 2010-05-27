/*!
	\class BAT::TestSuite
	\brief A test suite class for BAT
	
	Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
	All rights reserved.
	
	For the licensing terms see doc/COPYING.
*/ 

#ifndef BAT_TESTSUITE
#define BAT_TESTSUITE

#include <string>
#include <vector>

#include <include/PerfTest.h>

	class TestSuite
	{

	public:

		/** \name Constructors and destructors */
		/* @{ */

		/** The default constructor */
		TestSuite();
		
		/** The default destructor */
		~TestSuite();

		/* @} */
		/** \name Member functions (Set) */
		/* @{ */ 
		
		/* @} */ 
		/** \name Member functions (Get) */ 
		/* @{ */ 

		/** Get the number of tests which belong to this test. 
		 * @return the number of tests. */ 
		int GetNTests()
		{ return int(fTestContainer.size()); }; 
		
		/** Get the number of tests which belong to this test with
		 * specified status. 
		 * @param status the status code. 
		 * @return the number of tests. */
		int GetNTests(PerfSubTest::Status status); 

		/** Find a test by index
		 * @param index the index of the test.
		 * @return the subtest. */
		PerfTest * GetTest(int index)
		{ return fTestContainer.at(index); }; 

		/** Find a test by name
		 * @param name the name of the test.
		 * @return the test. */
		PerfTest * GetTest(std::string name);

		/* @} */
		/** \name Member functions (misc) */
		/* @{ */

		/** Define all tests here. */
		void DefineTests(); 

		/** Add a test. 
		 * @param test the test to be added.
		 * @return An error code. */
		int AddTest(PerfTest* test);

		/** Create the HTML output of the tests. 
		 * @param filename the name of the HTML file. */ 
		void PrintResultsHTML(std::string filename = "results.html"); 

		/** Create an ASCII output of the tests.
		 * @param filename the name of the ASCII file. */
		void PrintResultsASCII(std::string filename = "results.txt"); 

		/** Create an ASCII output of the tests.
		 * @param filename the name of the ASCII file. */
		void PrintResultsScreen(); 

		/** Run all tests.
		 * @return an error code. */
		int RunTests(); 
		
		/* @} */

	private:
			
		/** A container of tests which belong to the test suite. */
		std::vector <PerfTest *> fTestContainer; 
	};

#endif

