/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "include/TestSuite.h"
#include "include/PerfTest1DFunction.h" 

#include <iostream> 

//______________________________________________________________________________
TestSuite::TestSuite()
	: fTestContainer(std::vector<PerfTest *>(0)) 
{
	DefineTests(); 
}
	
//______________________________________________________________________________
TestSuite::~TestSuite()
{
}
	
//______________________________________________________________________________
int TestSuite::GetNTests(PerfSubTest::Status status)
{
	// get number of sub tests 
	int n = GetNTests(); 

	// initialize counter 
	int counter = 0; 

	// loop over all tests and compare status
	for (int i = 0; i < n; ++i) {
		if (fTestContainer.at(i) -> GetStatus() == status)
			counter++;
	}
		
	// return counter 
	return counter;
}
	
//______________________________________________________________________________
PerfTest * TestSuite::GetTest(std::string name)
{
	// get number of sub tests 
	int n = GetNTests(); 
		
	// loop over all subtests and compare status
	for (int i = 0; i < n; ++i) 
		{
			if (!name.compare(GetTest(i) -> GetName()))
				return GetTest(i); 
		}

	return 0;
}

//______________________________________________________________________________
int TestSuite::AddTest(PerfTest* test)
{
	// add test to container
	fTestContainer.push_back(test);

	// no error
	return 1;
}

//______________________________________________________________________________
int TestSuite::RunTests()
{
	// get number of sub tests 
	int n = GetNTests(); 
		
	// initialize error code
	int err = 1; 

	// loop over all tests, run and write the output to file 
	for (int i = 0; i < n; ++i) {
		err *= fTestContainer.at(i) -> Run(); 
		err *= fTestContainer.at(i) -> WriteResults(); 
	}
		
	// return error code
	return err;
}

//______________________________________________________________________________
void TestSuite::PrintResultsScreen()
{
	std::cout << std::endl; 
	std::cout << " Test results: " << std::endl << std::endl; 

	std::cout << " Overview: " << std::endl; 
	std::cout << " Number of tests               : " << GetNTests() << std::endl; 
	std::cout << " Number of successful tests    : " << GetNTests(PerfSubTest::kGood) << std::endl; 
	std::cout << " Number of flawed tests        : " << GetNTests(PerfSubTest::kFlawed) << std::endl; 
	std::cout << " Number of bad tests           : " << GetNTests(PerfSubTest::kBad) << std::endl; 
	std::cout << " Number of fatal tests         : " << GetNTests(PerfSubTest::kFatal) << std::endl; 
	std::cout << " Number of tests unkown status : " << GetNTests(PerfSubTest::kUnknown) << std::endl; 
	std::cout << std::endl; 

	// loop over tests
	std::cout << " Details: " << std::endl; 
	int n = GetNTests(); 
	for (int i = 0; i < n; ++i){
		std::cout << " Test \"" << (GetTest(i) -> GetName()).data() << "\" : " << GetTest(i) -> GetStatusString().data() << std::endl; 

		// loop over subtests
		int nsub = GetTest(i) -> GetNSubtests(); 
		for (int j = 0; j < nsub; ++j)
			std::cout << "  Subtest \"" << (GetTest(i) -> GetSubtest(j))->GetName() << "\" : " << (GetTest(i) -> GetSubtest(j))->GetStatusString().data() << std::endl; 
		std::cout << std::endl; 
	}
		
}

//______________________________________________________________________________
void TestSuite::DefineTests()
{
	// ...
}
	
//______________________________________________________________________________
	
