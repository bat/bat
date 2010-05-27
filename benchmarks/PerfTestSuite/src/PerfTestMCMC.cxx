/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "include/PerfTestMCMC.h"

namespace BAT {
	
	//______________________________________________________________________________
	PerfTestMCMC::PerfTestMCMC(std::string name) : PerfTest(name)
	{
		DefineSubtests(); 
	}
	
	//______________________________________________________________________________
	PerfTestMCMC::~PerfTestMCMC()
	{
	}
	
	//______________________________________________________________________________
	int PerfTestMCMC::Run()
	{
		// this function needs to be overloaded by the user 

		GetSubtest("chi2")->SetTestValue(0.3); 
		GetSubtest("mean")->SetTestValue(0.1); 
		GetSubtest("mode")->SetTestValue(1.4); 

		// no error 
		return 1; 
	}

	//______________________________________________________________________________
	void PerfTestMCMC::DefineSubtests()
	{
		PerfSubTest * subtest = new PerfSubTest("chi2"); 
		subtest->SetStatusRegion(PerfSubTest::kGood, 1.); 
		subtest->SetStatusRegion(PerfSubTest::kFlawed, 2.); 
		subtest->SetStatusRegion(PerfSubTest::kBad, 3.); 
		AddSubtest(subtest);

		subtest = new PerfSubTest("mean"); 
		subtest->SetStatusRegion(PerfSubTest::kGood, 0.); 
		subtest->SetStatusRegion(PerfSubTest::kFlawed, 1.); 
		subtest->SetStatusRegion(PerfSubTest::kBad, 2.); 
		AddSubtest(subtest);

		subtest = new PerfSubTest("mode"); 
		subtest->SetStatusRegion(PerfSubTest::kGood, 0.); 
		subtest->SetStatusRegion(PerfSubTest::kFlawed, 1.); 
		subtest->SetStatusRegion(PerfSubTest::kBad, 2.); 
		AddSubtest(subtest);
	}

	//______________________________________________________________________________
	
} // end namespace BAT
