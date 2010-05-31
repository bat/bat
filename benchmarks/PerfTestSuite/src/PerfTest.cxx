/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "include/PerfTest.h"

#include <TCanvas.h>
#include <TH1D.h>
#include <TPostScript.h>

#include <iostream> 
#include <fstream> 

//______________________________________________________________________________
PerfTest::PerfTest(std::string name) 
	: fSubtestContainer(std::vector<PerfSubTest *>(0))
	, fCanvasContainer(std::vector<TCanvas *>(0))
	, fTestType(PerfTest::kUnknown)
	, fRealTime(0.)
	, fCpuTime(0.)
{
	fName = name;
}
	
//______________________________________________________________________________
PerfTest::~PerfTest()
{
}
	
//______________________________________________________________________________
std::string PerfTest::ToString(PerfSubTest::Status status)
{
	PerfSubTest st; 
	return st.ToString(status); 
}

//______________________________________________________________________________
std::string PerfTest::ToStringHTML(PerfSubTest::Status status)
{
	PerfSubTest st; 
	return st.ToStringHTML(status); 
}

//______________________________________________________________________________
int PerfTest::GetNSubtests(PerfSubTest::Status status)
{
	// get number of sub tests 
	int n = GetNSubtests(); 

	// initialize counter 
	int counter = 0; 

	// loop over all subtests and compare status
	for (int i = 0; i < n; ++i) {
		if (fSubtestContainer.at(i)->GetStatus() == status)
			counter++;
	}

	// return counter 
	return counter;
}
	
//______________________________________________________________________________
PerfSubTest::Status PerfTest::GetStatus()
{
	// get number of active sub tests 
	int n = GetNSubtests() - GetNSubtests(PerfSubTest::kOff); 

	// get number of successful sub tests
	int ngood = GetNSubtests(PerfSubTest::kGood); 

	// get number of flawed sub tests
	int nflawed = GetNSubtests(PerfSubTest::kFlawed); 

	// get number of failed sub tests
	int nbad = GetNSubtests(PerfSubTest::kBad); 

	// get number of failed sub tests
	int nfatal = GetNSubtests(PerfSubTest::kFatal); 

	// get number of unkown
	int nunknown = GetNSubtests(PerfSubTest::kUnknown); 
		
	// calculate overall status
	if (n == ngood)
		return PerfSubTest::kGood; 

	else if (nunknown > 0)
		return PerfSubTest::kUnknown;

	else if (nfatal > 0)
		return PerfSubTest::kFatal;

	else if (nbad > 0)
		return PerfSubTest::kBad;

	else if (nflawed > 0)
		return PerfSubTest::kFlawed;

	return PerfSubTest::kUnknown;
}
	
//______________________________________________________________________________
PerfSubTest * PerfTest::GetSubtest(std::string name)
{
	// get number of sub tests 
	int n = GetNSubtests(); 

	// loop over all subtests and compare status
	for (int i = 0; i < n; ++i) 
		{
			if (!name.compare(GetSubtest(i)->GetName()))
				return GetSubtest(i); 
		}

	return 0; 
}

//______________________________________________________________________________
TCanvas* PerfTest::GetCanvas(int index)
{
	int ncanvases = GetNCanvases(); 

	// check index
	if (index < 0 || index >= ncanvases)
		return 0;
	
	// return canvas pointer
	return fCanvasContainer.at(index);
}

//______________________________________________________________________________
int PerfTest::ReadResults()
{
	/*
	// open file
	std::fstream file; 
	file.open((fName+std::string(".tst")).c_str(), std::fstream::in); 

	// check if file is open
	if (!file.is_open())
	{
	std::cout << "Could not open file." << std::endl;
	return 0; 
	}

	// read data from file 
	std::string dummy_string; 
	double dummy_double;
	bool dummy_bool;
	int dummy_int; 

	file >> dummy_string; 
	file >> dummy_int; 

	// check name
	if (fName.compare(dummy_string)) { 
	std::cout << "Test name and name in file do not agree." << std::endl;
	file.close();
	return 0; 
	}

	// check number of subtests
	if (GetNSubtests() != dummy_int) {
	std::cout << "Number of subtests no consistent." << std::endl;
	file.close();
	return 0; 
	}
		
	// loop over all subtests and read results
	for (int i = 0; i < dummy_int; ++i) {

	// get subtest
	PerfSubTest * subtest = GetSubtest(i); 

	file >> dummy_string; 

	// check name
	if (!(subtest->GetName().compare(dummy_string)==0)) { 
	std::cout << "Subtest name and name in file do not agree." << std::endl; 
	file.close();
	return 0; 
	}

	file >> dummy_double; 
	subtest->SetTestValue(dummy_double);
	file >> dummy_double; 
	subtest->SetStatusRegion(PerfSubTest::kGood, dummy_double);
	file >> dummy_double; 
	subtest->SetStatusRegion(PerfSubTest::kFlawed, dummy_double);
	file >> dummy_double; 
	subtest->SetStatusRegion(PerfSubTest::kBad, dummy_double);
	file >> dummy_bool; 
	subtest->SetStatusUnknown(dummy_bool);
	file >> dummy_bool; 
	subtest->SetStatusOff(dummy_bool);
	}	

	// close file
	file.close(); 
	*/

	// no error 
	return 1; 
}

//______________________________________________________________________________
int PerfTest::WriteResults()
{
	// open file
	std::fstream file; 
	file.open((fName.data()+std::string(".tst")).c_str(), std::fstream::out); 

	// check if file is open
	if (!file.is_open())
		{
			std::cout << "Could not open file." << std::endl;
			return 0; 
		}

	// write to file 
	file << fName.data() << std::endl;
	file << GetNSubtests() << std::endl; 
		
	// get number of active sub tests 
	int n = GetNSubtests() - GetNSubtests(PerfSubTest::kOff); 
		
	// loop over all subtests and write results
	for (int i = 0; i < n; ++i) {
		file << fSubtestContainer.at(i)->GetName().c_str() << std::endl;
		file << fSubtestContainer.at(i)->GetTestValue() << std::endl;
		file << fSubtestContainer.at(i)->GetTargetValue() << std::endl;
		file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kGood) << std::endl;
		file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kFlawed) << std::endl;
		file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kBad) << std::endl;
		file << fSubtestContainer.at(i)->GetStatusRegion(PerfSubTest::kFatal) << std::endl;
		file << fSubtestContainer.at(i)->GetStatusUnknown() << std::endl;
		file << fSubtestContainer.at(i)->GetStatusOff() << std::endl;
	}

	// close file 
	file.close(); 

	// create postscript
	TPostScript * ps = new TPostScript((fName.data()+std::string(".ps")).c_str());

	// get number of canvases
	int nhist = GetNCanvases(); 

	ps->NewPage();

	// loop over histograms
	for (int i = 0; i < nhist; ++i) {
		// get canvas
		TCanvas* c = GetCanvas(i);
		
		// update post script
		c->Update();
		if (i != nhist-1)
			ps->NewPage();
		c->cd();
	}

	// close ps
	ps->Close();

	// no error 
	return 1; 
}

//______________________________________________________________________________
	
