/*!
 * \class BAT::PerfTest
 * \brief A performance test class for BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_PERFTEST
#define BAT_PERFTEST

#include <string>
#include <vector>

#include <include/PerfSubTest.h>

class TH1D;
class TCanvas;

class PerfTest
{

 public:

	/** \name Enumerators  */
	/* @{ */
		 
	/** An enumerator for the test categories. */ 
	enum Status { kSampling }; 
		 
	/* @} */
	/** \name Constructors and destructors  */
	/* @{ */
		 
	/** The default constructor */
	PerfTest(std::string name = "unknown");

	/** The default destructor */
	~PerfTest();

	/* @} */
	/** \name Member functions (Set)  */
	/* @{ */

	/* @} */
	/** \name Member functions (Get)  */
	/* @{ */

	/** Return the name of the test. 
	 * @return the name of the test. */ 
	std::string GetName()
		{ return fName; }; 

	/** Get the number of subtests which belong to this test. 
	 * @return the number of subtests. */ 
	int GetNSubtests()
	{ return int(fSubtestContainer.size()); }; 

	/** Get the number of subtests which belong to this test with
	 * specified status. 
	 * @param status the status code. 
	 * @return the number of subtests. */ 		
	int GetNSubtests(PerfSubTest::Status status); 

	/** Get the number of canvases.
	 * @return the number of canvases. */
	int GetNCanvases() 
	{ return int(fCanvasContainer.size()); };
	 

	/** Calculate and return the overall status of the test. 
	 * @return a status code. */ 
	PerfSubTest::Status GetStatus(); 

	/** Return the current status as a string. 
	 * @return a string. */ 
	std::string GetStatusString()
		{ return ToString(GetStatus()); }; 

	/** Find a subtest by index
	 * @param index the index of the subtest.
	 * @return the subtest. */ 
	PerfSubTest * GetSubtest(int index)
	{ return fSubtestContainer.at(index); }; 

	/** Find a subtest by name
	 * @param name the name of the subtest.
	 * @return the subtest. */ 
	PerfSubTest * GetSubtest(std::string name);

	/** Return a canvas from the container.
	 * @param index the canvas index. 
	 * @return the canvas. */
	TCanvas* GetCanvas(int index); 

	/* @} */
	/** \name Member functions (misc)  */
	/* @{ */

	/** Return the status code as a string. 
	 * @param status the status code. 
	 * @return a string. */ 
	std::string ToString(PerfSubTest::Status status); 

	/** Add a subtest to the container. 
	 * @param a subtest. */ 
	void AddSubtest(PerfSubTest * test)
	{ fSubtestContainer.push_back(test); }; 

	/** Add a canvas to the container. 
	 * @param hist a canvas. */
	void AddCanvas(TCanvas* canvas)
	{ fCanvasContainer.push_back(canvas); }; 

	/** Read test results from file. 
	 * @return an error code. */ 
	int ReadResults(); 

	/** Writes the test to file. 
	 * @return an error code. */ 
	int WriteResults(); 

	/** Run the test. 
	 * @return an error code. */ 
	virtual int Run() = 0; 

	/** Defines the subtests. */ 
	virtual void DefineSubtests() = 0; 

	/* @} */

 private:

	/** A container of subtest which belong to the test. */ 
	std::vector <PerfSubTest *> fSubtestContainer; 

	/** A container of canvases for the test. */
	std::vector <TCanvas*> fCanvasContainer;
	
	/** The name of the test. */ 
	std::string fName; 
};

#endif

