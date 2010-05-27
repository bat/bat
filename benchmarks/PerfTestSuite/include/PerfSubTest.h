/*!
 * \class BAT::PerfSubTest
 * \brief A result class for performance tests in BAT
 */

/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef BAT_SUBTEST
#define BAT_SUBTEST

#include <string>
#include <vector>

class PerfSubTest
{

 public:

	/** \name Enumerators  */
	/* @{ */

	/** An enumerator for the status of a test. */ 
	enum Status { kGood, kFlawed, kBad, kFatal, kUnknown, kOff }; 
		 
	/* @} */
	/** \name Constructors and destructors  */
	/* @{ */

	/** The default constructor */
	PerfSubTest(std::string name = "unknown");

	/** The default destructor */
	~PerfSubTest();

	/* @} */
	/** \name Member functions (Set)  */
	/* @{ */

	/** Set the name of the subtest. 
	 * @param name the name of the subtest. */ 
	void SetName(std::string name)
	{ fName = name; }; 

	/** Set the test value range for the status. 
	 * @param status the status code.
	 * @param delta the allowed range around the target value. */
	void SetStatusRegion(PerfSubTest::Status status, double delta);

	/** Set the flag for status "unknown". 
	 * @param flag the flag. */ 
	void SetStatusUnknown(bool flag) 
	{ fStatusUnknown = flag; }; 

	/** Set the flag for status "off". 
	 * @param flag the flag. */ 
	void SetStatusOff(bool flag) 
	{ fStatusOff = flag; }; 

	/** Set the test value. 
	 * @param value the value returned from the test. */ 
	void SetTestValue(double value)
	{ fTestValue = value; }; 

	/** Set the aim value.
	 * @param value the value aimed at. */ 
	void SetTargetValue(double value) 
	{ fTargetValue = value; }; 

	/* @} */
	/** \name Member functions (Get)  */
	/* @{ */

	/** Return the name of the subtest. 
	 * @return a string. */
	std::string GetName()
		{ return fName; }; 

	/** Calculate and return the overall status of the subtest. 
	 * @return a status code. */ 
	PerfSubTest::Status GetStatus(); 

	/** Return the current status as a string. 
	 * @return a string. */ 
	std::string GetStatusString()
		{ return ToString(GetStatus()); }; 

	/** Return the test value range of the region for the status "good",
	 * "flawed" and "bad".
	 * @param status the status code.
	 * @return the minimum value. */ 
	double GetStatusRegion(PerfSubTest::Status status); 

	/** Return the current test value. 
	 * @return the test value. */ 
	double GetTestValue()
	{ return fTestValue; }; 

	/** Return the target value. 
	 * @return the target value. */ 
	double GetTargetValue()
	{ return fTargetValue; }; 

	/** Return flag for unknown status. 
	 * @return a flag. */ 
	bool GetStatusUnknown()
	{ return fStatusUnknown; }; 

	/** Return flag for off status. 
	 * @return a flag. */ 
	bool GetStatusOff()
	{ return fStatusOff; }; 

	/* @} */
	/** \name Member functions (misc)  */
	/* @{ */

	/** Return the status code as a string. 
	 * @param status the status code. 
	 * @return a string. */ 
	std::string ToString(PerfSubTest::Status status); 

	/* @} */

 private:

	/** The name of the subtest. */ 
	std::string fName; 

	/** The test value. */ 
	double fTestValue; 

	/** The target value. */ 
	double fTargetValue; 

	/** The test value range for the status. */
	std::vector <double> fStatusRegion; 

	/** A flag for status "unknown". */ 
	bool fStatusUnknown; 

	/** A flag for status "off". */ 
	bool fStatusOff; 
};

#endif
