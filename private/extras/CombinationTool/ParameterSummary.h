#ifndef __PARAMETERSUMMARY__H
#define __PARAMETERSUMMARY__H

/*!
 * \class ParameterSummary
 * \brief A class for summarizing a single parameter
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 04.2010
 * This class is a container for summary information of a single
 * parameter. The information from a BAT 1-D histogram can be copied
 * automatically.
 */

/*
 * Copyright (C) 2010, Kevin Kroeninger.
 * All rights reserved.
 */

// ---------------------------------------------------------

#include <string>

class BCH1D;

// ---------------------------------------------------------
class ParameterSummary 
{
 public:
	
	/** \name Constructors and destructors */
	/* @{ */
	
	/**
	 * Default constructor. 
	 * @param name The parameter name.
	 */
	ParameterSummary(const char* name);

	/**
	 * Default destructor.
	 */ 
	~ParameterSummary(); 

	/* @} */

	/** \name Public member functions (get) */
	/* @{ */

	/** 
	 * Return the name of the parameter. 
	 */ 
	std::string GetName()
	{ return fName; }; 

	/**
	 * Return mode.
	 */ 
	double GetMode()
	{ return fMode; }; 

	/**
	 * Return global mode.
	 */ 
	double GetGlobalMode()
	{ return fGlobalMode; }; 

	/**
	 * Return mean. 
	 */ 
	double GetMean()
	{ return fMean; };

	/**
	 * Return median. 
	 */ 
	double GetMedian()
	{ return fMedian; };

	/**
	 * Return 5% quantile.
	 */ 
	double GetQuantile5()
	{ return fQuantile5; }; 

	/**
	 * Return 10% quantile.
	 */ 
	double GetQuantile10()
	{ return fQuantile10; }; 

	/**
	 * Return 16% quantile.
	 */ 
	double GetQuantile16()
	{ return fQuantile16; }; 

	/**
	 * Return 84% quantile.
	 */ 
	double GetQuantile84()
	{ return fQuantile84; }; 

	/**
	 * Return 90% quantile.
	 */ 
	double GetQuantile90()
	{ return fQuantile90; }; 

	/**
	 * Return 95% quantile.
	 */ 
	double GetQuantile95()
	{ return fQuantile95; }; 

	/**
	 * Return standard deviation.
	 */ 
	double GetRMS()
	{ return fRMS; }; 

	/**
	 * Return arbitrary buffer
	 */
	double GetBuffer()
	{ return fBuffer; };

	/* @} */

	/** \name Public member functions (set) */
	/* @{ */

	/**
	 * Set the mode. 
	 * @param mode The mode.
	 */ 
	void SetMode(double mode)
	{ fMode = mode; }; 

	/**
	 * Set the global mode. 
	 * @param mode The global mode.
	 */ 
	void SetGlobalMode(double mode)
	{ fGlobalMode = mode; }; 


	/**
	 * Set the mean. 
	 * @param mean The mean.
	 */ 
	void SetMean(double mean)
	{ fMean = mean; }; 

	/**
	 * Set the median. 
	 * @param median The median.
	 */ 
	void SetMedian(double median)
	{ fMedian = median; }; 

	/**
	 * Set the 5% quantile. 
	 * @param quantile The 5% quantile.
	 */ 
	void SetQuantile5(double quantile)
	{ fQuantile5 = quantile; }; 

	/**
	 * Set the 10% quantile. 
	 * @param quantile The 10% quantile.
	 */ 
	void SetQuantile10(double quantile)
	{ fQuantile10 = quantile; }; 

	/**
	 * Set the 16% quantile. 
	 * @param quantile The 16% quantile.
	 */ 
	void SetQuantile16(double quantile)
	{ fQuantile16 = quantile; }; 

	/**
	 * Set the 84% quantile. 
	 * @param quantile The 84% quantile.
	 */ 
	void SetQuantile84(double quantile)
	{ fQuantile84 = quantile; }; 

	/**
	 * Set the 90% quantile. 
	 * @param quantile The 90% quantile.
	 */ 
	void SetQuantile90(double quantile)
	{ fQuantile90 = quantile; }; 

	/**
	 * Set the 95% quantile. 
	 * @param quantile The 95% quantile.
	 */ 
	void SetQuantile95(double quantile)
	{ fQuantile95 = quantile; }; 

	/**
	 * Set the standard deviation quantile. 
	 * @param rms The standard deviation.
	 */ 
	void SetRMS(double rms)
	{ fRMS = rms; }; 

	/**
	 * Set the arbitrary buffer. 
	 * @param val The buffer value.
	 */ 
	void SetBuffer(double val)
	{ fBuffer = val; }; 

	/* @} */

	/** \name Public member functions (misc) */
	/* @{ */
	
	/**
	 * Copy the information from a histogram.
	 * @param hist A BAT 1-D histogram.
	 */ 
	void Summarize(BCH1D* hist); 

	/* @} */

 private:

	/**
	 * The parameter name.
	 */ 
	std::string fName; 

	/**
	 * The mode.
	 */ 
	double fMode;

	/**
	 * The global mode.
	 */ 
	double fGlobalMode;

	/**
	 * The mean.
	 */ 
	double fMean;

	/**
	 * The median.
	 */ 
	double fMedian;

	/**
	 * The 5% quantile. 
	 */
	double fQuantile5;

	/**
	 * The 10% quantile. 
	 */
	double fQuantile10;

	/**
	 * The 16% quantile. 
	 */
	double fQuantile16;

	/**
	 * The 84% quantile. 
	 */
	double fQuantile84;

	/**
	 * The 90% quantile. 
	 */
	double fQuantile90;

	/**
	 * The 95% quantile. 
	 */
	double fQuantile95;

	/**
	 * The standard deviation.
	 */
	double fRMS; 

	/**
	 * A buffer for additional information. 
	 */ 
	double fBuffer;
};

// ---------------------------------------------------------

#endif

