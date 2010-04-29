#ifndef __PARAMETERSUMMARY__H
#define __PARAMETERSUMMARY__H

#include <string>

class BCH1D;

// ---------------------------------------------------------
class ParameterSummary 
{
 public:
	
	ParameterSummary(const char* name);

	~ParameterSummary(); 

	// getters 

	std::string GetName()
	{ return fName; }; 

	double GetMode()
	{ return fMode; }; 

	double GetMean()
	{ return fMean; };

	double GetMedian()
	{ return fMedian; };

	double GetQuantile5()
	{ return fQuantile5; }; 

	double GetQuantile10()
	{ return fQuantile10; }; 

	double GetQuantile16()
	{ return fQuantile16; }; 

	double GetQuantile84()
	{ return fQuantile84; }; 

	double GetQuantile90()
	{ return fQuantile90; }; 

	double GetQuantile95()
	{ return fQuantile95; }; 

	double GetRMS()
	{ return fRMS; }; 

	// setters

	void SetMode(double mode)
	{ fMode = mode; }; 

	void SetMean(double mean)
	{ fMean = mean; }; 

	void SetMedian(double median)
	{ fMedian = median; }; 

	void SetQuantile5(double quantile)
	{ fQuantile5 = quantile; }; 

	void SetQuantile10(double quantile)
	{ fQuantile10 = quantile; }; 

	void SetQuantile16(double quantile)
	{ fQuantile16 = quantile; }; 

	void SetQuantile84(double quantile)
	{ fQuantile84 = quantile; }; 

	void SetQuantile90(double quantile)
	{ fQuantile90 = quantile; }; 

	void SetQuantile95(double quantile)
	{ fQuantile95 = quantile; }; 

	void SetRMS(double rms)
	{ fRMS = rms; }; 

	// misc
	
	void Summarize(BCH1D* hist); 

 private:

	std::string fName; 
	double fMode;
	double fMean;
	double fMedian;
	double fQuantile5;
	double fQuantile10;
	double fQuantile16;
	double fQuantile84;
	double fQuantile90;
	double fQuantile95;
	double fRMS; 
};

// ---------------------------------------------------------

#endif

