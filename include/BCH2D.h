/*! \class BCH2D
 *  \brief A class for handling 2D distributions
 *
 * A class which contains a TH2D histogram and can be used for marginalized probabilties
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  D. Kollar, K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, 
 *          kevin.kroeninger *at* phys *dot* uni *minus* goettingen *dot* de 
 *
 * CREATED: 02.03.2007 
 * 
 * REVISION: 
 *
 * 02.03.2007  Kevin  * added comments and header\n
 * 22.05.2007  Kevin  * added nicer 2D plots including contours\n
 * 03.08.2007  Dano   * increase printout level of ROOT routines like Print()\n
 *
 * --------------------------------------------------------- 
 *
*/ 

// --------------------------------------------------------- 

#ifndef __BCH2D__H
#define __BCH2D__H

#include <vector.h> 

#include <TH1D.h> 
#include <TH2D.h> 

// --------------------------------------------------------- 

class BCH2D
{
  
 public:
  
	// constructors and destructor 

	/**
	 * The default constructor. 
	 */ 
	BCH2D(); 

	/** 
	 * The default destructor. 
	 */ 
	~BCH2D(); 

	// methods (get) 

	/**
	 * Returns the histogram. 
	 */ 
	TH2D * GetHistogram()
	{ return fHistogram; }; 

	/** 
	 * Returns the mean of the distribution 
	 */ 
	void GetMean(double& mean);

	/** 
	 * Returns the mode of the distribution
	 */ 
	void GetMode(double& mode);

	// methods (set) 

	/**
	 * Returns the histogram. 
	 */ 
	void SetHistogram(TH2D * hist)
	{ fHistogram = hist; };

	// methods 

	/** 
	 * Print 2-d histogram to file 
	 * @param filename The filename 
	 */ 
	void Print(char* filename, int options=0);

	/*
	 * Calculates the integral of the distribution as a function of the
	 * height.
	 */ 
	void CalculateIntegratedHistogram(); 

	/*
	 * Calculates the height below which the integrated probability has
	 * a certain value.
	 * @param p The integrated probability in the region below the height to be estimated. 
	 */ 
	double GetLevel(double p); 
  
 private: 

	/** 
	 * The 2-d histogram
	 */ 
	TH2D * fHistogram; 

	/**
	 * The integrated 2-d histogram 
	 */ 
	TH1D * fIntegratedHistogram; 

}; 

// --------------------------------------------------------- 

#endif 
