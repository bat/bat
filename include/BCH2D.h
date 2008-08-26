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

#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>

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
	 * Set the histogram.
	 */
	void SetHistogram(TH2D * hist)
		{ fHistogram = hist; };

	/**
	 * Set global mode.
	 */
	void SetGlobalMode(double mode[2])
		{ fMode[0]=mode[0]; fMode[1]=mode[1]; fModeFlag=1; };

	// methods 

	/**
	 * Print 2-d histogram to file
	 * @param filename The filename
	 * @param ww canvas size in pixels along X
	 * @param ww canvas size in pixels along Y
	 * If ww and wh are set to 0, default ROOT canvas size is used.
	 * For explanation of the parameter options see the Draw() method.
	 */
	void Print(const char * filename, int options=0, int ww=0, int wh=0);

	/**
	 * Draw 2-d distribution into the active canvas
	 * @param options explanation to come
	 */
	void Draw(int options=0);

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

	/**
	 *
	 */
	std::vector <int> GetNIntervalsY(TH2D * h, int &nfoundmax);

	/**
	 *
	 */
	TGraph * GetLowestBandGraph(TH2D * h, std::vector<int> nint);
	TGraph * GetLowestBandGraph(TH2D * h);


	std::vector <double> GetLevelBoundary(TH2D * h, double level);
	TGraph * GetBandGraph(TH2D * h , double level1, double level2);

	TGraph ** GetBandGraphs(TH2D * h);
	TGraph ** GetBandGraphs(TH2D * h, int &n);

 private:

	/**
	 * The 2-d histogram
	 */
	TH2D * fHistogram;

	/**
	 * The integrated 2-d histogram
	 */
	TH1D * fIntegratedHistogram;

	/**
	 * Global mode
	 */
	double fMode[2];

	/**
	 * "Is there a global mode?" flag
	 */
	int fModeFlag;

}; 

// --------------------------------------------------------- 

#endif 
