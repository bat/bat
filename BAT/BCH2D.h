#ifndef __BCH2D__H
#define __BCH2D__H

/*!
 * \class BCH2D
 * \brief  A class for handling 2D distributions.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Daniel Greenwald
 * \version 1.0
 * \date 08.2008
 * \detail This class contains a TH2D histogram and some additional
 * functions. It is used for marginalized distributions.
 */

/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>
#include <string>

// ROOT classes

#include "BCHistogramBase.h"

class TH2;
class TGraph;
class TObject;

// ---------------------------------------------------------

class BCH2D : public BCHistogramBase {
	
public:

	/** \name Enumerators */
	/** @{ */

	/**
	 * Enum for type of bands to draw on plot. */
	enum BCH2DBandType {
		kNoBands          = -1,
		kSmallestInterval = 0
		// kUpperLimitOfXGivenY = 1,
		// kUpperLimitOfYGivenX = 2,
		// kLowerLimitOfXGivenY = 3,
		// kLowerLimitOfYGivenX = 4
	};

	/**
	 * Enum for type of profile. */
	enum BCH2DProfileType {
		kProfileMean   = 0,
		kProfileMedian = 1,
		kProfileMode   = 2
	};

	/**
	 * Enum for axis to profile. */
	enum BCH2DProfileAxis {
		kProfileX = 0,
		kProfileY = 1
	};

	/** @} */

	/** \name Constructor and destructors */
	/** @{ */
	
	/**
	 * The complete constructor. */
	BCH2D(TH2 * h = 0);

	/**
	 * Copy constuctor. */
	BCH2D(const BCH2D & other);
	
	/**
	 * The default destructor. */
	~BCH2D();

	/** @} */
	/** \name Member functions (get)  */
	/** @{ */
	
	/**
	 * Return the TH2D histogram
	 * @return The TH2D histogram. */
	TH2 * GetHistogram()
	{ return (TH2*)fHistogram; }

	/**
	 * @return Band type. */
	BCH2DBandType GetBandType()
	{ return fBandType; }

	/**
	 * @return flag for plotting log on z axis. */
	bool GetLogz()
	{ return fLogz; }

	/**
	 * @return flag for drawing x profile. */
	bool GetDrawProfileX() 
	{ return fDrawProfileX; }

	/** 
	 * @return profile type for x profile. */
	BCH2DProfileType GetProfileXType()
	{ return fProfileXType; }

	/**
	 * @return line color for x profile. */
	int GetProfileXLineColor()
	{ return fProfileXLineColor; }

	/**
	 * @return line style for x profile. */
	int GetProfileXLineStyle()
	{ return fProfileXLineStyle; }

	/**
	 * @return flag for drawing y profile. */
	bool GetDrawProfileY() 
	{ return fDrawProfileY; }

	/** 
	 * @return profile type for y profile. */
	BCH2DProfileType GetProfileYType()
	{ return fProfileYType; }

	/**
	 * @return line color for y profile. */
	int GetProfileYLineColor()
	{ return fProfileYLineColor; }

	/**
	 * @return line style for y profile. */
	int GetProfileYLineStyle()
	{ return fProfileYLineStyle; }

	/** @} */


	/** \name Member functions (set)  */
	/** @{ */

	/**
	 * copy options from */
	void CopyOptions(const BCH2D & other);

	using BCHistogramBase::SetGlobalMode;

	/**
	 * Set global mode.
	 * @param x Global mode in x.
	 * @param y Global mode in y. */
	void SetGlobalMode(double x, double y)
	{ std::vector<double> m(1,x); m.push_back(y); SetGlobalMode(m); }

	/** Set band type. */
	void SetBandType(BCH2DBandType bt)
	{ fBandType = bt; }

	/**
	 * Sets drawing of z axis in log. */
	void SetLogz(bool flag=true)
	{ fLogz = true; }

	/**
	 * Set drawing of x profile. */
	void SetDrawProfileX(bool flag=true)
	{ fDrawProfileX = flag; }

	/**
	 * Set profile type of x profile. */
	void SetProfileXType(BCH2DProfileType pt)
	{ fProfileXType = pt; }

	/**
	 * Set line color of x profile. */
	void SetProfileXLineColor(int c)
	{ fProfileXLineColor = c; }
	
	/**
	 * Set line style of x profile. */
	void SetProfileXLineStyle(int s)
	{ fProfileXLineStyle = s; }

	/**
	 * Set drawing of y profile. */
	void SetDrawProfileY(bool flag=true)
	{ fDrawProfileY = flag; }

	/**
	 * Set profile type of y profile. */
	void SetProfileYType(BCH2DProfileType pt)
	{ fProfileYType = pt; }

	/**
	 * Set line color of y profile. */
	void SetProfileYLineColor(int c)
	{ fProfileYLineColor = c; }
	
	/**
	 * Set line style of y profile. */
	void SetProfileYLineStyle(int s)
	{ fProfileYLineStyle = s; }

	/** @} */


	/** \name Member functions (miscellaneous methods) */
	/** @{ */

	/**
	 * Draw distribution into the active canvas.
	 * @param options Drawing options: \n
	 * BTf : band type a filled area [default] \n
	 * BTc : band type is a contour \n
	 * mean : draw mean value and standard deviation [default] \n
	 * lmode : draw global mode [default] \n
	 * profilex : draw the profile line vs. x using the mode \n
	 * profiley : draw the profile line vs. y using the mode \n
	 * @param intervals the intervals for the bands
	 */
	void Draw(std::string options="BTfB3CS1meangmodelmode", std::vector<double> intervals=std::vector<double>(0));

	/**
	 * Calculates the integral of the distribution as a function of the
	 * height. */
	void CalculateIntegratedHistogram();

	/**
	 * Calculates the height below which the integrated probability has
	 * a certain value.
	 * @param p The integrated probability in the region below the height to be estimated. */
	double GetLevel(double p);

	/**
	 * Calculate the smallest area over which the integral of the function is p
	 * @param p The integrated probability in the region below the height to be estimated. */
	double GetArea(double p);

	/**
	 * Returns the number of intervals as a function of x
	 * @param h The histogram.
	 * @param nfoundmax The maximum number of intervals.
	 * @return A vector containing the number of intervals for all bins in x. */
	std::vector<int> GetNIntervalsY(TH2* h, int &nfoundmax);

	/**
	 * Return a graph of the profile along x or y. The profile is
	 * calculated by scanning through the one axis, e.g., y, and
	 * finding the mode, mean, median, etc. with respect to the
	 * other axis, e.g. x.
	 * @param axis x-axis (0) or y-axis (1)
	 * @param pt Type of profile to construct. */
	TGraph* CalculateProfileGraph(BCH2DProfileAxis axis, BCH2DProfileType pt=kProfileMean);

	/**
	 * Draw the profiles along x and y */
	void DrawProfileGraphs();

	/** @} */
	
protected:

	static bool Compare(std::pair<double,double> x, std::pair<double,double> y)
	{ return x.first < y.first; }

	/**
	 * The integrated 2D histogram */
	TH1D* fIntegratedHistogram;

	/**
	 * band type. */
	BCH2DBandType fBandType;

	/**
	 * flag for drawing log z. */
	bool fLogz;

	/**
	 * flag for drawing x profile. */
	bool fDrawProfileX;
	
	/**
	 * profile type of x profile. */
	BCH2DProfileType fProfileXType;

	/**
	 * color of x profile. */
	int fProfileXLineColor;

	/**
	 * line style of x profile. */
	int fProfileXLineStyle;

	/**
	 * flag for drawing y profile. */
	bool fDrawProfileY;

	/**
	 * profile type of y profile. */
	BCH2DProfileType fProfileYType;
	
	/**
	 * color of y profile. */
	int fProfileYLineColor;

	/**
	 * line style of y profile. */
	int fProfileYLineStyle;

};

// ---------------------------------------------------------

#endif
