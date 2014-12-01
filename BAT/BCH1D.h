#ifndef __BCH1D__H
#define __BCH1D__H

/*!
 * \class BCH1D
 * \brief A class for handling 1D distributions.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \author Daniel Greenwald
 * \version 1.0
 * \date 08.2008
 * \detail This class contains a TH1D histogram and some additional
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

#include "BCHistogramBase.h"

class TH1;

// ---------------------------------------------------------

class BCH1D : public BCHistogramBase {

public:

	/** \name Enumerators */
	/** @{ */

	/**
	 * Enum for type of bands to be drawn on plot. */
	enum BCH1DBandType {
		kNoBands          = -1,
		kCentralInterval  = 0,
		kSmallestInterval = 1,
		kUpperLimit       = 2,
		kLowerLimit       = 3,
		kUserSpecified    = 4,
	};

	/** @} */

  /** \name Constructors and destructors */
  /** @{ */

  /**
   * The default constructor. */
  BCH1D(TH1 * hist = 0);

	/**
	 * Copy constructor. */
	BCH1D(const BCH1D & other);

  /**
   * The default destructor. */
  virtual ~BCH1D();

  /** @} */

  /** \name Member functions (get)  */
  /** @{ */

  /**
   * @return The median of the distribution. */
  double GetMedian()
  { return this -> GetQuantile(0.5); };

  /**
   * Returns the quantile of the distribution.
   * @param probability The probability.
   * @return The quantile of the distribution for the probability.
   * @see GetLimit(double probability) */
  double GetQuantile(double probablity);

  /**
   * Return the quantile of the distribution
   * @param probability The probability.
   * @return The quantile of the distribution for the probability.
   * @see GetQuantile(double probablity) */
  double GetLimit(double probability)
  { return this->GetQuantile(probability); };

	/**
	 * @return Band type. */
	BCH1DBandType GetBandType()
	{ return fBandType; }

	/**
	 * @return Number of quantiles to draw. */
	unsigned GetNQuantiles()
	{ return fNQuantiles; }

	/**
	 * @return Quantile line color. */
	int GetQuantileLineColor()
	{ return fQuantileLineColor; }

	/**
	 *@return whether to draw median. */
	bool GetDrawMedian() 
	{ return fDrawMedian; }
	
	/**
	 * @return whether to draw central 68% interval. */
	bool GetDrawCentral68()
	{ return fDrawCentral68; }

  /** @} */

  /** \name Member functions (set)  */
  /** @{ */

	/**
	 * Copy options from. */
	void CopyOptions(const BCH1D & other);

  /**
   * Sets the color scheme.
   * @param scheme the scheme index \n
   * 0 : black and white
   * 1 : yellow-green-red
   * 2 : blueish colors
   * 2 : redish colors
   * 2 : blueish colors
   */
  void SetColorScheme(BCHColorScheme scheme);

	using BCHistogramBase::SetGlobalMode;

  /**
   * Set global mode */
  void SetGlobalMode(double mode)
	{ SetGlobalMode(std::vector<double>(1,mode)); }
	
	/**
	 * Set band type. */
	void SetBandType(BCH1DBandType bt)
	{ fBandType = bt; }

	/**
	 * Set draw quantiles. 
	 * @param n N divisions of quantiles to draw, set to zero or one to disable drawing of quantiles. */
	void SetDrawQuantiles(unsigned n)
	{ fNQuantiles = n; }

	/**
	 * Set quantile line color.
	 * @param c Quantile line color. */
	void SetQuantileLineColor(int c)
	{ fQuantileLineColor = c; }

	/**
	 * Set drawing of median.
	 * @param flag Toggles drawing of median.
	 * @param central68 Toggles drawing of arrows for central 68% interval. (Automatically suppressed if median is suppressed.)*/
	void SetDrawMedian(bool flag=true, bool central68=true)
	{ fDrawMedian = flag; fDrawCentral68 = central68;}

  /** @} */

  /** \name Member functions (miscellaneous methods) */
  /** @{ */

	using BCHistogramBase::Draw;

  /**
   * Draw distribution into the active canvas.
   * @param options ROOT drawing options
   * @param intervals the intervals */
  void Draw(std::string options="", std::vector<double> intervals=std::vector<double>(0));

  /**
   * Calculate the minimal interval of the distribution containing a given content.
   * @param min calculated minimum of the interval
   * @param max calculated maximum of the interval
   * @param content content of the interval [default is .68]
   * @return the content of the histogram between min and max */
  double GetSmallestInterval(double & min, double & max, double content=.68);

  /**
   * Create a histogram with the smallest intervals of the original
   * histogram containing a certain level. The histogram is yellow.
   * @param level the level or content of the histogram
   * @return the histogram.
   */
  TH1D* GetSmallestIntervalHistogram(double level);

  /**
   * Return a vector of vectors containing information about the set of smallest
   * intervals. \n
   * 0 : x_min \n
   * 1 : x_max \n
   * 2 : relative height
   * 3 : local mode \n
   * 4 : relative area.
   * @param content The content of the smallest interval
   * @return the vector.
   */
	std::vector<std::vector<double> > GetSmallestIntervals(double content = 0.68);

  /**
   * Calculate integral of the distribution between min and max.
   * @param min lower boundary of the integrated interval
   * @param max upper boundary of the integrated interval
   * @return integral calculated as sum of BinContent*BinWidth */
  double IntegralWidth(double min, double max);

  /**
   * Get histogram with bins outside min, max band being zero. The
   * new histogram can have 2 more bins than the original one as the
   * bins where min and max fall into will be split in two (except for the
   * case when min and/or max are equal to some of the original bin
   * boundaries.
   * @param min lower boundary of the non-zero interval
   * @param max upper boundary of the non-zero interval
	 * @param name Name for new histogram; empty string (default) appends "subhist" to histogram name.
	 * @param preserve_range If true, preserves original histograms range, setting bins outside subhistogram range to zero.
   * @return new histogram which is nonzero only between min and max */
  TH1D* GetSubHisto(double min, double max, std::string name="", bool preserve_range=false);

  /** @} */

protected:

	/**
	 * Band type */
	BCH1DBandType fBandType;

	/**
	 * Number of quantiles to draw. */
	unsigned fNQuantiles;

	/**
	 * Quantile line color. */
	int fQuantileLineColor;

	/**
	 * flag for drawing median. */
	bool fDrawMedian;

	/**
	 * flag for darwing central 68% interval arrows. */
	bool fDrawCentral68;

};

// ---------------------------------------------------------

#endif
