#ifndef __BCH1D__H
#define __BCH1D__H

/*!
 * \class BCH1D
 * \brief A class for handling 1D distributions.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class contains a TH1D histogram and some additional
 * functions. It is used for marginalized distributions.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>

#include <TH1.h>

// ---------------------------------------------------------

class BCH1D
{
 public:

  /** \name Constructors and destructors */
  /** @{ */

  /**
   * The default constructor. */
  BCH1D(TH1D * hist = 0);

  /**
   * The default destructor. */
  ~BCH1D();

  /** @} */

  /** \name Member functions (get)  */
  /** @{ */

  /**
   * @return The one-dimensional histogram. */
  TH1D* GetHistogram()
    { return fHistogram; };

  /**
   * @return The mean of the distribution. */
  double GetMean()
  { return fHistogram -> GetMean(); };

  /**
   * @return The mode of the distribution. */
  double GetMode();

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
   * @return The RMS of the distribution. */
  double GetRMS()
  { return fHistogram->GetRMS(); };

  /**
   * @return The standard deviation of the distribution. */
  double GetSTD()
  { return fHistogram->GetRMS(); };

  /**
   * @return The variance of the distribution. */
  double GetVariance()
  { return (GetSTD()*GetSTD()); };

  /**
   * @return The skew of the distribution. */
  double GetSkew()
  { return fHistogram->GetSkewness(); };

  /**
   * @return The STD of the distribution. */
  double GetKurtosis()
  { return fHistogram->GetKurtosis(); };

  /**
   * Returns the integral of distribution the between two values.
   * @param valuemin The value from which the intergration is done.
   * @param valuemax The value up to which the intergration is done.
   * @return The integral. */
  double GetIntegral(double valuemin, double valuemax);

  /**
   * Returns the p-value.
   * Returns the integral from 0 to the probability.
   * @param probability Upper limit of integration.
   * @return The p-value. */
  double GetPValue(double probability);

  /**
   * Returns a color of the current color scheme.
   * @param index the color index
   * @return the color number. */
  int GetColor(int index)
  { return fColors.at(index); };

  /** @} */

  /** \name Member functions (set)  */
  /** @{ */

  /**
   * Sets the color scheme.
   * @param scheme the scheme index \n
   * 0 : black and white
   * 1 : yellow-green-red
   * 2 : blueish colors
   * 2 : redish colors
   * 2 : blueish colors
   */
  void SetColorScheme(int scheme);

  /**
   * Sets the histogram. */
  void SetHistogram(TH1D * hist)
  { fHistogram = hist; };

  /**
   * Set default probability limits. Allowed values are between 68%
   * and 100%. The default value is 95%. */
  void SetDefaultCLLimit(double limit);

  /**
   * Set global mode */
  void SetGlobalMode(double mode)
  { fMode=mode;
    fModeFlag=1; };

  /** @} */

  /** \name Member functions (miscellaneous methods) */
  /** @{ */

  /**
   * Print distribution into a PostScript file.
   * @param filename Output filename
   * @param option the draw options (see Draw()), plus \n
   * logx : draw x-axis in log-scale \n
   * logy : draw y-axis in log-scale \n
   * R : rescale canvas to have a squared histogram
   * @param intervals the intervals for the bands
   * @param ww canvas size in pixels along X
   * @param ww canvas size in pixels along Y
   * If ww and wh are set to 0, default ROOT canvas size is used.
   * For explanation of parameters options and ovalue look at BCH1D::Draw()
   * method. */
  void Print(const char * filename, std::string options="BTsiB3CS1D0Lmeanmode", std::vector<double> intervals=std::vector<double>(0), int ww=0, int wh=0);

  /**
   *Print distribution into a PostScript file.
   * @param filename Output filename
   * @param option the draw options, @see Print(const char * filename, std::string options="BTsiB3CS1D0Lmeanmode", std::vector<double> intervals=std::vector<double>(0), int ww=0, int wh=0)
   * @param interval an upper or lower limit
   * @param ww canvas size in pixels along X
   * @param ww canvas size in pixels along Y
   * @see Print(const char * filename, std::string options="BTsiB3CS1D0Lmeanmode", std::vector<double> intervals=std::vector<double>(0), int ww=0, int wh=0)
   */
  void Print(const char * filename, std::string options, double interval, int ww=0, int wh=0);

  /**
   * Draw distribution into the active canvas.
   * @param options Drawing options: \n
   * BTci : band type is central interval [default] \n
   * BTsi : band type is/are smallest interval(s) \n
   * BTul : band type is upper limit \n
   * BTll : band type is lower limit \n
   * B1 : draw one band between values specified in intervals [default] \n
   * B2 : draw two bands between values specified in intervals \n
   * B3 : draw three bands between values specified in intervals \n
   * D0 : draw histogram [default] \n
   * D1 : draw smooth curve \n
   * CS0 : choose color scheme 0 (B&W) \n
   * CS1 : choose color scheme 1 (green/yellow/red) [default] \n
   * CS2 : choose color scheme 2 (blueish colors) \n
   * CS3 : choose color scheme 3 (redish colors) \n
   * smooth1 : use ROOT smoothing algorithm once \n
   * smooth3 : use ROOT smoothing algorithm three times \n
   * smooth5 : use ROOT smoothing algorithm five times \n
   * smooth10 : use ROOT smoothing algorithm ten times \n
   * median : draw median and central interval \n
   * mode : draw global mode and standard deviation \n
   * mean : draw mean value and standard deviation \n
   * quartiles : indicate quartiles \n
   * deciles : indicate deciles \n
   * percentiles : indicate percentiles \n
   * L : add a legend \n
   * same: add histogram on top of another histogram
   * @param intervals the intervals
   */
  void Draw(std::string options="BTsiB3CS1D0Lmeanmode", std::vector<double> intervals=std::vector<double>(0));

  /**
   *Draw distribution into the active canvas.
   * @param options Drawing options, @see Print(const char * filename, std::string options, double interval, int ww=0, int wh=0)
   * @param interval an upper or lower limit
   */
  void Draw(std::string options, double interval);

  /**
   * Draw the 1D marginal for a parameter fixed by a delta prior.
   * @param value The fixed value of the parameter. */
  void DrawDelta(double value) const;

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
   * Return a vector containing information about the set of smallest
   * intervals. \n
   * 0 : x_min \n
   * 1 : x_max \n
   * 2 : relative height
   * 3 : local mode \n
   * 4 : relative area.
   * @param content The content of the smallest interval
   * @return the vector.
   */
  std::vector<double> GetSmallestIntervals(double content = 0.68);

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
   * @return new histogram which is nonzero only between min and max */
  TH1D* GetSubHisto(double min, double max, const char * name);

  /** @} */

 private:

  /**
   * The 1D histogram */
  TH1D* fHistogram;

  /**
   * Default confidence level limit */
  double fDefaultCLLimit;

  /**
   * Global mode */
  double fMode;

  /**
   * "Is there a global mode?" flag */
  int fModeFlag;

  /**
   * The colors of the color scheme. */
  std::vector<int> fColors;

  /**
   * Storage for plot objects. */
  mutable std::vector<TObject*> fROOTObjects;

  /**
   * Helper method to get an unique number to be used in histogram naming */
  static unsigned int getNextIndex()
    { return ++fHCounter; }

  /**
   * Helper variable to get an unique number to be used in histogram naming */
  static unsigned int fHCounter;
};

// ---------------------------------------------------------

#endif
