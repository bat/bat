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

/**
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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
      BCH1D();

      /**
       * The default constructor. */
      BCH1D(TH1D * hist)
         { fHistogram = hist; };

      /**
       * The default destructor. */
      ~BCH1D();

      /** @} */

      /** \name Member functions (get)  */
      /** @{ */

      /**
       * @return The one-dimensional histogram. */
      TH1D * GetHistogram()
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
       * Return the quantily of the distribution
       * @param probability The probability.
       * @return The quantile of the distribution for the probability.
       * @see GetQuantile(double probablity) */
      double GetLimit(double probability)
         { return this -> GetQuantile(probability); };

      /**
       * @return The RMS of the distribution. */
      double GetRMS()
         { return fHistogram -> GetRMS(); };

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
      int GetColor(int index) {
	return fColors.at(index); };

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
         { fMode=mode; fModeFlag=1; };

      /** @} */

      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Print distribution into a PostScript file.
       * @param filename Output filename
       * @param ww canvas size in pixels along X
       * @param ww canvas size in pixels along Y
       * If ww and wh are set to 0, default ROOT canvas size is used.
       * For explanation of parameters options and ovalue look at BCH1D::Draw()
       * method. */
      void Print(const char * filename, int options=0, double ovalue=0., int ww=0, int wh=0);

     /**
       * Print distribution into a PostScript file.
       * @param filename Output filename
       * @param option the draw options
       * @param ww canvas size in pixels along X
       * @param ww canvas size in pixels along Y
       * If ww and wh are set to 0, default ROOT canvas size is used.
       * For explanation of parameters options and ovalue look at BCH1D::Draw()
       * method. */
      void myPrint(const char * filename, std::string options="BT0B3CS1D0pdf0L", std::vector<double> intervals=std::vector<double>(0), int ww=0, int wh=0);
      void myPrint(const char * filename, std::string options, double interval, int ww=0, int wh=0);

      /**
       * Draw distribution into the active canvas.
       * @param options Drawing options: \n 0 = band mode [default], \n
       *                1 = draw vertical line, \n
       *                2 = band mode with minimal interval
       *                3 = ??
       *                4 = Delta prior, i.e. fixed parameter value
       * @param ovalue Option specific value. For option 0, if ovalue is nonzero
       *    a limit is to be drawn rather than central band with ovalue being the
       *    per cent value of the limit. If negative, limit is drawn from minimum,
       *    if positive limit is drawn from maximum. Allowed values are
       *    68 < |limit| < 100. If mode is outside the band, the limit is
       *    drawn automatically. The default limit can be changed by
       *    BCH1D::SetDefaultCLLimit(int limit). \n
       *    For option 1 the ovalue defines
       *    where the line is drawn. \n
       *    For option 2 the ovalue sets the content of
       *    the minimal interval in per cent. If omitted a 68% minimal interval
       *    will be drawn.
       *    For option 3 ???
       *    For option 4 draw one bin representing the delta prior around ovalue. */
      void Draw(int options=0, double ovalue=0.);

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
       * CS2 : choose color scheme 1 (green/yellow/red) [default] \n
       * CS3 : choose color scheme 2 (blueish colors) \n
       * CS4 : choose color scheme 3 (redish colors) \n
       * pdf0 : draw pdf [default] \n
       * pdf1 : draw cumulative pdf \n
			 * median : draw median and central interval \n
			 * mode : draw mode and standard deviation \n
			 * quartiles : indicate quartiles \n
			 * deciles : indicate deciles \n
			 * percentiles : indicate percentiles \n
			 * L : add legend \n
       * @param intervals:
       */
      void myDraw(std::string options="BTciB3CS1D0pdf0L", std::vector<double> intervals=std::vector<double>(0));
      void myDraw(std::string options, double interval);

     /**
       * Draw the 1D marginal for a parameter fixed by a delta prior.
       * @param value The fixed value of the parameter. */
      void DrawDelta(double value) const;

      /**
       * Draw distribution with band between min and max and with marker at the mode.
       * Write the location of the mode with uncertainties. If limit is specified,
       * draw CL limit. Allowed values are 68 < |limit| < 100. */
      void DrawShadedLimits(double mode, double min, double max, double limit=0);

      /**
       * Draw distribution with bands so that the total shaded area is the
       * smallest possible containing and integral of "prob". Draw the location
       * of the mean and median if requested (default). */
      void DrawSmallest(double mode, double prob, bool drawmean=true);

      /**
       * Include a legend for the symbols of mean, mode, median and
       * confidence band used in 1D marginalized posterior distributions.
       * @param text the text used to name the legend entry for the confidence band
       */
      void DrawLegend(const char* text);

      /**
       * Calculate the minimal interval of the distribution containing a given content.
       * @param min calculated minimum of the interval
       * @param max calculated maximum of the interval
       * @param content content of the interval [default is .68]
       * @return the content of the histogram between min and max */
      double GetSmallestInterval(double & min, double & max, double content=.68);

      TH1D * GetSmallestIntervalHistogram(double level);

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
      TH1D * GetSubHisto(double min, double max, const char * name);

      /** @} */

   private:

      /**
       * The 1D histogram */
      TH1D * fHistogram;

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
       * The colors of the color scheme. */
      std::vector<TObject*> fROOTObjects;
};

// ---------------------------------------------------------

#endif
