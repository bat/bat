#ifndef __BCH2D__H
#define __BCH2D__H

/*!
 * \class BCH2D
 * \brief  A class for handling 2D distributions.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class contains a TH2D histogram and some additional
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
#include <string>

// ROOT classes

class TH1D;
class TH2D;
class TGraph;
class TObject;

// ---------------------------------------------------------

class BCH2D
{

public:

   /** \name Constructors and destructors */
   /** @{ */

   /**
    * The complete constructor. */
   BCH2D(TH2D* h = 0);

   /**
    * The default destructor. */
   ~BCH2D();

   /** @} */
   /** \name Member functions (get)  */
   /** @{ */

   /**
    * Return the TH2D histogram
    * @return The TH2D histogram. */
   TH2D* GetHistogram()
   { return fHistogram; };

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
    * Set the 2D histogram. */
   void SetHistogram(TH2D* hist);

   /**
    * Set global mode.
    * @param The global mode. */
   void SetGlobalMode(double mode[2])
   { fMode[0] = mode[0];
     fMode[1] = mode[1];
     fModeFlag =1; };

   /** @} */
   /** \name Member functions (miscellaneous methods) */
   /** @{ */

   /**
    * Print distribution into a PostScript file.
    * @param filename Output filename
    * @param option the draw options (see myDraw()) plus \n
    * @param intervals the intervals for the bands
    * logz : draw z-axis in log-scale \n
    * R : rescale canvas to have a squared histogram
    * @param ww canvas size in pixels along X
    * @param ww canvas size in pixels along Y
    * If ww and wh are set to 0, default ROOT canvas size is used. */
   void Print(const char* filename, std::string options="BTfB3CS1meangmode", std::vector<double> intervals=std::vector<double>(0), int ww=0, int wh=0);

   /**
    *Print distribution into a PostScript file.
    * @param filename Output filename
    * @param option the draw options, @see Print(const char* filename, std::string options="BTfB3CS1meangmode", std::vector<double> intervals=std::vector<double>(0), int ww=0, int wh=0)
    * @param interval an upper or lower limit
    * @param ww canvas size in pixels along X
    * @param ww canvas size in pixels along Y
    * @see Print(const char* filename, std::string options="BTfB3CS1meangmode", std::vector<double> intervals=std::vector<double>(0), int ww=0, int wh=0)
    */
   void Print(const char* filename, std::string options, double interval, int ww=0, int wh=0);


   /**
    * Draw distribution into the active canvas.
    * @param options Drawing options: \n
    * BTf : band type a filled area [default] \n
    * BTc : band type is a contour \n
    * B1 : draw one band corresponding to an integrated probability specified in intervals \n
    * B2 : draw two bands corresponding to an integrated probability  specified in intervals \n
    * B3 : draw three bands corresponding to an integrated probability  specified in intervals [default] \n
    * CS0 : choose color scheme 0 (B&W) \n
    * CS1 : choose color scheme 1 (green/yellow/red) [default] \n
    * CS2 : choose color scheme 2 (blueish colors) \n
    * CS3 : choose color scheme 3 (redish colors) \n
    * smooth1 : use ROOT smoothing algorithm once \n
    * smooth3 : use ROOT smoothing algorithm three times \n
    * smooth5 : use ROOT smoothing algorithm five times \n
    * smooth10 : use ROOT smoothing algorithm ten times \n
    * mean : draw mean value and standard deviation [default] \n
    * gmode : draw global mode [default] \n
    * lmode : draw global mode [default] \n
    * profilex : draw the profile line vs. x using the mode \n
    * profiley : draw the profile line vs. y using the mode \n
    * nL : remove legend \n
    * @param intervals the intervals for the bands
    */
   void Draw(std::string options="BTfB3CS1meangmodelmode", std::vector<double> intervals=std::vector<double>(0));

   /**
    * Draw distribution into the active canvas.
    * @param options Drawing options, @see Draw(std::string options="BTfB3CS1meangmodelmode", std::vector<double> intervals=std::vector<double>(0))
    * @param interval an upper or lower limit
    * @see Draw(std::string options="BTfB3CS1meangmodelmode", std::vector<double> intervals=std::vector<double>(0))
    */
   void Draw(std::string options, double interval);

   /**
    * Calculates the integral of the distribution as a function of the
    * height. */
   void CalculateIntegratedHistogram();

   /**
    * Print the integrated histogram.
    * @param filename the name of the file. */
   void PrintIntegratedHistogram(const char* filename);

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
   std::vector<int> GetNIntervalsY(TH2D* h, int &nfoundmax);

   /**
    * Return a graph of the profile along x or y. The profile is
    * calculated by scanning through the one axis, e.g., y, and
    * finding the mode, mean, median, etc. with respect to the
    * other axis, e.g. x.
    * @param axis x-axis (0) or y-axis (1)
    * @param options Options: \n
    * mode     : calculate the profile using the mode [default] \n
    * mean     : calculate the profile using the mean \n
    * median   : calculate the profile using the median \n
    */
   TGraph* CalculateProfileGraph(int axis, std::string options="mode");

   /**
    * Draw the profile along x or y
    * @param options The options of the calculation of the profile
    * @param drawoptions The drawing options: \n
    * black  : line color is black [default] \n
    * red    : line color is red \n
    * solid  : line is solid [default] \n
    * dashed : line is dashed \n
    * dotted : line is dotted \n
    * @return The profile graph
    */
   TGraph* DrawProfile(int axis, std::string options, std::string drawoptions="blacksolid");

   /** Draw the profile along x
    * @param options The options defined in DrawProfile
    * @param drawoptions The drawoptions defined DrawProfile
    * @return The profile graph
    */
   TGraph* DrawProfileX(std::string options, std::string drawoptions)
   { return DrawProfile(0, options, drawoptions); };

   /** Draw the profile along y
    * @param options The options defined in DrawProfile
    * @param drawoptions The drawoptions defined DrawProfile
    * @return The profile graph
    */
   TGraph* DrawProfileY(std::string options, std::string drawoptions)
   { return DrawProfile(1, options, drawoptions); };

   /** @} */

private:

   /**
    * The 2D histogram */
   TH2D* fHistogram;

   /**
    * The integrated 2D histogram */
   TH1D* fIntegratedHistogram;

   /**
    * Global mode */
   double fMode[2];

   /**
    * "Is there a global mode?" flag */
   int fModeFlag;

   /**
    * The colors of the color scheme. */
   std::vector<int> fColors;

   /**
    * The colors of the color scheme. */
   std::vector<TObject*> fROOTObjects;

   /** Helper method to get an unique number to be used in histogram naming */
   static unsigned int getNextIndex()
   { return ++fHCounter; }

   /** helper variable to get an unique number to be used in histogram naming */
   static unsigned int fHCounter;
};

// ---------------------------------------------------------

#endif
