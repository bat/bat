#ifndef __BCPARAMETER__H
#define __BCPARAMETER__H

/*!
 * \class BCParameter
 * \brief A class representing a parameter of a model.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a parameter of a model. It contains
 * information about the name and the range of the parameter.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include <string>

// ---------------------------------------------------------

class BCParameter
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      BCParameter();

      /**
       * A constructor.
       * @param name The name of the parameter.
       * @param lowerlimit The lower limit of the parameter values.
       * @param upperlimit The upper limit of the parameter values.
       * @param latexname The latex name of the parameter used in axis labeling.
        */
      BCParameter(const char* name, double lowerlimit, double upperlimit, const char* latexname = "");

      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The name of the parameter. */
      const std::string & GetName() const
         { return fName; }

      /**
       * Returns latex name if set, else identical to GetName().
       */
      const std::string & GetLatexName() const
         { return (fLatexName.empty()) ? fName : fLatexName; }

      /**
       * @return The lower limit of the parameter values. */
      double GetLowerLimit() const
         { return fLowerLimit; }

      /**
       * @return The upper limit of the parameter values. */
      double GetUpperLimit() const
         { return fUpperLimit; }

      /**
       * Returns the range width of the parameter values. It is
       * always a positive value.
       * @return The range width of the parameter values. */
      double GetRangeWidth() const
         { return (fUpperLimit > fLowerLimit) ? fUpperLimit - fLowerLimit : fLowerLimit - fUpperLimit; }

      bool FillHistograms() const
         { return fFillHistograms; }

      bool Fixed() const
         { return fFixed; }

      double GetFixedValue() const
         { return fFixedValue; }

      unsigned GetNbins() const
         { return fNbins; }

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

      /**
       * @param name The name of the parameter. */
      void SetName(const char * name)
         { fName = name; }

      void SetLatexName(const char * latex_name)
         { fLatexName = latex_name; }

      /**
       * Set the lower limit of the parameter values.
       * @param limit The lower limit of the parameter values. */
      void SetLowerLimit(double limit = 0)
         { fLowerLimit = limit; }

      /**
       * Set the upper limit of the parameter values.
       * @param limit The upper limit of the parameter values. */
      void SetUpperLimit(double limit = 1)
         { fUpperLimit = limit; }

      /**
       * Set the limits of the parameter values.
       * @param lowerlimit The lower limit of the parameter values.
       * @param upperlimit The upper limit of the parameter values. */
      void SetLimits(double lowerlimit = 0, double upperlimit = 1)
         { fLowerLimit = lowerlimit; fUpperLimit = upperlimit; }

      void Fix(double value)
      {
         fFixed = true;
         fFixedValue = value;
      }

      void Unfix()
         { fFixed = false; }

      void FillHistograms(bool flag)
         { fFillHistograms = flag; }

      void SetNbins(unsigned nbins)
         { fNbins = nbins; }
      /** @} */

      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Returns true if the value is at a parameter limit.
       * @return flag States if value is at parameter limit. */
      bool IsAtLimit(double value) const;

      bool IsValid(double value) const
      { return (fLowerLimit <= value) && (value <= fUpperLimit) ? true : false; }

      /**
       * Prints a parameter summary on the screen. */
      void PrintSummary() const;

      /** @} */

   private:
      /// The name of the parameter.
      std::string fName;

      ///  The lower limit of the parameter value.
      double fLowerLimit;

      /// The upper limit of the parameter value.
      double fUpperLimit;

      /// The latex name of the parameter.
      std::string fLatexName;

      /// Flag to fix parameter; useful for example, for integration.
      bool fFixed;

      /// The fixed value of the parameter.
      double fFixedValue;

      /// Flag to store MCMC samples in histograms.
      bool fFillHistograms;

      /// The number of equal-size bins used in histograms involving this parameter.
      unsigned fNbins;
};
#endif
