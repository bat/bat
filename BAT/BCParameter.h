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

/**
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <string>
#include <vector>

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
       * Returns the index of the parameter within the parameter
       * container of a BCModel.
       * @return The index of the parameter in the model. */
      int GetIndex() const
         { return fIndex; }

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
       * Set the index of the parameter within the parameter
       * container of a BCModel.
       * @param index The index of the parameter. */
      void SetIndex(int index)
         { fIndex = index; }

      /**
       * Set parameter to be nuisance.
       * @param nuisance false - not nuisance */
      void SetNuisance(bool nuisance = true)
         { fNuisance = nuisance; }

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

      /** @} */

      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Returns 1 if parameter is a nuisance parameter or 0 if not.
       * @return 1 - is nuisance paramete, 0 - is not nuisance parameter */
      bool IsNuisance() const
         { return fNuisance; }

      /**
       * Returns true if the value is at a parameter limit.
       * @return flag States if value is at parameter limit. */
      bool IsAtLimit(double value) const;

      /**
       * Prints a parameter summary on the screen. */
      void PrintSummary() const;

      /** @} */

   private:

      /**
       * The name of the parameter. */
      std::string fName;

      /**
       * The index of the parameter within the BCParameterSet of a BCModel. */
      int fIndex;

      /**
       * The lower limit of the parameter value. */
      double fLowerLimit;

      /**
       * The upper limit of the parameter value. */
      double fUpperLimit;

      /**
       * The latex name of the parameter. */
      std::string fLatexName;

      /**
       * Flag to specify whether to integrate over this parameter or not. */
      bool fNuisance;

};

// ---------------------------------------------------------

/**
 * \typedef
 * \brief A vector of pointer to BCParameter.*/
typedef std::vector<BCParameter*> BCParameterSet;

// ---------------------------------------------------------

#endif
