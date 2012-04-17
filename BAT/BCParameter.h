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
       * @param upperlimit The upper limit of the parameter values. */
      BCParameter(const char* name, double lowerlimit, double upperlimit);

      /**
       * The default copy constructor. */
      BCParameter(const BCParameter & parameter);

      /**
       * The default destructor. */
      ~BCParameter();

      /** @} */

      /** \name Assignment operators */
      /** @{ */

      /**
       * The defaut assignment operator. */
      BCParameter & operator = (const BCParameter & parameter);

      /** @} */

      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The name of the parameter. */
      std::string GetName()
         { return fName; };

      /**
       * Returns the index of the parameter within the parameter
       * container of a BCModel.
       * @return The index of the parameter in the model. */
      int GetIndex()
         { return fIndex; };

      /**
       * @return The lower limit of the parameter values. */
      double GetLowerLimit()
         { return fLowerLimit; };

      /**
       * @return The upper limit of the parameter values. */
      double GetUpperLimit()
         { return fUpperLimit; };

      /**
       * Returns the range width of the parameter values. It is
       * always a positive value.
       * @return The range width of the parameter values. */
      double GetRangeWidth()
         { return (fUpperLimit>fLowerLimit)?fUpperLimit-fLowerLimit:fLowerLimit-fUpperLimit; };

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

      /**
       * @param name The name of the parameter. */
      void SetName(const char * name)
         { fName = name; };

      /**
       * Set the index of the parameter within the parameter
       * container of a BCModel.
       * @param index The index of the parameter. */
      void SetIndex(int index)
         { fIndex = index; };

      /**
       * Set the lower limit of the parameter values.
       * @param limit The lower limit of the parameter values. */
      void SetLowerLimit(double limit = 0)
         { fLowerLimit = limit; };

      /**
       * Set the upper limit of the parameter values.
       * @param limit The upper limit of the parameter values. */
      void SetUpperLimit(double limit = 1)
         { fUpperLimit = limit; };

      /**
       * Set the limits of the parameter values.
       * @param lowerlimit The lower limit of the parameter values.
       * @param upperlimit The upper limit of the parameter values. */
      void SetLimits(double lowerlimit = 0, double upperlimit = 1)
         { fLowerLimit = lowerlimit; fUpperLimit = upperlimit; };

      /**
       * Set parameter to be nuisance.
       * @param nuisance 1 - nuisance, 0 - not nuisance */
      void SetNuisance(int nuisance = 1)
         { fNuisance = nuisance; };

      /** @} */

      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Returns 1 if parameter is a nuisance parameter or 0 if not.
       * @return 1 - is nuisance paramete, 0 - is not nuisance parameter */
      double IsNuisance()
         { return fNuisance; };

      /**
       * Returns true if the value is at a parameter limit.
       * @return flag States if value is at parameter limit. */
      bool IsAtLimit(double value);

      /**
       * Prints a parameter summary on the screen. */
      void PrintSummary();

      /** @} */

   private:

      /**
       * Copies this BCParameter into another one */
      void Copy(BCParameter & parameter) const;

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
       * Flag to specify whether to integrate over this parameter or not. */
      int fNuisance;

};

// ---------------------------------------------------------

/**
 * \typedef
 * \brief A vector of pointer to BCParameter.*/
typedef std::vector<BCParameter*> BCParameterSet;

// ---------------------------------------------------------

#endif
