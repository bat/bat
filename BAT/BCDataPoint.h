#ifndef __BCDATAPOINT__H
#define __BCDATAPOINT__H

/*!
 * \class BCDataPoint
 * \brief A class representing a data point.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a data point which is the basic unit
 * of information. A data point can be an event, a bin content, etc.
 * Each data point can store several variables of type  double.
 * The variables are organized in a vector.
 */

/**
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <vector>

// ---------------------------------------------------------

class BCDataPoint
{
   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * A constructor.
       * @param nvariables The number of variables stored in a data.
       * object */
      BCDataPoint(int nvariables);

      /**
       * A constructor.
       * @param x The vector containing the data. */
      BCDataPoint(std::vector<double> x);

      /**
       * The copy constructor. */
      BCDataPoint(const BCDataPoint & datapoint);

      /**
       * A destructor. */
      ~BCDataPoint();

      /** @} */
      /** \name Assignment operators */
      /** @{ */

      /**
       * Defaut assignment operator */
      BCDataPoint & operator = (const BCDataPoint & datapoint);

      /** @} */

      /** \name Member functions (get) */
      /** @{ */

      /**
       * @param index The index of the variable.
       * @return The value of the variable. */
      double GetValue(int index);

      /**
       * @return A vector of values. */
      std::vector<double> GetValues()
         { return fData; };

      /**
       * Returns the number of values. */
      unsigned int GetNValues()
         { return fData.size(); };

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

      /**
       * Set the value of a variable.
       * @param index The index of the variable
       * @param value The value of the variable */
      void SetValue(int index, double value);

      /**
       * Set the values of all variables.
       * @param values A vector of values */
      void SetValues(std::vector<double> values);

      /** @} */

   private:

      /**
       * The vector containing the values of the variables. */
      std::vector<double> fData;

};

// ---------------------------------------------------------

#endif

