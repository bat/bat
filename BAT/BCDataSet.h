#ifndef __BCDATASET__H
#define __BCDATASET__H

/*!
 * \class BCDataSet
 * \brief A class representing a set of data points.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 08.2008
 * \detail This class represents a data set containing a set of data
 * points. The data points are organized in a vector. The class
 * provides functions to read in data from a file.
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

// BAT classes
class BCDataPoint;

// ---------------------------------------------------------

typedef std::vector<BCDataPoint*> BCDataVector;

// ---------------------------------------------------------

class BCDataSet
{
   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * Default constructor */
      BCDataSet();

      /**
       * The copy constructor */
      BCDataSet(const BCDataSet & bcdataset);

      /**
       * Default destructor */
      virtual ~BCDataSet();

      /** @} */
      /** \name Assignment operators */
      /** @{ */

      /**
       * Defaut assignment operator */
      BCDataSet & operator = (const BCDataSet & bcdataset);

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The number of data points. */
      unsigned int GetNDataPoints();

      /**
       * @return number of values per data point (dimension of data). */
      unsigned int GetNValuesPerPoint();

      /**
       * @param index The index of the data point to be returned.
       * @return The data point at the index. */
      BCDataPoint * GetDataPoint(unsigned int index);

      /**
       * Viewing the data set as a table with one row per point,
       * this method returns a specified column.
       * @param index The index of the component to be returned.
       * @return The (index)th component of all data points */
      std::vector<double> GetDataComponents(int index);

      /** @} */

      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Reads data from a file. For a description see the following
       * member functions. */
      int ReadDataFromFile(const char * filename, const char * treename, const char * branchnames)
         { return this ->  ReadDataFromFileTree(filename, treename, branchnames); };

      int ReadDataFromFile(const char * filename, int nvariables)
         { return this -> ReadDataFromFileTxt(filename, nvariables); };

      /**
       * Reads a TTree from a .root file.  Opens a .root file and
       * gets a TTree. It creates data points containing the values
       * read from the file.
       * @param filename The name of the .root file.
       * @param treename The name of the TTree.
       * @param branchnames A list of the names of the branches
       * separated by a comma
       * @return An error code.
       * @see ReadDataFromFileTxt(char* filename, int nbranches);
       * @see ReadDataFromFileUser(const char * filename, std::vector<int> options_int, std::vector<double> options_double, const char * options_char);
       */
      int ReadDataFromFileTree(const char * filename, const char * treename, const char * branchnames);

      /**
       * Reads data from a .txt file.  Opens a .txt file and creates
       * data objects containing the values read from the file.
       * @param filename The name of the .txt file.
       * @param nvariables The number of variables.
       * @see ReadDataFromFileTree(char* filename, char* treename, std::vector<char*> branchnames)
       * @see ReadDataFromFileUser(const char * filename, std::vector<int> options_int, std::vector<double> options_double, const char * options_char);
       */
      int ReadDataFromFileTxt(const char * filename, int nvariables);

      /**
       * Adds a data point to the data set.
       * @param datapoint The data point to be added */
      void AddDataPoint(BCDataPoint * datapoint);

      /**
       * Resets the content of the data set */
      void Reset();

      /**
       * Dump the data to the standard output */
      void Dump();

      /** @} */

   private:
      // todo why a pointer? creates a huge amount of unnecessary checks
      /**
       * A vector containing the data points */
      BCDataVector * fBCDataVector;

};

// ---------------------------------------------------------

#endif

