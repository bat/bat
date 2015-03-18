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
#include "BCLog.h"
#include "BCDataPoint.h"

// ---------------------------------------------------------

class BCDataSet
{
public:

    /** \name Constructors and destructor */
    /** @{ */

    /**
     * Default constructor
     * @param n Dimensionality (Number of values inside) of data points. */
    BCDataSet(unsigned n = 0);

    /**
     * The copy constructor */
    BCDataSet(const BCDataSet& bcdataset);

    /**
     * Destructor */
    virtual ~BCDataSet();

    /** @} */
    /** \name Assignment operators */
    /** @{ */

    /**
     * Assignment operator */
    BCDataSet& operator = (const BCDataSet& bcdataset)
    { Copy(bcdataset); return *this; }

    /** @} */
    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return The number of data points. */
    unsigned GetNDataPoints() const
    { return fDataVector.size(); }

    /**
     * @return number of values per data point (dimension of data). */
    unsigned GetNValuesPerPoint() const
    { return fNValuesPerPoint; }

    /**
     * @param index The index of the data point to be returned.
     * @return The data point at the index. */
    BCDataPoint* GetDataPoint(unsigned index) const
    { return (index < fDataVector.size()) ? fDataVector[index] : 0; }

    /**
     * Viewing the data set as a table with one row per point,
     * this method returns a specified column.
     * @param index The index of the component to be returned.
     * @return The (index)th component of all data points */
    std::vector<double> GetDataComponents(unsigned index) const;

    /**
     * @return Whether bounds for the data set exist.*/
    bool BoundsExist() const;

    /**
     * @return BCDataPoint with values set to actual lower bounds of data. */
    const BCDataPoint& GetLowerBounds() const
    { return fLowerBounds; }

    /**
     * @return BCDataPoint with values set to actual upper bounds of data. */
    const BCDataPoint& GetUpperBounds() const
    { return fUpperBounds; }

    /**
     * @return BCDataPoint with values set to user-set lower bounds of data. */
    BCDataPoint& GetUserLowerBounds()
    { return fUserLowerBounds; }

    /**
     * @param check_size Flag for checking size of user-set bound object against data point size.
     * @return BCDataPoint with values set to user-set upper bounds of data. */
    BCDataPoint& GetUserUpperBounds()
    { return fUserUpperBounds; }

    /**
     * Return user-set lower bound on data, if set, otherwise actual lower bound.
     * @param index Index of data value to return lower bound of.
     * @return Lower bound on data values. */
    double GetLowerBound(unsigned index) const;

    /**
     * Return user-set upper bound on data, if set, otherwise actual upper bound.
     * @param index Index of data value to return upper bound of.
     * @return Upper bound on data values. */
    double GetUpperBound(unsigned index) const;

    /**
     * Return upper-bound minus lower-bound for data axis,
     * using user-set bounds, if provided, other actual bounds.
     * @return range width.*/
    double GetRangeWidth(unsigned index) const
    { return GetUpperBound(index) - GetLowerBound(index); }

    /**
     * Return wether data axis is fixed
     * @param index Index of axis to query
     * @return Whether data axis is fixed. */
    bool IsFixed(unsigned index) const
    { return (index < fFixed.size() and fFixed[index]); }

    /** @} */

    /** \name Setters */
    /** @{ */

    /**
     * Set number of values inside each data point.
     * If set to zero, then this will be set by first added data point.
     * @param n Size of data point. */
    void SetNValuesPerPoint(unsigned n);

    /**
     * Set bounds for data values
     * @param index Index of data axis to provide bounds for. */
    void SetBounds(unsigned index, double lower_bound, double upper_bound, bool fixed = false);

    /**
     * Set fixed flag of a data axis
     * @param i index of axis to fix
     * @param b whether to fix (true) or unfix (false) */
    void Fix(unsigned i, bool b = true)
    { if (i < fFixed.size()) fFixed[i] = b; }

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Copy a data set. */
    void Copy(const BCDataSet& other);

    /**
     * Reads data from a TTree in file.
     * @param filename Path to file containing TTree object.
     * @param treename Name of TTree inside file.
     * @param branchnames List of names of branches to be read from TTree
     * @param delim Character deliminating branch names (default: comma).
     * @return Success of action. */
    bool ReadDataFromFile(const char* filename, const char* treename, const char* branchnames, char delim = ',')
    { return ReadDataFromFileTree(filename, treename, branchnames, delim); };

    /**
     * Reads data from a file
     * @param filename Path to file containing data.
     * @param nvariables Number of variables (columns) in data file.
     * @return Success of action. */
    bool ReadDataFromFile(const char* filename, int nvariables)
    { return ReadDataFromFileTxt(filename, nvariables); };

    /**
     * Reads a TTree from a .root file.  Opens a .root file and
     * gets a TTree. It creates data points containing the values
     * read from the file.
     * @param filename The name of the .root file.
     * @param treename The name of the TTree.
     * @param branchnames A list of the names of the branches
     * @param delim Character deliminating between branch names (default: comma)
     * @return Success of action. */
    bool ReadDataFromFileTree(const char* filename, const char* treename, std::string branchnames, char delim = ',');

    /**
     * Reads data from a .txt file.  Opens a .txt file and creates
     * data objects containing the values read from the file.
     * @param filename The name of the .txt file.
     * @param nvariables The number of variables (columns) in file.
     * @return Success of action. */
    bool ReadDataFromFileTxt(const char* filename, int nvariables);

    /**
     * Adds a data point to the data set.
     * @param datapoint The data point to be added */
    bool AddDataPoint(BCDataPoint* datapoint);

    /**
     * Resets the content of the data set */
    void Reset();

    /**
     * Dump the data to the standard output
     * @param output Function of (const char *) to handle the output (default = BCLog::OutSummary). */
    void Dump(void (*output)(const char*) = BCLog::OutSummary) const;

    /** @} */

private:
    /**
     * A vector containing the data points */
    std::vector<BCDataPoint*> fDataVector;

    /** number of values per data point. */
    unsigned fNValuesPerPoint;

    /** pointer to BCDataPoint storing actual lower bounds of data set. */
    BCDataPoint fLowerBounds;

    /** pointer to BCDataPoint storing actual upper bounds of data set. */
    BCDataPoint fUpperBounds;

    /** pointer to BCDataPoint storing user-set lower bounds of data set. */
    BCDataPoint fUserLowerBounds;

    /** pointer to BCDataPoint storing user-set upper bounds of data set. */
    BCDataPoint fUserUpperBounds;

    /** vector of flags for whether data values are fixed. */
    std::vector<bool> fFixed;

};

// ---------------------------------------------------------

#endif

