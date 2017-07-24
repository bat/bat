#ifndef __BCDATASET__H
#define __BCDATASET__H

/**
 * @class BCDataSet
 * @brief A class representing a set of data points.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 08.2008
 * @details This class represents a data set containing a set of data
 * points. The data points are organized in a vector. The class
 * provides functions to read in data from a file.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <string>
#include <vector>

#include "BCDataPoint.h"
#include "BCLog.h"

class TGraph;
class TGraphErrors;
class TGraphAsymmErrors;
class TH2;

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
     * Destructor */
    virtual ~BCDataSet() {};

    /** @} */
    /** \name operators*/
    /** @{ */

    /**
     * Raw and fast access. */
    BCDataPoint& operator[](unsigned index)
    {	return fDataVector[index]; }

    /**
     * Raw and fast access. */
    const BCDataPoint& operator[](unsigned index) const
    {	return fDataVector[index]; }

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
     * Safer, but slower, access to data points
     * @param index The index of the data point to be returned.
     * @return The data point at the index. */
    BCDataPoint& GetDataPoint(unsigned index)
    { return fDataVector.at(index); }

    /**
     * Safer, but slower, access to data points
     * @param index The index of the data point to be returned.
     * @return The data point at the index. */
    const BCDataPoint& GetDataPoint(unsigned index) const
    { return fDataVector.at(index); }

    /**
     * Access to last added data point. */
    BCDataPoint& Back()
    { return fDataVector.back(); }

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
     * Reads data from a TTree in file.
     * @param filename Path to file containing TTree object.
     * @param treename Name of TTree inside file.
     * @param branchnames List of names of branches to be read from TTree
     * @param delim Character deliminating branch names (default: comma).
     * @return Success of action. */
    bool ReadDataFromFile(const std::string& filename, const std::string& treename, const std::string& branchnames, char delim = ',')
    { return ReadDataFromFileTree(filename, treename, branchnames, delim); };

    /**
     * Reads data from a file
     * @param filename Path to file containing data.
     * @param nvariables Number of variables (columns) in data file.
     * @return Success of action. */
    bool ReadDataFromFile(const std::string& filename, int nvariables)
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
    bool ReadDataFromFileTree(const std::string& filename, const std::string& treename, const std::string& branchnames, char delim = ',');

    /**
     * Reads data from a .txt file.  Opens a .txt file and creates
     * data objects containing the values read from the file.
     * @param filename The name of the .txt file.
     * @param nvariables The number of variables (columns) in file.
     * @return Success of action. */
    bool ReadDataFromFileTxt(const std::string& filename, int nvariables);

    /**
     * Adds a data point to the data set.
     * @param datapoint The data point to be added */
    bool AddDataPoint(const BCDataPoint& datapoint);

    /**
     * Recalculate a data axis bound accounting for uncertainties
     * specified by other data axes. If a second error index is
     * provided, the first is taken as error below the value, and the
     * second as error above the value.
     * @param i Index of the data axis to be recalculated
     * @param nSigma Multiples of the stored uncertainty to account for
     * @param i_err1 Index of the data axis containing uncertainty (below point, if next i_err2 also specified)
     * @param i_err2 Index of the data axis containing uncertainty above point*/
    void AdjustBoundForUncertainties(unsigned i, double nSigma, unsigned i_err1, int i_err2 = -1);

    /**
     * Resets the content of the data set */
    void Reset()
    { fDataVector.clear(); SetNValuesPerPoint(0); }

    /**
     * Print summary to string handler
     * @param output String handler (default = BCLog::OutSummary). */
    void PrintSummary(void (*output)(const std::string&) = BCLog::OutSummary) const;

    /**
     * Get data set as ROOT TGraph object,
     * @param x Index of data axis plotted as abscissa
     * @param y Index of data axis plotted as ordinate
     * @return pointer to filled ROOT TGraph */
    TGraph* GetGraph(unsigned x, unsigned y) const;

    /**
     * Get data set as ROOT TGraphErrors object.
     * Set error indices negative to leave errors unset.
     * @param x Index of data axis plotted as abscissa
     * @param y Index of data axis plotted as ordinate
     * @param ex Index of data axis for error on abscissa
     * @param ey Index of data axis for error on ordinate
     * @return pointer to filled ROOT TGraphErrors */
    TGraphErrors* GetGraph(unsigned x, unsigned y, int ex, int ey) const;

    /**
     * Get data set as ROOT TGraphAsymmErrors object.
     * Set error indices negative to leave errors unset.
     * @param x Index of data axis plotted as abscissa
     * @param y Index of data axis plotted as ordinate
     * @param ex_below Index of data axis for error on abscissa below data points
     * @param ex_above Index of data axis for error on abscissa below data points
     * @param ey_below Index of data axis for error on ordinate below data points
     * @param ey_above Index of data axis for error on ordinate below data points
     * @return pointer to filled ROOT TGraphAsymmErrors */
    TGraphAsymmErrors* GetGraph(unsigned x, unsigned y, int ex_below, int ex_above, int ey_below, int ey_above) const;

    /**
     * Get ROOT TH2 with ranges set to data bounds.
     * Padding is specified as fraction of boundary range.
     * @param x Index of data axis for abscissa
     * @param y Index of data axis for ordinate
     * @param nbins_x number of bins on abscissa (default 100)
     * @param nbins_y number of bins on ordinate (default 100)
     * @param x_padding Amount to pad on either side of abscissa boundaries (default = 0.10)
     * @param y_padding Amount to pad on either side of ordinate boundaries (default = 0.10)
     * @return pointer to empty ROOT TH2 */
    TH2* CreateH2(const char* name, const char* title, unsigned x, unsigned y, unsigned nbins_x = 100, unsigned nbins_y = 100, double x_padding = 0.10, double y_padding = 0.10) const;

    /** @} */

private:
    /**
     * A vector containing the data points */
    std::vector<BCDataPoint> fDataVector;

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
