/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCDataSet.h"

#include "BCDataPoint.h"
#include "BCLog.h"

#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

// ---------------------------------------------------------
BCDataSet::BCDataSet(unsigned n)
{
    SetNValuesPerPoint(n);
}

// ---------------------------------------------------------
std::vector<double> BCDataSet::GetDataComponents(unsigned index) const
{
    std::vector<double> components;

    if (index >= fNValuesPerPoint)
        return components;

    // reserve space
    components.reserve(fDataVector.size());

    // loop over data points
    for (unsigned i = 0; i < fDataVector.size(); ++i)
        components.push_back(fDataVector[i][index]);

    return components;
}

// ---------------------------------------------------------
bool BCDataSet::BoundsExist() const
{
    for (unsigned i = 0; i < GetNValuesPerPoint(); ++i)
        if (!std::isfinite(GetLowerBound(i)) or !std::isfinite(GetUpperBound(i)))
            return false;
    return true;
}

// ---------------------------------------------------------
double BCDataSet::GetLowerBound(unsigned index) const
{
    if (index > GetNValuesPerPoint()) {
        BCLog::OutError("BCDataSet::GetLowerBound : index out of range.");
        return std::numeric_limits<double>::infinity();
    }
    if (std::isfinite(fUserLowerBounds[index]))
        return fUserLowerBounds[index];
    return fLowerBounds[index];
}

// ---------------------------------------------------------
double BCDataSet::GetUpperBound(unsigned index) const
{
    if (index > GetNValuesPerPoint()) {
        BCLog::OutError("BCDataSet::GetUpperBound : index out of range.");
        return -std::numeric_limits<double>::infinity();
    }
    if (std::isfinite(fUserUpperBounds[index]))
        return fUserUpperBounds[index];
    return fUpperBounds[index];
}

// ---------------------------------------------------------
bool BCDataSet::ReadDataFromFileTree(const std::string& filename, const std::string& treename, const std::string& branchnames, char delim)
{
    // open root file
    TFile* file = TFile::Open(filename.data(), "READ");

    // check if file is open and warn if not.
    if (!file->IsOpen()) {
        BCLog::OutError("BCDataSet::ReadDataFromFileTree : Could not open file " + filename + ".");
        return false;
    }

    // get tree
    TTree* tree = (TTree*) file->Get(treename.data());

    // check if tree is there and warn if not.
    if (!tree) {
        BCLog::OutError("BCDataSet::ReadDataFromFileTree : Could not find TTree " + treename + ".");
        file->Close();
        return false;
    }

    // calculate maximum number of entries
    long nentries = tree->GetEntries();

    // check if there are any events in the tree and close file if not.
    if (nentries <= 0) {
        BCLog::OutError("BCDataSet::ReadDataFromFileTree : No events in TTree " + treename + ".");
        file->Close();
        return false;
    }

    // if data set contains data, clear data object container ...
    if (!fDataVector.empty()) {
        Reset();
        BCLog::OutDetail("BCDataSet::ReadDataFromFileTree : Overwrite existing data.");
    }

    // define a vector of std::strings which contain the tree names.
    std::vector<std::string> branches;
    // split branchnames string up into above vector
    std::stringstream branch_ss(branchnames);
    std::string branchname;
    while (std::getline(branch_ss, branchname, delim))
        if (!branchname.empty())
            branches.push_back(branchname);

    // create temporary vector with data and assign some zeros.
    std::vector<double> data(branches.size(), 0);

    // set the branch address.
    for (unsigned i = 0; i < branches.size(); ++i)
        tree->SetBranchAddress(branches[i].data(), &data[i]);

    // loop over entries
    for (long ientry = 0; ientry < nentries; ++ientry) {
        tree->GetEntry(ientry);
        AddDataPoint(BCDataPoint(data));
    }

    file->Close();

    // remove file pointer.
    if (file)
        delete file;

    return true;
}

// ---------------------------------------------------------
bool BCDataSet::ReadDataFromFileTxt(const std::string& filename, int nbranches)
{
    // open text file.
    std::fstream file;
    file.open(filename.data(), std::fstream::in);

    // check if file is open and warn if not.
    if (!file.is_open()) {
        BCLog::OutError("BCDataSet::ReadDataFromFileText : Could not open file " + filename + ".");
        return false;
    }

    // if data set contains data, clear data object container ...
    if (!fDataVector.empty()) {
        Reset();
        BCLog::OutDetail("BCDataSet::ReadDataFromFileTxt : Overwrite existing data.");
    }

    // create temporary vector with data and assign some zeros.
    std::vector<double> data(nbranches, 0);

    // reset counter
    int nentries = 0;

    // read data and create data points.
    while (!file.eof()) {

        // read data from file
        int i = 0;
        while (file >> data[i]) {
            if (i == nbranches - 1)
                break;
            i++;
        }

        // create data point.
        if (i == nbranches - 1) {
            AddDataPoint(BCDataPoint(data));
            ++nentries;
        }
    }

    // issue error if no entries were loaded
    if (nentries <= 0)
        BCLog::OutError("BCDataSet::ReadDataFromFileText : No events in the file " + filename + ".");

    file.close();

    return (nentries > 0);
}

// ---------------------------------------------------------
bool BCDataSet::AddDataPoint(const BCDataPoint& datapoint)
{
    if (fNValuesPerPoint == 0 and fDataVector.empty())
        SetNValuesPerPoint(datapoint.GetNValues());

    if (datapoint.GetNValues() != GetNValuesPerPoint())
        return false;

    fDataVector.push_back(datapoint);

    for (unsigned i = 0; i < GetNValuesPerPoint(); ++i) {
        // check lower bound
        if (fDataVector.back()[i] < fLowerBounds[i])
            fLowerBounds[i] = fDataVector.back()[i];
        // check upper bound
        if (fDataVector.back()[i] > fUpperBounds[i])
            fUpperBounds[i] = fDataVector.back()[i];
    }

    return true;
}

// ---------------------------------------------------------
void BCDataSet::AdjustBoundForUncertainties(unsigned i, double nSigma, unsigned i_err1, int i_err2)
{
    // check indices
    if (i >= GetNValuesPerPoint() or i_err1 >= GetNValuesPerPoint() or i_err2 >= (int)GetNValuesPerPoint())
        return;

    // if uncertainty above value is unassigned, use same data axis as for below.
    if (i_err2 < 0)
        i_err2 = i_err1;

    // recalculate bounds accounting for uncertainty
    for (unsigned j = 0; j < fDataVector.size(); ++j) {

        // check lower bound
        if (fDataVector[j][i] - nSigma * fDataVector[j][i_err1] < fLowerBounds[i])
            fLowerBounds[i] = fDataVector[j][i] - nSigma * fDataVector[j][i_err1];

        // check upper bound
        if (fDataVector[j][i] + nSigma * fDataVector[j][i_err2] > fUpperBounds[i])
            fUpperBounds[i] = fDataVector[j][i] + nSigma * fDataVector[j][i_err2];
    }
}

// ---------------------------------------------------------
void BCDataSet::SetNValuesPerPoint(unsigned n)
{
    fNValuesPerPoint = n;
    fLowerBounds.SetNValues(n, std::numeric_limits<double>::infinity());
    fUpperBounds.SetNValues(n, -std::numeric_limits<double>::infinity());
    fUserLowerBounds.SetNValues(n, std::numeric_limits<double>::infinity());
    fUserUpperBounds.SetNValues(n, -std::numeric_limits<double>::infinity());
    fFixed.assign(n, false);
}

// ---------------------------------------------------------
void BCDataSet::SetBounds(unsigned index, double lower_bound, double upper_bound, bool fixed)
{
    if (index >= GetNValuesPerPoint()) {
        BCLog::OutError("BCDataSet::SetBounds : index out of range.");
        return;
    }
    if (lower_bound >= upper_bound) {
        BCLog::OutWarning("BCDataSet::SetBounds : lower bound is greater than or equal to upper_bound.");
        return;
    }
    fUserLowerBounds[index] = lower_bound;
    fUserUpperBounds[index] = upper_bound;
    fFixed[index] = fixed;
}

// ---------------------------------------------------------
void BCDataSet::PrintSummary(void (*output)(const std::string&)) const
{
    output("Data set summary:");
    output(Form("Number of points           : %u", GetNDataPoints()));
    output(Form("Number of values per point : %u", GetNValuesPerPoint()));
    for (unsigned i = 0; i < fDataVector.size(); ++i) {
        output(Form("Data point %5u", i));
        fDataVector[i].PrintSummary(output);
    }
}

// ---------------------------------------------------------
TGraph* BCDataSet::GetGraph(unsigned x, unsigned y) const
{
    if (x >= GetNValuesPerPoint() or y >= GetNValuesPerPoint())
        return NULL;

    TGraph* G = new TGraph();

    // fill graph
    for (unsigned i = 0; i < fDataVector.size(); ++i)
        G->SetPoint(i, fDataVector[i][x], fDataVector[i][y]);

    return G;
}

// ---------------------------------------------------------
TGraphErrors* BCDataSet::GetGraph(unsigned x, unsigned y, int ex, int ey) const
{
    if (x >= GetNValuesPerPoint() or y >= GetNValuesPerPoint()
            or ex >= (int)GetNValuesPerPoint() or ey >= (int)GetNValuesPerPoint())
        return NULL;

    TGraphErrors* G = new TGraphErrors();

    // fill graph
    for (unsigned i = 0; i < fDataVector.size(); ++i) {
        G->SetPoint(i, fDataVector[i][x], fDataVector[i][y]);
        double EX = (ex >= 0) ? fDataVector[i][ex] : 0;
        double EY = (ey >= 0) ? fDataVector[i][ey] : 0;
        G->SetPointError(i, EX, EY);
    }

    return G;
}

// ---------------------------------------------------------
TGraphAsymmErrors* BCDataSet::GetGraph(unsigned x, unsigned y, int ex_below, int ex_above, int ey_below, int ey_above) const
{
    if (x >= GetNValuesPerPoint() or y >= GetNValuesPerPoint()
            or ex_below >= (int)GetNValuesPerPoint() or ex_above >= (int)GetNValuesPerPoint()
            or ey_below >= (int)GetNValuesPerPoint() or ey_above >= (int)GetNValuesPerPoint())
        return NULL;

    TGraphAsymmErrors* G = new TGraphAsymmErrors();

    // fill graph
    for (unsigned i = 0; i < fDataVector.size(); ++i) {
        G->SetPoint(i, fDataVector[i][x], fDataVector[i][y]);
        double EXb = (ex_below >= 0) ? fDataVector[i][ex_below] : 0;
        double EXa = (ex_above >= 0) ? fDataVector[i][ex_above] : 0;
        double EYb = (ey_below >= 0) ? fDataVector[i][ey_below] : 0;
        double EYa = (ey_above >= 0) ? fDataVector[i][ey_above] : 0;
        G->SetPointError(i, EXb, EXa, EYb, EYa);
    }

    return G;
}

// ---------------------------------------------------------
TH2* BCDataSet::CreateH2(const char* name, const char* title, unsigned x, unsigned y, unsigned nbins_x, unsigned nbins_y, double x_padding, double y_padding) const
{
    if (x >= GetNValuesPerPoint() or y >= GetNValuesPerPoint())
        return NULL;

    if (!BoundsExist())
        return NULL;

    double x_low = GetLowerBound(x);
    double x_high = GetUpperBound(x);
    double y_low = GetLowerBound(y);
    double y_high = GetUpperBound(y);

    if (x_padding > 0) {
        double dX = x_padding * (x_high - x_low);
        x_low  -= dX;
        x_high += dX;
    }
    if (y_padding > 0) {
        double dY = y_padding * (y_high - y_low);
        y_low  -= dY;
        y_high += dY;
    }

    return new TH2D(name, title, nbins_x, x_low, x_high, nbins_y, y_low, y_high);
}
