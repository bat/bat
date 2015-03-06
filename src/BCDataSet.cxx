/*
 * Copyright (C) 2007-2013, the BAT core developer team
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
#include <TTree.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <math.h>

// ---------------------------------------------------------
BCDataSet::BCDataSet(unsigned n) {
	SetNValuesPerPoint(n);
}

// ---------------------------------------------------------
BCDataSet::~BCDataSet() {
	for (int i=fDataVector.size()-1; i>=0; --i)
		if (fDataVector[i])
			delete fDataVector[i];
}

// ---------------------------------------------------------
BCDataSet::BCDataSet(const BCDataSet & bcdataset) {
	Copy(bcdataset);
}

// ---------------------------------------------------------
void BCDataSet::Copy(const BCDataSet & other) {
	Reset();
	fNValuesPerPoint = other.fNValuesPerPoint;
	for (unsigned i=0; i<other.fDataVector.size(); ++i)
		if (other.fDataVector[i]) // only use existing data points
			fDataVector.push_back(new BCDataPoint(*(other.fDataVector[i])));
	fLowerBounds = other.fLowerBounds;
	fUpperBounds = other.fUpperBounds;
	fUserLowerBounds = other.fUserLowerBounds;
	fUserUpperBounds = other.fUserUpperBounds;
	fFixed = other.fFixed;
}

// ---------------------------------------------------------
std::vector<double> BCDataSet::GetDataComponents(unsigned index) const {
   std::vector<double> components;

	 if (index>=fNValuesPerPoint)
		 return components;

	 // reserve space
	 components.reserve(fDataVector.size());

	 // loop over data points
	 for (unsigned i=0; i<fDataVector.size(); ++i)
		 if (fDataVector[i])				// only use existing data points
			 components.push_back(fDataVector[i]->GetValue(index));
	 return components;
}

// ---------------------------------------------------------
bool BCDataSet::BoundsExist() const {
	for (unsigned i=0; i<GetNValuesPerPoint(); ++i)
		if (!std::isfinite(GetLowerBound(i)) or !std::isfinite(GetUpperBound(i)))
			return false;
	return true;
}

// ---------------------------------------------------------
double BCDataSet::GetLowerBound(unsigned index) const {
	if (index>GetNValuesPerPoint()) {
		BCLog::OutError("BCDataSet::GetLowerBound : index out of range.");
		return std::numeric_limits<double>::infinity();
	}
	if (std::isfinite(fUserLowerBounds.GetValue(index)))
		return fUserLowerBounds.GetValue(index);
	return fLowerBounds.GetValue(index);
}

// ---------------------------------------------------------
double BCDataSet::GetUpperBound(unsigned index) const {
	if (index>GetNValuesPerPoint()) {
		BCLog::OutError("BCDataSet::GetUpperBound : index out of range.");
		return -std::numeric_limits<double>::infinity();
	}
	if (std::isfinite(fUserUpperBounds.GetValue(index)))
		return fUserUpperBounds.GetValue(index);
	return fUpperBounds.GetValue(index);
}
			
// ---------------------------------------------------------
bool BCDataSet::ReadDataFromFileTree(const char * filename, const char * treename, std::string branchnames, char delim) {
	// open root file
	TFile * file = TFile::Open(filename, "READ");
	
	// check if file is open and warn if not.
	if (!file->IsOpen()) {
		BCLog::OutError(Form("BCDataSet::ReadDataFromFileTree : Could not open file %s.", filename));
		return false;
	}
	
	// get tree
	TTree * tree = (TTree*) file->Get(treename);

	// check if tree is there and warn if not.
	if (!tree) {
		BCLog::OutError(Form("BCDataSet::ReadDataFromFileTree : Could not find TTree %s.", treename));
		file -> Close();
		return false;
	}

	// calculate maximum number of entries
	long nentries = tree -> GetEntries();

   // check if there are any events in the tree and close file if not.
	if (nentries <= 0) {
		BCLog::OutError(Form("BCDataSet::ReadDataFromFileTree : No events in TTree %s.", treename));
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
	while (std::getline(branch_ss,branchname,delim))
		if (!branchname.empty())
			branches.push_back(branchname);
	
	// create temporary vector with data and assign some zeros.
	std::vector<double> data(branches.size(),0);

	// set the branch address.
	for (unsigned i=0; i < branches.size(); ++i)
		tree -> SetBranchAddress(branches[i].data(),&data[i]);

	// loop over entries
	for (long ientry=0; ientry<nentries; ++ientry) {
		tree -> GetEntry(ientry);
		AddDataPoint(new BCDataPoint(data));
	}

	file -> Close();

	// remove file pointer.
	if (file)
		delete file;

	return true;
}

// ---------------------------------------------------------
bool BCDataSet::ReadDataFromFileTxt(const char * filename, int nbranches) {
	// open text file.
	std::fstream file;
	file.open(filename, std::fstream::in);
	
	// check if file is open and warn if not.
	if (!file.is_open()) {
		BCLog::OutError(Form("BCDataSet::ReadDataFromFileText : Could not open file %s.", filename));
		return false;
	}

   // if data set contains data, clear data object container ...
	if (!fDataVector.empty()) {
		Reset();
		BCLog::OutDetail("BCDataSet::ReadDataFromFileTxt : Overwrite existing data.");
	}

	// create temporary vector with data and assign some zeros.
	std::vector<double> data(nbranches,0);

	// reset counter
	int nentries = 0;

	// read data and create data points.
	while (!file.eof()) {

		// read data from file
		int i=0;
		while(file >> data[i]) {
			if (i==nbranches-1)
				break;
			i++;
		}

		// create data point.
		if(i == nbranches-1) {
			AddDataPoint(new BCDataPoint(data));
			++nentries;
		}
	}

	// issue error if no entries were loaded
	if (nentries <= 0)
		BCLog::OutError(Form("BCDataSet::ReadDataFromFileText : No events in the file %s.", filename));

	file.close();

	return (nentries>0);
}

// ---------------------------------------------------------
bool BCDataSet::AddDataPoint(BCDataPoint * datapoint) {
	if (!datapoint)
		return false;

	if (fNValuesPerPoint == 0 and fDataVector.empty())
		SetNValuesPerPoint(datapoint->GetNValues());

	if (datapoint->GetNValues() != GetNValuesPerPoint())
		return false;

	fDataVector.push_back(datapoint);

	for (unsigned i=0; i<GetNValuesPerPoint(); ++i) {
		// check lower bound
		if (fDataVector.back()->GetValue(i) < fLowerBounds.GetValue(i))
			fLowerBounds.SetValue(i,fDataVector.back()->GetValue(i));
		// check upper bound
		if (fDataVector.back()->GetValue(i) > fUpperBounds.GetValue(i))
			fUpperBounds.SetValue(i,fDataVector.back()->GetValue(i));
	}

	return true;
}

// ---------------------------------------------------------
void BCDataSet::SetNValuesPerPoint(unsigned n) {
	fNValuesPerPoint = n;
	fLowerBounds.SetNValues(n,std::numeric_limits<double>::infinity());
	fUpperBounds.SetNValues(n,-std::numeric_limits<double>::infinity());
	fUserLowerBounds.SetNValues(n,std::numeric_limits<double>::infinity());
	fUserUpperBounds.SetNValues(n,-std::numeric_limits<double>::infinity());
	fFixed.assign(n,false);
}

// ---------------------------------------------------------
void BCDataSet::SetBounds(unsigned index, double lower_bound, double upper_bound, bool fixed) {
	if (index >= GetNValuesPerPoint()) {
		BCLog::OutError("BCDataSet::SetBounds : index out of range.");
		return;
	}
	if (lower_bound >= upper_bound) {
		BCLog::OutWarning("BCDataSet::SetBounds : lower bound is greater than or equal to upper_bound.");
		return;
	}
	fUserLowerBounds.SetValue(index,lower_bound);
	fUserUpperBounds.SetValue(index,upper_bound);
	fFixed[index] = fixed;
}

// ---------------------------------------------------------
void BCDataSet::Reset() {
	for (int i=fDataVector.size()-1; i>=0; --i)
		if (fDataVector[i])
			delete fDataVector[i];
	fDataVector.clear();
	SetNValuesPerPoint(0);
}

// ---------------------------------------------------------
void BCDataSet::Dump(void (*output)(const char *)) const {
	output("Data set summary:");
	output(Form("Number of points           : %u", GetNDataPoints()));
	output(Form("Number of values per point : %u", GetNValuesPerPoint()));
	for (unsigned i=0; i<fDataVector.size(); ++i)
		if (fDataVector[i]) {
			output(Form("Data point %5u",i));
			fDataVector[i] -> Dump(output);
		}
}

