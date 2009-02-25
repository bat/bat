/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BAT/BCDataSet.h"

#include "BAT/BCDataPoint.h"
#include "BAT/BCLog.h"
#include "BAT/BCErrorCodes.h"

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <iostream>
#include <fstream>

// ---------------------------------------------------------

BCDataSet::BCDataSet()
{
	fBCDataVector = 0;
}

// ---------------------------------------------------------

BCDataSet::~BCDataSet()
{
	if (fBCDataVector)
		delete fBCDataVector;
}

// ---------------------------------------------------------

unsigned int BCDataSet::GetNDataPoints()
{
	// check if vector exists. Return number of data points if true ...
	if (fBCDataVector)
		return fBCDataVector -> size();

	// ... or give out warning and return 0 if not.
	BCLog::Out(BCLog::warning, BCLog::warning,"BCDataSet::GetNDataPoints : DataSet not yet created.");
	return 0;
}

// ---------------------------------------------------------

unsigned int BCDataSet::GetNValuesPerPoint()
{
	// check if vector exists and contains datapoints
	if (fBCDataVector && fBCDataVector -> size() > 0)
		return this -> GetDataPoint(0) -> GetNValues();

	BCLog::Out(BCLog::error, BCLog::error,
			"BCDataSet::GetNValuesPerPoint : Data set doesn't exist yet");
	return 0;
}

// ---------------------------------------------------------

BCDataPoint * BCDataSet::GetDataPoint(unsigned int index)
{
	if (!fBCDataVector || this -> GetNDataPoints()==0 )
	{
		BCLog::Out(BCLog::error, BCLog::error,"BCDataSet::GetDataPoint : Dataset is empty.");
		return 0;
	}

	// check if index is within range. Return the data point if true ...
	if(index >= 0 && index < this -> GetNDataPoints())
		return fBCDataVector -> at(index);

	// ... or give out warning and return 0 if not.
	BCLog::Out(BCLog::error, BCLog::error,"BCDataSet::GetDataPoint : Index out of range. Return 0.");
	return 0;
}

// ---------------------------------------------------------

int BCDataSet::ReadDataFromFileTree(const char * filename, const char * treename, const char * branchnames)
{
	// open root file
	TFile * file = new TFile(filename, "READ");

	// check if file is open and warn if not.
	if (!file -> IsOpen())
	{
		BCLog::Out(BCLog::error, BCLog::error,
				Form("BCDataSet::ReadDataFromFileTree : Could not open file %s.", filename));
		return ERROR_FILENOTFOUND;
	}

	// get tree
	TTree * tree = (TTree*) file -> Get(treename);

	// check if tree is there and warn if not.
	if (!tree)
	{
		BCLog::Out(BCLog::error, BCLog::error,
				Form("BCDataSet::ReadDataFromFileTree : Could not find TTree %s.", treename));

		// close file
		file -> Close();

		return ERROR_TREENOTFOUND;
	}

	// if data set contains data, clear data object container ...
	if (fBCDataVector != 0)
	{
		fBCDataVector -> clear();

		BCLog::Out(BCLog::detail, BCLog::detail,"BCDataSet::ReadDataFromFileTree : Overwrite existing data.");
	}

	// ... or allocate memory for the vector if not.
	else
		fBCDataVector = new BCDataVector();

	// get branch names.

	// first, copy the branchnames into a std::string.
	std::string branches(branchnames);

	// define a vector of std::strings which contain the tree names.
	std::vector<std::string> * branchnamevector = new std::vector<std::string>;

	// the names are supposed to be separated by commas. find first comma
	// entry in the string.
	int temp_index = branches.find_first_of(",");

	// reset number of branches
	int nbranches = 0;

	// repeat until the is nothing left in the string.
	while(branches.size() > 0)
	{
		// temporary string which contains the name of the current branch
		std::string branchname;

		// get current branch name

		// if there is no comma the current branchname corresponds to the whole string, ...
		if (temp_index == -1)
			branchname = branches;

		// ... if there is a comma, copy that part of the string into the current branchname.
		else
			branchname.assign(branches, 0, temp_index);

		// write branch name to a vector
		branchnamevector -> push_back(branchname);

		// increase the number of branches found
		nbranches++;

		// cut remaining string with branchnames

		// if there is no comma left empty the string, ...
		if (temp_index == -1)
				branches = "";

		// ... if there is a comma remove the current branchname from the string.
		else
			branches.erase(0, temp_index + 1);

		// find the next comma
		temp_index = branches.find_first_of(",");
	}

	// create temporary vector with data and assign some zeros.
	std::vector<double> data;
	data.assign(nbranches, 0.0);

	// set the branch address.
	for (int i = 0; i < nbranches; i++)
		tree -> SetBranchAddress(branchnamevector -> at(i).data(), &data.at(i));

	// calculate maximum number of entries
	int nentries = tree -> GetEntries();

	// check if there are any events in the tree and close file if not.
	if (nentries <= 0)
	{
		BCLog::Out(BCLog::error, BCLog::error,
				Form("BCDataSet::ReadDataFromFileTree : No events in TTree %s.", treename));

		// close file
		file -> Close();

		return ERROR_NOEVENTS;
	}

	// loop over entries
	for (int ientry = 0; ientry < nentries; ientry++)
	{
		// get entry
		tree -> GetEntry(ientry);

		// create data object
		BCDataPoint * datapoint = new BCDataPoint(nbranches);

		// copy data

		for (int i = 0; i < nbranches; i++)
			datapoint -> SetValue(i, data.at(i));

		// add data point to this data set.
		this -> AddDataPoint(datapoint);
	}

	// close file
	file -> Close();

	// remove file pointer.
	if (file)
		delete file;

	return 0;

}

// ---------------------------------------------------------

int BCDataSet::ReadDataFromFileTxt(const char * filename, int nbranches)
{
	// open text file.
	std::fstream file;
	file.open(filename, std::fstream::in);

	// check if file is open and warn if not.
	if (!file.is_open())
	{
		BCLog::Out(BCLog::error, BCLog::error,
				Form("BCDataSet::ReadDataFromFileText : Could not open file %s.", filename));

		return ERROR_FILENOTFOUND;
	}

	// if data set contains data, clear data object container ...
	if (fBCDataVector != 0)
	{
		fBCDataVector -> clear();

		BCLog::Out(BCLog::detail, BCLog::detail,"BCDataSet::ReadDataFromFileTxt : Overwrite existing data.");
	}

	// ... or allocate memory for the vector if not.
	else
		fBCDataVector = new BCDataVector();

	// create temporary vector with data and assign some zeros.
	std::vector<double> data;
	data.assign(nbranches, 0.0);

	// reset counter
	int nentries = 0;

	// read data and create data points.
	while (!file.eof())
	{
		// read data from file
		int i=0;
		while(file >> data[i])
		{
			if (i==nbranches-1)
				break;
			i++;
		}

		// create data point.
		if(i == nbranches-1)
		{
			BCDataPoint * datapoint = new BCDataPoint(nbranches);

			// copy data into data point
			for (int i = 0; i < nbranches; i++)
				datapoint -> SetValue(i, data.at(i));

			// add data point to this data set.
			this -> AddDataPoint(datapoint);

			// increase counter
			nentries++;
		}
	}

	// check if there are any events in the tree and close file if not.
	if (nentries <= 0)
	{
		BCLog::Out(BCLog::error, BCLog::error,
				Form("BCDataSet::ReadDataFromFileText : No events in the file %s.", filename));

		// close file
		file.close();

		return ERROR_NOEVENTS;
	}

	// close file
	file.close();

	return 0;

}

// ---------------------------------------------------------

void BCDataSet::AddDataPoint(BCDataPoint * datapoint)
{

	// check if memory for the vector has been allocated and
	// allocate if not.
	if (fBCDataVector == 0)
		fBCDataVector = new BCDataVector();

	// add data point to the data set.
	fBCDataVector -> push_back(datapoint);

}

// ---------------------------------------------------------

void BCDataSet::Reset()
{

	// if memory has been allocated to the data set
	// clear the content.
	if (fBCDataVector != 0)
		fBCDataVector -> clear();

}

// ---------------------------------------------------------

void BCDataSet::Dump()
{
	if (!fBCDataVector)
	{
		BCLog::Out(BCLog::error, BCLog::error, "BCDataSet::Dump : Data set is empty. Nothing to dump.");
		return;
	}

	std::cout << std::endl
		<< "Dumping dataset:" << std::endl
		<< "----------------" << std::endl
		<< " - number of points:            " << fBCDataVector -> size() << std::endl
		<< " - number of values per point:  " << this -> GetDataPoint(0) -> GetNValues() << std::endl
		<< " - values:" << std::endl;
	unsigned int n = this -> GetDataPoint(0) -> GetNValues();
	for (unsigned int i=0; i< fBCDataVector -> size(); i++)
	{
		std::cout << Form("%5d :  ", i);
		for (unsigned int j=0; j<n; j++)
			std::cout << Form("%12.5g", this -> GetDataPoint(i) -> GetValue(j));
		std::cout << std::endl;
	}
	std::cout << std::endl;

}


// ---------------------------------------------------------
