#include "BCDataSet.h" 
#include "BCLog.h" 
#include "BCErrorCodes.h" 

#include <TFile.h> 
#include <TTree.h> 

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

int BCDataSet::GetNDataPoints()
{

	int ndatapoints = 0; 

	// check if vector exists. Return number of data points if true ... 

	if (fBCDataVector) 
		ndatapoints = int(fBCDataVector -> size()); 

	// ... or give out warning and return 0 if not. 

	else 
		BCLog::Out(BCLog::warning, BCLog::warning,"BCDataSet::GetNDataPoints. Memory for vector fBCDataVector not yet allocated. Return 0."); 

	return ndatapoints; 

}

// --------------------------------------------------------- 

BCDataPoint * BCDataSet::GetDataPoint(int index)
{

	BCDataPoint * datapoint = 0; 

	// check if index is within range. Return the data point if true ... 

	if (fBCDataVector && index >= 0 && index < this -> GetNDataPoints())
		datapoint = fBCDataVector -> at(index); 

	// ... or give out warning and return 0 if not. 

	else
		BCLog::Out(BCLog::warning, BCLog::warning,"BCDataSet::GetDataPoint. Index out of range. Return 0."); 

	return datapoint; 

}

// --------------------------------------------------------- 

int BCDataSet::ReadDataFromFileTree(char * filename, char * treename, const char * branchnames) 
{

	// if data set contains data, clear data object container ... 

	if (fBCDataVector != 0)
		{
			fBCDataVector -> clear();

			BCLog::Out(BCLog::detail, BCLog::detail,"BCDataSet::ReadDataFromFileTree. Overwrite existing data."); 
		}

	// ... or allocate memory for the vector if not. 

	else
		fBCDataVector = new BCDataVector(); 

	// open root file 

	TFile * file = new TFile(filename, "READ"); 

	// check if file is open and warn if not. 

	if (!file -> IsOpen())
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 Form("BCDataSet::ReadDataFromFileTree. Could not open file %s.", filename)); 

			return ERROR_FILENOTFOUND; 
		}

	// get tree

	TTree * tree = (TTree*) file -> Get(treename); 

	// check if tree is there and warn if not. 

	if (!tree) 
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 Form("BCDataSet::ReadDataFromFileTree. Could not find TTree %s.", treename)); 

			// close file 

			file -> Close(); 

			return ERROR_TREENOTFOUND; 
		}

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
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 Form("BCDataSet::ReadDataFromFileTree. No events in TTree %s.", treename)); 
    
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

int BCDataSet::ReadDataFromFileTxt(char* filename, int nbranches) 
{

	// if data set contains data, clear data object container ... 

	if (fBCDataVector != 0)
		{
			fBCDataVector -> clear();

			BCLog::Out(BCLog::detail, BCLog::detail,"BCDataSet::ReadDataFromFileTxt. Overwrite existing data."); 
		}

	// ... or allocate memory for the vector if not. 

	else
		fBCDataVector = new BCDataVector(); 

	// open text file. 

	std::fstream file; 

	file.open(filename, std::fstream::in); 

	// check if file is open and warn if not. 

	if (!file.is_open())
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 Form("BCDataSet::ReadDataFromFileText. Could not open file %s.", filename)); 

			return ERROR_FILENOTFOUND; 
		}

	// create temporary vector with data and assign some zeros. 

	std::vector<double> data; 
	data.assign(nbranches, 0.0); 

	// reset counter 

	int nentries = 0; 

	// read data and create data points. 

	while (!file.eof())
		{
			// read data from file 

			for (int i = 0; i < nbranches; i++)
				file >> data[i]; 

			// create data point. 

			BCDataPoint * datapoint = new BCDataPoint(nbranches); 

			// copy data into data point 

			for (int i = 0; i < nbranches; i++)
				datapoint -> SetValue(i, data.at(i)); 

			// add data point to this data set. 

			this -> AddDataPoint(datapoint); 

			// increase counter 

			nentries++; 
		}

	// check if there are any events in the tree and close file if not. 

	if (nentries <= 0) 
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 Form("BCDataSet::ReadDataFromFileText. No events in the file %s.", filename)); 
    
			// close file 

			file.close(); 

			return ERROR_NOEVENTS; 
		}

	// close file 

	file.close(); 

	return 0; 

}  

// --------------------------------------------------------- 

int BCDataSet::ReadDataFromFileUser(char * filename, std::vector<int> options_int, std::vector<double> options_double, const char * options_char)
{

	// if this method is called without being overloaded give out warning. 

	BCLog::Out(BCLog::warning, BCLog::warning, 
						 "BCDataSet::ReadDataFromFileUser. Function needs to be overloaded by the user."); 

	return ERROR_METHODNOTOVERLOADED; 

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

}; 

// --------------------------------------------------------- 

void BCDataSet::Reset()
{ 

	// if memory has been allocated to the data set 
	// clear the content. 

	if (fBCDataVector != 0) 
		fBCDataVector -> clear(); 

};

// --------------------------------------------------------- 
