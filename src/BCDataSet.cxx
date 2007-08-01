#include "BCDataSet.h" 

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

  if (fBCDataVector) 
    ndatapoints = int(fBCDataVector -> size()); 

  return ndatapoints; 

}

// --------------------------------------------------------- 

BCDataPoint* BCDataSet::GetDataPoint(int index)
{
  
  BCDataPoint* datapoint = 0; 

  if (fBCDataVector && index >= 0 && index < this -> GetNDataPoints())
    datapoint = fBCDataVector -> at(index); 

  return datapoint; 

}

// --------------------------------------------------------- 

int BCDataSet::ReadDataFromFileTree(char* filename, char* treename, const char* branchnames) 
{

  // clear data object container 

  if (fBCDataVector != 0)
    fBCDataVector -> clear();

  else
    fBCDataVector = new BCDataVector(); 

  // open file 

  TFile* file = new TFile(filename, "READ"); 

  // check if file is open 
  
  if (!file -> IsOpen())
    {
      BCLog::Out(BCLog::warning, BCLog::warning, 
		 Form("BCDataSet::ReadDataFromFileTree. Could not open file %s.", filename)); 

      return ERROR_FILENOTFOUND; 
    }

  // get tree

  TTree* tree = (TTree*) file -> Get(treename); 

  // check if tree is there 

  if (!tree) 
    {
      BCLog::Out(BCLog::warning, BCLog::warning, 
		 Form("BCDataSet::ReadDataFromFileTree. Could not find TTree %s.", treename)); 

      // close file 

      file -> Close(); 

      return ERROR_TREENOTFOUND; 
    }

  // get branch names 

  TString branches = branchnames; 

  std::vector<TString> branchnamevector; 

  int temp_index = branches.First(","); 
  
  int nbranches = 0; 

  while(branches.Length() > 0)
    {
      // temporary string which contains the name of the current branch 

      TString branchname; 

      // get current branch name 

      if (temp_index == -1) 
	branchname = branches(0, branches.Length()); 
      else
	branchname = branches(0, temp_index); 

      branchnamevector.push_back(branchname); 

      // increase number of branches 

      nbranches++; 

      // cut remaining string with branchnames 
      
      if (temp_index == -1)
	branches = ""; 
      else
	branches = branches(temp_index + 1, branches.Length()); 

      temp_index = branches.First(","); 
    }     

  // create temporary vector with data 

  std::vector<double> data; 
  data.assign(nbranches, 0.0); 

  // set the branch address 

  for (int i = 0; i < nbranches; i++)
    {
      TString branchname = branchnamevector.at(i); 

      tree -> SetBranchAddress(branchname, &data.at(i)); 
    }

  // calculate maximum number of entries 

  int nentries = tree -> GetEntries(); 

  // check if there are any events in the tree 

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

      BCDataPoint* datapoint = new BCDataPoint(nbranches); 

      // copy data 

      for (int i = 0; i < nbranches; i++)
	datapoint -> SetValue(i, data.at(i)); 

      // add data object to container 

      this -> AddDataPoint(datapoint); 
    }

  // close file 

  file -> Close(); 

  return 0;  

}  

// --------------------------------------------------------- 

int BCDataSet::ReadDataFromFileTxt(char* filename, int nbranches) 
{

  // clear data object container 

  if (fBCDataVector != 0)
    fBCDataVector -> clear();

  else
    fBCDataVector = new BCDataVector(); 

  // open file 

  std::fstream file; 

  file.open(filename, std::fstream::in); 

  // check if file is open 

  if (!file.is_open())
    {
      BCLog::Out(BCLog::warning, BCLog::warning, 
		 Form("BCDataSet::ReadDataFromFileText. Could not open file %s.", filename)); 

      return ERROR_FILENOTFOUND; 
    }

  // create temporary vector with data 

  std::vector<double> data; 
  data.assign(nbranches, 0.0); 

  // reset counter 

  int nentries = 0; 

  // read data and create data objects 

  while (!file.eof())
    {
      // read data 

      for (int i = 0; i < nbranches; i++)
	file >> data[i]; 

      // create data object 

      BCDataPoint* datapoint = new BCDataPoint(nbranches); 

      // copy data 

      for (int i = 0; i < nbranches; i++)
	datapoint -> SetValue(i, data.at(i)); 

      // add data object to container 

      this -> AddDataPoint(datapoint); 

      // increase counter 

      nentries++; 
    }
  
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

int BCDataSet::ReadDataFromFileUser(char* filename, std::vector<int> options_int, std::vector<double> options_double, const char* options_char)
{

  BCLog::Out(BCLog::warning, BCLog::warning, 
	     "BCDataSet::ReadDataFromFileUser. Function needs to be overloaded by the user."); 
  
  return ERROR_METHODNOTOVERLOADED; 
  
}

// --------------------------------------------------------- 

void BCDataSet::AddDataPoint(BCDataPoint* datapoint)
{ 

  if (fBCDataVector == 0) 
    fBCDataVector = new BCDataVector();   

  fBCDataVector -> push_back(datapoint); 

}; 

// --------------------------------------------------------- 

void BCDataSet::Reset()
{ 

  if (fBCDataVector != 0) 
    fBCDataVector -> clear(); 

};

 // --------------------------------------------------------- 
