#include "BCModelManager.h" 
#include "BCLog.h" 

#include <TFile.h> 
#include <TH2D.h> 
#include <TTree.h> 
#include <TString.h> 

#include <fstream.h>

// --------------------------------------------------------- 

BCModelManager::BCModelManager()
{

  fModelContainer = new BCModelContainer(); 

  fDataSet = 0; 

}

// ---------------------------------------------------------

BCModelManager::~BCModelManager()
{

  delete fModelContainer; 

  if (fDataSet) 
    delete fDataSet; 

}

// ---------------------------------------------------------

void BCModelManager::SetDataSet(BCDataSet* dataset) 
{

  // set data set 

  fDataSet = dataset; 

  // set data set of all models in the manager 

  for (int i = 0; i < this -> GetNModels(); i++)
    this -> GetModel(i) -> SetDataSet(fDataSet); 

}

// ---------------------------------------------------------

void BCModelManager::AddModel(BCModel* model, double probability) 
{

  // create index 

  int index = int(fModelContainer -> size()); 

  // set index of new model 

  model -> SetIndex(index); 

  // set a priori probability of new model 

  model -> SetModelAPrioriProbability(probability); 

  // set data set 

  model -> SetDataSet(fDataSet); 

  // fill model into container 

  fModelContainer -> push_back(model); 

}   

// ---------------------------------------------------------

void BCModelManager::SetNIterationsMax(int niterations)
{
  // set maximum number of iterations of all models in the manager 

  for (int i = 0; i < this -> GetNModels(); i++)
    this -> GetModel(i) -> SetNIterationsMax(niterations); 

}

// --------------------------------------------------------- 

int BCModelManager::ReadDataFromFileTree(char* filename, char* treename, const char* branchnames)
{

  if (fModelContainer -> size() < 0) 
    {
      BCLog::Out(BCLog::warning, BCLog::warning, "BCModelManager::ReadDataFromFileTree. No model defined."); 
      return ERROR_NOMODELS;  
    }
    
  // create data set 

  if (!fDataSet)
    fDataSet = new BCDataSet(); 

  else 
    fDataSet -> Reset(); 

  // read data from tree 

  int read_file = fDataSet -> ReadDataFromFileTree(filename, treename, branchnames); 

  if (read_file >=0) 
    {
      this -> SetDataSet(fDataSet); 
      
      for (int i = 0; i < this -> GetNModels(); i++)
	fModelContainer -> at(i) -> SetDataSet(fDataSet); 
    }

  else if (read_file == ERROR_FILENOTFOUND) 
    {
      delete fDataSet; 
      
      return ERROR_FILENOTFOUND; 
    }
  
  return 0; 

}

// --------------------------------------------------------- 

int BCModelManager::ReadDataFromFileUser(char* filename, std::vector<int> options_int, std::vector<double> options_double) 
{
  
  if (fModelContainer -> size() < 0) 
    {
      BCLog::Out(BCLog::warning, BCLog::warning, "BCModelManager::ReadDataFromFileTree. No model defined."); 
      return ERROR_NOMODELS;  
    }
    
  // create data set 

  if (!fDataSet) 
    fDataSet = new BCDataSet(); 
  else
    fDataSet -> Reset(); 
  
  // read data from user specified file 

  int read_file = fDataSet -> ReadDataFromFileUser(filename, options_int, options_double, ""); 

  if (read_file >=0) 
    {
      this -> SetDataSet(fDataSet); 
      
      for (int i = 0; i < this -> GetNModels(); i++)
	fModelContainer -> at(i) -> SetDataSet(fDataSet); 
    }

  else
    {
      delete fDataSet; 
      
      return ERROR_FILENOTFOUND; 
    }

  return 0; 

}

// --------------------------------------------------------- 

int BCModelManager::ReadDataFromFileTxt(char* filename, int nbranches)
{

  if (fModelContainer -> size() < 0) 
    {
      BCLog::Out(BCLog::warning, BCLog::warning, "BCModelManager::ReadDataFromFileTree. No model defined."); 
      return ERROR_NOMODELS;  
    }
    
  // create data set 

  if (!fDataSet) 
    fDataSet = new BCDataSet(); 
  else
    fDataSet -> Reset(); 

  // read data from txt file 

  int read_file = fDataSet -> ReadDataFromFileTxt(filename, nbranches); 

  if (read_file >=0) 
    {
      this -> SetDataSet(fDataSet); 
      
      for (int i = 0; i < this -> GetNModels(); i++)
	fModelContainer -> at(i) -> SetDataSet(fDataSet); 
    }

  else
    {
      delete fDataSet; 
      
      return ERROR_FILENOTFOUND; 
    }

  return 0; 

}

// ---------------------------------------------------------

void BCModelManager::Initialize()
{
	// initialize likelihood norm
	double normalization = 0.0;

	for (int i = 0; i < this -> GetNModels(); i++)
	{
		// calculate model likelihood normalization
		BCLog::Out(BCLog::summary, BCLog::summary, Form("BCModelManager::Initialize. Normalize model %s.",
					this -> GetModel(i) -> GetName()));

		fModelContainer -> at(i) -> Normalize(); 

		// add to total normalization 
		normalization += (fModelContainer -> at(i) -> GetNormalization() * 
			fModelContainer -> at(i) -> GetModelAPrioriProbability()); 
	}

	// set model a posteriori probabilities 
	for (int i = 0; i < int(fModelContainer -> size()); i++)
		fModelContainer -> at(i) -> SetModelAPosterioriProbability(
					(fModelContainer -> at(i) -> GetNormalization() * 
					fModelContainer -> at(i) -> GetModelAPrioriProbability()) / 
					normalization); 
}

// ---------------------------------------------------------

void BCModelManager::PrintSummary()
{

  // model summary 

  int nmodels = int(fModelContainer -> size()); 

  std::cout << std::endl; 
  std::cout << " ===================================== " << std::endl; 
  std::cout << " Summary                               " << std::endl; 
  std::cout << " ===================================== " << std::endl; 
  std::cout << std::endl;   
  std::cout << " Number of models               : " << nmodels << std::endl; 
  std::cout << std::endl; 
  std::cout << " - Models: " << std::endl; 
  std::cout << std::endl; 

  for (int i = 0; i < nmodels; i++)
    fModelContainer -> at(i) -> PrintSummary(); 

  // data summary 

  std::cout << " - Data: " << std::endl; 
  std::cout << std::endl; 
  std::cout << "   Number of entries: " << fDataSet -> GetNDataPoints() << std::endl; 
  std::cout << std::endl; 
  
  // probability summary 

  std::cout << " - A priori probabilities: " << std::endl; 
  std::cout << std::endl; 
  
  for (int i = 0; i < nmodels; i++)
    std::cout << " p(" << fModelContainer -> at(i) -> GetName() 
	      << ") = " << fModelContainer -> at(i) -> GetModelAPrioriProbability() 
	      << std::endl; 
  std::cout << std::endl;       
  
  std::cout << " - A posteriori probabilities: " << std::endl; 
  std::cout << std::endl; 
  
  for (int i = 0; i < nmodels; i++)
    std::cout << " p(" << fModelContainer -> at(i) -> GetName() 
	      << "|data) = " << fModelContainer -> at(i) -> GetModelAPosterioriProbability() 
	      << std::endl; 
  std::cout << std::endl; 

}

// ---------------------------------------------------------

void BCModelManager::PrintSummaryToFile(char* filename)
{

  // open file 

  std::ofstream file; 

  file.open(filename); 

  // check if file is open 

  if (file.is_open() == false)
    std::cout << " File not open. " << std::endl; 

  // redirect 

  std::streambuf *old_buffer = std::cout.rdbuf(file.rdbuf());

  // print summary 

  this -> PrintSummary(); 

  // restore the old buffer

  std::cout.rdbuf(old_buffer);

}

// ---------------------------------------------------------
