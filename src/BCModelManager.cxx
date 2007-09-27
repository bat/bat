#include "BCModelManager.h" 
#include "BCLog.h" 
#include "BCErrorCodes.h"

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

void BCModelManager::SetDataSet(BCDataSet * dataset) 
{

	// set data set 

	fDataSet = dataset; 

	// set data set of all models in the manager 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetDataSet(fDataSet); 

}

// ---------------------------------------------------------

void BCModelManager::AddModel(BCModel * model, double probability) 
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

int BCModelManager::ReadDataFromFileTree(const char * filename, const char * treename, const char * branchnames)
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

int BCModelManager::ReadDataFromFileUser(const char * filename, std::vector<int> options_int, std::vector<double> options_double) 
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

int BCModelManager::ReadDataFromFileTxt(const char * filename, int nbranches)
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

	BCLog::Out(BCLog::summary, BCLog::summary, "Running normalization of all models.");

	for (int i = 0; i < this -> GetNModels(); i++)
	{
		// calculate model likelihood normalization
//		BCLog::Out(BCLog::summary, BCLog::summary, Form("BCModelManager::Initialize. Normalize model %s.",
//					this -> GetModel(i) -> GetName()));

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

void BCModelManager::PrintSummary(const char * file)
{

	ofstream out;
	std::streambuf * old_buffer = 0;

	if(file)
	{
		out.open(file);
		if (!out.is_open())
		{
			cerr<<"Couldn't open file "<<file<<endl;
			return;
		}
		old_buffer = std::cout.rdbuf(out.rdbuf());
	}

	// model summary
	int nmodels = int(fModelContainer -> size());
	cout
		<<endl
		<<"======================================"<<endl
		<<" Summary"<<endl
		<<"======================================"<<endl
		<<endl
		<<" Number of models               : "<<nmodels<<endl
		<<endl
		<<" - Models:"<<endl;

	for (int i = 0; i < nmodels; i++)
		fModelContainer -> at(i) -> PrintSummary();

	// data summary
	cout
		<<" - Data:"<<endl
		<<endl
		<<"     Number of entries: "<<fDataSet -> GetNDataPoints()<<endl
		<<endl;

	cout
		<<"======================================"<<endl
		<<" Model comparison"<<endl
		<<endl;

	// probability summary
	cout
		<<" - A priori probabilities:"<<endl
		<<endl;
  
	for (int i=0; i<nmodels; i++)
		cout
			<<"     p("<< fModelContainer -> at(i) -> GetName()
			<<") = "<< fModelContainer -> at(i) -> GetModelAPrioriProbability()
			<<endl;
	cout<<endl;

	cout
		<<" - A posteriori probabilities:"<<endl
		<<endl;

	for (int i = 0; i < nmodels; i++)
		cout
			<<"     p("<< fModelContainer -> at(i) -> GetName()
			<<"|data) = "<< fModelContainer -> at(i) -> GetModelAPosterioriProbability()
			<<endl;
	cout<<endl;

	cout
		<<"======================================"<<endl
		<<endl;

	cout.rdbuf(old_buffer);

}

// ---------------------------------------------------------
