#include "BCModelManager.h"
#include "BCLog.h"
#include "BCErrorCodes.h"

#include <fstream>

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

BCModelManager::BCModelManager(const BCModelManager & modelmanager)
{

	modelmanager.Copy(*this); 

}

// --------------------------------------------------------- 

BCModelManager & BCModelManager::operator = (const BCModelManager & modelmanager)
{

	if (this != &modelmanager) 
		modelmanager.Copy(* this); 

	return * this; 

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

void BCModelManager::SetSingleDataPoint(BCDataPoint * datapoint) 
{

	// create new data set consisting of a single data point 

	BCDataSet * dataset = new BCDataSet(); 

	// add the data point 

	dataset -> AddDataPoint(datapoint); 

	// set this new data set 

	this -> SetDataSet(dataset); 

}

// ---------------------------------------------------------

void BCModelManager::SetSingleDataPoint(BCDataSet * dataset, int index)
{

	if (index < 0 || index > dataset -> GetNDataPoints())
		return; 

	this -> SetSingleDataPoint(dataset -> GetDataPoint(index)); 

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
// DEBUG DELETE?  
/*
void BCModelManager::SetNIterationsMax(int niterations)
{

	// set maximum number of iterations of all models in the manager 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetNIterationsMax(niterations); 

}

// --------------------------------------------------------- 
*/ 

void BCModelManager::SetIntegrationMethod(BCIntegrate::BCIntegrationType method)
{

	// set integration method for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetIntegrationMethod(method); 

}

// --------------------------------------------------------- 

void BCModelManager::SetMarginalizationMethod(BCIntegrate::BCMarginalizationType method)
{ 

		// set marginalization method for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetMarginalizationMethod(method); 

};

// --------------------------------------------------------- 

void BCModelManager::SetModeFindingMethod(BCIntegrate::BCModeFindingType method)
{

	// set mode finding method for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetModeFindingMethod(method); 
	
}

// --------------------------------------------------------- 

void BCModelManager::SetNiterationsPerDimension(int niterations)
{

	// set number of iterations per dimension for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetNiterationsPerDimension(niterations); 

}

// --------------------------------------------------------- 

void BCModelManager::SetNSamplesPer2DBin(int n)
{

	// set samples per 2d bin for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetNSamplesPer2DBin(n); 

}

// --------------------------------------------------------- 

void BCModelManager::SetRelativePrecision(double relprecision) 
{

	// set relative precision for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetRelativePrecision(relprecision); 

}

// --------------------------------------------------------- 

void BCModelManager::SetNbins(int n)
{

	// set number of bins for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetNbins(n); 

}

// --------------------------------------------------------- 

void BCModelManager::SetFitFunctionIndexX(int index) 
{

	// set fit function x index for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetFitFunctionIndexX(index); 

}

// --------------------------------------------------------- 

void BCModelManager::SetFitFunctionIndexY(int index) 
{

	// set  fit function y index for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetFitFunctionIndexY(index); 

}

// --------------------------------------------------------- 

void BCModelManager::SetFitFunctionIndices(int indexx, int indexy)
{

	// set fit function indices for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetFitFunctionIndices(indexx, indexy); 

}

// --------------------------------------------------------- 

void BCModelManager::SetDataPointLowerBoundaries(BCDataPoint* datasetlowerboundaries)
{

	// set lower boundary point for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetDataPointLowerBoundaries(datasetlowerboundaries); 

}

// --------------------------------------------------------- 

void BCModelManager::SetDataPointUpperBoundaries(BCDataPoint* datasetupperboundaries)
{

	// set upper boundary point for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetDataPointUpperBoundaries(datasetupperboundaries); 

}

// --------------------------------------------------------- 

void BCModelManager::SetDataPointLowerBoundary(int index, double lowerboundary)
{

	// set lower bounday values for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetDataPointLowerBoundary(index, lowerboundary); 

}

// --------------------------------------------------------- 

void BCModelManager::SetDataPointUpperBoundary(int index, double upperboundary)
{

	// set upper boundary values for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetDataPointUpperBoundary(index, upperboundary); 

}

// --------------------------------------------------------- 

void BCModelManager::SetDataBoundaries(int index, double lowerboundary, double upperboundary) 
{
	
	// set lower and upper boundary values for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> SetDataBoundaries(index, lowerboundary, upperboundary); 

}

// --------------------------------------------------------- 

void BCModelManager::SetNChains(int n) 
{
	
	// set number of Markov chains for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
	  this -> GetModel(i) -> MCMCSetNChains(n); 

}

// --------------------------------------------------------- 

void BCModelManager::SetFlagPCA(bool flag)
{
	
	// sets the flag for PCA 

	for (int i = 0; i < this -> GetNModels(); i++)
	  this -> GetModel(i) -> MCMCSetFlagPCA(flag); 

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

void BCModelManager::Normalize()
{

	// initialize likelihood norm

	double normalization = 0.0;

	BCLog::Out(BCLog::summary, BCLog::summary, "Running normalization of all models.");

	for (int i = 0; i < this -> GetNModels(); i++)
	{
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

void BCModelManager::FindMode()
{

	// finds mode for all models registered 

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> FindMode(); 

}

// ---------------------------------------------------------

void BCModelManager::MarginalizeAll()
{

	// marginalizes all models registered

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> MarginalizeAll(); 

}

// ---------------------------------------------------------

void BCModelManager::WriteMarkovChain(bool flag)
{

	// marginalizes all models registered

	for (int i = 0; i < this -> GetNModels(); i++)
		this -> GetModel(i) -> WriteMarkovChain(flag); 

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
			std::cerr<<"Couldn't open file "<<file<<std::endl;
			return;
		}
		old_buffer = std::cout.rdbuf(out.rdbuf());
	}

	// model summary
	int nmodels = int(fModelContainer -> size());
	std::cout
		<<std::endl
		<<"======================================"<<std::endl
		<<" Summary"<<std::endl
		<<"======================================"<<std::endl
		<<std::endl
		<<" Number of models               : "<<nmodels<<std::endl
		<<std::endl
		<<" - Models:"<<std::endl;

	for (int i = 0; i < nmodels; i++)
		fModelContainer -> at(i) -> PrintSummary();

	// data summary
	std::cout
		<<" - Data:"<<std::endl
		<<std::endl
		<<"     Number of entries: "<<fDataSet -> GetNDataPoints()<<std::endl
		<<std::endl;

	std::cout
		<<"======================================"<<std::endl
		<<" Model comparison"<<std::endl
		<<std::endl;

	// probability summary
	std::cout
		<<" - A priori probabilities:"<<std::endl
		<<std::endl;
  
	for (int i=0; i<nmodels; i++)
		std::cout
			<<"     p("<< fModelContainer -> at(i) -> GetName()
			<<") = "<< fModelContainer -> at(i) -> GetModelAPrioriProbability()
			<<std::endl;
	std::cout<<std::endl;

	std::cout
		<<" - A posteriori probabilities:"<<std::endl
		<<std::endl;

	for (int i = 0; i < nmodels; i++)
		std::cout
			<<"     p("<< fModelContainer -> at(i) -> GetName()
			<<"|data) = "<< fModelContainer -> at(i) -> GetModelAPosterioriProbability()
			<<std::endl;
	std::cout<<std::endl;

	std::cout
		<<"======================================"<<std::endl
		<<std::endl;

	if (file)
		std::cout.rdbuf(old_buffer);

}

// --------------------------------------------------------- 

void BCModelManager::Copy(BCModelManager & modelmanager) const 
{

	// don't copy the content only the pointers 

	modelmanager.fModelContainer = this -> fModelContainer; 
	modelmanager.fDataSet        = this -> fDataSet; 

}

// ---------------------------------------------------------
