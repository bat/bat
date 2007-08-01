/**
 * A class which manages a set of models and a data set 
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, kroening *at* mppmu *dot* mppmu *dot* de 
 *
 * CREATED: 14.03.2007 by Kevin 
 * 
 * REVISION: 
 *
 * 16.03.2007 Kevin,  modified ROOT interface
 * 12.06.2007 Kevin,  renamed to BCModelManager 
 * 30.08.2007 Dano,   added method to set max number of iterations in all
 *                    attached models
 *
 * --------------------------------------------------------- 
 *
 *
 * The class manages a set of models and a data set 
 *
 */ 

// --------------------------------------------------------- 

#ifndef __BCMODELMANAGER__H
#define __BCMODELMANAGER__H

#include <TROOT.h>

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelManager
{

 public:

  // constructor and destructor 

  /**
   * A constructor. 
   */ 
  BCModelManager(); 

  /**
   * A destructor. 
   */ 
  virtual ~BCModelManager(); 

  // methods (get) 

  /** 
   * @return The number of models. 
   */ 
  int GetNModels() 
    { return int(fModelContainer -> size()); }; 

  /** 
   * @param index The index of the model 
   * @return The BCModel. 
   */ 
  BCModel* GetModel(int index) 
    { return fModelContainer -> at(index); }; 

  /**
   * @return The number of entries in the data container
   */ 
  int GetNDataPoints() 
    { return fDataSet -> GetNDataPoints(); }; 

  /**
   * @param index The index of the data point
   * @return The data point
   */ 
  BCDataPoint * GetDataPoint(int index) 
    { return fDataSet -> GetDataPoint(index); }; 

  /** 
   * @return The data set 
   */ 
  BCDataSet * GetDataSet()
    { return fDataSet; }; 

  // methods (set) 

   /**
    * @param dataset A data set 
    */ 
  void SetDataSet(BCDataSet * dataset); 

   /**
    * Set maximum number of iterations for all models added to the manager.
	 * Only works on models already added to the manager.
    * @param niterations
    */ 
  void SetNIterationsMax(int niterations);

  // methods 

  /** 
   * Adds a model to the container
   * @param model The model
   * @param probability The a priori probability 
   * @see AddModel(BCModel * model)
   * @see SetModelPrior(BCModel * model, double probability)
   */ 
  void AddModel(BCModel * model, double probability); 

  /** 
   * Adds a model to the container
   * @param model The model
   * @see AddModel(BCModel * model, double probability)
   * @see SetModelPrior(BCModel * model, double probability)
   */ 
  void AddModel(BCModel * model)
    { this -> AddModel(model, 0.0); }; 

  /**
   * Adds a data point to the data container. 
   * @param datapoint The data point
   */ 
  void AddDataPoint(BCDataPoint * datapoint) 
    { fDataSet -> AddDataPoint(datapoint); }; 

  /** 
   * Reads tree data from a ROOT file. 
   * Opens a ROOT file and gets a ROOT tree. It creates data set
   * containing the values read from the file. 
   * @param filename The filename of the ROOT file 
   * @param treename The name of the ROOT tree 
   * @param branchnames A vector of the names of the branches 
   * @see ReadDataFromFileHist(char * filename, char * histname, const char*  branchnames); 
   * @see ReadDataFromFileTxt(char * filename, int nbranches); 
   */ 
  int ReadDataFromFileTree(char* filename, char * treename, const char * branchnames); 

  /** 
   * Reads data from a txt file. 
   * Opens a txt file and creates data set
   * containing the values read from the file. 
   * @param filename The filename of the ROOT file 
   * @param nbranches The number of variables 
   * @see ReadDataFromFileTree(char * filename, char * treename, std::vector<char*> branchnames)
   * @see ReadDataFromFileHist(char * filename, char * histname, const char * branchnames); 
   */   
  int ReadDataFromFileTxt(char * filename, int nbranches); 

  /** 
   * Reads data from a txt file (user specifies). 
   * Opens a file and creates data set
   * containing the values read from the file (user specifies). 
   * @param filename The filename of the ROOT file 
   * @param nbranches The number of variables 
   * @see ReadDataFromFileTree(char * filename, char * treename, std::vector<char*> branchnames)
   * @see ReadDataFromFileHist(char * filename, char * histname, const char * branchnames); 
   * @see ReadDataFromFileTxt(char * filename, int nbranches); 
   */   
  int ReadDataFromFileUser(char*  filename, std::vector<int> options_int, std::vector<double> options_double); 
  
 /** 
   * Calculates the normalization of the likelihood for each model in the container. 
   */   
  void Initialize(); 

  /**
   * Resets the data set 
   */ 
  void ResetDataSet()
    { fDataSet -> Reset(); }; 

  /** 
   * Prints a summary on the screen. 
   */ 
  void PrintSummary(); 

  /** 
   * Prints a summary to a file. 
   */ 
  void PrintSummaryToFile(char * filename); 

 private: 

  /**
   * A model container. 
   */ 
  BCModelContainer * fModelContainer; 

  /** 
   * A data set 
   */ 
  BCDataSet * fDataSet; 

}; 

// --------------------------------------------------------- 

#endif 
