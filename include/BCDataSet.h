/*! \class BCDataSet
 *  \brief A class for handling a set of data points
 *
 *
 * The class defines a data set. It contains a vector of data points
 * and methods to read data from a file.
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, kroening *at* mppmu *dot* mppmu *dot* de 
 *
 * CREATED: 06.06.2007 
 * 
 * REVISION: 
 *
 * --------------------------------------------------------- 
 */ 

// --------------------------------------------------------- 

#ifndef __BCDATASET__H
#define __BCDATASET__H

#include <vector.h> 

#include "BCDataPoint.h" 
#include "BCErrorCodes.h" 

// --------------------------------------------------------- 

typedef std::vector<BCDataPoint*> BCDataVector; 

// --------------------------------------------------------- 

class BCDataSet
{

 public: 

  /* 
   * Default constructor 
   */ 
  BCDataSet(); 

  /* 
   * Default destructor 
   */ 
  virtual ~BCDataSet(); 

  // getter 

  /* 
   * @return The vector of data points 
   */ 
  BCDataVector* GetDataVector(); 

  /* 
   * @return The number of data points in the data set 
   */
  int GetNDataPoints(); 

  /* 
   * @return The index of the data point to be returned 
   */ 
  BCDataPoint* GetDataPoint(int index); 

  // methods 

  /** 
   * Reads a TTree from a .root file. 
   * Opens a .root file and gets a TTree. It creates data points 
   * containing the values read from the file. 
   * @param filename The name of the .root file 
   * @param treename The name of the TTree 
   * @param branchnames A list of the names of the branches separated by a comma
   * @return An error code 
   * @see ReadDataFromFileHist(char* filename, char* histname); 
   * @see ReadDataFromFileTxt(char* filename, int nbranches); 
   */ 
  int ReadDataFromFileTree(char* filename, char* treename, const char* branchnames); 

  /** 
   * Reads data from a .txt file. 
   * Opens a .txt file and creates data objects 
   * containing the values read from the file. 
   * @param filename The name of the .txt file 
   * @param nvariables The number of variables 
   * @see ReadDataFromFileTree(char* filename, char* treename, std::vector<char*> branchnames)
   * @see ReadDataFromFileHist(char* filename, char* histname); 
   */
  int ReadDataFromFileTxt(char* filename, int nvariables);

  /** 
   * Reads data from a user specified file. 
   * Opens a user specified file and creates data objects 
   * containing the values read from the file. 
   * @param filename The name of the file 
   * @param options_int A vector of options of type int 
   * @param options_double A vector of options of type double 
   * @param options_char A pointer of characters 
   * @see ReadDataFromFileTree(char* filename, char* treename, std::vector<char*> branchnames)
   * @see ReadDataFromFileHist(char* filename, char* histname, const char* branchnames); 
   * @see ReadDataFromFileTxt(char* filename, int nbranches); 
   */   
  virtual int ReadDataFromFileUser(char* filename, std::vector<int> options_int, std::vector<double> options_double, const char* options_char); 

  /**
   * Adds a data point to the data set 
   */ 
  void AddDataPoint(BCDataPoint* datapoint); 

  /**
   * Resets the content of the data set 
   */ 
  void Reset(); 
  
 private: 

  BCDataVector* fBCDataVector; 

}; 

// --------------------------------------------------------- 

#endif 

