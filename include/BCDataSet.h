/*! \class BCDataSet
 *  \brief A class for handling a set of data points
 *
 * The class defines a data set. It contains a vector of data \n 
 * points and methods to read data from a file. \n 
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, 
 *          kevin.kroeninger *at* phys *dot* uni *minus* goettingen *dot* de 
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
	 * Returns the vector of data points 
	 */ 
	BCDataVector* GetDataVector(); 

	/* 
	 * Returns the number of data points in the data set 
	 */
	int GetNDataPoints(); 

	/* 
	 * Returns the data point at index. 
	 * @param index The index of the data point to be returned 
	 */ 
	BCDataPoint* GetDataPoint(int index); 

	// methods 

	/** 
	 * Reads a TTree from a .root file. \n 
	 * Opens a .root file and gets a TTree. It creates \n 
	 * data points containing the values read from the file. \n 
	 * @param filename The name of the .root file 
	 * @param treename The name of the TTree 
	 * @param branchnames A list of the names of the branches separated by a comma
	 * @return An error code 
	 * @see ReadDataFromFileHist(char* filename, char* histname); 
	 * @see ReadDataFromFileTxt(char* filename, int nbranches); 
	 */ 
	int ReadDataFromFileTree(const char* filename, char* treename, const char* branchnames); 

	/** 
	 * Reads data from a .txt file. \n 
	 * Opens a .txt file and creates data objects \n 
	 * containing the values read from the file. 
	 * @param filename The name of the .txt file 
	 * @param nvariables The number of variables 
	 * @see ReadDataFromFileTree(char* filename, char* treename, std::vector<char*> branchnames)
	 * @see ReadDataFromFileHist(char* filename, char* histname); 
	 */
	int ReadDataFromFileTxt(const char* filename, int nvariables);

	/** 
	 * Reads data from a user specified file. \n 
	 * Opens a user specified file and creates data objects \n 
	 * containing the values read from the file. This method \n 
	 * needs to be overloaded. 
	 * @param filename The name of the file 
	 * @param options_int A vector of options of type int 
	 * @param options_double A vector of options of type double 
	 * @param options_char A pointer of characters 
	 * @see ReadDataFromFileTree(char* filename, char* treename, std::vector<char*> branchnames)
	 * @see ReadDataFromFileHist(char* filename, char* histname, const char* branchnames); 
	 * @see ReadDataFromFileTxt(char* filename, int nbranches); 
	 */   
	virtual int ReadDataFromFileUser(const char* filename, std::vector<int> options_int, std::vector<double> options_double, const char* options_char); 

	/**
	 * Adds a data point to the data set. 
	 * @param datapoint The data point to be added 
	 */ 
	void AddDataPoint(BCDataPoint* datapoint); 

	/**
	 * Resets the content of the data set 
	 */ 
	void Reset(); 

 private: 

	/* 
	 * A vector containing the data points 
	 */ 
	BCDataVector* fBCDataVector; 

}; 

// --------------------------------------------------------- 

#endif 

