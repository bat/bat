/*! \class BCDataPoint
 *  \brief A class for handling a data point. 
 *
 * The class defines a data point. A data object can be an event, \n
 * a bin content, etc. Each data point can store several variables \n 
 * of type double. \n 
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  D. Kollar, K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, 
 *          kevin.kroeninger *at* phys *dot* uni *minus* goettingen *dot* de 
 *
 * CREATED: 02.03.2007 
 * 
 * REVISION: 
 *
 * 02.03.2007 Kevin, added comments and header \n
 * 12.06.2007 Kevin, renamed to BCDataPoint \n
 *
 * --------------------------------------------------------- 
 */ 

// --------------------------------------------------------- 

#ifndef __BCDATAPOINT__H
#define __BCDATAPOINT__H

#include <vector.h> 

// --------------------------------------------------------- 

class BCDataPoint
{

 public:
	
	// constructor and destructor 

	/** 
	 * A constructor. 
	 * @param nvariables The number of variables stored in a data object 
	 */ 
	BCDataPoint(int nvariables); 

	/** 
	 * A constructor. 
	 * @param x The vector containing the data
	 */ 
	BCDataPoint(vector<double> x); 

	/** 
	 * A destructor. 
	 */ 
	~BCDataPoint(); 

	// methods (get) 

	/**
	 * Returns the value of a variable. 
	 * @param index The index of the variable 
	 * @return The value of the variable
	 */  
	double GetValue(int index); 

	/** 
	 * Returns the number of values. 
	 */ 
	int GetNValues()
	{ return int(fData.size()); }; 

	// methods (set) 

	/** 
	 * Set the value of a variable. 
	 * @param index The index of the variable
	 * @param value The value of the variable
	 */ 
	void SetValue(int index, double value); 

	/** 
	 * Set the values of all variables. 
	 * @param values A vector of values
	 */ 
	void SetValues(std::vector <double> values); 

 private: 

	/** 
	 * The values of the variables. 
	 */ 
	std::vector <double> fData; 

}; 

// --------------------------------------------------------- 

#endif 

