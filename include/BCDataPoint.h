/**
 * A class which defines a data object. 
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  D. Kollar, K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, kroening *at* mppmu *dot* mppmu *dot* de 
 *
 * CREATED: 02.03.2007 
 * 
 * REVISION: 
 *
 * 02.03.2007 Kevin, added comments and header 
 * 12.06.2007 Kevin, renamed to BCDataPoint
 *
 * --------------------------------------------------------- 
 *
 *
 * The class defines a data object. A data object can be an event, a bin content, etc. 
 * Each data object can store several variables. 
 *
*/ 

// --------------------------------------------------------- 

#ifndef __BCDATAPOINT__H
#define __BCDATAPOINT__H

#include <TROOT.h>

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
   * A destructor 
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
   * returns number of values. 
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

  /** Set the values of all variables 
   * @param values A vector of values
   */ 
  void SetValues(std::vector <double> values); 

 private: 

  /** 
   * The values of the variables
   */ 
  std::vector <double> fData; 

}; 

// --------------------------------------------------------- 

#endif 
