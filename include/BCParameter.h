/*! \class BCParameter
 *  \brief A model parameter class
 *
 * This class defines a parameter for a model as well as a container for
 * a parameter set.
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
 * 02.03.2007 Kevin, added comments and header \n
 * 16.05.2007 Dano, added nuisence flag\n
 * 12.06.2007 Kevin, renamed to BCParameter \n
 *
 * --------------------------------------------------------- 
 *
*/ 

// --------------------------------------------------------- 

#ifndef __BCPARAMETER__H
#define __BCPARAMETER__H

#include <TROOT.h>

#include <vector.h>

// --------------------------------------------------------- 

class BCParameter
{

 public:

  // constructors and destructor 

  /** 
   * The default constructor. 
   */ 
  BCParameter(); 

  /** 
   * A constructor. 
   * @param name The name of the parameter 
   * @param lowerlimit The lower limit of the parameter values
   * @param upperlimit The upper limit of the parameter values 
   */ 
  BCParameter(const char* name, double lowerlimit, double upperlimit); 

  /**
   * The default destructor 
   */ 
  ~BCParameter(); 

  // methods (get) 

  /**
   * Returns the name of the parameter. 
   * @return The name of the parameter
   */ 
  char* GetName()
    { return fName; }; 
  
  /** 
   * Returns the index of the parameter within the parameter container of a BCModel. 
   * @return The index of the parameter
   */ 
  int GetIndex()
    { return fIndex; }; 

  /** 
   * Returns the lower limit of the parameter values. 
   * @return The lower limit of the parameter values 
   */ 
  double GetLowerLimit()
    { return fLowerLimit; }; 

  /** 
   * Returns the upper limit of the parameter values. 
   * @return The upper limit of the parameter values 
   */ 
  double GetUpperLimit()
    { return fUpperLimit; }; 

  /** 
   * Returns 1 if parameter is a nuisence parameter or 0 if not.
   * @return 1 - is nuisence paramete, 0 - is not nuisence parameter
   */
  double IsNuisence()
    { return fNuisence; }; 

  // methods (set) 

  /** 
   * Set the name of the parameter. 
   * @param name The name of the parameter 
   */ 
  void SetName(char* name)
    { fName = name; }; 

  /** Set the index of the parameter within the parameter container of a BCModel. 
   * @param index The index of the parameter
   */ 
  void SetIndex(int index) 
    { fIndex = index; }; 

  /** 
   * Set the lower limit of the parameter values. 
   * @param limit The lower limit of the parameter values 
   */ 
  void SetLowerLimit(double limit) 
    { fLowerLimit = limit; }; 

  /** 
   * Set the upper limit of the parameter values. 
   * @param limit The upper limit of the parameter values 
   */ 
  void SetUpperLimit(double limit) 
    { fUpperLimit = limit; }; 

  /**
   * Set parameter to be nuisence.
   * @param nuisence 1 - nuisence, 0 - not nuisence
   */
  void SetNuisence(int nuisence=1)
    { fNuisence = nuisence; };

  // methods 

  /** 
   * Prints a summary on the screen. 
   */ 
  void PrintSummary(); 

 private: 
  
  /**
   * The name of the parameter. 
   */ 
  int fNameSize; 
  char* fName; //[fNameSize]

  /** 
   * The index of the parameter within the parameter container of a BCModel. 
   */ 
  int fIndex; 

  /** 
   * The lower limit of the parameter values. 
   */ 
  double fLowerLimit; 

  /** 
   * The lower limit of the parameter values. 
   */ 
  double fUpperLimit; 

  /**
   * Flag to specify whether to integarte over this parameter
   */
  int fNuisence;

}; 

// --------------------------------------------------------- 

typedef std::vector<BCParameter*> BCParameterSet; 

// --------------------------------------------------------- 

#endif 
