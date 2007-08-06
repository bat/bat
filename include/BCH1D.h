/**
 * A class which contains a TH1D histogram and can be used for marginalized probabilties 
 *
 * --------------------------------------------------------- 
 *
 * AUTHOR:  D. Kollar, K. Kroeninger 
 *
 * CONTACT: dkollar *at* mppmu *dot* mppmu *dot* de, kroening *at* mppmu *dot* mppmu *dot* de 
 *
 * CREATED: 12.06.2007 
 * 
 * REVISION:
 *
 * 03.08.2007  Dano  * increase printout level of ROOT routines like Print()
 * 06.08.2007  Dano  * change Print "if" structure to "case" structure
 *
 *
 * --------------------------------------------------------- 
 *
 *
*/ 

// --------------------------------------------------------- 

#ifndef __BCH1D__H
#define __BCH1D__H

#include <vector.h> 

#include <TH1D.h> 

#include <TError.h>

#include "BCLog.h"

// --------------------------------------------------------- 

class BCH1D
{
  
 public:
  
  // constructors and destructor 
  
  /**
   * The default constructor. 
   */ 
  BCH1D(); 

  /** 
   * The default destructor. 
   */ 
  ~BCH1D(); 
  
  // methods (get) 

  /**
   * @return The histogram of the marginalized probability 
   */ 
  TH1D* GetHistogram()
    { return fHistogram; }; 

  /** 
   * @return The mean of the distribution 
   */ 
  double GetMean()
    { return fHistogram -> GetMean(); }; 

  /** 
   * @return The mode of the distribution 
   */ 
  double GetMode(); 

  /** 
   * @return The median of the distribution 
   */ 
  double GetMedian()
    { return this -> GetQuantile(0.5); }; 

  /** 
   * @param probabilitysum The probability sum 
   * @return The quantile of the distribution for the sum  
   */ 
  double GetQuantile(double probabilitysum); 

  /** 
   * @param valuemin The value from which the intergration is done 
   * @param valuemax The value up to which the intergration is done 
   * @return The integral  
   */ 
  double GetIntegral(double valuemin, double valuemax); 

  /** 
   * Returns the value below which 90% of the marginalized probability lie. 
   * @return The 90% limit 
   */ 
  double GetLimit90()
    { return this -> GetQuantile(0.9); }; 

  /** 
   * Returns the value below which 95% of the marginalized probability lie. 
   * @return The 95% limit 
   */ 
  double GetLimit95()
    { return this -> GetQuantile(0.95); }; 

  /** 
   * Returns the value below which 99% of the marginalized probability lie. 
   * @return The 99% limit 
   */ 
  double GetLimit99()
    { return this -> GetQuantile(0.99); }; 

  /** 
   * Returns the positive uncertainty from the mode to the 84% value. 
   * Value at which the probability integral reaches 84% minus the mode. 
   * @return The positive uncertainty. 
   */ 
  double GetUncertaintyPlusFromMode()
    { return this -> GetQuantile(0.84) - this -> GetMode(); }; 

  /** 
   * Returns the negative uncertainty from the mode to the 16% probability 
   * The mode - the value at which the probability integral reaches 16%. 
   * @return The negative uncertainty. 
   */ 
  double GetUncertaintyMinusFromMode()
    { return this -> GetMode() - this -> GetQuantile(0.16); }; 

  /** 
   * Returns the positive uncertainty from the mean to the 84% value. 
   * Value at which the probability integral reaches 84% minus the mean. 
   * @return The positive uncertainty. 
   */ 
  double GetUncertaintyPlusFromMean()
    { return this -> GetQuantile(0.84) - this -> GetMean(); }; 

  /** 
   * Returns the negative uncertainty from the mean to the 16% probability 
   * The mean - the value at which the probability integral reaches 16%. 
   * @return The negative uncertainty. 
   */ 
  double GetUncertaintyMinusFromMean()
    { return this -> GetMean() - this -> GetQuantile(0.16); }; 

  /** 
   * Returns the p-value. 
   * Returns the intergral from probability to the maximum. 
   * @param probability Lower limit of integration 
   * @return The p-value 
   */ 
  double GetPValue(double probability);     

  // methods (set) 

  /**
   * @param hist The histogram 
   */ 
  void SetHistogram(TH1D* hist) 
    { fHistogram = hist; }; 

  // methods 

  /** 
   * @param filename The filename 
   * @param options 0 = no lines, 1 = lines at mode, 16% and 84% probability, 2 = line at value 
   * @param value Line position for option = 2 
   */ 
  void Print(char* filename, int options, double value); 

  void Print(char* filename, int options); 

  void Print(char* filename); 
    
 private: 

  /** 
   * The histogram
   */ 
  TH1D * fHistogram; 

};

// --------------------------------------------------------- 

#endif 
