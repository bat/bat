/**
 * A class which defines a marginalized probability. 
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
 * 22.05.2007 Kevin, added nicer 2-d plots including contours 
 *
 * --------------------------------------------------------- 
 *
 *
 * The class defines a marginalized probability. 
 *
*/ 

// --------------------------------------------------------- 

#ifndef __BCH2D__H
#define __BCH2D__H

#include <vector.h> 

#include <TH1D.h> 
#include <TH2D.h> 

// --------------------------------------------------------- 

class BCH2D
{
  
 public:
  
  // constructors and destructor 
  
  /**
   * The default constructor. 
   */ 
  BCH2D(); 

  /** 
   * The default destructor. 
   */ 
  ~BCH2D(); 
  
  // methods (get) 

  /**
   * @return The histogram
   */ 
  TH2D * GetHistogram()
    { return fHistogram; }; 

  /** 
   * @param mean The mean of the distribution 
   */ 
  void GetMean(double& mean);

  /** 
   * @param The mode of the distribution
   */ 
  void GetMode(double& mode);

  // methods (set) 

  /**
   * @param hist The histogram
   */ 
  void SetHistogram(TH2D * hist)
    { fHistogram = hist; };

  // methods 

  /** 
   * Print 2-d histogram to file 
   * @param filename The filename 
   */ 
  void Print(char* filename, int options); 

  void Print(char* filename); 

  void CalculateIntegratedHistogram(); 

  double GetLevel(double p); 
    
 private: 

  /** 
   * The 2-d histogram
   */ 
  TH2D * fHistogram; 

  /**
   * The integrated 2-d histogram 
   */ 
  TH1D * fIntegratedHistogram; 

}; 

// --------------------------------------------------------- 

#endif 
