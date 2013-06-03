#ifndef __BCMVCOBSERVABLE__H
#define __BCMVCOBSERVABLE__H

#include <string>

// ---------------------------------------------------------
class BCMVCObservable  
{
 public:

  // Constructor
  BCMVCObservable();

  // Destructor
  ~BCMVCObservable();

  // setters

  // set the name
  void SetName(std::string name)
  { fName = name; };

  // set the minimum and maximum
  void SetMinMax(double min, double max)
  { fMin = min; fMax = max; };

  // getters

  // return the name
  std::string GetName() 
    { return fName; };

  // return the minimum value
  double GetMinimum()
  { return fMin; };

  // return the maximum value
  double GetMaximum()
  { return fMax; };

 private:

  // name of the observable
  std::string fName;

  // the minimum value
  double fMin;

  // the maximum value
  double fMax;

};
// ---------------------------------------------------------

#endif

