#ifndef __BCFBUNORMSYSTEMATIC__H
#define __BCFBUNORMSYSTEMATIC__H

#include <vector>
#include <string>
#include <BCFBUSystematic.h>
#include <algorithm>

class TH1D;
class TF1;
class TH1;

class BCFBUNormSystematic:public BCFBUSystematic
{
 public:

   /** \name Constructors and destructors */
   /** @{ */

         /**
    * The default constructor. 
    * @param name The name of the background process. */

  BCFBUNormSystematic(std::string name, double uncertainty);

  void SetUncertainty(double uncertainty)
  { fNormUncertainty = uncertainty;};


  double GetUncertainty()
  { return fNormUncertainty;};
 
  void AddBackground(std::string sampleName)
  { fBackgroundsAffected.push_back(sampleName); };
 
  bool ApplySystematic(std::string sampleName)
  {
    bool res = false;
    
    if (std::find(fBackgroundsAffected.begin(), fBackgroundsAffected.end(), sampleName) != fBackgroundsAffected.end())
      res = true;

    return res;
  }

  /**
   * The default destructor. */
  ~BCFBUNormSystematic();
 

 private:
  double fNormUncertainty;
  
  std::vector<std::string> fBackgroundsAffected;

};

#endif
