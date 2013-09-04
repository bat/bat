#ifndef __BCFBUBKGSYSTEMATIC__H
#define __BCFBUBKGSYSTEMATIC__H

#include <string>
#include <map>
#include <BCFBUSystematic.h>
#include <vector>
#include <algorithm>

class TH1D;
class TH2D;
class TF1;
class TH1;

class BCFBUBkgSystematic:public BCFBUSystematic
{
 public:

   /** \name Constructors and destructors */
   /** @{ */

         /**
    * The default constructor. 
    * @param name The name of the background process. */

  BCFBUBkgSystematic(std::string systname);

  BCFBUBkgSystematic(std::string systname, std::string samplename, TH1 *h_up, TH1 *h_down, TH2D *responseup, TH2D *responsedown);
  
  BCFBUBkgSystematic(std::string systname, std::string samplename, TH1 *h_up, TH1 *h_down);
  


  TH2D *GetResponseUp() {return fResponseUp;};
  TH2D *GetResponseDown() {return fResponseDown;};

  void SetResponseUp(TH2D *h_responseUp) { fResponseUp = h_responseUp;};
  void SetResponseDown(TH2D *h_responseDown) { fResponseDown = h_responseDown;};

  void SetHistoUp(std::string backgroundName, TH1 *h_up)
  { fHistogramUp[backgroundName]=h_up;};

  void SetHistoDown(std::string backgroundName, TH1 *h_down)
  { fHistogramDown[backgroundName]=h_down;};


  TH1* GetHistoUp(std::string backgroundName)
    { return fHistogramUp[backgroundName];};

  TH1* GetHistoDown(std::string backgroundName)
    { return fHistogramDown[backgroundName];};

  
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
  ~BCFBUBkgSystematic();
 

 private:
  TH2D *fResponseUp;
  TH2D *fResponseDown;
  std::map<std::string, TH1*> fHistogramUp;
  std::map<std::string, TH1*> fHistogramDown;

  std::vector<std::string> fBackgroundsAffected;
};

#endif
