#ifndef __BCMTFSYSTEMATIC__H
#define __BCSYSTEMAITC__H

#include <string>

// ---------------------------------------------------------
class BCMTFSystematic
{
 public:

   // Constructors and destructor
   BCMTFSystematic(const char * name);
   ~BCMTFSystematic();

   // setters

   void SetFlagSystematicActive(bool flag)
      { fFlagSystematicActive = flag; };

   // getters
   std::string GetName()
      { return fSystematicName; };

   // return flag
   bool GetFlagSystematicActive()
      { return fFlagSystematicActive; };

 private:

      // name of the systematic source
      std::string fSystematicName;

      // flag: systematic is used (true) or not (false) in fit
      bool fFlagSystematicActive;
};
// ---------------------------------------------------------

#endif

