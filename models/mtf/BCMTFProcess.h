#ifndef __BCMTFPROCESS__H
#define __BCMTFPROCESS__H

#include <string>

// ---------------------------------------------------------
class BCMTFProcess
{
   public:

      // Constructors and destructor
      BCMTFProcess(const char * name);
      ~BCMTFProcess();

      // setters

      // set name
      void SetName(const char * name)
         { fName = name; };

      // getters
      std::string GetName()
         { return fName; };

 private:

      // name of the channel
      std::string fName;

};
// ---------------------------------------------------------

#endif

