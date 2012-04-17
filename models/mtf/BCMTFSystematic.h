#ifndef __BCMTFSYSTEMATIC__H
#define __BCSYSTEMAITC__H

/*!
 * \class BCMTFSystematic
 * \brief A class for ...
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 04.2012
 * \detail
 *
 *
 *
 *
 */

/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

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

