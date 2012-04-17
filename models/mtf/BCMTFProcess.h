#ifndef __BCMTFPROCESS__H
#define __BCMTFPROCESS__H

/*!
 * \class BCMTFProcess
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

