#ifndef __BCMTFPROCESS__H
#define __BCMTFPROCESS__H

/*!
 * \class BCMTFProcess
 * \brief A class describing a process.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.1
 * \date 06.2012
 * \detail This class describes a process.
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

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. 
       * name The name of the process. */
      BCMTFProcess(const char * name);

      /**
       * The default destructor. */
      ~BCMTFProcess();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

			/**
       * @return The name of the process. */
      std::string GetName()
         { return fName; };

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

			/** 
       * Set the name of the process.
			 * @param name The name of the process. */
      void SetName(const char * name)
         { fName = name; };

      /** @} */

 private:

			/**
			 * The name of the process. */
      std::string fName;

};
// ---------------------------------------------------------

#endif

