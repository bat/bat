#ifndef __BCMTFSYSTEMATIC__H
#define __BCMTFSYSTEMATIC__H

/*!
 * \class BCMTFSystematic
 * \brief A class desribing a systematic uncertainty.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.1
 * \date 06.2012
 * \detail This class describes a systematic uncertainty.
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <string>

// ---------------------------------------------------------
class BCMTFSystematic
{
public:

        /** \name Constructors and destructors */
        /** @{ */

        /**
         * The default constructor.
         * @param name The name of the source of systematic uncertainty. */
        BCMTFSystematic(const char * name);

        /**
         * The default destructor. */
        ~BCMTFSystematic();

        /** @} */
        /** \name Member functions (get) */
        /** @{ */

        /**
         * @return The name of the systematic uncertainty. */
        std::string GetName()
        { return fSystematicName; };

        /**
         * @return A flag defining if this uncertainty is active or not. */
        bool GetFlagSystematicActive()
        { return fFlagSystematicActive; };

        /** @} */
        /** \name Member functions (get) */
        /** @{ */

        /**
         * Set a flag defining if this uncertainty is active or not.
         * @param flag The flag. */
        void SetFlagSystematicActive(bool flag)
        { fFlagSystematicActive = flag; };

        /** @} */

private:

        /**
         * The name of the source of the systematic uncertainty. */
        std::string fSystematicName;

        /**
         * A flag defining if this uncertainty is active or not. */
        bool fFlagSystematicActive;

};
// ---------------------------------------------------------

#endif

