#ifndef __BCMTFPROCESS__H
#define __BCMTFPROCESS__H

/**
 * @class BCMTFProcess
 * @brief A class describing a process.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.1
 * @date 06.2012
 * @details This class describes a process.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "../../BAT/BCAux.h"

#include <string>

// ---------------------------------------------------------
class BCMTFProcess
{
public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * The default constructor.
     * @param name The name of the process. */
    BCMTFProcess(const std::string& name);

    /**
     * The default destructor. */
    ~BCMTFProcess();

    /** @} */
    /** \name Member functions (get) */
    /** @{ */

    /**
     * @return The name of the process. */
    const std::string& GetName()
    { return fName; };

    /**
     * @return The safe name of the process. */
    const std::string& GetSafeName()
    { return fSafeName; };

    /**
     * @return The histogram color. */
    int GetHistogramColor()
    { return fHistogramColor; };

    /**
     * @return The histogram fill style. */
    int GetHistogramFillStyle()
    { return fHistogramFillStyle; };

    /**
     * @return the Histogram line style. */
    int GetHistogramLineStyle()
    { return fHistogramLineStyle; };

    /** @} */

    /** \name Member functions (set) */
    /** @{ */

    /**
     * Set the name of the process.
     * @param name The name of the process. */
    void SetName(const std::string& name)
    { fName = name; fSafeName = BCAux::SafeName(name); };

    /**
     * Set the histogram color.
     * @param color The color. */
    void SetHistogramColor(int color)
    { fHistogramColor = color; };

    /**
     * Set the histogram fill style.
     * @param color The fill style. */
    void SetHistogramFillStyle(int style)
    { fHistogramFillStyle = style; };

    /**
     * Set the histogram line style.
     * @param color The line style. */
    void SetHistogramLineStyle(int style)
    { fHistogramLineStyle = style; };

    /** @} */

private:

    /**
     * The name of the process. */
    std::string fName;

    /**
     * The safe name of the process. */
    std::string fSafeName;

    /**
     * The histogram color */
    int fHistogramColor;

    /**
     * The fill style. */
    int fHistogramFillStyle;

    /**
     * The line style. */
    int fHistogramLineStyle;

};
// ---------------------------------------------------------

#endif
