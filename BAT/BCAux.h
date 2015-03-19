#ifndef __BCAUX__H
#define __BCAUX__H

/*!
 * \namespace BCAux
 * \brief Some functions not fitting anywhere else
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 01.2009
 * \detail A namespace which encapsulates auxiliary functions
 * necessary for BAT.
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

class TH2;

// ---------------------------------------------------------

namespace BCAux
{

/**
 * Sets the default BAT style for drawing plots. */
void SetStyle();

/**
 * Force file extension to be .pdf if not already .pdf or .ps
 * @param filename Filename to be altered */
void DefaultToPDF(std::string& filename);

/**
 * Transpose a TH2.
 * @param h 2D histogram to transpose
 * @param name Name to give to transposed histogram, if empty, "_tr" is added to original hist's name.
 * @return Transposed independent histogram. */
TH2* Transpose(TH2 const* const h, std::string name = "");

/** Range types. */
enum BCRange {
    kFiniteRange           = 0, //!< lower and upper limits finite
    kNegativeInfiniteRange = 1, //!< lower limit infinite, upper limit finite
    kPositiveInfiniteRange = 2, //!< lower limit finite, upper limit infinite
    kInfiniteRange         = 3, //!< lower and upper limits infinite
    kEmptyRange            = 4  //!< lower limit == upper limit
};

/**
 * Return type of range as a BCAux::BCRange enum.
 * @param xmin lower limit of range
 * @param xmax upper limit of range
 * @return range type. */
BCAux::BCRange RangeType(double xmin, double xmax);

/**
 * Make an infinite range finite by setting inf values to max.
 * @param xmin lower limit to coerce
 * @param xmax upper limit to coerce */
void MakeFinite(double& xmin, double& xmax);

}

// ---------------------------------------------------------

#endif
