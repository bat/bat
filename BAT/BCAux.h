#ifndef __BCAUX__H
#define __BCAUX__H

/**
 * @namespace BCAux
 * @brief Some functions not fitting anywhere else
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @version 1.0
 * @date 01.2009
 * @details A namespace which encapsulates auxiliary functions
 * necessary for BAT.
 */

/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <string>

#include "BCH1D.h"
#include "BCH2D.h"

#include <TROOT.h>

class TH2;

// ---------------------------------------------------------

namespace BCAux
{

/**
 * A trash to keep heap-allocated objects of type T alive until the
 * trash goes out of scope. Then they are deleted. Ownership is not
 * transferred in copies, only during a swap. */
template <typename T>
class BCTrash
{
public:
    //! @nowarn
    BCTrash() {}
    BCTrash(const BCTrash<T>&) {}
    BCTrash<T>& operator=(const BCTrash<T>&)
    {
        return *this;
    }

    ~BCTrash()
    {
        for (unsigned i = 0; i < fStorage.size(); ++i)
            delete fStorage[i];
    }

    void Put(T* object)
    {
        fStorage.push_back(object);
    }
    //! @endnowarn

private:
    std::vector<T*> fStorage;
};

/**
 * @deprecated This function is empty but for a warning message. Please do not call it. */
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
TH2* Transpose(const TH2* const h, const std::string& name = "");

/** Range types. */
enum BCRange {
    kFiniteRange           = 1, //!< lower < upper, lower and upper limits finite
    kNegativeInfiniteRange = 2, //!< lower < upper, lower limit infinite, upper limit finite
    kPositiveInfiniteRange = 3, //!< lower < upper, lower limit finite, upper limit infinite
    kInfiniteRange         = 4, //!< lower < upper, lower and upper limits infinite
    kEmptyRange            = 5, //!< lower limit == upper limit
    kReverseRange          = 6  //!< lower > upper
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

/**
 * Convert a name into a safe name for use in ROOT object naming. */
std::string SafeName(const std::string& name);

/**
 * @return Whether character is allowed in a safe name. */
bool AllowedCharacter(char c);

/** An enumerator for the knowledge update drawing style presets. */
enum BCKnowledgeUpdateDrawingStyle {
    kKnowledgeUpdateDefaultStyle      = 0, ///< Simple line-drawn histograms
    kKnowledgeUpdateDetailedPosterior = 1, ///< Posterior drawn with detailed info, prior drawn as overlaid line
    kKnowledgeUpdateDetailedPrior     = 2  ///< Prior drawn with detailed info, posterior drawn as overlaid line
};

/**
 * Use pre-made drawing options for knowledge update plots.
 * @param prior Prior histogram container
 * @param posterior Posterior histogram container
 * @param style Style option */
void SetKnowledgeUpdateDrawingStyle(BCH1D& prior, BCH1D& posterior, BCAux::BCKnowledgeUpdateDrawingStyle style = BCAux::kKnowledgeUpdateDefaultStyle);

/**
 * Use pre-made drawing options for knowledge update plots.
 * @param prior Prior histogram container
 * @param posterior Posterior histogram container
 * @param style Style option*/
void SetKnowledgeUpdateDrawingStyle(BCH2D& prior, BCH2D& posterior, BCAux::BCKnowledgeUpdateDrawingStyle style = BCAux::kKnowledgeUpdateDefaultStyle);

/**
 * Draw knowledge update plot into current TPad
 * @param prior BCHistogramBase containing prior
 * @param posterior BCHistogramBase containing posterior
 * @param draw_prior_first Flag for deciding drawing order.
 * @param trash Where to keep references of plot objects for later deletion. */
void DrawKnowledgeUpdate(BCHistogramBase& prior, BCHistogramBase& posterior, bool draw_prior_first, BCTrash<TObject>& trash);

/**
 * Print plots
 * @param h1 Vector of 1D histograms to plot
 * @param h2 Vector of 2D histograms to plot
 * @param filename Path to file to print to
 * @param hdiv Number of columns of plots per page
 * @param vdiv Number of rows of plots per page
 * @return Number of plots printed */
unsigned PrintPlots(std::vector<BCH1D>& h1, std::vector<BCH2D>& h2, const std::string& filename, unsigned hdiv = 1, unsigned vdiv = 1);

/**
 * A guard object to prevent ROOT from taking over ownership of `TNamed` objects.
 *
 * The guard should be created in the scope in which `TNamed` objects are created, for example:
 *
 * ```{.cxx}
 * TH1D* hist;
 * {
 *     RootSideEffectGuard g;
 *     hist = new TH1D(...);
 * }
 * ```
 *
 * */
class RootSideEffectGuard
{
public:
    RootSideEffectGuard();
    ~RootSideEffectGuard();

private:
    // not implemented
    RootSideEffectGuard(const RootSideEffectGuard&);
    RootSideEffectGuard& operator=(const RootSideEffectGuard&);

    TDirectory* fDirectory;
};

/**
 * Create a clone of the input but avoid registering the object with ROOT so it
 * cannot be deleted twice */
template <class T>
T* OwnClone(const T* o)
{
    return static_cast<T*>(gROOT->CloneObject(o, false));
}

/**
 * Create a clone of the input, change the name but avoid registering the object
 * with ROOT so it cannot be deleted twice */
template <class T>
T* OwnClone(const T* o, const std::string& name)
{
    T* res = OwnClone(o);
    res->SetName(name.c_str());
    return res;
}

}

// ---------------------------------------------------------

#endif
