#ifndef __BCMTFCOMPARISONTOOL__H
#define __BCMTFCOMPARISONTOOL__H

/*!
 * \class BCMTFComparisonTool
 * \brief A helper class for BCMTFAnalysisFacility storing information.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.1
 * \date 06.2012
 * \detail This is a helper class for BCMTFAnalysisFacility storing information.
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
#include <vector>

#include <TH1D.h>

// ---------------------------------------------------------
class BCMTFComparisonTool
{

public:

        /** \name Constructors and destructors */
        /** @{ */

        /**
         * The default constructor.
         * @param name The name of the class. */
        BCMTFComparisonTool(const char * name);

        /**
         * The defaul destructor. */
        ~BCMTFComparisonTool();

        /** @} */
        /** \name Member functions (get) */
        /** @{ */

        /**
         * @return The name of the class. */
        std::string GetName()
        { return fName; };

        /**
         * @return The number of contributions. */
        int GetNContributions()
        { return (int) fHistogramContainer.size(); };

        /** @} */
        /** \name Member functions (miscellaneous methods) */
        /** @{ */

        /**
         * Add a constribution.
         * @param name The name of the contribution.
         * @param hist The histogram. */
        void AddContribution(const char * name, TH1D hist);

        /**
         * Add a constribution.
         * @param name The name of the contribution.
         * @param centralvalue The central value.
         * @param uncertainty The uncertainty. */
        void AddContribution(const char * name, double centralvalue, double uncertainty);

        /**
         * Draw an overview. */
        void DrawOverview();

        /** @} */
        /** \name Member functions (output methods) */
        /** @{ */

        /**
         * Print all histograms to a file.
         * @param filename The name of the file. */
        void PrintHistograms(const char * filename);

        /**
         * Print an overview to a file.
         * @param filename The name of the file. */
        void PrintOverview(const char * filename);

        /** @} */

private:

        /**
         * The name of the class. */
        std::string fName;

        /**
         * The names of the contributions. */
        std::vector<std::string> fNameContainer;

        /**
         * A container of TH1D histograms. */
        std::vector<TH1D *> fHistogramContainer;

        /**
         * A container of central values. */
        std::vector<double> fCentralValueContainer;

        /**
         * A container of uncertainties. */
        std::vector<double> fUncertaintyContainer;

};
// ---------------------------------------------------------

#endif

