#ifndef __BCMTFTEMPLATE__H
#define __BCMTFTEMPLATE__H

/*!
 * \class BCMTFTemplate
 * \brief A class describing a template
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.1
 * \date 06.2012
 * \detail This class describes a template.
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

#include <TH1D.h>
#include <TRandom3.h>

class TF1;

// ---------------------------------------------------------
class BCMTFTemplate
{
public:

        /** \name Constructors and destructors */
        /** @{ */

        /**
         * The default constructor.
         * @param channelname The name of the channel.
         * @param process name The name of the process. */
        BCMTFTemplate(const char * channelname, const char * processname);

        /**
         * The default destructor. */
        ~BCMTFTemplate();

        /** @} */

        /** \name Member functions (get) */
        /** @{ */

        /**
         * @return The name of the channel. */
        std::string GetChannelName()
        { return fChannelName; };

        /**
         * @return The name of the process. */
        std::string GetProcessName()
        { return fProcessName; };

        /**
         * @return The efficiency. */
        double GetEfficiency()
        { return fEfficiency; };

        /**
         * @return The TH1D histogram. */
        TH1D * GetHistogram()
        { return fHistogram; };

        /**
         * @return The normalization. */
        double GetNorm()
        { return fNormalization; };

        /**
         * @return The original normalization. */
        double GetOriginalNorm()
        { return fOriginalNormalization; };

        /**
         * Fluctuate the original template histogram by the uncertainty on the bin content.
         * @param options A set of options.
         * "P" : use Poisson model, the expectation value parameter is the bin content and also defines the uncertainties.
         * "G" [default] : use a Gaussian mode, the expectation value paramer mu is the bin content, the uncertainty sigma is defined by the uncertainty in the histogram.
         * "Z" [default] : make sure that the bin content is positive.
         * @param norm The target normalization.
         * @return A histogram with each bin fluctuated by the uncertainty on the bin content. \n
         */
        TH1D FluctuateHistogram(std::string options = "GZ", double norm = 1);

        /**
         * @return The function container. */
        std::vector<TF1 *> * GetFunctionContainer()
        { return fFunctionContainer; };

        /**
         * @return The number of bins. */
        int GetNBins()
        { return fNBins; };

        /** @} */

        /** \name Member functions (set) */
        /** @{ */

        /**
         * Set the efficiency.
         * @param eff The efficiency. */
        void SetEfficiency(double eff)
        { fEfficiency = eff; };

        /**
         * Set the histogram.
         * @param hist The TH1D histogram.
         * @param norm The target normalization. */

        void SetHistogram(TH1D * hist, double norm = 1);

        /**
         * Set the original normalization.
         * @param norm The normalization. */
        void SetOrignialNormalization(double norm) {
                fOriginalNormalization = norm; };

        /** Set a function container
         * funccont The function container
         * nbins The number of bins (and functions) */
        void SetFunctionContainer(std::vector<TF1 *> * funccont, int nbins);

        /** @} */

private:

        /**
         * The efficiency of the contribution. */
        double fEfficiency;

        /**
         * The TH1D histogram. */
        TH1D * fHistogram;

        /**
         * A histogram alternative for templates: a vector of TF1 functions. */
        std::vector<TF1 *> * fFunctionContainer;

        /**
         * The number of bins in the histogram. */
        int fNBins;

        /**
         * The normalization. */
        double fNormalization;

        /**
         * The original normalization. */
        double fOriginalNormalization;

        /**
         * The name of the channel. */
        std::string fChannelName;

        /**
         * The name of the process. */
        std::string fProcessName;

        /**
         * A random number generator. */
        TRandom3* fRandom;


};
// ---------------------------------------------------------

#endif

