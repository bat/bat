#ifndef __BCMTFSYSTEMATICVARIATION__H
#define __BCMTFSYSTEMATICVARIATION__H

/*!
 * \class BCMTFSystematicVariation
 * \brief A class describing a systematic variation.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.1
 * \date 06.2012
 * \detail This class describes the impact of a systematic
 * uncertainty.
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

class TH1D;

// ---------------------------------------------------------
class BCMTFSystematicVariation
{
public:

        /** \name Constructors and destructors */
        /** @{ */

        /**
         * The default constructor.
         * @param channelname The name of the channel.
         * @param systematicname The name of the systematic.
         * @param nprocesses The number of processes. */
        BCMTFSystematicVariation(const char * channelname, const char * systematicname, int nprocesses);

        /**
         * The default destructor. */
        ~BCMTFSystematicVariation();

        /** @} */
        /** \name Member functions (get) */
        /** @{ */

        /**
         * Returns the histogram correponding to the up-scale variation of
         * the systematic.
         * @param index The process index.
         * @return The histogram. */
        TH1D * GetHistogramUp(int index)
        { return fHistogramUpContainer.at(index); };

        /**
         * Returns the histogram correponding to the down-scale variation
         * of the systematic.
         * @param index The process index.
         * @return The histogram. */
        TH1D * GetHistogramDown(int index)
        { return fHistogramDownContainer.at(index); };

        /** @} */
        /** \name Member functions (set) */
        /** @{ */

        /**
         * Set the histogram correponding to the up-scale variation of the
         * systematic.
         * @param index The process index.
         * @param hist The histogram.
         * @see SetHistogramDown(int index, TH1D * hist)
         * @see SetHistograms(int index, TH1D * hist_up, TH1D * hist_down)*/
        void SetHistogramUp(int index, TH1D * hist)
        { fHistogramUpContainer[index] = hist; };

        /**
         * Set the histogram correponding to the down-scale variation of
         * the systematic.
         * @param index The process index.
         * @param hist The histogram.
         * @see SetHistogramUp(int index, TH1D * hist)
         * @see SetHistograms(int index, TH1D * hist_up, TH1D * hist_down)*/
        void SetHistogramDown(int index, TH1D * hist)
        { fHistogramDownContainer[index] = hist; };

        /**
         * Set the histograms correponding to the up- and down-scale
         * variations of the systematic.
         * @param index The process index.
         * @param hist_up The up-scale histogram.
         * @param hist_down The down-scale histogram.
         * @see SetHistogramUp(int index, TH1D * hist)
         * @see SetHistogramDown(int index, TH1D * hist) */
        void SetHistograms(int index, TH1D * hist_up, TH1D * hist_down)
        { fHistogramUpContainer[index] = hist_up;
                fHistogramDownContainer[index] = hist_down; };

        /** @} */
        /** \name Member functions (miscellaneous methods) */
        /** @{ */

        /**
         * Add a histogram for up-scale variations.
         * @param hist The histogram.
         * @see AddHistogramDown(TH1D * hist)
         * @see AddHistograms(TH1D * hist_up, TH1D * hist_down)*/
        void AddHistogramUp(TH1D * hist)
        { fHistogramUpContainer.push_back(hist); };

        /**
         * Add a histogram for down-scale variations.
         * @param hist The histogram.
         * @see AddHistogramUp(TH1D * hist)
         * @see AddHistograms(TH1D * hist_up, TH1D * hist_down)*/
        void AddHistogramDown(TH1D * hist)
        { fHistogramDownContainer.push_back(hist); };

        /**
         * Add a histograms for up- and down-scale variations.
         * @param hist_up The up-scale histogram.
         * @param hist_down The down-scale histogram.
         * @see AddHistogramUp(TH1D * hist)
         * @see AddHistogramDown(TH1D * hist) */
        void AddHistograms(TH1D * hist_up, TH1D * hist_down)
        { fHistogramUpContainer.push_back(hist_up);
                fHistogramDownContainer.push_back(hist_down); };

        /** @} */

private:

        /**
         * A container of histograms. */
        std::vector<TH1D *> fHistogramUpContainer;

        /**
         * A container of histograms. */
        std::vector<TH1D *> fHistogramDownContainer;

        /**
         * The name of the corresponding channel. */
        std::string fChannelName;

        /**
         * The name of the corresponding source of systematic
         * uncertainty. */
        std::string fSystematicName;

};
// ---------------------------------------------------------

#endif

