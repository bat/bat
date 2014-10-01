#ifndef __BCMTFCHANNEL__H
#define __BCMTFCHANNEL__H

/*!
 * \class BCMTFChannel
 * \brief A class describing a physics channel.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.1
 * \date 06.2012
 * \detail This class describes a physics channel.
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

class BCMTFTemplate;
class BCMTFSystematicVariation;
class TH2D;

// ---------------------------------------------------------
class BCMTFChannel
{

public:

        /** \name Constructors and destructors */
        /** @{ */

        /**
         * The default constructor.
         * @param name The name of the channel. */
        BCMTFChannel(const char * name);

        /**
         * The default destructor. */
        ~BCMTFChannel();

        /** @} */
        /** \name Member functions (get) */
        /** @{ */

        /**
         * @return The name of the channel. */
        std::string GetName()
        { return fName; };

        /**
         * @return The data. */
        BCMTFTemplate * GetData()
        { return fData; };

        /**
         * Return a template
         * @param index The template index.
         * @return The template. */
        BCMTFTemplate * GetTemplate(int index)
        { return fTemplateContainer.at(index); };

        /**
         * Return a systematic variation
         * @param index The systematic index.
         * @return The systematic variation. */
        BCMTFSystematicVariation * GetSystematicVariation(int index)
        { return fSystematicVariationContainer.at(index); };

        /**
         * @return Flag defining if the channel is active or not. */
        bool GetFlagChannelActive()
        { return fFlagChannelActive; };

        /**
         * Return a histogram ued for the calculation of the error band
         * of the expectation.
         * @return The histogram. */
        TH2D* GetHistUncertaintyBandExpectation()
        { return fHistUncertaintyBandExpectation; };

        /**
         * Return a histogram used for the calculation of the Poisson fluctuations.
         * @return The histogram. */
        TH2D* GetHistUncertaintyBandPoisson()
        { return fHistUncertaintyBandPoisson; };

        /**
         * @return The minimal y-range for printing. */
        double GetRangeYMin()
        { return fRangeYMin; };

        /**
         * @return The maximal y-range for printing. */
        double GetRangeYMax()
        { return fRangeYMax; };

        /** @} */
        /** \name Member functions (set) */
        /** @{ */

        /**
         * Set the name of the channel.
         * @param name The name of the channel. */
        void SetName(const char * name)
        { fName = name; };

        /**
         * Set the data set.
         * @param bctemplate The data set. */
        void SetData(BCMTFTemplate * bctemplate)
        { fData = bctemplate; };

        /**
         * Set a histogram ued for the calculation of the error band of
         * the expectation.
         * @param hist The histogram. */
        void SetHistUncertaintyBandExpectation(TH2D* hist) {
                fHistUncertaintyBandExpectation = hist; }

        /**
         * Set a histogram used for the calculation of the Poisson fluctuations.
         * @param The histogram. */
        void SetHistUncertaintyBandPoisson(TH2D* hist) {
                fHistUncertaintyBandPoisson = hist; }

        /**
         * Set flag to define if the channel is active or not.
         * @param flag The flag. */
        void SetFlagChannelActive(bool flag)
        { fFlagChannelActive = flag; };

        /**
         * Set the y-ranges for printing.
         * @param min The minimum range.
         * @param max The maximum range. */
        void SetRangeY(double min, double max) {
                fRangeYMin = min; fRangeYMax = max; };

        /** @} */

        /** \name Member functions (miscellaneous methods) */
        /** @{ */

        /**
         * Add a template.
         * @param bctemplate The template. */
        void AddTemplate(BCMTFTemplate * bctemplate)
        { fTemplateContainer.push_back(bctemplate); };

        /**
         * Add a systematic variation.
         * @param variation The variation. */
        void AddSystematicVariation(BCMTFSystematicVariation * variation)
        { fSystematicVariationContainer.push_back(variation); };

        /**
         * Calculate histogram for uncertainty band calculation. */
        void CalculateHistUncertaintyBandPoisson();

        /**
         * Calculate histogram for uncertainty band calculation and
         * return a TH1D.
         * @param minimum The minimum value on the expectation.
         * @param maximum The maximum value on the expectation.
         * @param color The color scheme.
         * @return A TH1D histogram. */
        TH1D* CalculateUncertaintyBandPoisson(double minimum, double maximumm, int color);

        /** @} */

        /** \name Member functions (output methods) */
        /** @{ */

        /**
         * Print the templates in this channel.
         * @param filename The name of the file. */
        void PrintTemplates(std::string filename);

        /**
         * Print a particular template with systematics.
         * @param index The template index.
         * @param filename The name of the file. */
        void PrintTemplate(int index, const char * filename);

        /**
         * Print histogram for uncertainty band calculation.
         * @param filename The name of the file. */
        void PrintHistUncertaintyBandExpectation(const char* filename);

        /**
         * Print histogram for uncertainty band calculation.
         * @param filename The name of the file. */
        void PrintHistUncertaintyBandPoisson(const char* filename, const char* options="COLZ");

        /**
         * Print cumulative histogram for uncertainty band calculation.
         * @param filename The name of the file. */
        void PrintHistCumulativeUncertaintyBandPoisson(const char* filename);

        /**
         * Print uncertainty band.
         * @param filename The name of the file. */
        void PrintUncertaintyBandPoisson(const char* filename, double minimum, double maximum, int color);

        /** @} */

private:

        /**
         * The name of the channel. */
        std::string fName;

        /**
         * The data set. */
        BCMTFTemplate * fData;

        /**
         * The minimal y-range for printing. */
        double fRangeYMin;

        /**
         * The maximal y-range for printing. */
        double fRangeYMax;

        /**
         * A container of templates. */
        std::vector<BCMTFTemplate *> fTemplateContainer;

        /**
         * A container of systematics. */
        std::vector<BCMTFSystematicVariation *> fSystematicVariationContainer;

        /**
         * Flag defining if the channel is active or not. */
        bool fFlagChannelActive;

        /**
         * A histogram for the calculation of uncertainty bands. */
        TH2D* fHistUncertaintyBandExpectation;

        /**
         * A histogram for the calculation of uncertainty bands. */
        TH2D* fHistUncertaintyBandPoisson;

};
// ---------------------------------------------------------

#endif

