#ifndef __BCMTFCHANNEL__H
#define __BCMTFCHANNEL__H

/*!
 * \class BCMTFChannel
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
#include <vector>

class BCMTFTemplate;
class BCMTFSystematicVariation;
class TH2D; 

// ---------------------------------------------------------
class BCMTFChannel
{

   public:

      // Constructors and destructor
      BCMTFChannel(const char * name);
      ~BCMTFChannel();

      // setters

      // set name
      void SetName(const char * name)
         { fName = name; };

      // set data
      void SetData(BCMTFTemplate * bctemplate)
         { fData = bctemplate; };

			// set histogram for uncertainty band on expectation
			void SetHistUncertaintyBandExpectation(TH2D* hist) {
				fHistUncertaintyBandExpectation = hist; }
			
			// set histogram for uncertainty band on observed number of events
			void SetHistUncertaintyBandPoisson(TH2D* hist) {
				fHistUncertaintyBandPoisson = hist; }
			
      // set a flag for using this channel or not
      void SetFlagChannelActive(bool flag)
         { fFlagChannelActive = flag; };

			// set y-ranges for printing
			void SetRangeY(double min, double max) {
				fRangeYMin = min; fRangeYMax = max; };

      // getters
      std::string GetName()
         { return fName; };

      // return a template
      BCMTFTemplate * GetTemplate(int index)
         { return fTemplateContainer.at(index); };

      // return a systematicvariation
      BCMTFSystematicVariation * GetSystematicVariation(int index)
         { return fSystematicVariationContainer.at(index); };

      // return the data
      BCMTFTemplate * GetData()
         { return fData; };

      // return flag
      bool GetFlagChannelActive()
         { return fFlagChannelActive; };

			// return error band for expectation
			TH2D* GetHistUncertaintyBandExpectation()
			{ return fHistUncertaintyBandExpectation; }; 

			// return error band for number of events
			TH2D* GetHistUncertaintyBandPoisson()
			{ return fHistUncertaintyBandPoisson; }; 

			// return minimal y-range for printing
			double GetRangeYMin()
			{ return fRangeYMin; }; 

			// return maximal y-range for printing
			double GetRangeYMax()
			{ return fRangeYMax; }; 

      // misc

      // add a template
      void AddTemplate(BCMTFTemplate * bctemplate)
         { fTemplateContainer.push_back(bctemplate); };

      // add a systematic variation
      void AddSystematicVariation(BCMTFSystematicVariation * variation)
         { fSystematicVariationContainer.push_back(variation); };

      // print templates
      void PrintTemplates(const char * filename);

      // print template with systematics
      void PrintTemplate(int index, const char * filename);

			// print histogram for uncertainty band calculation 
			void PrintHistUncertaintyBandExpectation(const char* filename); 

			// print histogram for uncertainty band calculation 
			void PrintHistUncertaintyBandPoisson(const char* filename); 

 private:

      // name of the channel
      std::string fName;

      // the data set
      BCMTFTemplate * fData;

			// minimal y-range for printing
			double fRangeYMin;

			// maximal y-range for printing
			double fRangeYMax;

      // a container of templates
      std::vector<BCMTFTemplate *> fTemplateContainer;

      // a container of systematics
      std::vector<BCMTFSystematicVariation *> fSystematicVariationContainer;

      // flag: channel is used (true) or not (false) in fit
      bool fFlagChannelActive;

			// a histogram for the calculation of uncertainty bands
			TH2D* fHistUncertaintyBandExpectation; 

			// a histogram for the calculation of uncertainty bands
			TH2D* fHistUncertaintyBandPoisson;

};
// ---------------------------------------------------------

#endif

