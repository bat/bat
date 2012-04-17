#ifndef __BCMTFTEMPLATE__H
#define __BCMTFTEMPLATE__H

/*!
 * \class BCMTFTemplate
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

class TH1D;
class TF1;

// ---------------------------------------------------------
class BCMTFTemplate
{
 public:

   // Constructors and destructor
   BCMTFTemplate(const char * channelname, const char * processname);
   ~BCMTFTemplate();

   // setters

   // set efficiency
   void SetEfficiency(double eff)
      { fEfficiency = eff; };

   // set histogram
   void SetHistogram(TH1D * hist);

   // set function container
   void SetFunctionContainer(std::vector<TF1 *> * funccont, int nbins);

   // getters

   // return the name of the channel
   std::string GetChannelName()
      { return fChannelName; };

   // return the name of the process
   std::string GetProcessName()
      { return fProcessName; };

   // return efficiency
   double GetEfficiency()
      { return fEfficiency; };

   // return histogram
   TH1D * GetHistogram()
      { return fHistogram; };

   std::vector<TF1 *> * GetFunctionContainer()
      { return fFunctionContainer; };

   // return the number of bins
   int GetNBins()
      { return fNBins; };

 private:

      // the efficiency of the contribution
      double fEfficiency;

      // the template histogram
      TH1D * fHistogram;

      // a histogram alternative for templates: a vector of TF1 functions
      std::vector<TF1 *> * fFunctionContainer;

      // number of bins in the histogram
      int fNBins;

      // channel name
      std::string fChannelName;

      // process name
      std::string fProcessName;

};
// ---------------------------------------------------------

#endif

