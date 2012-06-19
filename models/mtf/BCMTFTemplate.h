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
    * @param hist The TH1D histogram. */
   void SetHistogram(TH1D * hist);

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
    * The name of the channel. */
   std::string fChannelName;

   /** 
    * The name of the process. */
   std::string fProcessName;

};
// ---------------------------------------------------------

#endif

