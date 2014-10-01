/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <TH1D.h>
#include <TRandom.h>

#include <iostream>

#include "BCMTFTemplate.h"

// ---------------------------------------------------------
BCMTFTemplate::BCMTFTemplate(const char * channelname, const char * processname)
   : fEfficiency(0)
   , fHistogram(0)
   , fNBins(0)
   , fNormalization(0)
   , fOriginalNormalization(0)
{
   fChannelName = channelname;
   fProcessName = processname;
   fFunctionContainer = new std::vector<TF1 *>(0);
   fRandom = new TRandom3(0);
}

// ---------------------------------------------------------
BCMTFTemplate::~BCMTFTemplate()
{
}

// ---------------------------------------------------------
void BCMTFTemplate::SetHistogram(TH1D * hist, double norm)
{
   // set histogram
   fHistogram = hist;

   // check if histogram exists
   if (!hist)
      return;

   // get number of bins
   fNBins = fHistogram->GetNbinsX();

   // set original normalization
   double orignorm = fHistogram->Integral();
   SetOrignialNormalization(orignorm);

   // normalize histogram
   if (orignorm && norm)
      fHistogram->Scale(norm / orignorm);

   // set normalization
   if (norm)
      fNormalization = norm;
}

// ---------------------------------------------------------
void BCMTFTemplate::SetFunctionContainer(std::vector<TF1 *> * funccont, int nbins)
{
   fFunctionContainer = funccont;
   fNBins = nbins;
}

// ---------------------------------------------------------
TH1D BCMTFTemplate::FluctuateHistogram(std::string options, double norm)
{
   // option flags
   bool flag_p = false;
   bool flag_g = false;
   bool flag_z = false;

   // check content of options string
   if (options.find("P") < options.size()) {
      flag_p = true;
   }

   if (options.find("G") < options.size()) {
      flag_g = true;
   }

   if (options.find("Z") < options.size()) {
      flag_z = true;
   }

   if (flag_p && flag_g) {
      flag_g = false;
   }

   TH1D hist_temp = TH1D(*fHistogram);

   for (int i = 1; i <= fNBins; ++i) {
      double expectation = fOriginalNormalization * hist_temp.GetBinContent(i);
      double error = fOriginalNormalization * hist_temp.GetBinError(i);
      double n = 0;

      // throw random number according to Poisson distribution
      if (flag_p) {
         n = (double) fRandom->Poisson(expectation);
      }

      // throw random number according to Gauss distribution
      else if (flag_g) {
         double dn = fRandom->Gaus(expectation, error);

         // make it a truncated Gaussian
         if (flag_z) {
            while (n + dn < 0)
               dn = fRandom->Gaus(expectation, error);
         }
         n += dn;
      }

      // set the number of events in the template
      hist_temp.SetBinContent(i, n);
   }

   // normalize histogram
   double orignorm = hist_temp.Integral();

   if (orignorm)
      hist_temp.Scale(norm / orignorm);

   return hist_temp;
}

// ---------------------------------------------------------
