#ifndef __BCMTFANALYSISFACILITY__H
#define __BCMTFANALYSISFACILITY__H

/*!
 * \class BCMTFAnalysisFacility
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

#include <TH1D.h>

class BCMTF;
class TTree;
class TRandom3;

// ---------------------------------------------------------
class BCMTFAnalysisFacility
{

   public:

      // Constructors and destructor
      BCMTFAnalysisFacility(BCMTF * mtf);
      ~BCMTFAnalysisFacility();

      // setters

      // set BCMTF
      void SetBCMTF(BCMTF * mtf)
         { fMTF = mtf; };

      // set flag for using MCMC (true) or not (false)
      void SetFlagMCMC(bool flag)
         { fFlagMCMC = flag; };

      // getters

      // return BCMTF
      BCMTF * GetBCMTF()
         { return fMTF; };

      // misc

      // perform the full set of single channel analyses and the
      // combination
      int PerformSingleChannelAnalyses(const char * dirname, const char * options = "");

      // perform the analysis using one systematic at a time, without
      // systematic and with all systematics
      // flag_nuisance: use nuisance parameters (true) or delta method (false)
      int PerformSingleSystematicAnalyses(const char * dirname, const char * options = "");

      // perform calibration curve
      int PerformCalibrationAnalysis(const char * dirname, const std::vector<double> & default_parameters, int index, const std::vector<double> & parametervalues, int nensembles = 1000);

      // perform full ensemble test
      int PerformEnsembleTest(const std::vector<double> & parameters);

      // build a single ensemble based on a single set of parameters
      std::vector<TH1D> BuildEnsemble(const std::vector<double> & parameters);

      // build ensembles based on a single set of parameters
      TTree * BuildEnsembles(const std::vector<double> & parameters, int nensembles);

      // build ensembles based on a varying sets of parameters, e.g., using the prior or posterior
      TTree * BuildEnsembles(TTree * tree, int nensembles);

      // perform ensemble test based on one set of parameters
      TTree * PerformEnsembleTest(const std::vector<double> & parameters, int nensembles);

      // perform ensemble test based on varying sets of parameters
      TTree * PerformEnsembleTest(TTree * tree, int nensembles, int start = 0);

      // transform a matrix to a set of histograms
      std::vector<TH1D> MatrixToHistograms(const std::vector< std::vector<double> > & matrix);

 private:

      // the multi template fitter
      BCMTF * fMTF;

      // a random number generator
      TRandom3 * fRandom;

      // flag: use MCMC for analysis
      bool fFlagMCMC;
};
// ---------------------------------------------------------

#endif

