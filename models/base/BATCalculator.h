/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * Original authors: Gregory Schott and Stefan A. Schmitz with inspiration from
 * Kyle Cranmer, Lorenzo Moneta, Gregory Schott,  and Wouter Verkerke.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#ifndef ROOSTATS_BATCalculator
#define ROOSTATS_BATCalculator

#include "BCRooInterface.h"

#include <TH1D.h>
#include <TNamed.h>

#include <RooStats/IntervalCalculator.h>
#include <RooStats/MCMCInterval.h>
#include <RooStats/ModelConfig.h>
#include <RooStats/SimpleInterval.h>

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooArgSet.h>
#include <RooPlot.h>

#include <map>
#include <vector>

namespace RooStats
{

class BATCalculator : public IntervalCalculator, public TNamed
{

public:

    /// constructor
    BATCalculator( );

    BATCalculator( RooAbsData& data,
                   RooAbsPdf&   pdf,
                   RooArgSet&   POI,
                   RooAbsPdf&   prior,
                   RooArgSet*   params = 0,
                   bool fillChain = false );

    BATCalculator( RooAbsData&   data,
                   ModelConfig& model,
                   bool fillChain = false);

    /// destructor
    virtual ~BATCalculator();

    RooPlot* GetPosteriorPlot1D() const;

    /// return posterior pdf
    RooAbsPdf* GetPosteriorPdf1D() const;
    RooAbsPdf* GetPosteriorPdf1D(const char* POIname) const;

    /// return SimpleInterval object: central interval (one poi only)
    virtual SimpleInterval* GetInterval1D() const;
    virtual SimpleInterval* GetInterval1D(const char* POIname) const;

    /// return SimpleInterval object: shortest interval (not necessarily connected, one poi only)
    SimpleInterval* GetShortestInterval1D() const;
    SimpleInterval* GetShortestInterval1D(const char* POIname, bool& checkConnected) const;

    /// temporary solution ?
    Double_t GetOneSidedUperLim();

    virtual void SetData( RooAbsData& data )
    { fData = &data; ClearAll(); }

    virtual void SetModel(const ModelConfig&) {}

    /// set the size of the test (rate of Type I error) ( Eg. 0.05 for a 95% Confidence Interval)
    virtual void SetTestSize( Double_t size )
    { fSize = size; fValidInterval = false; }

    /// set left side tail fraction (only for 1D interval, not meaningful for shortest interval)
    void SetLeftSideTailFraction(Double_t leftSideFraction );

    /// set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
    virtual void SetConfidenceLevel( Double_t cl )
    { SetTestSize(1. - cl); }

    /// Get the size of the test
    virtual Double_t Size() const
    { return fSize; }
    /// Get left side tail fraction (only for 1D interval, not meaningful for shortest interval)
    double GetLeftSideTailFraction()
    { return fLeftSideFraction; }

    /// Get the Confidence level for the test
    virtual Double_t ConfidenceLevel() const
    { return 1. - fSize; }

    void SetBrfPrecision( double precision )
    { fBrfPrecision = precision; }

    double GetBrfPrecision()
    { return fBrfPrecision; }

    ///set number of iterations per chain
    void SetnMCMC(int nMCMC)
    { _nMCMC = nMCMC; }

    ///return number of iterations per chain
    int  GetnMCMC()
    { return _nMCMC; }

    BCRooInterface* GetBCRooInterface() const
    { return _myRooInterface; }

    RooStats::MarkovChain* GetRooStatsMarkovChain() const
    { return _myRooInterface->GetRooStatsMarkovChain();}

    /// returns MCMCInterval object
    virtual MCMCInterval* GetInterval() const;

    ///returns if last calculated shortest interval is connected (1 poi only)
    bool GetConnected()
    { return fConnectedInterval; }

    /// returns  interval borders of shortest interval (1 poi only)
    std::vector<double> GetIntervalBorders1D()
    { return _intervalBorders1D; }

    ///set the number of histogram bins for a specific parameter
    void SetNumBins(const char* parname, int nbins);

    ///set the number of histogram bins for all parameters
    void SetNumBins(int nbins);

    void CleanCalculatorForNewData()
    { ClearAll(); }

protected:

    void ClearAll() const;

private:
    RooAbsData* fData;
    RooAbsPdf* fPdf;
    const RooArgSet fPOI;
    RooAbsPdf* fPrior;
    const RooArgSet* fparams;
    BCRooInterface* _myRooInterface;
    mutable TH1D* _posteriorTH1D;


    mutable RooAbsPdf* fProductPdf;
    mutable RooAbsReal* fLogLike;
    mutable RooAbsReal* fLikelihood;
    mutable RooAbsReal* fIntegratedLikelihood;
    mutable RooAbsPdf* fPosteriorPdf;
    mutable Double_t fLower;
    mutable Double_t fUpper;
    double fBrfPrecision;
    mutable Bool_t fValidInterval;
    mutable Bool_t fValidMCMCInterval;
    mutable bool fConnectedInterval;

    int _nMCMC; /// number of chain elements per Markov Chain
    double fSize;  /// size used for the interval object
    double fLeftSideFraction; //
    mutable std::vector<double> _intervalBorders1D;

protected:

    ClassDef(BATCalculator, 1) // BATCalculator class

};
}

#endif
