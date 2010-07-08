// 
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke, Stefan A. Schmitz

#ifndef ROOSTATS_BATCalculator
#define ROOSTATS_BATCalculator

#include <TNamed.h>

#include <RooStats/IntervalCalculator.h>
#include <RooStats/SimpleInterval.h>

#include <RooArgSet.h>
#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooPlot.h>
#include <RooAbsReal.h>

#include <BCRooInterface.h>


namespace RooStats
{

   class BATCalculator : public IntervalCalculator, public TNamed
	{

   public:

      // constructor
      BATCalculator( );

      BATCalculator( RooAbsData & data,
                     RooAbsPdf & pdf,
                     RooArgSet & POI,
                     RooAbsPdf & priorPOI,
                     RooAbsPdf * PriorNuisance = 0,
                     RooArgSet * params = 0 );


      BATCalculator( RooAbsData & data,
                     ModelConfig & model ); 

      // destructor
      virtual ~BATCalculator();

      /*RooPlot * GetPosteriorPlot() const; */

      // return posterior pdf (object is managed by the BayesianCalculator class)
      RooAbsPdf * GetPosteriorPdf() const ; 

      virtual SimpleInterval * GetInterval() const ; 

      virtual void SetData( RooAbsData & data ) {
         fData = &data;
         ClearAll();
      }

      virtual void SetModel( const ModelConfig& model ); 

      // set the size of the test (rate of Type I error) ( Eg. 0.05 for a 95% Confidence Interval)
      virtual void SetTestSize( Double_t size ) {
         fSize = size;
         fValidInterval = false; 
      }

      // set the confidence level for the interval (eg. 0.95 for a 95% Confidence Interval)
      virtual void SetConfidenceLevel( Double_t cl )
		   { SetTestSize(1.-cl); }

      // Get the size of the test (eg. rate of Type I error)
      virtual Double_t Size() const
		   { return fSize; }

      // Get the Confidence level for the test
      virtual Double_t ConfidenceLevel() const
		   { return 1.-fSize; }

      void SetBrfPrecision( double precision ) { fBrfPrecision = precision; }

      void CleanCalculatorForNewData() { ClearAll();}

   protected:

      void ClearAll() const; 
   
   private:

      // compute the most probable value: move to public once implemented
      // returns a RooArgSet
      RooArgSet * GetMode( RooArgSet * parameters ) const;
      // plan to replace the above: return a SimpleInterval integrating 
      // over all other parameters except the one specified as argument
      //virtual SimpleInterval* GetInterval( RooRealVar* parameter  ) const { return 0; }
    
      RooAbsData * fData;
      RooAbsPdf  * fPdf;
      mutable RooArgSet fPOI;
      RooAbsPdf  * fPriorPOI;
      RooAbsPdf  * fPriorNuisance;
      RooArgSet  * fparams;
      BCRooInterface * _myRooInterface;


      mutable RooAbsPdf  * fProductPdf; 
      mutable RooAbsReal * fLogLike; 
      mutable RooAbsReal * fLikelihood; 
      mutable RooAbsReal * fIntegratedLikelihood; 
      mutable RooAbsPdf  * fPosteriorPdf; 
      mutable Double_t  fLower; 
      mutable Double_t  fUpper; 
      double  fBrfPrecision;
      mutable Bool_t    fValidInterval;

      double fSize;  // size used for getting the interval

   protected:

      ClassDef(BATCalculator,1)  // BATCalculator class

   };
}

#endif
