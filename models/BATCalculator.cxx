// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke, Stefan A. Schmitz


// include other header files

#include <TAxis.h>
#include <TFile.h>
#include <TCanvas.h>

#include <RooAbsFunc.h>
//#include <RooAbsReal.h>
#include <RooRealVar.h>
//#include <RooArgSet.h>
#include <RooBrentRootFinder.h>
#include <RooFormulaVar.h>
#include <RooGenericPdf.h>
//#include <RooPlot.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>
//#include <RooAbsPdf.h>
#include <RooHistPdf.h>

// include header file of this class 
//#include <RooStats/ModelConfig.h>

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>

#include "BCRooInterface.h"
#include "BATCalculator.h"

ClassImp(RooStats::BATCalculator)

namespace RooStats { 


BATCalculator::BATCalculator() :
  fData(0),
  fPdf(0),
  fPriorPOI(0),
  fProductPdf (0), fLogLike(0), fLikelihood (0), fIntegratedLikelihood (0), fPosteriorPdf(0), 
  fLower(0), fUpper(0),
  fBrfPrecision(0.00005),
  fValidInterval(false),
  fSize(0.05)
{
   // default constructor
   _myRooInterface = new BCRooInterface();
}

BATCalculator::BATCalculator( /* const char * name, const char * title, */						   
                             RooAbsData & data,
                             RooAbsPdf & pdf,
                             RooArgSet & POI,
                             RooAbsPdf & priorPOI,
                             RooAbsPdf * PriorNuisance,
                             RooArgSet * params) :
   //TNamed( TString(name), TString(title) ),
  fData(&data),
  fPdf(&pdf),
  fPOI(POI),
  fPriorPOI(&priorPOI),
  fPriorNuisance(PriorNuisance),
  fparams(params),
  fProductPdf (0), fLogLike(0), fLikelihood (0), fIntegratedLikelihood (0), fPosteriorPdf(0),
  fLower(0), fUpper(0),
  fBrfPrecision(0.00005),
  fValidInterval(false),  
  fSize(0.05)
{
   // constructor
   //if (nuisanceParameters) fNuisanceParameters.add(*nuisanceParameters); 
   _myRooInterface = new BCRooInterface();
}

/*BATCalculator::BATCalculator( RooAbsData & data,
                       ModelConfig & model) : 
   fData(&data), 
   fPdf(model.GetPdf()),
   fPriorPOI( model.GetPriorPdf()),
   fProductPdf (0), fLogLike(0), fLikelihood (0), fIntegratedLikelihood (0), fPosteriorPdf(0),
   fLower(0), fUpper(0),
   fBrfPrecision(0.00005),
   fValidInterval(false),
   fSize(0.05)
{
   // constructor from Model Config
   SetModel(model);
   _myRooInterface = new BCRooInterface();
}*/


BATCalculator::~BATCalculator()
{
   // destructor
   ClearAll(); 
   delete _myRooInterface;
}

void BATCalculator::ClearAll() const
{ 
   // clear cached pdf objects
   if (fProductPdf) delete fProductPdf; 
   if (fLogLike) delete fLogLike; 
   if (fLikelihood) delete fLikelihood; 
   if (fIntegratedLikelihood) delete fIntegratedLikelihood; 
   if (fPosteriorPdf) delete fPosteriorPdf;      
   fPosteriorPdf = 0; 
   fProductPdf = 0;
   fLogLike = 0; 
   fLikelihood = 0; 
   fIntegratedLikelihood = 0; 
   fLower = 0;
   fUpper = 0;
   fValidInterval = false;
}


void BATCalculator::SetModel(const ModelConfig & model)
{
   // set the model
   //fPdf = model.GetPdf();
   //fPriorPOI =  model.GetPriorPdf(); 
   // assignment operator = does not do a real copy the sets (must use add method) 
   //fPOI.removeAll();
   //*fparams.removeAll();
   //if (model.GetParametersOfInterest()) fPOI.add( *(model.GetParametersOfInterest()) );
   //if (model.GetNuisanceParameters())  fNuisanceParameters.add( *(model.GetNuisanceParameters() ) );

   // invalidate the cached pointers
   //ClearAll(); 
}


RooArgSet * BATCalculator::GetMode(RooArgSet * /* parameters */) const
{
  /// Returns the value of the parameters for the point in
  /// parameter-space that is the most likely.
  // Should cover multi-dimensional cases...
  // How do we do if there are points that are equi-probable?

  return 0;
}


RooAbsPdf * BATCalculator::GetPosteriorPdf() const
{
  /// build and return the posterior PDF

   // set nice style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file with default level of logging
   BCLog::OpenLog("results/bat_log.txt");
   BCLog::SetLogLevel(BCLog::detail);


   // run some sanity checks
   if (!fPdf ) {
	  BCLog::OutError("BATCalculator::GetPosteriorPdf - missing pdf model");
     //std::cerr << "BATCalculator::GetPosteriorPdf - missing pdf model" << std::endl;
     return 0;
   }
   if (!fPriorPOI) {
      BCLog::OutWarning("BATCalculator::GetPosteriorPdf - missing prior pdf");
      //std::cerr << "BATCalculator::GetPosteriorPdf - missing prior pdf" << std::endl;
   }
   if (fPOI.getSize() == 0) {
      BCLog::OutError("BATCalculator::GetPosteriorPdf - missing parameter of interest");
      //std::cerr << "BATCalculator::GetPosteriorPdf - missing parameter of interest" << std::endl;
      return 0;
   }
   if (fPOI.getSize() > 1) { 
      BCLog::OutError("BATCalculator::GetPosteriorPdf - current implementation works only on 1D intervals");
      //std::cerr << "BATCalculator::GetPosteriorPdf - current implementation works only on 1D intervals" << std::endl;
      return 0; 
   }


   _myRooInterface->Initialize(*fData,*fPdf,*fPriorPOI,fPriorNuisance,fparams,fPOI);
   int nMCMC = 30000;  //remove this hardcoded stuff later
   std::cout << "set interations to: " << nMCMC << endl;

   _myRooInterface->MCMCSetNIterationsRun(nMCMC);
   _myRooInterface->MarginalizeAll();
   _myRooInterface->FindMode();

   _myRooInterface->PrintAllMarginalized1D("results/bat_plots.ps");
   //_myRooInterface -> PrintAllMarginalized("bat_plots.ps");
   _myRooInterface -> PrintResults("results/bat_results.txt");

   const char * POIname = fPOI.first()->GetName();
   BCParameter * myPOI = _myRooInterface->GetParameter(POIname);
   BCH1D * myPosterior =_myRooInterface->GetMarginalized(myPOI);
   TH1D * posteriorTH1D = myPosterior->GetHistogram(); //this seems to be healthy

   RooDataHist * posteriorRooDataHist = new RooDataHist("posteriorhist","", fPOI,posteriorTH1D);
   fPosteriorPdf = new RooHistPdf("posteriorPdf","",fPOI,*posteriorRooDataHist);
   
   /*
   //plots for debugging 
   TFile* debugfile = new TFile( "debug_posterior.root" ,"RECREATE");  
   if ( debugfile->IsOpen() ) cout << "File debug_posterior opened successfully" << endl;  
   debugfile->cd();
   //posteriorTH1D->Write("posteriorTH1D");
   //posteriorRooDataHist->Write("posteriorRooDataHist");
   fPosteriorPdf->Write("fPosteriorPdf");
   TCanvas myCanvas("canvasforfposterior","canvas for fposterior");
   myCanvas.cd();
   //fPosteriorPdf->Draw();
   //myCanvas.Write();
   RooRealVar nSig("nSig","",0.0001,200.); 
   RooPlot* pl = nSig.frame(); 
   fPosteriorPdf->plotOn(pl); 
   pl->Draw();
   pl->Write();
   myCanvas.Write();
   debugfile->Close();

   //just for testing: create plots with BAT
   char* outputFile = "bat_plots_debugging.ps"; 
   _myRooInterface -> PrintAllMarginalized(outputFile);
   //*/


/*
   // create a unique name for the product pdf 
   TString prodName = TString("product_") + TString(fPdf->GetName()) + TString("_") + TString(fPriorPOI->GetName() );   
   fProductPdf = new RooProdPdf(prodName,"",RooArgList(*fPdf,*fPriorPOI));

   RooArgSet* constrainedParams = fProductPdf->getParameters(*fData);

   // use RooFit::Constrain() to make product of likelihood with prior pdf
   fLogLike = fProductPdf->createNLL(*fData, RooFit::Constrain(*constrainedParams) );

   TString likeName = TString("likelihood_") + TString(fProductPdf->GetName());   
   fLikelihood = new RooFormulaVar(likeName,"exp(-@0)",RooArgList(*fLogLike));
   RooAbsReal * plike = fLikelihood; 
   if (fNuisanceParameters.getSize() > 0) { 
      fIntegratedLikelihood = fLikelihood->createIntegral(fNuisanceParameters);
      plike = fIntegratedLikelihood; 
   }

   // create a unique name on the posterior from the names of the components
   TString posteriorName = this->GetName() + TString("_posteriorPdf_") + plike->GetName(); 
   fPosteriorPdf = new RooGenericPdf(posteriorName,"@0",*plike);

   delete constrainedParams;*/

   return fPosteriorPdf; //is of type RooAbsPdf*
}

 
/*RooPlot* BATCalculator::GetPosteriorPlot() const
{ 
  /// return a RooPlot with the posterior PDF and the credibility region

  if (!fPosteriorPdf) GetPosteriorPdf(); 
  if (!fValidInterval) GetInterval(); 

  RooAbsRealLValue* poi = dynamic_cast<RooAbsRealLValue*>( fPOI.first() );
  assert(poi);

   RooPlot* plot = poi->frame();

   plot->SetTitle(TString("Posterior probability of parameter \"")+TString(poi->GetName())+TString("\""));  
   fPosteriorPdf->plotOn(plot,RooFit::Range(fLower,fUpper,kFALSE),RooFit::VLines(),RooFit::DrawOption("F"),RooFit::MoveToBack(),RooFit::FillColor(kGray));
   fPosteriorPdf->plotOn(plot);
   plot->GetYaxis()->SetTitle("posterior probability");
   
   return plot; 
}*/


SimpleInterval * BATCalculator::GetInterval() const
{
  /// returns a SimpleInterval with lower and upper bounds on the
  /// parameter of interest. Applies the central ordering rule to
  /// compute the credibility interval. Covers only the case with one
  /// single parameter of interest

   if (fValidInterval) 
	   BCLog::OutWarning("BATCalculator::GetInterval : recomputing interval for the same CL and same model");
      //std::cout << "BATCalculator::GetInterval:" 
      //          << "Warning : recomputing interval for the same CL and same model" << std::endl;

   RooRealVar* poi = dynamic_cast<RooRealVar*>( fPOI.first() ); 
   assert(poi);

   //cout << "Hei-ho 1" << endl;
   if (!fPosteriorPdf) fPosteriorPdf = (RooAbsPdf*) GetPosteriorPdf();
   //cout << "Hei-ho 2" << endl;
   RooAbsReal* cdf = fPosteriorPdf->createCdf(fPOI,RooFit::ScanParameters(100,2));
   //RooAbsReal* cdf = fPosteriorPdf->createCdf(fPOI,RooFit::ScanNoCdf());
   //cout << "Hei-ho 3" << endl;
   RooAbsFunc* cdf_bind = cdf->bindVars(fPOI,&fPOI);
   RooBrentRootFinder brf(*cdf_bind);
   brf.setTol(fBrfPrecision); // set the brf precision

   double tmpVal = poi->getVal();  // patch used because findRoot changes the value of poi

   double y = fSize/2;
   brf.findRoot(fLower,poi->getMin(),poi->getMax(),y);

   y=1-fSize/2;
   bool ret = brf.findRoot(fUpper,poi->getMin(),poi->getMax(),y);
   if (!ret)
	   BCLog::OutWarning("BATCalculator::GetInterval : Error returned from Root finder, estimated interval is not fully correct");
	   //std::cout << "BATCalculator::GetInterval: Warning:"
      //          << "Error returned from Root finder, estimated interval is not fully correct" 
      //          << std::endl;

   poi->setVal(tmpVal); // patch: restore the original value of poi

   delete cdf_bind;
   delete cdf;
   fValidInterval = true; 

   TString interval_name = TString("BayesianInterval_a") + TString(this->GetName());
   SimpleInterval * interval = new SimpleInterval(interval_name,*poi,fLower,fUpper,ConfidenceLevel());
   interval->SetTitle("SimpleInterval from BATCalculator");

   return interval;
}

} // end namespace RooStats

