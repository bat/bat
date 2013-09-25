// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke, Stefan A. Schmitz

#include <TAxis.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TVector.h>

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

#include <RooStats/MarkovChain.h>

#include "../../BAT/BCLog.h"
#include "../../BAT/BCAux.h"
#include "../../BAT/BCH1D.h"

#include "BCRooInterface.h"
#include "BATCalculator.h"

#include <algorithm>
#include <cassert>


ClassImp(RooStats::BATCalculator)

namespace RooStats {

// ---------------------------------------------------------
BATCalculator::BATCalculator()
   : fData(0)
   , fPdf(0)
   , fPrior(0)
   , fparams(0)
   , _posteriorTH1D(0)
   , fProductPdf(0)
   , fLogLike(0)
   , fLikelihood(0)
   , fIntegratedLikelihood(0)
   , fPosteriorPdf(0)
   , fLower(0)
   , fUpper(0)
   , fBrfPrecision(0.00005)
   , fValidInterval(false)
   , fValidMCMCInterval(false)
   , fConnectedInterval(false)
   , _nMCMC(0)
   , fSize(0.05)
   , fLeftSideFraction(0.5)
{
   // default constructor
   _myRooInterface = new BCRooInterface();
}

// ---------------------------------------------------------
BATCalculator::BATCalculator( /* const char * name,  const char * title, */
                              RooAbsData & data,
                              RooAbsPdf  & pdf,
                              RooArgSet  & POI,
                              RooAbsPdf  & prior,
                              RooArgSet  * params,
                              bool fillChain)
   : fData(&data)
   , fPdf(&pdf)
   , fPOI(POI)
   , fPrior(&prior)
   , fparams(params)
   , _posteriorTH1D(0)
   , fProductPdf(0)
   , fLogLike(0)
   , fLikelihood(0)
   , fIntegratedLikelihood(0)
   , fPosteriorPdf(0)
   , fLower(0)
   , fUpper(0)
   , fBrfPrecision(0.00005)
   , fValidInterval(false)
   , fValidMCMCInterval(false)
   , fConnectedInterval(false)
   , _nMCMC(1000000)
   , fSize(0.05)
   , fLeftSideFraction(0.5)
   //, TNamed( TString(name), TString(title) )
{
   // constructor
   //if (nuisanceParameters)
   //   fNuisanceParameters.add(*nuisanceParameters);
   _myRooInterface = new BCRooInterface("BCRooInterfaceForBAT",fillChain);
}

// ---------------------------------------------------------
BATCalculator::BATCalculator( RooAbsData & data, ModelConfig & model, bool fillChain)
   : fData(&data)
   , fPdf(model.GetPdf())
   , fPOI(*model.GetParametersOfInterest())
   , fPrior(model.GetPriorPdf())
   , fparams(model.GetNuisanceParameters())
   , _posteriorTH1D(0)
   , fProductPdf(0)
   , fLogLike(0)
   , fLikelihood(0)
   , fIntegratedLikelihood(0)
   , fPosteriorPdf(0)
   , fLower(0)
   , fUpper(0)
   , fBrfPrecision(0.00005)
   , fValidInterval(false)
   , fValidMCMCInterval(false)
   , fConnectedInterval(false)
   , _nMCMC(1000000)
   , fSize(0.05)
   , fLeftSideFraction(0.5)
{
   std::cout << "BATCalculator calling constructor ..." << std::endl;
   // constructor from Model Config
   //  SetModel(model);
   _myRooInterface = new BCRooInterface("BCRooInterfaceForBAT",fillChain);
   std::cout << "BATCalculator constructed" << std::endl;
}

// ---------------------------------------------------------
BATCalculator::~BATCalculator()
{
   std::cout << "BATCalculator calling destructor ..." << std::endl;
   // destructor
   ClearAll();
   delete _myRooInterface;
}

// ---------------------------------------------------------
void BATCalculator::ClearAll() const
{
   // clear cached pdf objects
   if (fProductPdf)
      delete fProductPdf;
   if (fLogLike)
      delete fLogLike;
   if (fLikelihood)
      delete fLikelihood;
   if (fIntegratedLikelihood)
      delete fIntegratedLikelihood;
   if (fPosteriorPdf)
      delete fPosteriorPdf;
   fPosteriorPdf = 0;
   fProductPdf = 0;
   fLogLike = 0;
   fLikelihood = 0;
   fIntegratedLikelihood = 0;
   fLower = 0;
   fUpper = 0;
   fValidInterval = false;
   fValidMCMCInterval = false;
}

// ---------------------------------------------------------
void BATCalculator::SetModel(const ModelConfig & /* model */)
{
   // set the model
   //fPdf = model.GetPdf();
   //fPriorPOI =  model.GetPriorPdf();
   // assignment operator = does not do a real copy the sets (must use add method)
   //fPOI.removeAll();
   //*fparams.removeAll();
   //if (model.GetParametersOfInterest())
   //   fPOI.add( *(model.GetParametersOfInterest()) );
   //if (model.GetNuisanceParameters())
   //   fNuisanceParameters.add( *(model.GetNuisanceParameters() ) );

   // invalidate the cached pointers
   //ClearAll();
}

// ---------------------------------------------------------
RooArgSet * BATCalculator::GetMode(RooArgSet * /* parameters */) const
{
  /// Returns the value of the parameters for the point in
  /// parameter-space that is the most likely.
  // Should cover multi-dimensional cases...
  // How do we do if there are points that are equi-probable?

  return 0;
}

// ---------------------------------------------------------
//returns posterior with respect to first entry in the POI ArgSet
RooAbsPdf * BATCalculator::GetPosteriorPdf1D() const
{
   const char * POIname = fPOI.first()->GetName();
   return GetPosteriorPdf1D(POIname);
}

// ---------------------------------------------------------
// returns posterior with respect to provided name in the POI ArgSet
RooAbsPdf * BATCalculator::GetPosteriorPdf1D(const char * POIname) const
{

   // run some sanity checks
   if (!fPdf ) {
     std::cerr << "BATCalculator::GetPosteriorPdf - missing pdf model" << std::endl;
     return 0;
   }

   if (!fPrior) {
      std::cerr << "BATCalculator::GetPosteriorPdf - missing prior pdf" << std::endl;
   }

   if (fPOI.getSize() == 0) {
     std::cerr << "BATCalculator::GetPosteriorPdf - missing parameter of interest" << std::endl;
     return 0;
   }

   if (fPOI.getSize() > 1) {
      std::cerr << "BATCalculator::GetPosteriorPdf - current implementation works only on 1D intervals" << std::endl;
      return 0;
   }

    //initialize RooInterface object
   _myRooInterface->Initialize(*fData,*fPdf,*fPrior,fparams,fPOI);
   _myRooInterface->MCMCSetNIterationsRun(_nMCMC);
   _myRooInterface->MarginalizeAll();
   _myRooInterface->FindMode();
   const BCParameter * myPOI = _myRooInterface->GetParameter(POIname);
   BCH1D * myPosterior =_myRooInterface->GetMarginalized(myPOI);
   TH1D * posteriorTH1D = myPosterior->GetHistogram();
   _posteriorTH1D = static_cast<TH1D *>(posteriorTH1D->Clone("_posteriorTH1D"));
   RooDataHist * posteriorRooDataHist = new RooDataHist("posteriorhist","", fPOI,posteriorTH1D);
   fPosteriorPdf = new RooHistPdf("posteriorPdf","",fPOI,*posteriorRooDataHist);

   /*
   // plots for debugging
   TFile * debugfile = TFile::Open( "debug_posterior.root" ,"RECREATE");
   if ( debugfile->IsOpen() )
      std::cout << "File debug_posterior opened successfully" << std::endl;
   debugfile->cd();
   //posteriorTH1D->Write("posteriorTH1D");
   //posteriorRooDataHist->Write("posteriorRooDataHist");
   fPosteriorPdf->Write("fPosteriorPdf");
   TCanvas myCanvas("canvasforfposterior","canvas for fposterior");
   myCanvas.cd();
   //fPosteriorPdf->Draw();
   //myCanvas.Write();
   RooRealVar nSig("nSig","",0.0001,200.);
   RooPlot * pl = nSig.frame();
   fPosteriorPdf->plotOn(pl);
   pl->Draw();
   pl->Write();
   myCanvas.Write();
   debugfile->Close();

   // just for testing: create plots with BAT
   char * outputFile = "bat_plots_debugging.ps";
   _myRooInterface -> PrintAllMarginalized(outputFile);
   //*/

   return fPosteriorPdf; // is of type RooAbsPdf*
}

// ---------------------------------------------------------
// return a RooPlot with the posterior PDF and the credibility region
RooPlot * BATCalculator::GetPosteriorPlot1D() const
{


   if (!fPosteriorPdf){
      std::cout << "BATCalculator::GetPosteriorPlot1D:"
                << "Warning : posterior not available" << std::endl;
      GetPosteriorPdf1D();
   }
   if ((!fValidInterval) && (!fValidMCMCInterval)){
      std::cout << "BATCalculator::GetPosteriorPlot1D:"
                << "Warning : interval not available" << std::endl;
      GetInterval1D();
   }

   RooAbsRealLValue * poi = dynamic_cast<RooAbsRealLValue *>( fPOI.first() );
   assert(poi);

   RooPlot* plot = poi->frame();

   plot->SetTitle(TString("Posterior probability of parameter \"")+TString(poi->GetName())+TString("\""));
   fPosteriorPdf->plotOn(plot,RooFit::Range(fLower,fUpper,kFALSE),RooFit::VLines(),RooFit::DrawOption("F"),RooFit::MoveToBack(),RooFit::FillColor(kGray));
   fPosteriorPdf->plotOn(plot);
   plot->GetYaxis()->SetTitle("posterior probability");

   return plot;
}

// ---------------------------------------------------------
// returns central interval for first POI in the POI ArgSet
SimpleInterval * BATCalculator::GetInterval1D() const
{
   const char * POIname = fPOI.first()->GetName();
   return GetInterval1D(POIname); //is of type RooAbsPdf *
}


// ---------------------------------------------------------
// returns central interval for the defined POI in the POI ArgSet -> test code because orginal version is not working anymore
SimpleInterval * BATCalculator::GetInterval1D(const char * POIname) const
{
   //const char * POIname = fPOI.first()->GetName();

   if (fValidInterval)
   std::cout << "BATCalculator::GetShortestInterval1D:"
             << "Warning : recomputing interval for the same CL and same model" << std::endl;

   // get pointer to selected parameter of interest
   RooRealVar * poi = dynamic_cast<RooRealVar*>( fPOI.find(POIname) );
   assert(poi);

   // get pointer to poserior pdf
   if (!fPosteriorPdf)
      fPosteriorPdf = (RooAbsPdf*) GetPosteriorPdf1D();
   //RooAbsReal * cdf = fPosteriorPdf->createCdf(fPOI,RooFit::ScanParameters(100,2));

   // range of poi
   Double_t minpoi = poi->getMin();
   Double_t maxpoi = poi->getMax();

   // bin number of histogram for posterior
   Int_t stepnumber = _posteriorTH1D->GetNbinsX();
   std::cout << "stepnumber is: " << stepnumber << std::endl;

   // width of one bin in histogram of posterior
   Double_t stepsize = (maxpoi-minpoi)/stepnumber;
   std::cout << "stepsize is: " << stepsize << std::endl;

   // pair: first entry: bin number , second entry: value of posterior
   std::vector< std::pair< Int_t,Double_t > > posteriorVector;

   // for normalization
   Double_t histoIntegral = 0;
   // give posteriorVector the right length
   posteriorVector.resize(stepnumber);

   // see in BayesianCalculator for details about this "feature"
   // FB: comment in because 'unused variable'
   //   Double_t tmpVal = poi->getVal();

  // initialize elements of posteriorVector
   int i = 0;
   std::vector< std::pair< Int_t,Double_t > >::iterator vecit = posteriorVector.begin();
   std::vector< std::pair< Int_t,Double_t > >::iterator vecit_end = posteriorVector.end();
   for( ; vecit != vecit_end ; ++vecit) {
      poi->setVal(poi->getMin()+i*stepsize);
      posteriorVector[i] = std::make_pair(i, _posteriorTH1D->GetBinContent(i+1) ); // hope this is working, +1 necessary, because GetBinContent(0) returns the underflow bin
      histoIntegral+=_posteriorTH1D->GetBinContent(i); // better to get this directly from the histogram ?!
      //std::cout << "pair with binnumber: " << i << " and postriorprob: " << _posteriorTH1D->GetBinContent(i+1) << std::endl;
      i++;
   }

   std::cout << "histoIntegral is: " << histoIntegral << std::endl;

   double lowerProbLim = (1.-ConfidenceLevel())*fLeftSideFraction;
   double upperProbLim = 1. -( (1.-ConfidenceLevel())*(1.-fLeftSideFraction) );
   std::cout << "lowerProbLim is: " << lowerProbLim  << "upperProbLim is: " << histoIntegral << std::endl;

   fLower = -1.;
   fUpper = -1.;
   bool lowerlimset = false;
   bool upperlimset = false;


   if(fLeftSideFraction == 0){
      fLower = minpoi;
      lowerlimset = true;
   }

   if(fLeftSideFraction == 1){
      fUpper = maxpoi;
      upperlimset = true;
   }

   // keep track of integrated posterior in the interval
   Double_t integratedposterior = 0.;

   i = 0;
   vecit = posteriorVector.begin();
   for( ; vecit != vecit_end ; ++vecit) {
      integratedposterior += posteriorVector[i].second;
      if( (lowerlimset != true) && ((integratedposterior/histoIntegral) >= lowerProbLim) ){
         fLower = poi->getMin()+i*stepsize;
         lowerlimset = true;
      }
      if( (upperlimset != true) && ((integratedposterior/histoIntegral) >= upperProbLim) ){
         fUpper = poi->getMin()+i*stepsize;
         upperlimset = true;
         break;
      }
      i++;
   }

   TString interval_name = TString("CentralBayesianInterval_a") + TString(this->GetName());
   SimpleInterval * interval = new SimpleInterval(interval_name,*poi,fLower,fUpper,ConfidenceLevel());
   interval->SetTitle("SimpleInterval from BATCalculator");

   return interval;

}

// ---------------------------------------------------------
// returns central interval for requested POI -> no longer working since root v5_32
// SimpleInterval * BATCalculator::GetInterval1Dv0(const char * POIname) const
// {
//   /// returns a SimpleInterval with lower and upper bounds on the
//   /// parameter of interest. Applies the central ordering rule to
//   /// compute the credibility interval. Covers only the case with one
//   /// single parameter of interest
//
//    if (fValidInterval)
//       std::cout << "BATCalculator::GetInterval1D:"
//                 << "Warning : recomputing interval for the same CL and same model" << std::endl;
//
//    RooRealVar * poi = dynamic_cast<RooRealVar *>( fPOI.find(POIname) );
//    assert(poi);
//
//    if (!fPosteriorPdf)
// 	   fPosteriorPdf = (RooAbsPdf*) GetPosteriorPdf1D();
//
//    Double_t* myarr = new Double_t[1];
//    myarr[0]= 0.001;
//    std::cout << "fPosteriorPdf: " << fPosteriorPdf->maxVal(1) << std::endl;
//    std::cout << "fPosteriorPdf: " << fPosteriorPdf->getVal() << std::endl;
// //   std::cout << "fPosteriorPdf: " << fPosteriorPdf->evaluate() << std::endl;
//
//    RooAbsReal * cdf = fPosteriorPdf->createCdf(fPOI,RooFit::ScanParameters(100,2));
//    //RooAbsReal* cdf = fPosteriorPdf->createCdf(fPOI,RooFit::ScanNoCdf());
//
//    RooAbsFunc * cdf_bind = cdf->bindVars(fPOI,&fPOI);
//    std::cout << "cdf_bind: " << std::endl;
//    std::cout << "cdf_bind: " << cdf_bind->isValid() << std::endl;
//    std::cout << "cdf_bind: " << cdf_bind->getDimension() << std::endl;
//    //Double_t test = (*_function)(&a) - value;
//    std::cout << "cdf_bind: " << (*cdf_bind)(myarr) << std::endl;
//    std::cout << "cdf_bind: " << cdf_bind->getMaxLimit(1) << std::endl;
//
//    RooBrentRootFinder brf(*cdf_bind);
//    brf.setTol(fBrfPrecision); // set the brf precision
//
//    double tmpVal = poi->getVal();  // patch used because findRoot changes the value of poi
//
//    double y = fSize*fLeftSideFraction; //adjust lower tail prob. according to fLeftSideFraction
//
//    brf.findRoot(fLower,poi->getMin(),poi->getMax(),y);
//
//    y = 1.-(fSize*(1.-fLeftSideFraction) ); //adjust upper tail prob. according to fLeftSideFraction
//
//    bool ret = brf.findRoot(fUpper,poi->getMin(),poi->getMax(),y);
//    std::cout << " brf.findRoot: " << fUpper << " " << poi->getMin() << " " << poi->getMax() << " " << y << std::endl;
//    if (!ret)
//       std::cout << "BATCalculator::GetInterval1D: Warning:"
//                 << "Error returned from Root finder, estimated interval is not fully correct"
//                 << std::endl;
//
//    poi->setVal(tmpVal); // patch: restore the original value of poi
//
//    delete cdf_bind;
//    delete cdf;
//    fValidInterval = true;
//    fConnectedInterval = true;
//
//    TString interval_name = TString("CentralBayesianInterval_a") + TString(this->GetName());
//    SimpleInterval * interval = new SimpleInterval(interval_name,*poi,fLower,fUpper,ConfidenceLevel());
//    interval->SetTitle("SimpleInterval from BATCalculator");
//
//    return interval;
// }

// ---------------------------------------------------------
SimpleInterval * BATCalculator::GetShortestInterval1D() const
{
   const char * POIname = fPOI.first()->GetName();
   bool checkConnected = true;
   return GetShortestInterval1D(POIname, checkConnected);
}

// ---------------------------------------------------------
// returns a SimpleInterval with lower and upper bounds on the
// parameter of interest. Applies the shortest interval rule to
// compute the credibility interval. The resulting interval is not necessarily connected.
// Methods for constructing the shortest region with respect given CL are now available in
// different places (here; MCMCCalculator, BAT)-> this might require cleaning at some point.
// Result is approximate as CL is not reached exactly (due to finite bin number/bin size in posterior histogram)-> make CL/interval a bit smaller or bigger ?
SimpleInterval * BATCalculator::GetShortestInterval1D(const char * POIname, bool & checkConnected) const
{

   if (fValidInterval)
      std::cout << "BATCalculator::GetShortestInterval1D:"
                << "Warning : recomputing interval for the same CL and same model" << std::endl;

   // get pointer to selected parameter of interest
   RooRealVar * poi = dynamic_cast<RooRealVar*>( fPOI.find(POIname) );
   assert(poi);

   // get pointer to poserior pdf
   if (!fPosteriorPdf)
      fPosteriorPdf = (RooAbsPdf*) GetPosteriorPdf1D();
   //RooAbsReal * cdf = fPosteriorPdf->createCdf(fPOI,RooFit::ScanParameters(100,2));

   // range of poi
   Double_t minpoi = poi->getMin();
   Double_t maxpoi = poi->getMax();

   // bin number of histogram for posterior
   Int_t stepnumber = _posteriorTH1D->GetNbinsX();
   std::cout << "stepnumber is: " << stepnumber << std::endl;

   // width of one bin in histogram of posterior
   Double_t stepsize = (maxpoi-minpoi)/stepnumber;
   std::cout << "stepsize is: " << stepsize << std::endl;

   // pair: first entry: bin number , second entry: value of posterior
   std::vector< std::pair< Int_t,Double_t > > posteriorVector;

   // for normalization
   Double_t histoIntegral = 0;
   // give posteriorVector the right length
   posteriorVector.resize(stepnumber);

   // see in BayesianCalculator for details about this "feature"
   Double_t tmpVal = poi->getVal();

   // initialize elements of posteriorVector
   int i = 0;
   std::vector< std::pair< Int_t,Double_t > >::iterator vecit = posteriorVector.begin();
   std::vector< std::pair< Int_t,Double_t > >::iterator vecit_end = posteriorVector.end();
   for( ; vecit != vecit_end ; ++vecit) {
      poi->setVal(poi->getMin()+i*stepsize);
      posteriorVector[i] = std::make_pair(i, _posteriorTH1D->GetBinContent(i+1) ); // hope this is working, +1 necessary, because GetBinContent(0) returns the underflow bin
      histoIntegral+=_posteriorTH1D->GetBinContent(i); // better to get this directly from the histogram ?!
      //std::cout << "pair with binnumber: " << i << " and postriorprob: " << _posteriorTH1D->GetBinContent(i+1) << std::endl;
      i++;
   }

   std::cout << "histoIntegral is: " << histoIntegral << std::endl;

   // sort posteriorVector with respect to posterior pdf
   std::sort(posteriorVector.begin(), posteriorVector.end(), sortbyposterior);

   // keep track of integrated posterior in the interval
   Double_t integratedposterior = 0.;

   // keep track of lowerLim and upperLim
   Double_t lowerLim=posteriorVector.size();
   Double_t upperLim=0;

   // store the bins in the intervall
   std::vector<bool> inInterval;
   inInterval.resize(posteriorVector.size());

   // set all values in inInterval to false
   for (unsigned int k = 0; k < inInterval.size(); k++)
      inInterval[k] = false;

   unsigned int j = 0;
   // add bins to interval while CL not reached
   //std::cout << "integratedposterior: " << integratedposterior << "histoIntegral: " << histoIntegral << std::endl;
   //std::cout << " integratedposterior/histoIntegral : " << integratedposterior/histoIntegral << std::endl;
   //std::cout << " 1-fSize : " << 1-fSize << std::endl;
   //std::cout << " posteriorVector.size(): " << posteriorVector.size() << std::endl;
   while(((integratedposterior/histoIntegral) < (1-fSize)) && (j < posteriorVector.size())) {
      //std::cout << "j is: " << j << " , int prob: " << integratedposterior/histoIntegral << std::endl;
      integratedposterior+=posteriorVector[j].second;

      //std::cout << "bin number: " << posteriorVector[j].first << " with posterior prob.: " << posteriorVector[j].second << std::endl;

      // update vector with bins included in the interval
      inInterval[posteriorVector[j].first] = true;

      if(posteriorVector[j].first < lowerLim) {
         lowerLim = posteriorVector[j].first; // update lowerLim
         std::cout << "updating lower lim to: " << lowerLim << std::endl;
      }
      if(posteriorVector[j].first > upperLim) {
         upperLim = posteriorVector[j].first; // update upperLim
         std::cout << "updating upper lim to: " << upperLim << std::endl;
      }

      fLower = lowerLim * stepsize; //
      fUpper = upperLim * stepsize;
      j++;
   }

   // initialize vector with interval borders

   bool runInside = false;
   for (unsigned int l = 0; l < inInterval.size(); l++) {
      if ( (runInside == false) && (inInterval[l] == true) ) {
         _intervalBorders1D.push_back(static_cast<double>(l)* stepsize);
         runInside = true;
      }
      if ( ( runInside == true) && (l < (inInterval.size()-1) ) && (inInterval[l+1] == false) ) {
         _intervalBorders1D.push_back(static_cast<double>(l)* stepsize);
         runInside = false;
      }
      if ( ( runInside == true) && (l == (inInterval.size()-1)) ) {
         _intervalBorders1D.push_back(static_cast<double>(l)* stepsize);
      }
   }

   // check if the intervall is connected
   if(checkConnected) {
      if (_intervalBorders1D.size() > 2) {
         fConnectedInterval = false;
      }
      else {
         fConnectedInterval = true;
      }
   }

   poi->setVal(tmpVal); // patch: restore the original value of poi

   fValidInterval = true;

   TString interval_name = TString("ShortestBayesianInterval_a") + TString(this->GetName());
   SimpleInterval * interval = new SimpleInterval(interval_name,*poi,fLower,fUpper,ConfidenceLevel());
   interval->SetTitle("Shortest SimpleInterval from BATCalculator");

   //if(assert)
   return interval;
}

// ---------------------------------------------------------
MCMCInterval* BATCalculator::GetInterval() const{

   // run some sanity checks
   if (!fPdf ) {
     std::cerr << "BATCalculator::GetInterval - missing pdf model" << std::endl;
     return 0;
   }

   if (!fPrior) {
      std::cerr << "BATCalculator::GetInterval - missing prior pdf" << std::endl;
   }

   if (fPOI.getSize() == 0) {
     std::cerr << "BATCalculator::GetInterval - missing parameter of interest" << std::endl;
     return 0;
   }

   if (!fPosteriorPdf){
      //initialize RooInterface object
      _myRooInterface->Initialize(*fData,*fPdf,*fPrior,fparams,fPOI);
      _myRooInterface->MCMCSetNIterationsRun(_nMCMC);
      _myRooInterface->MarginalizeAll();
      _myRooInterface->FindMode();
   }

   MarkovChain * roostats_chain = GetBCRooInterface()->GetRooStatsMarkovChain();
   MCMCInterval * mcmcInterval = new MCMCInterval("roostatsmcmcinterval", fPOI , *roostats_chain);
   //MCMCInterval * mcmcInterval = new MCMCInterval("roostatsmcmcinterval", *(GetBCRooInterface()->GetArgSetForMarkovChain()) , *roostats_chain);
   mcmcInterval->SetUseKeys(false);
   mcmcInterval->SetConfidenceLevel(1.-fSize);
   fValidMCMCInterval = true;
   return mcmcInterval;
}

// ---------------------------------------------------------
void BATCalculator::SetNumBins(const char * parname, int nbins)
{
   _myRooInterface->SetNumBins(parname, nbins);
}

void BATCalculator::SetNumBins(int nbins)
{
   _myRooInterface->SetNumBins(nbins);
}

void BATCalculator::SetLeftSideTailFraction(Double_t leftSideFraction ){
   if( (leftSideFraction >= 0.) && (leftSideFraction <= 1.) ){
      fLeftSideFraction = leftSideFraction;
   }
   else{
      std::cout << "BATCalculator::SetLeftSideTailFraction(Double_t leftSideFraction ) - value needs to be in the interval [0.,1.] to be meaningful, your value: " << leftSideFraction <<  " ,left side tail fraction has not been changed!" << std::endl;
   }

}


// ---------------------------------------------------------
Double_t BATCalculator::GetOneSidedUperLim()
{
//   double safeVal = fSize;
//   fSize = safeVal/2.;
   std::cout << "calculating " << (1.-fSize/2) << "upper limit" << std::endl;
   return GetInterval1D()->UpperLimit();
}

/*
// ---------------------------------------------------------
int BATCalculator::GetNbins(const char * parname)
{
   ;
}
*/

// ---------------------------------------------------------
bool sortbyposterior(std::pair< Int_t,Double_t > pair1, std::pair< Int_t,Double_t > pair2)
{
   return (pair1.second > pair2.second);
}



} // end namespace RooStats
