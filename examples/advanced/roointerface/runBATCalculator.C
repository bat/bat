/////////////////////////////////////////////////////////////////////////
//
// 'Using BATCalculator ' macro No. 1
// author: Stefan A. Schmitz
// date: July 2010
// last update: October 2012
/////////////////////////////////////////////////////////////////////////

// this macro shows how to use the RooStats Interface (BATCalculator) to BAT
//
// about the statistics model in this tutorial:
// This tutorial addresses the task of evaluating a Bayesian counting experiment with
// background and systematic uncertainties. The observable is the measured number of events.
// Here the bayesian 90% interval (upper limit corresponds to one-sided 95% upper limit) is
// evaluated as a function of the observed number of events in a hypothetical experiment.

// The macro will need some time to run. You can adjust the number of Markov chain elements and
// the range of tested numbers of observed events to make the macro terminate faster. The result
// is a plot showing the 90% (central) interval // as a function of the observed number of events

#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"

#include "RooStats/SimpleInterval.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

void runBATCalculator()
{
   // Definiton of a RooWorkspace containing the statistics model. Later the
   // information for BATCalculator is retrieved from the workspace. This is
   // certainly a bit of overhead but better from an educative point of view.
   cout << "preparing the RooWorkspace object" << endl;

   RooWorkspace * myWS = new RooWorkspace("myWS",true);

   // combined prior for signal contribution
   myWS->factory("Product::signal({sigma_s[0,20],L[5,15],epsilon[0,1]})");
   myWS->factory("N_bkg[0,3]");
   // define prior functions
   // uniform prior for signal crosssection
   myWS->factory("Uniform::prior_sigma_s(sigma_s)");
   // (truncated) prior for efficiency
   myWS->factory("Gaussian::prior_epsilon(epsilon,0.51,0.0765)");
   // (truncated) Gaussian prior for luminosity
   myWS->factory("Gaussian::prior_L(L,10,1)");
   // (truncated) Gaussian prior for bkg crosssection
   myWS->factory("Gaussian::prior_N_bkg(N_bkg,0.52,0.156)");

   // Poisson distribution with mean signal+bkg
   myWS->factory("Poisson::model(n[0,300],sum(signal,N_bkg))");

   // define the global prior function
   myWS->factory("PROD::prior(prior_sigma_s,prior_epsilon,prior_L,prior_N_bkg)");

   // Definition of observables and parameters of interest
   myWS->defineSet("obsSet","n");
   myWS->defineSet("poiSet","sigma_s");
   myWS->defineSet("nuisanceSet","N_bkg,L,epsilon");

   // ->model complete (Additional information can be found in the
   // RooStats manual)

   //  feel free to vary the parameters, but don't forget to choose reasonable ranges for the
   // variables. Currently the Bayesian methods will often not work well if the variable ranges
   // are either too short (for obvious reasons) or too large (for technical reasons).

   // A ModelConfig object is used to associate parts of your workspace with their statistical
   // meaning (it is also possible to initialize BATCalculator directly with elements from the
   // workspace but if you are sharing your workspace with others or if you want to use several
   // different methods the use of ModelConfig will most often turn out to be the better choice.)

   // setup the ModelConfig object
   cout << "preparing the ModelConfig object" << endl;

   ModelConfig modelconfig("modelconfig","ModelConfig for this example");
   modelconfig.SetWorkspace(*myWS);

   modelconfig.SetPdf(*(myWS->pdf("model")));
   modelconfig.SetParametersOfInterest(*(myWS->set("poiSet")));
   modelconfig.SetPriorPdf(*(myWS->pdf("prior")));
   modelconfig.SetNuisanceParameters(*(myWS->set("nuisanceSet")));
   modelconfig.SetObservables(*(myWS->set("obsSet")));


   // use BATCalculator to the derive confidence intervals as a function of the observed number of
   // events in the hypothetical experiment

   // define vector with tested numbers of events
   TVectorD obsEvents;
   // define vectors which will be filled with the lower and upper limits for each tested number
   // of observed events
   TVectorD BATul;
   TVectorD BATll;

   // fix upper limit of tested observed number of events
   int obslimit = 10;

   obsEvents.ResizeTo(obslimit);
   BATul.ResizeTo(obslimit);
   BATll.ResizeTo(obslimit);


   cout << "starting the calculation of Bayesian confidence intervals with BATCalculator" << endl;
   // loop over observed number of events in the hypothetical experiment
   for (int obs = 1; obs<=obslimit; obs++){

      obsEvents[obs-1] = (static_cast<double>(obs));

      // prepare data input for the the observed number of events
      // adjust number of observed events in the workspace. This is communicated to ModelConfig!
      myWS->var("n")->setVal(obs);
      // create data
      RooDataSet data("data","",*(modelconfig.GetObservables()));
      data.add( *(modelconfig.GetObservables()));

      // prepare BATCalulator
      BATCalculator batcalc(data, modelconfig);

      // give the BATCalculator a unique name (always a good idea in ROOT)
      TString namestring = "mybatc_";
      namestring += obs;
      batcalc.SetName(namestring);

      // fix amount of posterior probability in the calculated interval.
      batcalc.SetConfidenceLevel(0.90);
      // fix number of Markov chain elements. (in general: the longer the Markov chain the more
      // precise will be the results)
      batcalc.SetnMCMC(300000);

      // retrieve SimpleInterval object containing the information about the interval (this
      // triggers the actual calculations)
      SimpleInterval * interval = batcalc.GetInterval1D("sigma_s");

      std::cout << "BATCalculator: 90% CL interval: [ " << interval->LowerLimit() << " - " << interval->UpperLimit() << " ] or 95% CL upper limit\n";

      // add the interval borders for the current number of observed events to the vectors
      // containing the lower and upper limits
      BATll[obs-1]= interval->LowerLimit();
      BATul[obs-1]= interval->UpperLimit();

      // clean up for next loop element
      batcalc.CleanCalculatorForNewData();
      delete interval;
   }
   cout << "all limits calculated" << endl;

   // summarize the results in a plot

   TGraph *grBATll = new TGraph(obsEvents,BATll);
   grBATll->SetLineColor(kGreen);
   grBATll->SetLineWidth(200);
   grBATll->SetFillStyle(3001);
   grBATll->SetFillColor(kGreen);

   TGraph *grBATul = new TGraph(obsEvents,BATul);
   grBATul->SetLineColor(kGreen);
   grBATul->SetLineWidth(-200);
   grBATul->SetFillStyle(3001);
   grBATul->SetFillColor(kGreen);

   // create and draw multigraph
   TMultiGraph *mg = new TMultiGraph("BayesianLimitsBATCalculator","BayesianLimitsBATCalculator");
   mg->SetTitle("example of Bayesian confidence intervals derived with BATCAlculator ");

   mg->Add(grBATll);
   mg->Add(grBATul);

   mg->Draw("AC");

   mg->GetXaxis()->SetTitle ("# observed events");
   mg->GetYaxis()->SetTitle("limits on signal S (size of test: 0.1)");

   mg->Draw("AC");

}
