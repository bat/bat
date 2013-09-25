// root -q -x -l 'PrepareWorkspace_GaussOverFlat_withSystematics.C()'

void PrepareWorkspace_GaussOverFlat_withSystematics( TString fileName = "WS_GaussOverFlat_withSystematics.root" )
{
  // In this macro a PDF model is built assuming signal has a Gaussian
  // PDF and the background a flat PDF.  The parameter of interest is
  // the signal yield and we assume for it a flat prior.  It is shown
  // how two types of systematics uncertainties can be expressed;
  // those are a sytematic uncertainty on the background yield and
  // another on one of the parameters (sigma) of the signal shape.
  // All needed objects are stored in a ROOT file (within a
  // RooWorkspace container); this ROOT file can then be fed as input
  // to various statistical methods.

  using namespace RooFit;
  using namespace RooStats;

  // use an observable for this shape-based analysis
  RooRealVar* mass = new RooRealVar("mass","mass",0,500,"GeV/c^{2}");
  mass->setBins(100);
  RooArgSet* observables = new RooArgSet(*mass,"observables");

  // signal (Gaussian) and background (flat) PDFs
  RooRealVar* sigSigma = new RooRealVar("sigSigma","sigma in signal PDF",0,100);
  RooAbsPdf* sigPdf = new RooGaussian("sigPdf","signal PDF",*mass,RooConst(200),*sigSigma);
  RooAbsPdf* bkgPdf = new RooPolynomial("bkgPdf","background PDF",*mass,RooFit::RooConst(0));
  
  // S+B model: the sum of both shapes weighted with the yields
  RooRealVar* S = new RooRealVar("S","signal yield",10,0,50);
  RooRealVar* B = new RooRealVar("B","background yield",10,0,50);
  RooAbsPdf* model = new RooAddPdf("model","S+B PDF",RooArgList(*sigPdf,*bkgPdf),RooArgList(*S,*B));
  
  // B-only model: the same as with a signal yield fixed to 0
  RooAbsPdf* modelBkg = new RooExtendPdf("modelBkg","B-only PDF",*bkgPdf,*B);

  // assume a Gaussian uncertainty on the detector resolution affecting the signal width (of 10%)
  // another nuisance parameter is the background yield (apply an uncertainty of 20%)
  RooAbsPdf* prior_sigSigma = new RooGaussian("prior_sigSigma","prior probability on sigSigma",*sigSigma,RooConst(sigSigma->getVal()),RooConst(sigSigma->getVal()*0.10));
  RooAbsPdf* prior_B = new RooGaussian("prior_B","prior probability on B",*B,RooConst(B->getVal()),RooConst(B->getVal()*0.20));
  RooAbsPdf* priorNuisance = new RooProdPdf("priorNuisance","prior on the nuisance parameters",*prior_B,*prior_sigSigma);
  RooArgSet* parameters = new RooArgSet(*B,*sigSigma,"parameters");

  // assume a flat prior on our parameter of interest (POI) which is the signal yield
  RooAbsPdf* priorPOI = new RooPolynomial("priorPOI","flat prior on the POI",*S,RooFit::RooConst(0));
  RooArgSet* POI = new RooArgSet(*S,"POI");
  
  // different options are shown for the data generation from the model
  
  // unbinned data with Poisson fluctuations
//   RooAbsData* data = (RooDataSet*) model->generate(*observables,RooFit::Extended(),Name("data"));

  // binned data with Poisson fluctuations
//   RooAbsData* data = (RooDataHist*) model->generateBinned(*observables,Extended(),Name("data"));
  
  // binned without any fluctuations (average case)
  RooAbsData* data = (RooDataHist*) model->generateBinned(*observables,Name("data"),ExpectedData());

  // control plot of the generated data
//   RooPlot* plot = mass->frame();
//   data->plotOn(plot);
//   plot->Draw();

  // use a RooWorkspace to store the pdf models, prior informations, list of parameters,...
  RooWorkspace myWS("myWS");
  myWS.import(*data,Rename("data"));
  myWS.import(*model,RecycleConflictNodes());
  myWS.import(*modelBkg,RecycleConflictNodes());
  myWS.import(*priorPOI,RecycleConflictNodes());
  myWS.import(*priorNuisance,RecycleConflictNodes());  
  myWS.defineSet("observables",*observables,kTRUE);
  myWS.defineSet("parameters",*parameters,kTRUE);
  myWS.defineSet("POI",*POI,kTRUE);

  // store the workspace in a ROOT file  
  TFile file(fileName,"RECREATE");
  file.cd();
  myWS.Write();
  file.Write();
  file.Close();
  
  std::cout << "\nRooFit model initialized and stored in " << fileName << std::endl;
}
