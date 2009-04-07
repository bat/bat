// root -q -x -l PrepareWorkspace_GaussOverFlat.C

void PrepareWorkspace_GaussOverFlat()
{

	// use a observable for this shape-based analysis
	RooRealVar* x = new RooRealVar("x","discriminating variable for event count",100,50,150);
	RooArgSet* fObservables = new RooArgSet(*x);

	// signal and background PDFs
	RooGaussian* sig = new RooGaussian("sig","sig pdf",*x,RooFit::RooConst(100),RooFit::RooConst(10));
	RooPolynomial* bkg = new RooPolynomial("bkg","bkg pdf",*x,RooFit::RooConst(0));

	// signal yield (with a range used for the integration and posterior plotting)
	RooRealVar* nsig = new RooRealVar("S","sig yield and range",0,60);

	// background yield (with a default value)
	RooRealVar* nbkg = new RooRealVar("B","bkg yield and range",10,0,200);

	// total model for the signal+background PDF
	RooAddPdf* fModel = new RooAddPdf("fModel","pdf for sig+bkg model",RooArgList(*sig,*bkg),RooArgList(*nsig,*nbkg));

	// flat prior on the signal yield
	RooPolynomial* fPrior = new RooPolynomial("fPrior","flat prior",*nsig,RooFit::RooConst(0));

	// list of parameters of interest
	RooArgList* fParams = new RooArgList(*nsig);

	// generate a toy data sample assuming the background-only hypothesis (without Poisson fluctuation on the number of observed events but with fluctuations in the distribution of the events over the 'x' spectrum)
	double fNObserved = nbkg->getVal();
	RooDataSet* fData = fModel->generate(*fObservables,fNObserved);
	std::cout << "Generated a RooFit dataset of " << fData->numEntries() << " events\n";

	// set the names of the objects (hard coded in the interface so far)
	fData->SetName("fData");
	fObservables->setName("fObservables");
	fParams->setName("fParams");
	fPrior->SetName("fPrior");

	// store the pdf and prior informations in a RooWorkspace
	RooWorkspace bat_ws("bat_ws");
	bat_ws.import(*fData);
	bat_ws.import(*fModel);
	bat_ws.import(*fPrior);

	// store the workspace in a ROOT file
	TString fileName = "roodata.root";
	TFile file(fileName,"RECREATE");
	file.cd();
	bat_ws.Write();
	// temporarily needed (patch)
	fObservables->Write();
	fParams->Write();
	// end of patch
	file.Write();
	file.Close();

	std::cout << "\nRooFit model initialized and stored in " << fileName << std::endl;
}
