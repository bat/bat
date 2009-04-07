// root -q -x -l PrepareWorkspace_Poisson.C

void PrepareWorkspace_Poisson()
{

	// use a dummy observable for this counting-type analysis
	RooRealVar* x = new RooRealVar("x","dummy discriminating variable for event count",0,0,1);
	RooArgSet* fObservables = new RooArgSet(*x);

	// signal and background PDF are flat (does not depend on the actual value of the observable)
	RooPolynomial* sig = new RooPolynomial("sig","sig pdf",*x,RooFit::RooConst(0));
	RooPolynomial* bkg = new RooPolynomial("bkg","bkg pdf",*x,RooFit::RooConst(0));

	// signal yield (with a range used for the integration and posterior plotting)
	RooRealVar* nsig = new RooRealVar("S","sig yield and range",0,60);

	// background yield (with a default value)
	RooRealVar* nbkg = new RooRealVar("B","bkg yield and range",10);

	// total model for the signal+background PDF
	RooAddPdf* fModel = new RooAddPdf("fModel","pdf for sig+bkg model",RooArgList(*sig,*bkg),RooArgList(*nsig,*nbkg));

	// flat prior on the signal yield
	RooPolynomial* fPrior = new RooPolynomial("fPrior","flat prior",*nsig,RooFit::RooConst(0));

	// list of parameters of interest
	RooArgList* fParams = new RooArgList(*nsig);

	// create a toy data sample assuming the measurement will be exactly at the expected background value
	// only one entry, with a proper weight, is created to fasten the calculation in this case
	double fNObserved = nbkg->getVal();
	std::cout << "Creating a RooFit dataset of one event with a weight of " << fNObserved << " events\n";
	RooDataSet* fData = new RooDataSet();
	fData->add(RooArgSet(*x));
	RooRealVar* weightVar = new RooRealVar("weightVar","weight variable",fNObserved);
	fData->addColumn(*weightVar);
	fData->setWeightVar(*weightVar);

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
