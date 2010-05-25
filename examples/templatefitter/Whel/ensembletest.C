int ensembletest()
{
	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	// remember old directory
	TDirectory* f = gDirectory;

	// open file
	TFile * file = new TFile("templates.root", "READ");

	// check if file is open
	if (!file->IsOpen()) {
		std::cout << "Could not open file. Exit." << std::endl;
		std::cerr << "Run macro CreateHistograms.C in Root to create the file." << std::endl;
		return 1;
	}

	// go back to old directory for memory handling
	f->cd();

	// get histograms
	TH1D hist_process1 = *((TH1D*) file->Get("hist_bkg"));
	TH1D hist_process2 = *((TH1D*) file->Get("hist_sgn_h0"));
	TH1D hist_process3 = *((TH1D*) file->Get("hist_sgn_hL"));
	TH1D hist_process4 = *((TH1D*) file->Get("hist_sgn_hR"));
	TH1D hist_sum = *((TH1D*) file->Get("hist_sum"));

	// close file
	file->Close();

	// delete file
	delete file;

	// ----------------------------------------------------
	// configure BAT
	// ----------------------------------------------------

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

 	// ----------------------------------------------------
	// Create new model
	// ----------------------------------------------------

	// create new BCTemplateFitter object
	BCTemplateFitter * model = new BCTemplateFitter("model");

	// set options
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000); 
	model->SetFlagPhysicalLimits(true);

	// set data histogram
	model->SetData(hist_sum);

	// add template histograms
	model->AddTemplate(hist_process1, "background",    0.0, 2000000.0);
	model->AddTemplate(hist_process2, "signal (h= 0)", 0.0,    8000.0);
	model->AddTemplate(hist_process3, "signal (h=-1)", 0.0,    8000.0);
	model->AddTemplate(hist_process4, "signal (h=+1)", 0.0,    8000.0);

	// set efficiencies
	model->SetTemplateEfficiency("background", 0.001, 0.0005);
	model->SetTemplateEfficiency("signal (h= 0)", 0.20, 0.05);
	model->SetTemplateEfficiency("signal (h=-1)", 0.20, 0.05);
	model->SetTemplateEfficiency("signal (h=+1)", 0.20, 0.05);

	// set priors 
	model->SetTemplatePrior("background", 1300000.0, 2000.0);

	// set constraints
	std::vector <int> indices; 
	indices.push_back(1);
	indices.push_back(2);
	indices.push_back(3);
	model->ConstrainSum(indices, 7000.0, 100); 

	// set up printing of fractions
	model->CalculateRatio(1, indices); 
	model->CalculateRatio(2, indices); 
	model->CalculateRatio(3, indices); 

	// ----------------------------------------------------
	// create ensemble test tool
	// ----------------------------------------------------

	// create ensemble test tool
	BCTemplateEnsembleTest* tet = new BCTemplateEnsembleTest(); 

	// stetings
	tet->SetTemplateFitter(model);
	tet->SetEnsembleTemplate(hist_sum);
	tet->SetNEnsembles(1000); 
	tet->SetEnsembleExpectation(2700); 
	tet->SetFlagMCMC(true); 

	// perform ensemble tests
	tet->PerformEnsembleTest(); 

	// write results to file
	tet->Write("ensemble.root"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete model
 	delete model;

	return 0;
}

