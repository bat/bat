int templatefit()
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
		return 1;
	}

	// go back to old directory for memory handling
	f->cd();

	// templates
	TH1D hist_bkg = *((TH1D*) file->Get("hist_bkg"));
	TH1D hist_sgn_h0 = *((TH1D*) file->Get("hist_sgn_h0"));
	TH1D hist_sgn_hL = *((TH1D*) file->Get("hist_sgn_hL"));
	TH1D hist_sgn_hR = *((TH1D*) file->Get("hist_sgn_hR"));

	// efficiency parameterization
	TH1D hist_efficiency_bkg = *((TH1D*) file->Get("hist_efficiency_bkg")); 
	TH1D hist_efficiency_sgn_h0 = *((TH1D*) file->Get("hist_efficiency_sgn_h0")); 
	TH1D hist_efficiency_sgn_hL = *((TH1D*) file->Get("hist_efficiency_sgn_hL")); 
	TH1D hist_efficiency_sgn_hR = *((TH1D*) file->Get("hist_efficiency_sgn_hR")); 

	// efficiency uncertainty parameterization
	TH1D hist_efferror_bkg = *((TH1D*) file->Get("hist_efferror_bkg")); 
	TH1D hist_efferror_sgn_h0 = *((TH1D*) file->Get("hist_efferror_sgn_h0")); 
	TH1D hist_efferror_sgn_hL = *((TH1D*) file->Get("hist_efferror_sgn_hL")); 
	TH1D hist_efferror_sgn_hR = *((TH1D*) file->Get("hist_efferror_sgn_hR")); 

	// systematic uncertainty 1 for all four contributions
	TH1D hist_systerror1_bkg = *((TH1D*) file->Get("hist_systerror1_bkg"));
	TH1D hist_systerror1_sgn_h0 = *((TH1D*) file->Get("hist_systerror1_sgn_h0"));
	TH1D hist_systerror1_sgn_hL = *((TH1D*) file->Get("hist_systerror1_sgn_hL"));
	TH1D hist_systerror1_sgn_hR = *((TH1D*) file->Get("hist_systerror1_sgn_hR"));

	// systematic uncertainty 2 for all four contributions
	TH1D hist_systerror2_bkg = *((TH1D*) file->Get("hist_systerror2_bkg"));
	TH1D hist_systerror2_sgn_h0 = *((TH1D*) file->Get("hist_systerror2_sgn_h0"));
	TH1D hist_systerror2_sgn_hL = *((TH1D*) file->Get("hist_systerror2_sgn_hL"));
	TH1D hist_systerror2_sgn_hR = *((TH1D*) file->Get("hist_systerror2_sgn_hR"));

	// systematic uncertainty 3 for all four contributions
	TH1D hist_systerror3_bkg = *((TH1D*) file->Get("hist_systerror3_bkg"));
	TH1D hist_systerror3_sgn_h0 = *((TH1D*) file->Get("hist_systerror3_sgn_h0"));
	TH1D hist_systerror3_sgn_hL = *((TH1D*) file->Get("hist_systerror3_sgn_hL"));
	TH1D hist_systerror3_sgn_hR = *((TH1D*) file->Get("hist_systerror3_sgn_hR"));

	// data
	TH1D hist_data = *((TH1D*) file->Get("hist_data"));

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
	// create model
	// ----------------------------------------------------

	// create new TemplateModel object
	BCTemplateFitter * model = new BCTemplateFitter("model");

	// set MCMC options 
	model->MCMCSetNLag(10); 
	model->MCMCSetNChains(10); 
	//	model->MCMCSetNIterationsRun(1000000); // use more iterations
	//	model->MCMCSetNIterationsRun(10000000); // use even more iterations

	// set template fitting options
	model->SetFlagPhysicalLimits(true);

	// set data histogram
	model->SetData(hist_data);

	// add template histograms
	model->AddTemplate(hist_bkg, "background", 1300000.-10000., 1300000.+10000.);
	model->AddTemplate(hist_sgn_h0, "signal (h= 0)", 0.0, 15000.0);
	model->AddTemplate(hist_sgn_hL, "signal (h=-1)", 0.0, 15000.0);
	model->AddTemplate(hist_sgn_hR, "signal (h=+1)", 0.0, 15000.0);

	// set efficiencies and uncertainties
	model->SetTemplateEfficiency("background", hist_efficiency_bkg, hist_efferror_bkg);
	model->SetTemplateEfficiency("signal (h= 0)", hist_efficiency_sgn_h0, hist_efferror_sgn_h0);
	model->SetTemplateEfficiency("signal (h=-1)", hist_efficiency_sgn_hL, hist_efferror_sgn_hL);
	model->SetTemplateEfficiency("signal (h=+1)", hist_efficiency_sgn_hR, hist_efferror_sgn_hR);

	// set prior on background
	model->SetTemplatePrior("background", 1300000., 2000.);

	// set prior on signal contribution
	std::vector<int> indices; 
	indices.push_back(1);
	indices.push_back(2);
	indices.push_back(3);
	model->ConstrainSum(indices, 10000., 100.); 

	// set up printing of fractions
	model->CalculateRatio(3, indices, 0.0, 1.0); 

	// add systematic uncertainty 1
	model->AddSystError("syst. error 1", "gauss");
	model->SetTemplateSystError("syst. error 1", "background", hist_systerror1_bkg);
	model->SetTemplateSystError("syst. error 1", "signal (h= 0)", hist_systerror1_sgn_h0);
	model->SetTemplateSystError("syst. error 1", "signal (h=-1)", hist_systerror1_sgn_hL);
	model->SetTemplateSystError("syst. error 1", "signal (h=+1)", hist_systerror1_sgn_hR);

	// add systematic uncertainty 2
	model->AddSystError("syst. error 2", "gauss");
	model->SetTemplateSystError("syst. error 2", "background", hist_systerror2_bkg);
	model->SetTemplateSystError("syst. error 2", "signal (h= 0)", hist_systerror2_sgn_h0);
	model->SetTemplateSystError("syst. error 2", "signal (h=-1)", hist_systerror2_sgn_hL);
	model->SetTemplateSystError("syst. error 2", "signal (h=+1)", hist_systerror2_sgn_hR);

	// add systematic uncertainty 3
	model->AddSystError("syst. error 3", "gauss");
	model->SetTemplateSystError("syst. error 3", "background", hist_systerror3_bkg);
	model->SetTemplateSystError("syst. error 3", "signal (h= 0)", hist_systerror3_sgn_h0);
	model->SetTemplateSystError("syst. error 3", "signal (h=-1)", hist_systerror3_sgn_hL);
	model->SetTemplateSystError("syst. error 3", "signal (h=+1)", hist_systerror3_sgn_hR);

	// ----------------------------------------------------
	// perform analysis
	// ----------------------------------------------------

	// initialize model
	model->Initialize();

	// run MCMC
	model->MarginalizeAll(); 
	
	// find global mode
	model->FindMode();

	// ----------------------------------------------------
	// print
	// ----------------------------------------------------

	// print chi2 to screen 
	std::cout << " Chi2 / ndf = " << model->CalculateChi2() << " / " << model->GetNDF() << " (" << model->CalculateChi2Prob() << ")" << std::endl;
	
	// print KS test results to screen
	std::cout << " KS probability = " << model->CalculateKSProb() << std::endl;

	// create summary tool
	BCSummaryTool* st = new BCSummaryTool(model); 

	// print data
 	TCanvas c1("c1");
	c1.cd();
	hist_data.Draw();
	c1.Print("data.ps");

	// print results
	model->PrintAllMarginalized("model_marginalized.eps"); 
	model->PrintStack("model_stack.eps");
	model->PrintRatios("model_fraction.ps", 0, -95.);
	model->PrintResults("model_results.txt"); 

	// print templates 
	model->PrintTemplate("background", "background.eps");
	model->PrintTemplate("signal (h= 0)", "signal_h0.eps");
	model->PrintTemplate("signal (h=-1)", "signal_hL.eps");
	model->PrintTemplate("signal (h=+1)", "signal_hR.eps");

	// print summary plots
	st->PrintParameterPlot("model_parameters.eps"); 
	st->PrintCorrelationPlot("model_correlation.eps"); 
	st->PrintKnowlegdeUpdatePlot("model_update.eps"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete model
	delete model;

	// delete summary tool
	delete st; 

	// no error
	return 0;
}

