void CreateHistograms()
{
	// settings 

	// number of expected events to be produced
	int nev_bkg    =  1300000; 
	int nev_sgn_h0 =  0.5 * 10000.0; 
	int nev_sgn_hL =  0.3 * 10000.0;
	int nev_sgn_hR =  0.2 * 10000.0; 

	// histogram settings
	int nbins = 25;
	double xmin = -1.0; 
	double xmax =  1.0; 

	// histograms

	// templates
	TH1D* hist_bkg = new TH1D("hist_bkg", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D* hist_sgn_h0 = new TH1D("hist_sgn_h0", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D* hist_sgn_hL = new TH1D("hist_sgn_hL", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D* hist_sgn_hR = new TH1D("hist_sgn_hR", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);

	// sum of all contributions
	TH1D* hist_sum = new TH1D("hist_sum", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);

	// prior histograms 
	TH1D* hist_prior_bkg = new TH1D("hist_prior_bkg", ";x;1/N dN/dx", 100, 3000.0, 10000.0); 

	// data 
	TH1D* hist_data = new TH1D("hist_data", ";cos #theta^{*};dN/dcos #theta^{*}", nbins,xmin, xmax);

	// efficiency parameterization
	TH1D* hist_efficiency_bkg = new TH1D("hist_efficiency_bkg", ";;", nbins, xmin, xmax); 
	TH1D* hist_efficiency_sgn_h0 = new TH1D("hist_efficiency_sgn_h0", ";;", nbins, xmin, xmax); 
	TH1D* hist_efficiency_sgn_hL = new TH1D("hist_efficiency_sgn_hL", ";;", nbins, xmin, xmax); 
	TH1D* hist_efficiency_sgn_hR = new TH1D("hist_efficiency_sgn_hR", ";;", nbins, xmin, xmax); 

	// efficiency uncertainty parameterization
	TH1D* hist_efferror_bkg = new TH1D("hist_efferror_bkg", ";;", nbins, xmin, xmax); 
	TH1D* hist_efferror_sgn_h0 = new TH1D("hist_efferror_sgn_h0", ";;", nbins, xmin, xmax); 
	TH1D* hist_efferror_sgn_hL = new TH1D("hist_efferror_sgn_hL", ";;", nbins, xmin, xmax); 
	TH1D* hist_efferror_sgn_hR = new TH1D("hist_efferror_sgn_hR", ";;", nbins, xmin, xmax); 

	// systematic uncertainty 1 for all four contributions
	TH1D* hist_systerror1_bkg = new TH1D("hist_systerror1_bkg", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror1_sgn_h0 = new TH1D("hist_systerror1_sgn_h0", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror1_sgn_hL = new TH1D("hist_systerror1_sgn_hL", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror1_sgn_hR = new TH1D("hist_systerror1_sgn_hR", ";;", nbins, xmin, xmax); 

	// systematic uncertainty 2 for all four contributions
	TH1D* hist_systerror2_bkg = new TH1D("hist_systerror2_bkg", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror2_sgn_h0 = new TH1D("hist_systerror2_sgn_h0", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror2_sgn_hL = new TH1D("hist_systerror2_sgn_hL", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror2_sgn_hR = new TH1D("hist_systerror2_sgn_hR", ";;", nbins, xmin, xmax); 

	// systematic uncertainty 3 for all four contributions
	TH1D* hist_systerror3_bkg = new TH1D("hist_systerror3_bkg", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror3_sgn_h0 = new TH1D("hist_systerror3_sgn_h0", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror3_sgn_hL = new TH1D("hist_systerror3_sgn_hL", ";;", nbins, xmin, xmax); 
	TH1D* hist_systerror3_sgn_hR = new TH1D("hist_systerror3_sgn_hR", ";;", nbins, xmin, xmax); 

	// fill templates 
	for (int i = 1; i <= nbins; ++i) {
		double x = hist_data->GetBinCenter(i); 
		double y1 = 1; 
		double y2 = 1-x*x; 
		double y3 = (1-x)*(1-x); 
		double y4 = (1+x)*(1+x); 

		hist_bkg->SetBinContent(i, y1); 
		hist_sgn_h0->SetBinContent(i, y2); 
		hist_sgn_hL->SetBinContent(i, y3); 
		hist_sgn_hR->SetBinContent(i, y4); 
	}

	// fill efficiencies and uncertainties
	for (int i = 1; i <= nbins; ++i) {
		double x = hist_data->GetBinCenter(i); 
		double eff = 0.15 + 0.5*(x + 1.)*0.1; 
		double efferror = 0.5*(x + 1.)*0.05; 
		hist_efficiency_bkg->SetBinContent(i, 0.001);
		hist_efferror_bkg->SetBinContent(i, 0.0005);
		hist_efficiency_sgn_h0->SetBinContent(i, eff);
		hist_efferror_sgn_h0->SetBinContent(i, efferror);
		hist_efficiency_sgn_hL->SetBinContent(i, eff);
		hist_efferror_sgn_hL->SetBinContent(i, efferror);
		hist_efficiency_sgn_hR->SetBinContent(i, eff);
		hist_efferror_sgn_hR->SetBinContent(i, efferror);
	}	

	// fill systematic uncertainties for uncertainty 1
	for (int i = 1; i <= nbins; ++i) {
		hist_systerror1_bkg->SetBinContent(i, 0.0001); 
		hist_systerror1_sgn_h0->SetBinContent(i, 0.01); 
		hist_systerror1_sgn_hL->SetBinContent(i, 0.01); 
		hist_systerror1_sgn_hR->SetBinContent(i, 0.01); 
	}	

	// fill systematic uncertainties for uncertainty 2
	for (int i = 1; i <= nbins; ++i) {
		double x = hist_data->GetBinCenter(i); 
		double syst = 0.5*(x + 1.) * 0.005;
		hist_systerror2_bkg->SetBinContent(i, syst/10.); 
		hist_systerror2_sgn_h0->SetBinContent(i, syst); 
		hist_systerror2_sgn_hL->SetBinContent(i, syst); 
		hist_systerror2_sgn_hR->SetBinContent(i, syst); 
	}	

	// fill systematic uncertainties for uncertainty 3
	for (int i = 1; i <= nbins; ++i) {
		double x = hist_data->GetBinCenter(i); 
		double syst = 0.01 - 0.5*(x + 1.) * 0.5*(x + 1.) * 0.005;
		hist_systerror3_bkg->SetBinContent(i, syst/10.); 
		hist_systerror3_sgn_h0->SetBinContent(i, syst+0.01); 
		hist_systerror3_sgn_hL->SetBinContent(i, syst+0.01); 
		hist_systerror3_sgn_hR->SetBinContent(i, syst+0.01); 
	}	

	// fill priors
	for (int i = 1; i <= 100; ++i) {
		double x = hist_prior_bkg->GetBinCenter(i); 
		hist_prior_bkg->SetBinContent(i, TMath::Gaus(x, 6000.0, 800.0)); 
	}

	// scale histograms
	hist_bkg->Scale(1.0/hist_bkg->Integral());
	hist_sgn_h0->Scale(1.0/hist_sgn_h0->Integral());
	hist_sgn_hL->Scale(1.0/hist_sgn_hL->Integral());
	hist_sgn_hR->Scale(1.0/hist_sgn_hR->Integral());
	hist_prior_bkg->Scale(1.0/hist_prior_bkg->Integral());

	// fill data histogram
	gRandom = new TRandom3(1000); 

	for (int i = 1; i <= nbins; ++i) {
		// calculate expectation for each contribution
		double exp1 = nev_bkg * hist_bkg->GetBinContent(i) * hist_efficiency_bkg->GetBinContent(i); 
		double exp2 = nev_sgn_h0 * hist_sgn_h0->GetBinContent(i) * hist_efficiency_sgn_h0->GetBinContent(i); 
		double exp3 = nev_sgn_hL * hist_sgn_hL->GetBinContent(i) * hist_efficiency_sgn_hL->GetBinContent(i); 
		double exp4 = nev_sgn_hR * hist_sgn_hR->GetBinContent(i) * hist_efficiency_sgn_hR->GetBinContent(i); 

		// total expectation
		double exptotal = exp1 + exp2 + exp3 + exp4; 

		// fill data and sum histograms
		hist_data->SetBinContent(i, gRandom->Poisson(exptotal)); 
		hist_sum->SetBinContent(i, exptotal);
	}
	
	// write histograms to file 
 	TFile * file = new TFile("templates.root", "RECREATE");  
 	file -> cd();  

 	hist_bkg->Write();  
 	hist_sgn_h0->Write();  
 	hist_sgn_hL->Write();  
 	hist_sgn_hR->Write();  
	hist_sum->Write();
 	hist_data->Write();  
 	hist_prior_bkg->Write();  
	hist_efficiency_bkg->Write();
	hist_efficiency_sgn_h0->Write();
	hist_efficiency_sgn_hL->Write();
	hist_efficiency_sgn_hR->Write();
	hist_efferror_bkg->Write();
	hist_efferror_sgn_h0->Write();
	hist_efferror_sgn_hL->Write();
	hist_efferror_sgn_hR->Write();
	hist_systerror1_bkg->Write();
	hist_systerror1_sgn_h0->Write();
	hist_systerror1_sgn_hL->Write();
	hist_systerror1_sgn_hR->Write();
	hist_systerror2_bkg->Write();
	hist_systerror2_sgn_h0->Write();
	hist_systerror2_sgn_hL->Write();
	hist_systerror2_sgn_hR->Write();
	hist_systerror3_bkg->Write();
	hist_systerror3_sgn_h0->Write();
	hist_systerror3_sgn_hL->Write();
	hist_systerror3_sgn_hR->Write();

	// print .eps file
	TCanvas * c1 = new TCanvas("c1", "", 1000, 1000); 
	c1->Divide(3, 2); 
	c1->cd(1); 
	hist_bkg->Draw(); 
	c1->cd(2); 
	hist_sgn_h0->Draw(); 
	c1->cd(3); 
	hist_sgn_hL->Draw(); 
	c1->cd(4); 
	hist_sgn_hR->Draw(); 
	c1->cd(5); 
	hist_sum->Draw(); 
	c1->cd(6); 
	hist_data->Draw(); 
	c1->Print("hist.eps"); 

	// print integral
	cout << "Number of events in data: " << hist_data->Integral() << endl;

	// close file 
	file -> Close(); 

	// free memory
	delete c1; 
	delete hist_bkg;  
 	delete hist_sgn_h0;  
 	delete hist_sgn_hL;  
 	delete hist_sgn_hR;  
	delete hist_sum;
 	delete hist_data;  
 	delete hist_prior_bkg;  
	delete hist_efficiency_bkg;
	delete hist_efficiency_sgn_h0;
	delete hist_efficiency_sgn_hL;
	delete hist_efficiency_sgn_hR;
	delete hist_efferror_bkg;
	delete hist_efferror_sgn_h0;
	delete hist_efferror_sgn_hL;
	delete hist_efferror_sgn_hR;
	delete hist_systerror1_bkg;
	delete hist_systerror1_sgn_h0;
	delete hist_systerror1_sgn_hL;
	delete hist_systerror1_sgn_hR;
	delete hist_systerror2_bkg;
	delete hist_systerror2_sgn_h0;
	delete hist_systerror2_sgn_hL;
	delete hist_systerror2_sgn_hR;
	delete hist_systerror3_bkg;
	delete hist_systerror3_sgn_h0;
	delete hist_systerror3_sgn_hL;
	delete hist_systerror3_sgn_hR;
}
