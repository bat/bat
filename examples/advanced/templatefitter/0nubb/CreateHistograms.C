void CreateHistograms()
{
	// settings 

	// number of expected events to be produced
	int nev_bkg =  300; 
	int nev_sgn =  100;

	// histogram settings
	int nbins = 100;
	double xmin = 2039. - 50.; 
	double xmax = 2039. + 50.; 

	// histograms

	// templates
	TH1D* hist_bkg = new TH1D("hist_bkg", ";E [keV];p(E)", nbins, xmin, xmax);
	TH1D* hist_sgn = new TH1D("hist_sgn", ";E [keV];p(E)", nbins, xmin, xmax);

	// sum of all contributions
	TH1D* hist_sum = new TH1D("hist_sum", ";E [keV];p(E)", nbins, xmin, xmax);

	// data 
	TH1D* hist_data = new TH1D("hist_data", ";E [keV];p(E)", nbins,xmin, xmax);

	// fill templates 
	for (int i = 1; i <= nbins; ++i) {
		double x = hist_data->GetBinCenter(i); 
		double y1 = 1; 
		double y2 = TMath::Gaus(x, 2039.0, 5.0); 

		hist_bkg->SetBinContent(i, y1); 
		hist_sgn->SetBinContent(i, y2); 
	}

	// scale histograms
	hist_bkg->Scale(1.0/hist_bkg->Integral());
	hist_sgn->Scale(1.0/hist_sgn->Integral());

	// fill data histogram
	gRandom = new TRandom3(1000); 

	for (int i = 1; i <= nbins; ++i) {
		// calculate expectation for each contribution
		double exp1 = nev_bkg * hist_bkg->GetBinContent(i); 
		double exp2 = nev_sgn * hist_sgn->GetBinContent(i); 

		// total expectation
		double exptotal = exp1 + exp2; 

		// fill data and sum histograms
		hist_data->SetBinContent(i, gRandom->Poisson(exptotal)); 
		hist_sum->SetBinContent(i, exptotal);
	}
	
	// write histograms to file 
 	TFile * file = new TFile("templates.root", "RECREATE");  
 	file -> cd();  

 	hist_bkg->Write();  
 	hist_sgn->Write();  
	hist_sum->Write();
 	hist_data->Write();  

	// print .eps file
	TCanvas * c1 = new TCanvas("c1", "", 1000, 1000); 
	c1->Divide(2, 2); 
	c1->cd(1); 
	hist_bkg->Draw(); 
	c1->cd(2); 
	hist_sgn->Draw(); 
	c1->cd(3); 
	hist_sum->Draw(); 
	c1->cd(4); 
	hist_data->Draw(); 
	c1->Print("hist.eps"); 

	// print integral
	cout << "Number of events in data: " << hist_data->Integral() << endl;

	// close file 
	file -> Close(); 

	// free memory
	delete c1; 
	delete hist_bkg;  
 	delete hist_sgn;  
	delete hist_sum;
 	delete hist_data;  
}
