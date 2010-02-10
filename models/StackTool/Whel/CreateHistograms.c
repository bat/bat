void CreateHistograms()
{
	// settings 

	int nev_process1 =  1300000; // W + jets @ 100 pb-1
	int nev_process2 =  0.7 * 7000.0;
	int nev_process3 =  0.3 * 7000.0;
	int nev_process4 =  0.0 * 7000.0; 
	int nev_process5 =  0.0;

	double eff_process1 = 0.001;
	double eff_process2 = 0.20; 
	double eff_process3 = 0.20; 
	double eff_process4 = 0.20; 
	double eff_process5 = 0.0; 

	int nbins = 100;
	double xmin = -1.0; 
	double xmax =  1.0; 

	// histograms 
	TH1D hist_process1("hist_process1", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D hist_process2("hist_process2", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D hist_process3("hist_process3", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D hist_process4("hist_process4", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D hist_process5("hist_process5", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D hist_sum("hist_sum", ";cos #theta^{*};1/N dN/dcos #theta^{*}", nbins, xmin, xmax);
	TH1D hist_prior_process1("hist_prior_process1", ";x;1/N dN/dx", 100, 3000.0, 10000.0); 
	TH1D hist_prior_process2("hist_prior_process2", ";x;1/N dN/dx", 100,    0.0, 10000.0);
	TH1D hist_prior_process3("hist_prior_process3", ";x;1/N dN/dx", 100,    0.0, 10000.0);
	TH1D hist_prior_process4("hist_prior_process4", ";x;1/N dN/dx", 100,    0.0, 10000.0);
	TH1D hist_prior_process5("hist_prior_process5", ";x;1/N dN/dx", 100,    0.0,  1000.0);
	TH1D hist_data("hist_data", ";cos #theta^{*};dN/dcos #theta^{*}", nbins,xmin, xmax);

	// fill templates 
	for (int i = 1; i <= nbins; ++i) {
		double x = hist_data.GetBinCenter(i); 
		
		double y1 = 1; 
		double y2 = 1-x*x; 
		double y3 = (1-x)*(1-x); 
		double y4 = (1+x)*(1+x); 
		double y5 = TMath::Gaus(x, 0.3, 0.1); 
		
		hist_process1.SetBinContent(i, y1); 
		hist_process2.SetBinContent(i, y2); 
		hist_process3.SetBinContent(i, y3); 
		hist_process4.SetBinContent(i, y4); 
		hist_process5.SetBinContent(i, y5); 
	}

	// fill priors
	for (int i = 1; i <= 100; ++i) {
		double x = hist_prior_process1.GetBinCenter(i); 
		hist_prior_process1.SetBinContent(i, TMath::Gaus(x, 6000.0, 800.0)); 
	}
	for (int i = 1; i <= 100; ++i) {
		double x = hist_prior_process2.GetBinCenter(i); 
		hist_prior_process2.SetBinContent(i, TMath::Gaus(x, 5500.0, 100.0)); 
	}
	for (int i = 1; i <= 100; ++i) {
		double x = hist_prior_process3.GetBinCenter(i); 
		hist_prior_process3.SetBinContent(i, 1.0); 
	}
	for (int i = 1; i <= 100; ++i) {
		double x = hist_prior_process4.GetBinCenter(i); 
		hist_prior_process4.SetBinContent(i, TMath::Exp(-x/5000.0)); 
	}
	for (int i = 1; i <= 100; ++i) {
		double x = hist_prior_process5.GetBinCenter(i); 
		hist_prior_process5.SetBinContent(i, TMath::Gaus(x, 100.0, 100.0)); 
	}

	// scale histograms
	hist_process1.Scale(1.0/hist_process1.Integral());
	hist_process2.Scale(1.0/hist_process2.Integral());
	hist_process3.Scale(1.0/hist_process3.Integral());
	hist_process4.Scale(1.0/hist_process4.Integral());
	hist_process5.Scale(1.0/hist_process5.Integral());

	hist_prior_process1.Scale(1.0/hist_prior_process1.Integral());
	hist_prior_process2.Scale(1.0/hist_prior_process2.Integral());
	hist_prior_process3.Scale(1.0/hist_prior_process3.Integral());
	hist_prior_process4.Scale(1.0/hist_prior_process4.Integral());
	hist_prior_process5.Scale(1.0/hist_prior_process5.Integral());

	// fill data histogram
	gRandom = new TRandom3(1000); 

	for (int i = 1; i <= nbins; ++i) {
		double exp1 = nev_process1 * eff_process1 * hist_process1.GetBinContent(i); 
		double exp2 = nev_process2 * eff_process2 * hist_process2.GetBinContent(i); 
		double exp3 = nev_process3 * eff_process3 * hist_process3.GetBinContent(i); 
		double exp4 = nev_process4 * eff_process4 * hist_process4.GetBinContent(i); 
		double exp5 = nev_process5 * eff_process5 * hist_process5.GetBinContent(i); 

		double exptotal = exp1 + exp2 + exp3 + exp4 + exp5; 

		hist_data.SetBinContent(i, gRandom->Poisson(exptotal)); 
		hist_sum.SetBinContent(i, exptotal);
	}
	
	// write histograms to file 
 	TFile * file = new TFile("templates.root", "RECREATE");  
 	file -> cd();  

 	hist_process1.Write();  
 	hist_process2.Write();  
 	hist_process3.Write();  
 	hist_process4.Write();  
 	hist_process5.Write();  
	hist_sum.Write();
 	hist_data.Write();  
 	hist_prior_process1.Write();  
 	hist_prior_process2.Write();  
 	hist_prior_process3.Write();  
 	hist_prior_process4.Write();  
 	hist_prior_process5.Write();  
	
	// print .eps file
	TCanvas * c1 = new TCanvas("c1"); 
	c1->Divide(3, 2); 
	c1->cd(1); 
	hist_process1.Draw(); 
	c1->cd(2); 
	hist_process2.Draw(); 
	c1->cd(3); 
	hist_process3.Draw(); 
	c1->cd(4); 
	hist_process4.Draw(); 
	c1->cd(5); 
	hist_sum.Draw(); 
	c1->cd(6); 
	hist_data.Draw(); 
	c1 -> Print("hist.eps"); 

	cout << hist_data->Integral() << endl;

	// close file 
	file -> Close(); 

}
