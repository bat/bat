void CreateHistograms()
{
	// settings

	// resolution
	double resolution = 10.;

	// histogram settings
	int nbins_reco = 10;
	double xmin_reco = 20.;
	double xmax_reco = 250.;

	int nbins_truth = 20;
	double xmin_truth = 0.;
	double xmax_truth = 200.;

	int nevents_migration = 100000;

	int nevents_sgn = 1000;
	int nevents_bkg = 1000;

	// initialize TRandom3 object
	gRandom = new TRandom3(1000);

	// histograms

	// distribution of events used for creation of migration matrix
	TH1D* hist_migration_truth = new TH1D("hist_migration_truth", ";x;N_{ev}",
																			 nbins_truth, xmin_truth, xmax_truth);
	hist_migration_truth->SetStats(kFALSE);

	TH1D* hist_migration_reco = new TH1D("hist_migration_reco", ";x;",
																			nbins_reco, xmin_reco, xmax_reco);
	hist_migration_reco->SetStats(kFALSE);

	// migration matrix
	TH2D* hist_migration = new TH2D("hist_migration", ";truth;reco",
																 nbins_truth, xmin_truth, xmax_truth,
																 nbins_reco, xmin_reco, xmax_reco);
	hist_migration->SetStats(kFALSE);

	// truth distribution
	TH1D* hist_truth = new TH1D("hist_truth", ";x;N_{ev}",
															nbins_truth, xmin_truth, xmax_truth);
	hist_truth->SetStats(kFALSE);

	// reco distribution
	TH1D* hist_reco = new TH1D("hist_reco", ";x;N_{ev}",
														 nbins_reco, xmin_reco, xmax_reco);
	hist_reco->SetStats(kFALSE);

	// bkg distribution
	TH1D* hist_bkg = new TH1D("hist_bkg", ";x;N_{ev}",
														nbins_reco, xmin_reco, xmax_reco);
	hist_bkg->SetStats(kFALSE);

	// data distribution
	TH1D* hist_data = new TH1D("hist_data", ";x;N_{ev}",
														nbins_reco, xmin_reco, xmax_reco);
	hist_data->SetStats(kFALSE);

	// fill migration matrix
	for (int i = 0; i < nevents_migration; ++i) {
		double truth = gRandom->Landau(xmin_truth+0.2*(xmax_truth-xmin_truth), 0.1*(xmax_truth-xmin_truth));
		double reco  = gRandom->Gaus(truth, resolution);

		hist_truth->Fill(truth);
		hist_reco->Fill(reco);
		if ((truth > xmin_truth) && (truth < xmax_truth) && (reco > xmin_reco) && (reco < xmax_reco))
			hist_migration->Fill(truth, reco);
	}

	// fill truth and reco distribution
	for (int i = 0; i < nevents_sgn; ++i) {
		double truth = gRandom->Landau(xmin_truth+0.2*(xmax_truth-xmin_truth), 0.1*(xmax_truth-xmin_truth));
		double reco  = gRandom->Gaus(truth, resolution);

		hist_data->Fill(reco);
	}

	// fill background histograms
	for (int i = 0; i < nevents_bkg; ++i) {
		double reco = xmin_reco + (xmax_reco-xmin_reco)*gRandom->Uniform();

		hist_bkg->Fill(reco);
		hist_data->Fill(reco);
	}

	// write histograms to file
	TFile * file = new TFile("histograms.root", "RECREATE");
 	file -> cd();

 	hist_migration->Write();
	hist_truth->Write();
	hist_reco->Write();
	hist_bkg->Write();
	hist_data->Write();

	// print migration matrix
	TCanvas* c1 = new TCanvas("c1", "c1", 900, 300);
	c1->Divide(3, 1);
	c1->cd(1);
	hist_truth->Draw();
	c1->cd(2);
	hist_reco->Draw();
	c1->cd(3);
	hist_migration->Draw("COLZ");
	c1->Print("migration.pdf");

	// print truth, reco and background distributions
	TCanvas* c2 = new TCanvas("c2", "c2", 900, 300);
	c2->Divide(3, 1);
	c2->cd(1);
	hist_reco->Draw();
	c2->cd(2);
	hist_bkg->Draw();
	c2->cd(3);
	hist_data->Draw();
	c2->Print("reco.pdf");

	// close file
	file->Close();

	// clean up
	//	delete hist_migration;
	//	delete file;
}
