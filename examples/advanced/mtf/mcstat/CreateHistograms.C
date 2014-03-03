
void CreateHistograms()
{
   // settings

   // number of expected events to be produced
   int nev_bkg =  300;
   int nev_sgn =  100;

   // efficiencies
   double eff_sgn1 = 1.0;
   double eff_bkg1 = 1.0;

   // histogram settings
   int nbins = 20;
   double xmin = 1000. - 50.;
   double xmax = 1000. + 50.;
	 double mcstats = 200;  // the "MC" statistics of the templates

	 // random numbers
   gRandom = new TRandom3(1000);

   // histograms

   // templates
   TH1D * hist_bkg = new TH1D("hist_bkg", ";E [a.u.];p(E)", nbins, xmin, xmax);
   TH1D * hist_sgn = new TH1D("hist_sgn", ";E [a.u.];p(E)", nbins, xmin, xmax);

   // data
   TH1D * hist_data = new TH1D("hist_data", ";E [a.u.];p(E)", nbins,xmin, xmax);

   // fill templates
   for (int i = 1; i <= nbins; ++i) {
      double x = hist_data->GetBinCenter(i);
      double y1 = 1;
      double y2 = TMath::Gaus(x, 1000, 5);

      hist_bkg->SetBinContent(i, y1);
      hist_sgn->SetBinContent(i, y2);
   }

   // scale histograms
   hist_bkg->Scale(mcstats/hist_bkg->Integral());
   hist_sgn->Scale(mcstats/hist_sgn->Integral());

   // fluctuate templates
	 for (int i = 1; i <= nbins; ++i) {
		 hist_bkg->SetBinContent(i, gRandom->Poisson(hist_bkg->GetBinContent(i)));
		 hist_sgn->SetBinContent(i, gRandom->Poisson(hist_sgn->GetBinContent(i)));
	 }

   // fill data histogram
   for (int i = 1; i <= nbins; ++i) {
      // calculate expectation for each contribution
      double exp_bkg = nev_bkg * hist_bkg->GetBinContent(i) / mcstats;
      double exp_sgn = nev_sgn * hist_sgn->GetBinContent(i) / mcstats;

      // total expectation
      double exptotal1 = eff_bkg1 *exp_bkg + eff_sgn1 *exp_sgn;

      // fill data and sum histograms
      hist_data->SetBinContent(i, gRandom->Poisson(exptotal1));
   }

   // write histograms to file
   TFile * file = TFile::Open("templates.root", "RECREATE");
   file->cd();

   hist_bkg->Write();
   hist_sgn->Write();
   hist_data->Write();

   // print .pdf file
   TCanvas * c1 = new TCanvas("c1", "", 1500, 500);
   c1->Divide(3, 1);
   c1->cd(1);
   hist_bkg->Draw();
   c1->cd(2);
   hist_sgn->Draw();
   c1->cd(3);
   hist_data->Draw();
   c1->Print("hist.pdf");

   // close file
   file->Close();

   // free memory
//   delete c1;
//   delete hist_bkg;
//   delete hist_sgn;
//   delete hist_data;
}
