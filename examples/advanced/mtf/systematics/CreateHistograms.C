void CreateHistograms()
{
   // settings

   // number of expected events to be produced
   int nev_bkg =  300;
   int nev_sgn =  100;

   // efficiencies
   double eff_sgn1 = 1.0;
   double eff_bkg1 = 1.0;

   double eff_sgn2 = 0.5;
   double eff_bkg2 = 1.0;

   // histogram settings
   int nbins = 100;
   double xmin = 2039. - 50.;
   double xmax = 2039. + 50.;

   // histograms

   // templates
   TH1D * hist_bkg = new TH1D("hist_bkg", ";E [keV];p(E)", nbins, xmin, xmax);
   TH1D * hist_sgn = new TH1D("hist_sgn", ";E [keV];p(E)", nbins, xmin, xmax);

   // data
   TH1D * hist_data1 = new TH1D("hist_data1", ";E [keV];p(E)", nbins,xmin, xmax);
   TH1D * hist_data2 = new TH1D("hist_data2", ";E [keV];p(E)", nbins,xmin, xmax);

   // systematics
   TH1D * hist_syst1_bkg_1 = new TH1D("hist_syst1_bkg_1", ";E [keV];variation", nbins, xmin, xmax);
   TH1D * hist_syst1_sgn_1 = new TH1D("hist_syst1_sgn_1", ";E [keV];variation", nbins, xmin, xmax);
   TH1D * hist_syst1_bkg_2 = new TH1D("hist_syst1_bkg_2", ";E [keV];variation", nbins, xmin, xmax);
   TH1D * hist_syst1_sgn_2 = new TH1D("hist_syst1_sgn_2", ";E [keV];variation", nbins, xmin, xmax);

   TH1D * hist_syst2_bkg_1 = new TH1D("hist_syst2_bkg_1", ";E [keV];variation", nbins, xmin, xmax);
   TH1D * hist_syst2_sgn_1 = new TH1D("hist_syst2_sgn_1", ";E [keV];variation", nbins, xmin, xmax);
   TH1D * hist_syst2_bkg_2 = new TH1D("hist_syst2_bkg_2", ";E [keV];variation", nbins, xmin, xmax);
   TH1D * hist_syst2_sgn_2 = new TH1D("hist_syst2_sgn_2", ";E [keV];variation", nbins, xmin, xmax);

   // fill templates
   for (int i = 1; i <= nbins; ++i) {
      double x = hist_data1->GetBinCenter(i);
      double y1 = 1;
      double y2 = TMath::Gaus(x, 2039.0, 5.0);

      hist_bkg->SetBinContent(i, y1);
      hist_sgn->SetBinContent(i, y2);
   }

   // scale histograms
   hist_bkg->Scale(1.0/hist_bkg->Integral());
   hist_sgn->Scale(1.0/hist_sgn->Integral());

   // fill systematic histograms
   for (int i = 1; i <= nbins; ++i) {
      double x = hist_data1->GetBinCenter(i);

      hist_syst1_bkg_1->SetBinContent(i, 0.10);
      hist_syst1_sgn_1->SetBinContent(i, 0.05);
      hist_syst1_bkg_2->SetBinContent(i, 0.10);
      hist_syst1_sgn_2->SetBinContent(i, 0.05);

      hist_syst2_bkg_1->SetBinContent(i, 0.10);
      hist_syst2_sgn_1->SetBinContent(i, 0.20);
      hist_syst2_bkg_2->SetBinContent(i, 0.10);
      hist_syst2_sgn_2->SetBinContent(i, 0.20);
   }

   // fill data histogram
   gRandom = new TRandom3(1000);

   for (int i = 1; i <= nbins; ++i) {
      // calculate expectation for each contribution
      double exp_bkg = nev_bkg * hist_bkg->GetBinContent(i);
      double exp_sgn = nev_sgn * hist_sgn->GetBinContent(i);

      // total expectation
      double exptotal1 = eff_bkg1 *exp_bkg + eff_sgn1 *exp_sgn;
      double exptotal2 = eff_bkg2 *exp_bkg + eff_sgn2 *exp_sgn;

      // fill data and sum histograms
      hist_data1->SetBinContent(i, gRandom->Poisson(exptotal1));
      hist_data2->SetBinContent(i, gRandom->Poisson(exptotal2));
   }

   // write histograms to file
   TFile * file = TFile::Open("templates.root", "RECREATE");
   file -> cd();

   hist_bkg->Write();
   hist_sgn->Write();
   hist_data1->Write();
   hist_data2->Write();
   hist_syst1_bkg_1->Write();
   hist_syst1_sgn_1->Write();
   hist_syst1_bkg_2->Write();
   hist_syst1_sgn_2->Write();
   hist_syst2_bkg_1->Write();
   hist_syst2_sgn_1->Write();
   hist_syst2_bkg_2->Write();
   hist_syst2_sgn_2->Write();

   // print .pdf file
   TCanvas * c1 = new TCanvas("c1", "", 1000, 1000);
   c1->Divide(2, 2);
   c1->cd(1);
   hist_bkg->Draw();
   c1->cd(2);
   hist_sgn->Draw();
   c1->cd(3);
   hist_data1->Draw();
   c1->cd(4);
   hist_data2->Draw();
   c1->Print("hist.pdf");

   // close file
   file -> Close();

   // free memory
   delete c1;
   delete hist_bkg;
   delete hist_sgn;
   delete hist_data1;
   delete hist_data2;
   delete hist_syst1_bkg_1;
   delete hist_syst1_sgn_1;
   delete hist_syst1_bkg_2;
   delete hist_syst1_sgn_2;
   delete hist_syst2_bkg_1;
   delete hist_syst2_sgn_1;
   delete hist_syst2_bkg_2;
   delete hist_syst2_sgn_2;
}
