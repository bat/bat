
void CreateHistograms()
{
   // settings

   // number of expected events to be produced
   int nev_bkg1 =  800;
   int nev_bkg2 =  500;
   int nev_sgn =  200;

   // efficiencies
   double eff_sgn1 = 0.5;
   double eff_bkg1 = 1.0;

   double eff_sgn2 = 1.0;
   double eff_bkg2 = 1.0;

   // histogram settings
   int nbins = 100;
   double xmin = 1000. - 50.;
   double xmax = 1000. + 50.;

   // histograms

   // templates
   TH1D * hist_bkg1 = new TH1D("hist_bkg1", ";E [a.u.];p(E)", nbins, xmin, xmax);
   TH1D * hist_bkg2 = new TH1D("hist_bkg2", ";E [a.u.];p(E)", nbins, xmin, xmax);
   TH1D * hist_sgn1 = new TH1D("hist_sgn1", ";E [a.u.];p(E)", nbins, xmin, xmax);
   TH1D * hist_sgn2 = new TH1D("hist_sgn2", ";E [a.u.];p(E)", nbins, xmin, xmax);

   // data
   TH1D * hist_data1 = new TH1D("hist_data1", ";E [a.u.];p(E)", nbins,xmin, xmax);
   TH1D * hist_data2 = new TH1D("hist_data2", ";E [a.u.];p(E)", nbins,xmin, xmax);

   // fill templates
   for (int i = 1; i <= nbins; ++i) {
      double x = hist_data1->GetBinCenter(i);
      double ybkg1 = exp(-0.01*x);
      double ybkg2 = 1+0.001*x;
      double ysgn1 = TMath::Gaus(x, 1000.,  5.);
      double ysgn2 = TMath::Gaus(x, 1000., 10.);

      hist_bkg1->SetBinContent(i, ybkg1);
      hist_bkg2->SetBinContent(i, ybkg2);
      hist_sgn1->SetBinContent(i, ysgn1);
      hist_sgn2->SetBinContent(i, ysgn2);
   }

   // scale histograms
   hist_bkg1->Scale(1./hist_bkg1->Integral());
   hist_bkg2->Scale(1./hist_bkg2->Integral());
   hist_sgn1->Scale(1./hist_sgn1->Integral());
   hist_sgn2->Scale(1./hist_sgn2->Integral());

   // fill data histogram
   gRandom = new TRandom3(1000);

   for (int i = 1; i <= nbins; ++i) {
      // calculate expectation for each contribution
      double exp_bkg1 = nev_bkg1 * hist_bkg1->GetBinContent(i);
      double exp_bkg2 = nev_bkg2 * hist_bkg2->GetBinContent(i);
      double exp_sgn1 = nev_sgn * hist_sgn1->GetBinContent(i);
      double exp_sgn2 = nev_sgn * hist_sgn2->GetBinContent(i);

      // total expectation
      double exptotal1 = eff_bkg1 *exp_bkg1 + eff_sgn1 *exp_sgn1;
      double exptotal2 = eff_bkg2 *exp_bkg2 + eff_sgn2 *exp_sgn2;

      // fill data and sum histograms
      hist_data1->SetBinContent(i, gRandom->Poisson(exptotal1));
      hist_data2->SetBinContent(i, gRandom->Poisson(exptotal2));
   }

   // write histograms to file
   TFile * file = TFile::Open("templates.root", "RECREATE");
   file->cd();

   hist_bkg1->Write();
   hist_bkg2->Write();
   hist_sgn1->Write();
   hist_sgn2->Write();
   hist_data1->Write();
   hist_data2->Write();

   // print .pdf file
   TCanvas * c1 = new TCanvas("c1", "", 1000, 1000);
   c1->Divide(2, 2);
   c1->cd(1);
   hist_bkg1->Draw();
   hist_bkg1->GetYaxis()->SetRangeUser(0., 1.1*hist_bkg1->GetMaximum());
   c1->cd(2);
   hist_sgn1->Draw();
   hist_sgn1->GetYaxis()->SetRangeUser(0., 1.1*hist_sgn1->GetMaximum());
   c1->cd(3);
   hist_bkg2->Draw();
   hist_bkg2->GetYaxis()->SetRangeUser(0., 1.1*hist_bkg2->GetMaximum());
   c1->cd(4);
   hist_sgn2->Draw();
   hist_sgn2->GetYaxis()->SetRangeUser(0., 1.1*hist_sgn2->GetMaximum());
   c1->Print("hist1.pdf");

   TCanvas * c2 = new TCanvas("c2", "", 500, 500);
   c2->Divide(2, 1);
   c2->cd(1);
   hist_data1->Draw();
   c2->cd(2);
   hist_data2->Draw();
   c2->Print("hist2.pdf");

   std::cout << "Number of data events in channel 1: " << hist_data1->Integral() << std::endl;
   std::cout << "Number of data events in channel 2: " << hist_data2->Integral() << std::endl;

   // close file
   file->Close();

   // free memory
//   delete c1;
//   delete c2;
//   delete hist_bkg1;
//   delete hist_bkg2;
//   delete hist_sgn1;
//   delete hist_sgn2;
//   delete hist_data1;
//   delete hist_data2;
}
