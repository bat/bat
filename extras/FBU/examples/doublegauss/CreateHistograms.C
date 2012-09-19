void CreateHistograms()
{
	// settings 


  // convention followed throughout: 1st index: reco, 2nd index: truth

  // migration matrix transposed when it's filled 
  // at the end, since convention in FBU code
  // is 1st index: truth, 2nd index: reco
        

  int nbin = 10;

  TH1F *hist_truth = new TH1F("hist_truth","hist_truth",nbin,0,1);
  
  TF1 *f1 = new TF1("f1", "gaus(0) + gaus(3)", 0., 1.);
  
  f1->SetParameters(500., 0.25., 0.10, 500, 0.75, 0.10);

  hist_truth->FillRandom("f1",100000);

  hist_truth->Draw();

  TMatrix response(nbin,nbin);

  for (int j=1;j<=nbin;j++)
    {
      TF1 *f2 = new TF1("f2", "gaus(0)", 0., 1.);
  
      f2->SetParameters(1.,(j-1)*(1.0/nbin)+0.5/nbin, (1.0/nbin));
      
      TH1F *h_tmp = new TH1F("h_tmp","h_tmp",nbin,0,1);
      
      h_tmp->FillRandom("f2",100000);

      h_tmp->Scale(hist_truth->GetBinContent(j)/h_tmp->Integral());

      for (int i=1;i<=nbin;i++)
	response(i-1,j-1) = h_tmp->GetBinContent(i); 
      // sum over first index will give truth number of events in bin
      // corresponding to efficiency 1
    }

  /*
    for (int j=1;j<=nbin;j++)
     for (int i=1;i<=nbin;i++)
	response(i-1,j-1) = 0; 

   for (int j=1;j<=nbin;j++)
     response(j-1,j-1) = hist_truth->GetBinContent(j);

  */

  for (int i=1;i<=nbin;i++)
    {
      for (int j=1;j<=nbin;j++)
	cout << response(i-1,j-1) << " ";
      cout << endl;
    }
  
    TMatrix normresponse(nbin,nbin);
       
    TH1F *hist_eff = new TH1F("hist_eff","hist_eff",nbin,0,1);
    
    for (int j=1;j<=nbin;j++)
      {
	double eff = 0;
	
	for (int i=1;i<=nbin;i++)
	  {
	    eff+= response(i-1,j-1);
	  }
	
	cout << "efficiency " << j << " = " << eff/hist_truth->GetBinContent(j) << endl;

	hist_eff->SetBinContent(j, eff/hist_truth->GetBinContent(j));

	
	for (int i=1; i<=nbin;i++)
	  normresponse(i-1,j-1) = response(i-1,j-1)/eff;
	
      }


    // check that efficiency = sum_{truth ind} normresponse(reco,truth)


      for (int j=1;j<=nbin;j++)
      {
	
	double tmp = 0;


	for (int i=1;i<=nbin;i++)	  
	  {
	    tmp += normresponse(i-1,j-1);
	  }
	cout << "efficiency for bin " << j << " " << hist_eff->GetBinContent(j) << 
	  " calculated from normalised response " << tmp << endl;
	
      }
      

    TH1F *hist_reco = new TH1F("hist_reco","hist_reco",nbin,0,1);

    for (int i=1;i<=nbin;i++)
      {
	double temp = 0;

	for (int j=1;j<=nbin;j++)
	  {
	    temp += normresponse(i-1,j-1)*hist_truth->GetBinContent(j);
	  }
	
	hist_reco->SetBinContent(i,temp);
      }

    hist_reco->SetLineStyle(2);

    hist_reco->Draw("same");

  
    TH2D *hist_migration = new TH2D("hist_migration","hist_migration",nbin,0,1,nbin,0,1);

    for (int i=1;i<=nbin;i++)
      {
	cout << "reco index " << i << " ";
	for (int j=1;j<=nbin;j++)
	  {
	    double tmp = response(i-1,j-1);
	    cout << tmp << " ";
	    hist_migration->SetBinContent(j,i,tmp);
	  }
	cout << endl;
      }

    TRandom3 *r = new TRandom3();

    r->SetSeed(1000);
    
    TH1F *hist_data = new TH1F("hist_data","hist_data",nbin,0,1);
    
    for (int i=1;i<=nbin;i++)
      {
       	double newentry = r->PoissonD(hist_reco->GetBinContent(i));
	
	// double newentry = hist_reco->GetBinContent(i);

	hist_data->SetBinContent(i,newentry);
	
      }


    for (int i=1;i<=nbin;i++)
      {
	sum = 0;

	for (int j=1;j<=nbin;j++)
	  {
	    sum+ = hist_truth->GetBinContent(j)*hist_eff->GetBinContent(j)*normresponse(i-1,j-1);

	    std::cout << " bin " << i << " truth contrib " << sum << " truth " << hist_truth->GetBinContent(j) << " eff " << hist_eff->GetBinContent(j) << " response " << normresponse(i-1,j-1) << std::endl;
	  }

	cout << " data " << hist_data->GetBinContent(i) << " sum " << sum << " truth " << hist_truth->GetBinContent(i) << endl;
	  
      }
    

    TH1F *hist_bkg = new TH1F("hist_bkg","hist_bkg",nbin,0,1);
    
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
    gStyle->SetPalette(1);
    hist_migration->Draw("COLZ");
    c1->Print("migration.eps");
    
    // print truth, reco and background distributions
    TCanvas* c2 = new TCanvas("c2", "c2", 900, 300);
    c2->Divide(3, 1);
    c2->cd(1);
    hist_reco->Draw();
    c2->cd(2);
    hist_bkg->Draw();
    c2->cd(3);
    hist_data->Draw();
    c2->Print("reco.eps");
    
    // close file
    file->Close();
    
    // clean up
    //	delete hist_migration;
    //	delete file;
}
