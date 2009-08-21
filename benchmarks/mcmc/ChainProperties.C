// a simple ROOT script to study properties of Markov Chain
// Jing Liu 2009-08-21
//
{
	//Reset ROOT and connect tree file
	gROOT->Reset();
	TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcmcgaus1d.root");
	if (!f) {
		f = new TFile("mcmcgaus1d.root");
	}
	TTree *MarkovChainTree_0 = (TTree*)gDirectory->Get("MarkovChainTree_4");

	//Declaration of leaves types
	Int_t           fIteration;
	Int_t           fNParameters;
	Double_t        fLogLikelihood;
	Double_t        fParameter0;

	// Set branch addresses.
	MarkovChainTree_0->SetBranchAddress("fIteration",&fIteration);
	MarkovChainTree_0->SetBranchAddress("fNParameters",&fNParameters);
	MarkovChainTree_0->SetBranchAddress("fLogLikelihood",&fLogLikelihood);
	MarkovChainTree_0->SetBranchAddress("fParameter0",&fParameter0);

	const Long64_t nentries = MarkovChainTree_0->GetEntries();
	//const Long64_t nentries = 5000;
	cout<<"Nentries = "<<nentries<<endl;
	Double_t iter[nentries]={0};
	Double_t mean[nentries]={0};
	Double_t variance[nentries]={0};
	Double_t skewness[nentries]={0};

	Double_t sum=0, sum2=0, sum3=0;

	// book histograms
	const Int_t Nl=20;
	TH1D *distr[Nl];
    for (Int_t j=0; j<Nl; j++) {
		distr[j] = new TH1D(Form("distr%d",j+1),Form("lag=%d",j+1), 40,-4,16);
	}

	for (Long64_t i=0; i<nentries; i++) {
		MarkovChainTree_0->GetEntry(i);
		if (i%10000==0) 
			cout<<"now event = "<<i<<endl;

		iter[i]=i+1;

		sum += fParameter0;
		mean[i] = sum/float(i+1);

		sum2 += (fParameter0-mean[i])*(fParameter0-mean[i]);
		variance[i] = sqrt(sum2/float(i+1));

		sum3 += (fParameter0-mean[i])*(fParameter0-mean[i])*(fParameter0-mean[i]);
		if (variance[i]!=0) skewness[i] = (sum3/float(i+1))/(variance[i]*variance[i]*variance[i]);

		distr[0]->Fill(fParameter0);
		for (Int_t j=2; j<=Nl; j++) {
			if (i%j==0) distr[j-1]->Fill(fParameter0);
		}
	}

//	TGraph *gmean = new TGraph(nentries,iter,mean);
//	TGraph *gvari = new TGraph(nentries,iter,variance);
//	TGraph *gskew = new TGraph(nentries,iter,skewness);
//
//	gmean->SetMarkerColor(kBlue);
//	gvari->SetMarkerColor(kRed);
//	gskew->SetMarkerColor(kGreen);
//
//	TMultiGraph *mg = new TMultiGraph();
//	mg->Add(gmean,"p");
//	mg->Add(gvari,"p");
//	mg->Add(gskew,"p");

	TF1 *f1 = new TF1("f1","gaus",-4,16);
	f1->FixParameter(1,3.0);
	f1->FixParameter(2,4.0);

	TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
//	c1->Divide(2,3,0.00001,0.00001);
//	
//	c1->cd(1);
////	mg->Draw("a");
//	distrLag1->Fit(f1,"b","e");
//
//	c1->cd(2);
//	distrLag10->Fit(f1,"b","e");
//	c1->cd(3);
//	distrLag20->Fit(f1,"b","e");
//	c1->cd(4);
//	distrLag40->Fit(f1,"b","e");
//	c1->cd(5);
//	distrLag60->Fit(f1,"b","e");
//	c1->cd(6);
//	distrLag80->Fit(f1,"b","e");

	Float_t chi[Nl] = {0};
	Float_t lag[Nl] = {0};
	for (Int_t k=0; k<Nl; k++) {
		distr[k]->Fit(f1,"b","e");
		chi[k]=f1->GetChisquare()/39.;
		lag[k]=k+1;
	}
	TGraph *glag = new TGraph(Nl,lag,chi);
	glag->SetMarkerStyle(20);
	glag->SetTitle("");
	glag->Draw("ap");
	glag->GetXaxis()->SetTitle("Lag");
	glag->GetYaxis()->SetTitle("#chi^{2}/NDF");
	c1->Update();
}
