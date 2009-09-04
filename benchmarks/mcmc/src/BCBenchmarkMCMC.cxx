/*
 * Copyright (C) 2009, 
 * Daniel Kollar, Kevin Kroeninger and Jing Liu.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

//=============================================================================

#include <TFile.h>
#include <TTree.h>

#include "BAT/BCLog.h"

#include "BCBenchmarkMCMC.h"

//=============================================================================

BCBenchmarkMCMC::BCBenchmarkMCMC(
		TF1* testFunction,
		const char* outputFile,
		const char* modelName)
:BCModel(modelName), BCModelOutput(), fNbinx(100)
{
	BCLog::Out(BCLog::summary,BCLog::summary," setup test function ...");
	if (testFunction==NULL) {
		BCLog::Out(BCLog::error,BCLog::error,
				" The test function doesn't exist!");
		abort();
	}

	fTestFunction = testFunction;
	fXmin=fTestFunction->GetXmin(); 
	fXmax=fTestFunction->GetXmax();


	BCLog::Out(BCLog::summary,BCLog::summary," setup fitting function ...");
	fFitFunction = new TF1("fFitFunction",
			fTestFunction->GetExpFormula(),fXmin,fXmax);
	for (int i=1; i<fTestFunction->GetNpar(); i++)
		fFitFunction->FixParameter(i,fTestFunction->GetParameter(i));


	BCLog::Out(BCLog::summary,BCLog::summary," add parameter x ...");
	this->AddParameter("x",fXmin,fXmax);


	BCLog::Out(BCLog::summary,BCLog::summary," setup trees for chains ...");
	SetModel(this);
	BCModelOutput::WriteMarkovChain(true);
	SetFile(outputFile);
	

	BCLog::Out(BCLog::summary,BCLog::summary," book histograms for analysis ...");
	for (int i=0; i<fMaxChains; i++) {
		for (int j=1; j<fMaxLags; j++)
			fHistXLags[i][j] = new TH1F(Form("fHistXChain%dLag%d",i,j),
					Form("lag of %d",j), fNbinx,fXmin,fXmax);

		for (int j=1; j<fMax10thOfIters; j++)
			fHistXIter[i][j] = new TH1F(Form("fHistXChain%dIter%d0",i,j),
					Form("after %d%% of total iterations",j*10),
					fNbinx,fXmin,fXmax);

		fHChi2vsLags[i] = new TH1F(Form("HChi2vsLagsChain%d",i),"",
				fMaxLags-1,1,fMaxLags);
		fHChi2vsIter[i] = new TH1F(Form("HChi2vsIterChain%d",i),"",
				fMax10thOfIters-1,1,fMax10thOfIters);
	}
}

//=============================================================================

BCBenchmarkMCMC::~BCBenchmarkMCMC()
{
	if (fFitFunction) delete fFitFunction;

	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		for (int j=1; j<fMaxLags; j++)
			if (fHistXLags[i][j]) delete fHistXLags[i][j];

		for (int j=1; j<fMax10thOfIters; j++)
			if (fHistXIter[i][j]) delete fHistXIter[i][j];

		if (fHChi2vsLags[i]) delete fHChi2vsLags[i];
		if (fHChi2vsIter[i]) delete fHChi2vsIter[i];
	}
}

//=============================================================================

void BCBenchmarkMCMC::ProcessMCTrees()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) ProcessMCTree(i);
}

//=============================================================================

void BCBenchmarkMCMC::ProcessMCTree(int chainID)
{
	TTree *chain = this->MCMCGetMarkovChainTree(chainID);

	Double_t fMean,fVariance,fSkewness;
	chain->Branch("fMean",&fMean,"fMean/D");
	chain->Branch("fVariance",&fVariance,"fVariance/D");
	chain->Branch("fSkewness",&fSkewness,"fSkewness/D");

	Double_t par0;
	chain->SetBranchAddress("fParameter0",&par0);

	Double_t sum1=0, sum2=0, sum3=0;
	int Niters = chain->GetEntries();
	for (int i=0; i<Niters; i++) {
		chain->GetEntry(i);

		sum1 += par0;
		fMean = sum1/Double_t(i+1);

		sum2 += (par0-fMean)*(par0-fMean);
		fVariance = sqrt(sum2/Double_t(i+1));

		sum3 += (par0-fMean)*(par0-fMean)*(par0-fMean);
		if (fVariance!=0) 
			fSkewness = (sum3/Double_t(i+1))/(fVariance*fVariance*fVariance);

		chain->Fill();
		
		for (int j=1; j<fMaxLags; j++) {
			if (i%j==0) fHistXLags[chainID][j]->Fill(par0);
		}

		for (int j=1; j<fMax10thOfIters; j++) {
			if (i<=j*Niters/10) 
				fHistXIter[chainID][j]->Fill(par0);
		}
	}
}

//=============================================================================

void BCBenchmarkMCMC::PerformLagsTest()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		Chi2vsLagsOfChain(i);
	}
}

//=============================================================================

void BCBenchmarkMCMC::Chi2vsLagsOfChain(int chainID)
{
	for (int i=1; i<fMaxLags; i++) {
		fHistXLags[chainID][i]->Fit(fFitFunction,"bq");
		fHChi2vsLags[chainID]->SetBinContent(i,
				fFitFunction->GetChisquare()/fFitFunction->GetNDF());
	}
}

//=============================================================================

void BCBenchmarkMCMC::PerformIterationsTest()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		Chi2vsIterOfChain(i);
	}
}

//=============================================================================

void BCBenchmarkMCMC::Chi2vsIterOfChain(int chainID)
{
	for (int i=1; i<fMax10thOfIters; i++) {
		fHistXIter[chainID][i]->Fit(fFitFunction,"bq");
		fHChi2vsIter[chainID]->SetBinContent(i,
				fFitFunction->GetChisquare()/fFitFunction->GetNDF());
	}
}

//=============================================================================

void BCBenchmarkMCMC::WriteResults()
{
	TFile *file = this->GetFile();

	file->cd();
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		file->mkdir(Form("HistsForMCTree%d",i),
				Form("Histograms for Markov chain tree %d",i));
		file->cd(Form("HistsForMCTree%d",i));

		for (int j=1; j<fMaxLags; j++) {
			fHistXLags[i][j]->SetXTitle("x");
			fHistXLags[i][j]->SetYTitle("Entries");
			fHistXLags[i][j]->Write();
		}

		for (int j=1; j<fMax10thOfIters; j++) {
			fHistXIter[i][j]->SetXTitle("x");
			fHistXIter[i][j]->SetYTitle("Entries");
			fHistXIter[i][j]->Write();
		}

		fHChi2vsLags[i]->SetXTitle("Lags");
		fHChi2vsLags[i]->SetYTitle("#chi^{2}/NDF");
		fHChi2vsLags[i]->Write();

		fHChi2vsIter[i]->SetXTitle("10% of total iteration");
		fHChi2vsIter[i]->SetYTitle("#chi^{2}/NDF");
		fHChi2vsIter[i]->Write();
	}
	file->cd();

	this->Close();
}
