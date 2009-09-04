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

#include "BCBenchmarkMCMC2D.h"

//=============================================================================

BCBenchmarkMCMC2D::BCBenchmarkMCMC2D(
		TF2* testFunction,
		const char* outputFile,
		const char* modelName)
:BCModel(modelName), BCModelOutput(), fNbinx(100), fNbiny(100)
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
	fYmin=fTestFunction->GetYmin(); 
	fYmax=fTestFunction->GetYmax();


	BCLog::Out(BCLog::summary,BCLog::summary," setup fitting function ...");
	fFitFunction = new TF2("fFitFunction",
			fTestFunction->GetExpFormula(),fXmin,fXmax,fYmin,fYmax);
	for (int i=1; i<fTestFunction->GetNpar(); i++)
		fFitFunction->FixParameter(i,fTestFunction->GetParameter(i));


	BCLog::Out(BCLog::summary,BCLog::summary," add parameter x ...");
	this->AddParameter("x",fXmin,fXmax);
	BCLog::Out(BCLog::summary,BCLog::summary," add parameter y ...");
	this->AddParameter("y",fYmin,fYmax);


	BCLog::Out(BCLog::summary,BCLog::summary," setup trees for chains ...");
	SetModel(this);
	BCModelOutput::WriteMarkovChain(true);
	SetFile(outputFile);
	

	BCLog::Out(BCLog::summary,BCLog::summary," book histograms for analysis ...");
	for (int i=0; i<fMaxChains; i++) {
		for (int j=1; j<fMaxLags; j++)
			fHistXYLags[i][j] = new TH2F(Form("fHistXYChain%dLag%d",i,j),
					Form("lag of %d",j),
					fNbinx,fXmin,fXmax,fNbiny,fYmin,fYmax);

		for (int j=1; j<fMax10thOfIters; j++)
			fHistXYIter[i][j] = new TH2F(Form("fHistXYChain%dIter%d0",i,j),
					Form("after %d%% of total iterations",j*10),
					fNbinx,fXmin,fXmax,fNbiny,fYmin,fYmax);

		fHChi2vsLags[i] = new TH1F(Form("HChi2vsLagsChain%d",i),"",
				fMaxLags-1,1,fMaxLags);
		fHChi2vsIter[i] = new TH1F(Form("HChi2vsIterChain%d",i),"",
				fMax10thOfIters-1,1,fMax10thOfIters);
	}
}

//=============================================================================

BCBenchmarkMCMC2D::~BCBenchmarkMCMC2D()
{
	if (fFitFunction) delete fFitFunction;

	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		for (int j=1; j<fMaxLags; j++)
			if (fHistXYLags[i][j]) delete fHistXYLags[i][j];

		for (int j=1; j<fMax10thOfIters; j++)
			if (fHistXYIter[i][j]) delete fHistXYIter[i][j];

		if (fHChi2vsLags[i]) delete fHChi2vsLags[i];
		if (fHChi2vsIter[i]) delete fHChi2vsIter[i];
	}
}

//=============================================================================

void BCBenchmarkMCMC2D::ProcessMCTrees()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) ProcessMCTree(i);
}

//=============================================================================

void BCBenchmarkMCMC2D::ProcessMCTree(int chainID)
{
	TTree *chain = this->MCMCGetMarkovChainTree(chainID);

	Double_t par0,par1;
	chain->SetBranchAddress("fParameter0",&par0);
	chain->SetBranchAddress("fParameter1",&par1);

	int Niters = chain->GetEntries();
	for (int i=0; i<Niters; i++) {
		chain->GetEntry(i);

		for (int j=1; j<fMaxLags; j++) {
			if (i%j==0) fHistXYLags[chainID][j]->Fill(par0,par1);
		}

		for (int j=1; j<fMax10thOfIters; j++) {
			if (i<=j*Niters/10) 
				fHistXYIter[chainID][j]->Fill(par0,par1);
		}
	}
}

//=============================================================================

void BCBenchmarkMCMC2D::PerformLagsTest()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		Chi2vsLagsOfChain(i);
	}
}

//=============================================================================

void BCBenchmarkMCMC2D::Chi2vsLagsOfChain(int chainID)
{
	for (int i=1; i<fMaxLags; i++) {
		fHistXYLags[chainID][i]->Fit(fFitFunction,"bq");
		fHChi2vsLags[chainID]->SetBinContent(i,
				fFitFunction->GetChisquare()/fFitFunction->GetNDF());
	}
}

//=============================================================================

void BCBenchmarkMCMC2D::PerformIterationsTest()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		Chi2vsIterOfChain(i);
	}
}

//=============================================================================

void BCBenchmarkMCMC2D::Chi2vsIterOfChain(int chainID)
{
	for (int i=1; i<fMax10thOfIters; i++) {
		fHistXYIter[chainID][i]->Fit(fFitFunction,"bq");
		fHChi2vsIter[chainID]->SetBinContent(i,
				fFitFunction->GetChisquare()/fFitFunction->GetNDF());
	}
}

//=============================================================================

void BCBenchmarkMCMC2D::WriteResults()
{
	TFile *file = this->GetFile();

	file->cd();
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		file->mkdir(Form("HistsForMCTree%d",i),
				Form("Histograms for Markov chain tree %d",i));
		file->cd(Form("HistsForMCTree%d",i));

		for (int j=1; j<fMaxLags; j++) {
			fHistXYLags[i][j]->SetXTitle("x");
			fHistXYLags[i][j]->SetYTitle("y");
			fHistXYLags[i][j]->Write();
		}

		for (int j=1; j<fMax10thOfIters; j++) {
			fHistXYIter[i][j]->SetXTitle("x");
			fHistXYIter[i][j]->SetYTitle("y");
			fHistXYIter[i][j]->Write();
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
