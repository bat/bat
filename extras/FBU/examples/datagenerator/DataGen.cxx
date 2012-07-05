#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>

#include <TROOT.h>
#include <TFile.h>

#include <string>

#include "DataGen.h"

// ---------------------------------------------------------
DataGen::DataGen(std::string name, int nbinstruth, double truthmin, double truthmax, int nbinsreco, double recomin, double recomax)
	: fName(name)
	, fFuncTruthSignal(0)
	, fFuncBackground(0)
	, fHistData(0)
	, fHistTruthData(0)
	, fHistTruthSignal(0)
	, fHistBackground(0)
	, fHistMigrationMatrix(0)
	, fHistEfficiency(0)
	, fNbinsTruth(nbinstruth)
	, fTruthMin(truthmin)
	, fTruthMax(truthmax)
	, fNbinsReco(nbinsreco)
	, fRecoMin(recomin)
	, fRecoMax(recomax)
{
	CreateHistograms();
}

// ---------------------------------------------------------
DataGen::~DataGen()
{
	if (fHistData)
		delete fHistData;

	if (fHistTruthSignal)
		delete fHistTruthSignal;

	if (fHistTruthData)
		delete fHistTruthData;

	if (fHistBackground)
		delete fHistBackground;

	if (fHistMigrationMatrix)
		delete fHistMigrationMatrix;

	if (fHistEfficiency)
		delete fHistEfficiency;
}

// ---------------------------------------------------------
int DataGen::FillHistograms(int nevents, int nsignal, int nbackground)
{
	// check if signal function exists
	if (!fFuncTruthSignal) {
		std::cout << "Warning: signal function does not exist. Migration matrix not filled" << std::endl;
		return 0;
	}

	// check if background function exists
	if (nbackground && !fFuncBackground) {
		std::cout << "Warning: background function does not exist. Migration matrix not filled" << std::endl;
		return 0;
	}
	// check if resolution function exists
	if (!fFuncResolution) {
		std::cout << "Warning: resolution function does not exist. Migration matrix not filled" << std::endl;
		return 0;
		}

	// check if migration matrix exists
	if (!fHistMigrationMatrix) {
		std::cout << "Warning: migration matrix histogram does not exist. Migration matrix not filled" << std::endl;
		return 0;
	}

	// check if data histogram exists
	if (!fHistData) {
		std::cout << "Warning: data histogram does not exist. Migration matrix not filled" << std::endl;
		return 0;
	}

	// check if data truth histogram exists
	if (!fHistTruthData) {
		std::cout << "Warning: data truth histogram does not exist. Migration matrix not filled" << std::endl;
		return 0;
	}

	// check if signal histogram exists
	if (!fHistTruthSignal) {
		std::cout << "Warning: signal histogram does not exist. Migration matrix not filled" << std::endl;
		return 0;
	}

	// check if background histogram exists
	if (!fHistBackground) {
		std::cout << "Warning: background histogram does not exist. Migration matrix not filled" << std::endl;
		return 0;
	}

	// loop over all events
	for (int i = 0; i < nevents; ++i) {
		// generate truth value
		double t = fFuncTruthSignal->GetRandom();
		
		// get smearing
		double s = fFuncResolution->GetRandom();

		// calculate reco value
		double r = t+s;

		// get background value
		double b = 0; 
		if (fFuncBackground)
			b = fFuncBackground->GetRandom();
		
		// fill histograms
		fHistMigrationMatrix->Fill(t, r);
		fHistTruthSignal->Fill(t);
		if (fFuncBackground)
			fHistBackground->Fill(b);
	}

	// normalize histograms
	fHistTruthSignal->Scale(1./double(nevents));
	fHistMigrationMatrix->Scale(1./double(nevents));

	// loop over all truth bins
	for (int i = 1; i <= fNbinsTruth; ++i) {
		double sum = 0;

		for (int j = 1; j <= fNbinsReco; ++j) 
			sum+= fHistMigrationMatrix->GetBinContent(i, j);
		
		fHistEfficiency->SetBinContent(i, sum/fHistTruthSignal->GetBinContent(i));
	}
	
	// loop over all signal events
	for (int i = 0; i < nsignal; ++i) {
		// generate truth value
		double t = fFuncTruthSignal->GetRandom();
		
		// get smearing
		double s = fFuncResolution->GetRandom();

		// calculate reco value
		double r = t+s;

		// fill histograms
		fHistData->Fill(r);
		fHistTruthData->Fill(t);
	}

	// loop over all background events
	for (int i = 0; i < nbackground; ++i) {
		// generate background value
		double b = fFuncBackground->GetRandom();

		// fill migration matrix
		fHistData->Fill(b);
	}


	// no error
	return 1;	
}

// ---------------------------------------------------------
int DataGen::Write(std::string filename)
{
	// create new file
	TFile* f = new TFile(filename.c_str(), "RECREATE");
	f->cd();

	if (fHistData)
		fHistData->Write();

	if (fHistTruthData)
		fHistTruthData->Write();

	if (fHistTruthSignal)
		fHistTruthSignal->Write();

	if (fHistBackground)
		fHistBackground->Write();

	if (fHistMigrationMatrix)
		fHistMigrationMatrix->Write();

	if (fHistEfficiency)
		fHistEfficiency->Write();

	// close file
	f->Close();

	// delete file
	delete f;

	// no error
	return 1;	
}

// ---------------------------------------------------------
int DataGen::CreateHistograms()
{
	// create data histogram
	fHistData = new TH1D("hist_data", ";;", fNbinsReco, fRecoMin, fRecoMax );
	fHistData->SetStats(kFALSE);

	// create histogram for truth data signal
	fHistTruthData = new TH1D("hist_truthdata", ";;", fNbinsTruth, fTruthMin, fTruthMax);
	fHistTruthData->SetStats(kFALSE);

	// create histogram for truth signal pdf
	fHistTruthSignal = new TH1D("hist_truthsignal", ";;", fNbinsTruth, fTruthMin, fTruthMax);
	fHistTruthSignal->SetStats(kFALSE);

	// create histogram for background pdf
	fHistBackground = new TH1D("hist_background", ";;", fNbinsReco, fRecoMin, fRecoMax);
	fHistBackground->SetStats(kFALSE);

	// create histogram for efficiency
	fHistEfficiency = new TH1D("hist_efficiency", ";;", fNbinsTruth, fTruthMin, fTruthMax);
	fHistEfficiency->SetStats(kFALSE);

	// create migration matrix
	fHistMigrationMatrix = new TH2D("hist_migrationmatrix", ";;", fNbinsTruth, fTruthMin, fTruthMax, fNbinsReco, fRecoMin, fRecoMax);
	fHistMigrationMatrix->SetStats(kFALSE);

	// no error
	return 1;		
}

// ---------------------------------------------------------
