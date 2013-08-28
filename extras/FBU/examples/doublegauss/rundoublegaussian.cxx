#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BCFBU.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "TCanvas.h"

int main()
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

	// ---- read histograms from a file ---- //

	// open file
	TFile * file = TFile::Open("histograms.root", "READ");

	// check if file is open
	if (!file->IsOpen()) {
		std::cout << "Could not open file. Exit." << std::endl;
		std::cerr << "Run macro CreateHistograms.C in Root to create the file." << std::endl;
		return 1;
	}

	// read histograms
	TH2D* hist_migration = (TH2D*) file->Get("hist_migration");
	TH1D* hist_truth     = (TH1D*) file->Get("hist_truth");
	TH1D* hist_bkg       = (TH1D*) file->Get("hist_bkg");
	TH1D* hist_data      = (TH1D*) file->Get("hist_data");

	// check if histograms are really there
	if (!hist_migration) {
		std::cout << "Could not find migration matrix. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	if (!hist_truth) {
		std::cout << "Could not find truth matrix. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	if (!hist_bkg) {
		std::cout << "Could not find bkg matrix. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	if (!hist_data) {
		std::cout << "Could not find data matrix. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	// create new BCFBU object
	BCFBU* m = new BCFBU();

	// prepare response matrix and calculate efficiency
	m->PrepareResponseMatrix(hist_migration, hist_truth, hist_bkg);

	// set data
	m->SetDataHistogram(hist_data);

	// print
	m->PrintResponseMatrix("response_matrix.pdf");
	m->PrintMigrationMatrix("migration_matrix.pdf");
	m->PrintEfficiencyHistogram("efficiency.pdf");
	m->PrintTruthHistogram("truth.pdf");
	m->PrintDataHistogram("data.pdf");

	m->MarginalizeAll();

	TH1D *unfolded = (TH1D*) m->GetUnfoldedResult();

	TCanvas *bla = new TCanvas("bla","bla",600,600);

	//	bla->SetLogy();


	unfolded->SetMarkerColor(kBlue);
	unfolded->SetLineColor(kBlue);
	unfolded->Draw("E");

	hist_truth->Draw("same");

	bla->SaveAs("unfolded.pdf");


	// close log file
	BCLog::CloseLog();

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	// clean up
	delete hist_migration;
	delete hist_truth;
	delete hist_bkg;
	delete hist_data;
	delete m;
  delete file;

	// no error
	return 0;

}

