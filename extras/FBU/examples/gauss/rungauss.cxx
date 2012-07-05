#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BCFBU.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

int main()
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

	// ---- read histograms from a file ---- // 

	// open file
	TFile * file = new TFile("histograms.root", "READ");

	// check if file is open
	if (!file->IsOpen()) {
		std::cout << "Could not open file. Exit." << std::endl;
		std::cerr << "Run macro CreateHistograms.C in Root to create the file." << std::endl;
		return 1;
	}

	// read histograms
	TH2D* hist_migration = (TH2D*) file->Get("hist_migrationmatrix"); 
	TH1D* hist_truth     = (TH1D*) file->Get("hist_truthsignal"); 
	TH1D* hist_bkg       = (TH1D*) file->Get("hist_background");
	TH1D* hist_data      = (TH1D*) file->Get("hist_data");
	TH1D* hist_truthdata = (TH1D*) file->Get("hist_truthdata");

	// check if histograms are really there 
	if (!hist_migration) {
		std::cout << "Could not find migration matrix. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	if (!hist_truth) {
		std::cout << "Could not find truth histogram. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	if (!hist_bkg) {
		std::cout << "Could not find bkg histogram. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	if (!hist_data) {
		std::cout << "Could not find data histogram. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	if (!hist_truthdata) {
		std::cout << "Could not find truth data histogram. Run CreateHistograms before running the example. Exit." << std::endl;
		return(1);
	}

	// create new BCFBU object
	BCFBU* m = new BCFBU();

	// define parameter limits
	std::vector<double> paramin = std::vector<double>(10);
	std::vector<double> paramax = std::vector<double>(10);

	paramin[0]= 0.;
	paramax[0]= 100;

	paramin[1]= 0;
	paramax[1]= 100;

	paramin[2]= 0;
	paramax[2]= 100;

	paramin[3]= 0;
	paramax[3]= 500;

	paramin[4]= 0;
	paramax[4]= 800;

	paramin[5]= 0;
	paramax[5]= 800;

	paramin[6]= 0;
	paramax[6]= 500;

	paramin[7]= 0;
	paramax[7]= 100;

	paramin[8]= 0;
	paramax[8]= 100;

	paramin[9]= 0;
	paramax[9]= 100;

	// prepare response matrix and calculate efficiency
  m->PrepareResponseMatrix(hist_migration, hist_truth, hist_bkg, paramin, paramax);
 
	// set data
	m->SetDataHistogram(hist_data);
	
	// prepare summary tool
	BCSummaryTool* st = new BCSummaryTool(m);

	// perform unfolding
	m->MarginalizeAll();

	// print
	m->PrintResponseMatrix("response_matrix.eps");
	m->PrintMigrationMatrix("migration_matrix.eps");
	m->PrintEfficiencyHistogram("efficiency.eps");
	m->PrintTruthHistogram("truth.eps");
	m->PrintDataHistogram("data.eps");

	m->PrintAllMarginalized("plots.ps");
	m->PrintResults("results.txt");

	st->PrintParameterPlot("summary_parameters.ps");
	st->PrintCorrelationPlot("summary_correlation.ps");
	st->PrintCorrelationMatrix("summary_matrix.ps");


	TH1D *unfolded = (TH1D*) m->GetUnfoldedResult();

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);	
	c1->cd();
	unfolded->Draw();
	hist_truthdata->Draw("SAME");
	unfolded->SetMarkerColor(kBlue);
	unfolded->SetLineColor(kBlue);
	unfolded->Draw("E same");
	c1->SaveAs("unfolded.eps"); 

	// close log file
	BCLog::CloseLog();

	// clean up
	delete hist_migration;
	delete hist_truth;
	delete hist_bkg;
	delete hist_data;
	delete m;
	delete st;
  delete file;
	delete c1; 

	// no error 
	return 0;

}

