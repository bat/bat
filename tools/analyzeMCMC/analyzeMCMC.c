#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>

#include <math.h>
#include <iostream.h>

int main()
{

	// open file containing the MCMC tree

	TFile * file = new TFile("./MCMC.root", "READ");

	// get the MCMC tree from the file

	TTree * tree = (TTree*) file -> Get("MarkovChainTree");

	// get the number of parameters

	int nparameters;

	tree -> SetBranchAddress("fNParameters", &nparameters);

	tree -> GetEntry(0);

	// get probability density

	double logpdf;

	tree -> SetBranchAddress("fLogLikelihood", &logpdf);

	// get the parameters

	std::vector<double> parameter;

	for (int i = 0; i < nparameters; ++i)
		parameter.push_back(0.0);

	// loop over parameters

	for (int i = 0; i < nparameters; ++i)
		tree -> SetBranchAddress(Form("fParameter%i", i), &(parameter.at(i)));

	// loop over iterations and calculate limits

	int niterations = tree -> GetEntries();

	// debug
	niterations = 100001;

	double minlogpdf;
	double maxlogpdf;

	double meanlogpdf = 0;

	std::vector<double> minparameter;
	std::vector<double> maxparameter;

	std::vector<double> meanparameter;

	for (int i = 0; i < nparameters; ++i)
		{
			minparameter.push_back(0.0);
			maxparameter.push_back(0.0);
			meanparameter.push_back(0.0);
		}

	for (int iiteration = 0; iiteration < niterations; ++iiteration)
		{
			tree -> GetEntry(iiteration);

			if (logpdf < minlogpdf || iiteration == 0)
				minlogpdf = logpdf;

			if (logpdf > maxlogpdf || iiteration == 0)
				maxlogpdf = logpdf;

			for (int i = 0; i < nparameters; ++i)
				{
					if (parameter.at(i) < minparameter.at(i) ||  iiteration == 0)
						minparameter[i] = parameter.at(i);

					if (parameter.at(i) > maxparameter.at(i) ||  iiteration == 0)
						maxparameter[i] = parameter.at(i);

					meanparameter[i] += parameter.at(i);
				}

			meanlogpdf += logpdf;
		}

	// calculate mean values

	meanlogpdf = meanlogpdf / double(niterations);

	for (int i = 0; i < nparameters; ++i)
		meanparameter[i] = meanparameter.at(i) / double(niterations);

	// define histograms

	double minlimit = minlogpdf;
	double maxlimit = maxlogpdf;

	if (minlimit < 0)
		minlimit = 1.1 * minlimit;

	if (minlimit > 0)
		minlimit = 0.9 * minlimit;

	if (maxlimit < 0)
		maxlimit = 0.9 * maxlimit;

	if (maxlimit > 0)
		maxlimit = 1.1 * maxlimit;

	TH2D * hist_pdf_vs_iteration = new TH2D("hist_pdf_vs_iteration", "",
																					100, 0.0, double(niterations),
																					100, minlimit, maxlimit);
	hist_pdf_vs_iteration -> SetStats(kFALSE);
	hist_pdf_vs_iteration -> SetXTitle("Iteration");
	hist_pdf_vs_iteration -> SetYTitle("log(pdf)");

	// loop over iterations

	int order = 0;

	for (int iiteration = 0; iiteration < niterations; ++iiteration)
		{
			tree -> GetEntry(iiteration);

			// fill histograms

			hist_pdf_vs_iteration -> Fill(iiteration, logpdf);
		}

	// define autocorrelation functions

	std::vector<double> aclogpdf;

	int norders = int(floor(log10(niterations)));

	TGraph * graph_aclogpdf = new TGraph(norders);
	graph_aclogpdf -> SetMarkerStyle(20);

	std::vector <TGraph *> graph_acparameter;

	for (int i = 0; i < nparameters; ++i)
		{
			TGraph * g = new TGraph();
			g -> SetMarkerStyle(20+i+1);
			graph_acparameter.push_back(g);
		}

	std::vector <double> Aparameter;
	std::vector <double> Bparameter;
	std::vector <double> tempAparameter;

	for (int i = 0; i < nparameters; ++i)
		{
			Aparameter.push_back(0.0);
			Bparameter.push_back(0.0);
			tempAparameter.push_back(0.0);
		}

	for (int iorder = 0; iorder < norders; ++iorder)
		{
			double Apdf = 0;
			double Bpdf = 0;

			int k = int(pow(10, double(iorder)));

			for (int i = 0; i < niterations - k; ++i)
				{
					tree -> GetEntry(i);

					double tempApdf = logpdf - meanlogpdf;

					for (int iparameter = 0; iparameter < nparameters; ++iparameter)
						tempAparameter[iparameter] = parameter.at(iparameter) - meanparameter.at(iparameter);

					tree -> GetEntry(i + k);

					Apdf += tempApdf * (logpdf - meanlogpdf);

					for (int iparameter = 0; iparameter < nparameters; ++iparameter)
						Aparameter[iparameter] += tempAparameter.at(iparameter) *
							(parameter.at(iparameter) - meanparameter.at(iparameter));

				}

			for (int i = 0; i < niterations; ++i)
				{
					tree -> GetEntry(i);

					Bpdf += ((logpdf - meanlogpdf) *
									 (logpdf - meanlogpdf));

					for (int iparameter = 0; iparameter < nparameters; ++iparameter)
						Bparameter[iparameter] += ((parameter.at(iparameter) - meanparameter.at(iparameter)) *
																			 (parameter.at(iparameter) - meanparameter.at(iparameter)));
				}

			graph_aclogpdf -> SetPoint(iorder, double(iorder), Apdf/Bpdf);

			for (int iparameter = 0; iparameter < nparameters; ++iparameter)
				{
					TGraph * g = graph_acparameter.at(iparameter);

					g -> SetPoint(iorder,
												double(iorder),
												Aparameter.at(iparameter) / Bparameter.at(iparameter));

				}
		}

	// print to screen

	TCanvas * canvas1 = new TCanvas("canvas1");
	canvas1 -> cd();

	hist_pdf_vs_iteration -> Draw("COLZ");

	canvas1 -> Print("canvas1.pdf");


	TCanvas * canvas_aclogpdf = new TCanvas("canvas_aclogpdf");
	canvas_aclogpdf -> cd();

	TH2D * hist_acaxes = new TH2D("hist_acaxes", "",
																1, -0.5, double(norders)+0.5,
																1, -0.1, 1.1);
	hist_acaxes -> SetStats(kFALSE);
	hist_acaxes -> Draw();

	graph_aclogpdf -> Draw("SAMEPL");

	for (int iparameter = 0; iparameter < nparameters; ++iparameter)
		graph_acparameter.at(iparameter) -> Draw("SAMEPL");

	canvas_aclogpdf -> Print("canvas_aclogpdf.pdf");

	// close root file

	//	file -> Close();

	return 0;

}

// check acceptance (should be around 25%)
// trial distribution for dx. Now flat, choose Gaussian, Cauchy, arbitrary function, take correlations into account
// autocorrelation function p(l) = 1/N sum y(i) y(i-l)
// -> p(l) = exp(-l/lambda), eta= 1/(1*2*lambda), need lambda-1

