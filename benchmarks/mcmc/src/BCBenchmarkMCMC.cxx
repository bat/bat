#include "BCBenchmarkMCMC.h" 

#include <BCMath.h> 

// --------------------------------------------------------- 

BCBenchmarkMCMC::BCBenchmarkMCMC(const char* name) : BCModel(name)
{

}

// --------------------------------------------------------- 

double BCBenchmarkMCMC::LogLikelihood(std::vector <double> parameters)
{

	double loglikelihood = 0.; 

	double x = parameters.at(0); 

	loglikelihood = log (fTestFunction -> Eval(x)); 

	return loglikelihood; 

}

// --------------------------------------------------------- 

void BCBenchmarkMCMC::PerformTest(std::vector<double> parameters, int index, BCH1D * hist, double * chi2, bool flag_print, const char * filename)
{

	// get histogram from BCH1D 
	TH1D * hist_temp = hist -> GetHistogram();

	// create new histogram
	TH1D * hist_prob = (TH1D*) hist_temp -> Clone(); 

	// initialize chi2
	chi2 = new double; 
	*chi2 = 0; 

	// get number of bins 
	int nbins = hist_temp -> GetNbinsX(); 

	// get number of entries
	double nentries = hist_temp -> GetEntries(); 

	// loop over bins 
	for (int i = 1; i <= nbins; ++i)
		{
			// get the parameter value which is looped over 
			double x = hist_temp -> GetBinCenter(i); 

			// set parameter 
			parameters[index] = x; 

			// calculate posterior probability 
			double prob = this -> Eval(parameters); 

			hist_prob -> SetBinContent(i, prob * hist_prob -> GetBinWidth(i)); 
		}

	// normalize
	if (hist_prob -> Integral() > 0.)
		hist_prob -> Scale(1.0 / hist_prob -> Integral() / hist_prob -> GetBinWidth(1)); 
	else
		{
			BCLog::Out(BCLog::warning, BCLog::warning,"BCModel::Calculate1DLikelihood(). Histogram has zero integral.");
			return; 
		}

	// loop over bins 
	for (int i = 1; i <= nbins; ++i)
		{
			// get bin content from BCH1D and Poisson uncertainty 
			//			double x = hist_temp -> GetBinContent(i) * nentries * hist_prob -> GetBinWidth(i);
			double x = hist_temp -> GetBinContent(i); 
			double sigma2 = x; 

			double x0 = hist_prob -> GetBinContent(i) * nentries * hist_prob -> GetBinWidth(i);

			if (sigma2 > 0.)
				*chi2 += (x-x0)*(x-x0)/2./sigma2; 

		}

	// devide by the d.o.f. 
	*chi2 /= double(nbins); 	

	// print to screen 
	std::cout << " Test results : " << std::endl; 
	std::cout << " chi2 = " << *chi2 << " (" << int(nbins) << " d.o.f.)" << std::endl; 
	std::cout << std::endl; 
	std::cout << " Print results to " << filename << std::endl; 

	// print to file 
	if (flag_print)
		{
			TCanvas * can = new TCanvas("can"); 
			can -> cd(); 
			hist -> Draw(); 
			hist_prob -> SetLineColor(kRed);
			hist_prob -> SetLineStyle(0); 
			hist_prob -> Draw("SAME"); 
			can -> Print(filename); 
		}

	// free memory 
	delete hist_prob; 

	return; 

}

// --------------------------------------------------------- 
