#include "BCModelEfficiency.h" 

#include <BCMath.h> 

// debug
#include <TMath.h> 

// --------------------------------------------------------- 

BCModelEfficiency::BCModelEfficiency() : BCModel()
{
	// define parameters 
	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

BCModelEfficiency::BCModelEfficiency(const char* name) : BCModel(name)
{
	// define parameters 
	this -> DefineParameters(); 
}

// --------------------------------------------------------- 

void BCModelEfficiency::DefineParameters()
{
	this -> AddParameter("mu1", 0.0,  1.0); // index 0 
	this -> AddParameter("mu2", 0.0, 1.0); // index 0 

	hist_bestfit = new TH1D("hist_bestfit", ";efficiency;N", 100, 0.0, 1.0); 

	hist_efficiency = new TH1D("hist_efficiency", ";efficiency;N", 100, 0.0, 1.0); 

}

// --------------------------------------------------------- 

double BCModelEfficiency::LogAPrioriProbability(std::vector <double> parameters)
{
	// get parameter ranges 
	double eff_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double eff_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	// calculate probabilities 
	double logprob = 0.;

	logprob -= log(eff_upper - eff_lower);

	return logprob;
}

// --------------------------------------------------------- 

double BCModelEfficiency::LogLikelihood(std::vector <double> parameters)
{
	// get parameters
	double eff = parameters.at(0);

	// get data values
	int k = BCMath::Nint(this -> GetDataPoint(0) -> GetValue(0));
	int n = BCMath::Nint(this -> GetDataPoint(0) -> GetValue(1));

	//	return BCMath::LogApproxBinomial(n, k, eff);

	return log(2.0 * TMath::Gaus(parameters.at(0), 0.1, 0.1, true) + TMath::Gaus(parameters.at(0), 0.6, 0.1, true)) + 
		log(2.0 * TMath::Gaus(parameters.at(1), 0.5, 0.1, true) + TMath::Gaus(parameters.at(1), 0.8, 0.1, true)); 
}

// --------------------------------------------------------- 

void BCModelEfficiency::MCMCUserInterface()
{

  for (int i = 0; i < fMCMCNChains; ++i)
      for (int j = 0; j < fMCMCNParameters; ++j)
	{
	  //	  cout << i << " " << j << this -> MCMCGetx(i, j) << " " << this -> MCMCGetLogProbx(i) << endl;  
	} 

  for (int i = 0; i < fMCMCNChains; ++i)
    hist_efficiency -> Fill(fMCMCx[i]); 

  if (fMCMCNIterations[0] == fMCMCNIterationsMax - 1)
    for (int i = 0; i < fMCMCNChains; ++i)
      hist_bestfit -> Fill((this -> MCMCGetMinimumPoints()).at(i)); 
}

// --------------------------------------------------------- 
