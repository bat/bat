#ifndef __COMBINATIONXSEC__H
#define __COMBINATIONXSEC__H

#include <CombinationModel.h>

#include "ParameterSummary.h" 

class BCH1D; 

// ---------------------------------------------------------
class CombinationXSec : public CombinationModel
{
 public:
	
	// Constructors and destructor
	CombinationXSec(const char * name, double xmin, double xmax);
	~CombinationXSec();
	
	int AddChannel(const char* channelname); 

	int GetParIndexChannelEfficiency(const char* channelname); 
	int GetParIndexChannelLuminosity(const char* channelname); 

	int SetChannelEfficiencyPrior(const char* channelname, TF1* prior); 
	int SetChannelEfficiencyPriorGauss(const char* channelname, double mean, double sigma);

	int SetChannelLuminosityPrior(const char* channelname, TF1* prior); 
	int SetChannelLuminosityPriorGauss(const char* channelname, double mean, double sigma); 

	int SetChannelBR(const char* channelname, double BR);

	ParameterSummary PerformSingleChannelAnalysis(const char* channelname, bool flag_syst = true);

	void DefineParameters();
	double LogAPrioriProbability(std::vector <double> parameters);
	double LogLikelihood(std::vector <double> parameters);
	void MCMCUserIterationInterface();

 private:

	std::vector<TF1*> fChannelEfficiencyPriorContainer;
	std::vector<int> fParIndexChannelEfficiency;

	std::vector<TF1*> fChannelLuminosityPriorContainer;
	std::vector<int> fParIndexChannelLuminosity;
	
	std::vector<double> fChannelBR; 
};
// ---------------------------------------------------------

#endif

