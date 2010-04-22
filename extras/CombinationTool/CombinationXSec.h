#ifndef __COMBINATIONXSEC__H
#define __COMBINATIONXSEC__H

#include <CombinationModel.h>

class BCH1D; 

// ---------------------------------------------------------
class CombinationXSec : public CombinationModel
{
 public:
	
	// Constructors and destructor
	CombinationXSec(const char * name, double xmin, double xmax);
	~CombinationXSec();
	
	int AddChannel(const char* channelname); 

	int SetChannelEfficiency(const char* channelname, double efficiency); 
	int SetChannelEfficiencyPrior(const char* channelname, TF1* prior); 
	int SetChannelEfficiencyPriorGauss(const char* channelname, double mean, double sigma);

	int SetChannelLuminosity(const char* channelname, double efficiency); 
	int SetChannelLuminosityPrior(const char* channelname, TF1* prior); 
	int SetChannelLuminosityPriorGauss(const char* channelname, double mean, double sigma); 

	int SetChannelBR(const char* channelname, double BR);

	void DefineParameters();
	double LogAPrioriProbability(std::vector <double> parameters);
	double LogLikelihood(std::vector <double> parameters);
	void MCMCUserIterationInterface();

	void PrintChannelOverview(const char* filename); 
	void PrintChannels(const char* filename); 

 private:

	std::vector<double> fChannelEfficiency;
	std::vector<TF1*> fChannelEfficiencyPriorContainer;

	std::vector<double> fChannelLuminosity;
	std::vector<TF1*> fChannelLuminosityPriorContainer;
	
	std::vector<double> fChannelBR; 

	std::vector<BCH1D*> fChannelSignal; 

};
// ---------------------------------------------------------

#endif

