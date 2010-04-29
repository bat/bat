#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

#include <BAT/BCModel.h>

#include "ParameterSummary.h"

class TF1;

// ---------------------------------------------------------
class CombinationModel : public BCModel
{
	public:

		// Constructors and destructor
		CombinationModel(const char * name, double xmin, double xmax);
		~CombinationModel();

		int PerformAnalysis(); 
		int PerformFullAnalysis(); 

		virtual ParameterSummary PerformSingleChannelAnalysis(const char* channelname, bool flag_syst) = 0;

		int SetSignalPrior(TF1* prior);

		int AddChannel(const char* channelname); 
		int AddChannelBackground(const char* channelname, const char* backgroundname, double xmin, double xmax);

		int AddSystError(const char* systerrorname); 
		int SetSystErrorChannelSignal(const char* channelname, double sigmadown, double sigmaup);
		int SetSystErrorChannelBackground(const char* systerrorname, const char* channelname, const char* backgroundname, double sigmadown, double sigmaup); 

		int GetNSystErrors()
		{ return int(fSystErrorNameContainer.size()); };

		int GetNChannels()
		{ return int(fChannelNameContainer.size()); }; 

		int GetNChannelBackgrounds(int channelindex)
		{ return int(fChannelBackgroundNameContainer.at(channelindex).size()); }; 

		int GetContIndexChannel(const char* channelname);
		int GetContIndexChannelBackground(const char* channelname, const char* backgroundname);
		int GetContIndexSystError(const char* systerrorname);

		int GetParIndexChannel(const char* channelname);
		int GetParIndexChannelBackground(const char* channelname, const char* backgroundname);
		int GetParIndexChannelBackground(int channelindex, int backgroundindex);
		int GetParIndexSystError(int systerrorindex);

		int SetChannelObservation(const char* channelname, double observation);
		int SetChannelUncertainty(const char* channelname, TF1* uncertainty);

		int SetChannelSignalPrior(const char* channelname, TF1* prior); 
		int SetChannelSignalPriorGauss(const char* channelname, double mean, double sigma); 
		int SetChannelSignalPriorGauss(const char* channelname, double mean, double sigmadown, double sigmaup); 

		int SetChannelBackgroundPrior(const char* channelname, const char* backgroundname, TF1* prior);
		int SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigma);
		int SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigmadown, double sigmaup);

		int PrintChannelOverview(const char* filename); 
		int PrintChannelSummary(const char* filename); 

		void SetFlagSystErrors(bool flag) 
		{ fFlagSystErrors = flag; }; 

		bool GetFlagSystErrors()
		{ return fFlagSystErrors; }; 

		// Methods to overload, see file CombinationModel.cxx
		double LogAPrioriProbability(std::vector <double> parameters) = 0;
		double LogLikelihood(std::vector <double> parameters) = 0;

 protected:

		// flags
		bool fFlagSystErrors; 

		// name container
		std::vector<std::string> fChannelNameContainer; 
		std::vector< std::vector<std::string> > fChannelBackgroundNameContainer;
		std::vector<std::string> fSystErrorNameContainer; 

		// observation container
		std::vector<double> fChannelObservation;

		// efficiency container
		std::vector<double> fChannelEfficiencyContainer; 

		// uncertainty container
		std::vector<TF1*> fChannelUncertaintyContainer;
		std::vector<TF1*> fChannelEfficiencyPriorContainer;
		std::vector< std::vector< std::vector<double> > > fSystErrorSigmaUpContainer;
		std::vector< std::vector< std::vector<double> > > fSystErrorSigmaDownContainer;

		// prior container
		std::vector<TF1*> fChannelSignalPriorContainer;
		std::vector< std::vector<TF1*> > fChannelBackgroundPriorContainer;

		// index container
		int fParIndexSignal; 
		std::vector< std::vector<int> > fParIndexChannelBackground; 
		std::vector<int> fParIndexChannelObservationContainer;
		std::vector<int> fParIndexChannelEfficiencyContainer;
		std::vector<int> fParIndexSystErrorContainer;

		// function container
		std::vector<TF1*> fFunctionContainer; 

		// analysis summaries
		ParameterSummary* fSummaryCombinationNoSyst; 
		ParameterSummary* fSummaryCombinationSyst; 
	
		std::vector<ParameterSummary*> fSummaryChannelNoSyst; 
		std::vector<ParameterSummary*> fSummaryChannelSyst; 

};
// ---------------------------------------------------------

#endif

