// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project CombinationTool
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

#include <BAT/BCModel.h>

class TF1;

// This is a CombinationModel header file.
// Model source code is located in file CombinationTool/CombinationModel.cxx

// ---------------------------------------------------------
class CombinationModel : public BCModel
{
	public:

		// Constructors and destructor
		CombinationModel(const char * name, double xmin, double xmax);
		~CombinationModel();

		int SetSignalPrior(TF1* prior);

		int AddChannel(const char* channelname); 
		int AddChannelBackground(const char* channelname, const char* backgroundname, double xmin, double xmax);

		int GetNChannels()
		{ return int(fChannelNameContainer.size()); }; 

		int GetNChannelBackgrounds(int channelindex)
		{ return int(fChannelBackgroundNameContainer.at(channelindex).size()); }; 

		int GetContIndexChannel(const char* channelname);
		int GetContIndexChannelBackground(const char* channelname, const char* backgroundname);

		int GetParIndexChannel(const char* channelname);
		int GetParIndexChannelBackground(const char* channelname, const char* backgroundname);
		int GetParIndexChannelBackground(int channelindex, int backgroundindex);

		int SetChannelObservation(const char* channelname, double observation);
		int SetChannelUncertainty(const char* channelname, TF1* uncertainty);

		int SetChannelSignalPrior(const char* channelname, TF1* prior); 
		int SetChannelSignalPriorGauss(const char* channelname, double mean, double sigma); 
		int SetChannelSignalPriorGauss(const char* channelname, double mean, double sigmadown, double sigmaup); 
		int SetChannelBackgroundPrior(const char* channelname, const char* backgroundname, TF1* prior);
		int SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigma);
		int SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigmadown, double sigmaup);

		//		int SetChannelEfficiency(const char* channelname, double efficiency); 
		//		int SetChannelEfficiencyPrior(const char* channelname, TF1* prior);
		//		int SetChannelBackgroundEfficiency(const char* channelname, double efficiency); 
		//		int SetChannelBackgroundEfficiencyPrior(const char* channelname, TF1* prior);

		
		// set lumi and uncertainty
		// set efficiency and uncertainty
		// ...

		// Methods to overload, see file CombinationModel.cxx
		double LogAPrioriProbability(std::vector <double> parameters);
		double LogLikelihood(std::vector <double> parameters);

 protected:

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

		// prior container
		std::vector<TF1*> fChannelSignalPriorContainer;
		std::vector< std::vector<TF1*> > fChannelBackgroundPriorContainer;

		// index container
		int fParIndexSignal; 
		std::vector< std::vector<int> > fParIndexChannelBackground; 
		std::vector<int> fParIndexChannelObservationContainer;
		std::vector<int> fParIndexChannelEfficiencyContainer;

		// function container
		std::vector<TF1*> fFunctionContainer; 
};
// ---------------------------------------------------------

#endif

