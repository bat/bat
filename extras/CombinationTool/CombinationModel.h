#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

/*!
 * \class CombinationModel
 * \brief A class for combining different measurements.
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 04.2010
 * This class is a base class for combining different measurements.
 */

#include <BAT/BCModel.h>

class ParameterSummary;
class TF1;

// ---------------------------------------------------------
class CombinationModel : public BCModel
{
 public:
	
	/** \name Constructors and destructors */
	/* @{ */
	
	/**
	 * Default constructor. 
	 * @param name The name of the signal to be combined.
	 * @param xmin The minimum of the signal.
	 * @param xmax The maximum of the signal.
	 */ 
	CombinationModel(const char * name, double xmin, double xmax);

	/**
	 * Default destructor. 
	 */ 
	~CombinationModel();

	/* @} */

	/** \name Public member functions (get) */
	/* @{ */

	/**
	 * Return the flag for using systematics. 
	 */ 
	bool GetFlagSystErrors()
	{ return fFlagSystErrors; }; 

	/**
	 * Return the number of systematic uncertainties
	 */ 
	int GetNSystErrors()
	{ return int(fSystErrorNameContainer.size()); };

	/**
	 * Return the number of channels.
	 */ 
	int GetNChannels()
	{ return int(fChannelNameContainer.size()); }; 

	/**
	 * Return the number of background sources for a channel. 
	 * @param channelindex The index of the channel. 
	 */ 
	int GetNChannelBackgrounds(int channelindex)
	{ return int(fChannelBackgroundNameContainer.at(channelindex).size()); }; 

	/* @} */

	/** \name Public member functions (set) */
	/* @{ */

	/**
	 * Set flag for using systematics (1) or not (0). 
	 * @param flag The flag.
	 */ 
	void SetFlagSystErrors(bool flag) 
	{ fFlagSystErrors = flag; }; 

	/**
	 * Set the number of observed events in one channel. 
	 * @param channelname The name of the channel. 
	 * @param observation The number of observed events.
	 * @return An error code. 
	 */ 
	int SetChannelObservation(const char* channelname, double observation);

	/**
	 * Set a prior for the signal to be combined. 
	 * @param prior The prior TF1 function
	 * @return An error code. 
	 */ 
	int SetSignalPrior(TF1* prior);

	/**
	 * Set the signal prior for one channel. *
	 * @param channelname The name of the channel. 
	 * @param prior The prior TF1 function.
	 * @return An error code. 
	 */ 
	int SetChannelSignalPrior(const char* channelname, TF1* prior); 

	/**
	 * Set a Gaussian signal prior for one channel. *
	 * @param channelname The name of the channel. 
	 * @param mean The mean value of the Gauss function.
	 * @param sigma The standard deviation of the Gauss function.
	 * @return An error code. 
	 */ 
	int SetChannelSignalPriorGauss(const char* channelname, double mean, double sigma); 

	/**
	 * Set a Gaussian signal prior for one channel. The Gaussian has
	 * different standard deviations for values left and right of the
	 * mean.
	 * @param channelname The name of the channel. 
	 * @param mean The mean value of the Gauss function.
	 * @param sigmadown The (left) standard deviation of the Gauss function.
	 * @param sigmaup The (right) standard deviation of the Gauss function.
	 * @return An error code. 
	 */ 
	int SetChannelSignalPriorGauss(const char* channelname, double mean, double sigmadown, double sigmaup); 

	/**
	 * Set the background prior for one channel. *
	 * @param channelname The name of the channel. 
	 * @param backgroundname The name of the background. 
	 * @param prior The prior TF1 function.
	 * @return An error code. 
	 */ 
	int SetChannelBackgroundPrior(const char* channelname, const char* backgroundname, TF1* prior);

	/**
	 * Set a Gaussian background prior for one channel. *
	 * @param channelname The name of the channel. 
	 * @param backgroundname The name of the background. 
	 * @param mean The mean value of the Gauss function.
	 * @param sigma The standard deviation of the Gauss function.
	 * @return An error code. 
	 */ 
	int SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigma);

	/**
	 * Set a Gaussian background prior for one channel. The Gaussian has
	 * different standard deviations for values left and right of the
	 * mean.
	 * @param channelname The name of the channel. 
	 * @param backgroundname The name of the background. 
	 * @param mean The mean value of the Gauss function.
	 * @param sigmadown The (left) standard deviation of the Gauss function.
	 * @param sigmaup The (right) standard deviation of the Gauss function.
	 * @return An error code. 
	 */ 
	int SetChannelBackgroundPriorGauss(const char* channelname, const char* backgroundname, double mean, double sigmadown, double sigmaup);

	/**
	 * Set the systematic uncertainty on the signal for one channel. 
	 * @param systerrorname The name of the systematic uncertainty. 
	 * @param channelname The name of the channel.
	 * @param sigmadown The down-scale uncertainty. 
	 * @param sigmaup The up-scale uncertainty. 
	 * @return An error code. 
	 */ 
	int SetSystErrorChannelSignal(const char* systerrorname, const char* channelname, double sigmadown, double sigmaup);

	/**
	 * Set the systematic uncertainty on the background for one channel. 
	 * @param systerrorname The name of the systematic uncertainty. 
	 * @param channelname The name of the channel.
	 * @param backgroundname The name of the background. 
	 * @param sigmadown The down-scale uncertainty. 
	 * @param sigmaup The up-scale uncertainty. 
	 * @return An error code. 
	 */ 
	int SetSystErrorChannelBackground(const char* systerrorname, const char* channelname, const char* backgroundname, double sigmadown, double sigmaup); 

	/* @} */

	/** \name Public member functions (run) */
	/* @{ */

	/**
	 * Run the combination with current settings.
	 * @return An error code.
	 */ 
	int PerformAnalysis(); 

	/**
	 * Run the combination with and without systematics. Perform
	 * analysis for each channel separately.
	 */ 
	int PerformFullAnalysis(); 

	/**
	 * Perform analysis on one channel.
	 * @param channelname The name of the channel
	 * @param flag_syst Flag: systematics on (1) or off (0).
	 * @return A summary object of the signal parameter.
	 */ 
	virtual ParameterSummary PerformSingleChannelAnalysis(const char* channelname, bool flag_syst) = 0;

	/* @} */

	/** \name Public member functions (misc) */
	/* @{ */

	/**
	 * Add a channel to be combined. 
	 * @param channelname The name of the channel. 
	 * @return An error code. 
	 */ 
	virtual int AddChannel(const char* channelname); 

	/**
	 * Add a background source to a channel. 
	 * @param channelname The name of the channel.
	 * @param backgroundname The name of the background
	 * @param xmin The minimum value of the background contribution. 
	 * @param xmax The maximum value of the background contribution. 
	 * @return An error code. 
	 */ 
	int AddChannelBackground(const char* channelname, const char* backgroundname, double xmin, double xmax);

	/**
	 * Add a source of systematic uncertainty. 
	 * @param systerrorname The name of the systematic uncertainty. 
	 * @return An error code. 
	 */ 
	int AddSystError(const char* systerrorname); 

	/**
	 * Print an overview on all channels and the combination.
	 * @param filename The name of the file. 
	 * @return An error code. 
	 */ 
	int PrintChannelOverview(const char* filename); 

	/**
	 * Write the summary of the analysis to a file.
	 * @param filename The name of the file. 
	 * @return An error code. 
	 */ 
	int PrintChannelSummary(const char* filename); 

	/* @} */

	/** \name Public member functions (overloaded from BCModel) */
	/* @{ */

	/**
	 * Calculate the log of the prior probability. 
	 * @param parameters A vector of parameters. 
	 * @return The log of the prior probability.
	 */
	double LogAPrioriProbability(std::vector <double> parameters) = 0;

	/**
	 * Calculate the log of the likelihood
	 * @param parameters A vector of parameters. 
	 * @return The log of the likelihood.
	 */
	double LogLikelihood(std::vector <double> parameters) = 0;

 protected:

	/** \name Protected member functions (get) */
	/* @{ */

	/**
	 * Return the container index for a channel. 
	 * @param channelname The name of the channel.
	 */ 
	int GetContIndexChannel(const char* channelname);

	/**
	 * Return the container index for a background. 
	 * @param channelname The name of the channel.
	 * @param backgroundname The name of the background.
	 */ 
	int GetContIndexChannelBackground(const char* channelname, const char* backgroundname);

	/**
	 * Return the container index for a systematic uncertainty. 
	 * @param systerrorname The name of the systematic uncertainty.
	 */ 
	int GetContIndexSystError(const char* systerrorname);

	/**
	 * Return the parameter index for a channel. 
	 * @param channelname The name of the channel.
	 */ 
	int GetParIndexChannel(const char* channelname);

	/**
	 * Return the parameter index for a background. 
	 * @param channelname The name of the channel.
	 * @param backgroundname The name of the background.
	 */ 
	int GetParIndexChannelBackground(const char* channelname, const char* backgroundname);

	/**
	 * Return the parameter index for a background. 
	 * @param channelindex The container index of the channel.
	 * @param backgroundindex The container index of the background.
	 */ 
	int GetParIndexChannelBackground(int channelindex, int backgroundindex);

	/**
	 * Return the container index for a systematic uncertainty. 
	 * @param systerrorindex The container index of the systematic uncertainty.
	 */ 
	int GetParIndexSystError(int systerrorindex);

	/* @} */

	/**
	 * A flag for including the systematic uncertainties (1) or not
	 * (0). 
	 */ 
	bool fFlagSystErrors;

	/** 
	 * The container of channel names. Index is channel container index. 
	 */
	std::vector<std::string> fChannelNameContainer; 

	/** 
	 * The container of background names. Indeces are channel container
	 * index and background container index.
	 */ 
	std::vector< std::vector<std::string> > fChannelBackgroundNameContainer;

	/**
	 * The container of systematic uncertainy names. Index is error
	 * container index. 
	 */ 
	std::vector<std::string> fSystErrorNameContainer; 

	/**
	 * The container of observations. 
	 */ 
	std::vector<double> fChannelObservation;

	/**
	 * The container of up-scale signal systematic
	 * uncertainties. Indeces are error container index and channel
	 * container index.
	 */ 
	std::vector< std::vector<double> > fSystErrorChannelSigmaUpContainer;

	/**
	 * The container of down-scale signal systematic
	 * uncertainties. Indeces are error container index and channel
	 * container index.
	 */ 
	std::vector< std::vector<double> > fSystErrorChannelSigmaDownContainer;

	/**
	 * The container of up-scale background systematic
	 * uncertainties. Indeces are error container index, channel
	 * container index and background container index.
	 */ 
	std::vector< std::vector< std::vector<double> > > fSystErrorSigmaUpContainer;

	/**
	 * The container of down-scale background systematic
	 * uncertainties. Indeces are error container index, channel
	 * container index and background container index.
	 */ 
	std::vector< std::vector< std::vector<double> > > fSystErrorSigmaDownContainer;

	/**
	 * The container of signal priors. The index is the channel container index. 
	 */ 
	std::vector<TF1*> fChannelSignalPriorContainer;

	/**
	 * The container of background priors. The indeces are channel
	container index and background container index. 
	*/ 
	std::vector< std::vector<TF1*> > fChannelBackgroundPriorContainer;

	/**
	 * The parameter index of the signal. 
	 */ 
	int fParIndexSignal; 

	/**
	 * The parameter index container of the background. Indeces are the
	 * channel container index and the background container index.
	 */ 
	std::vector< std::vector<int> > fParIndexChannelBackground; 

	/**
	 * The parameter index container for systematic errors. Index is the
	 * error container index.
	 */ 
	std::vector<int> fParIndexSystErrorContainer;

	/**
	 * A container of function for memory handling. 
	 */ 
	std::vector<TF1*> fFunctionContainer; 

	/**
	 * A summary of the signal parameter for the combination without
	 * systematic errors.
	 */ 
	ParameterSummary* fSummaryCombinationNoSyst; 

	/**
	 * A summary of the signal parameter for the combination with
	 * systematic errors.
	 */ 
	ParameterSummary* fSummaryCombinationSyst; 
	
	/**
	 * The container of summaries of the signal parameter for the single
	 * channel analyses without systematic errors.
	 */
 	std::vector<ParameterSummary*> fSummaryChannelNoSyst; 

	/**
	 * The container of summaries of the signal parameter for the single
	 * channel analyses with systematic errors.
	 */
	std::vector<ParameterSummary*> fSummaryChannelSyst; 

};
// ---------------------------------------------------------

#endif

