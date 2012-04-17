#ifndef __BCMTF__H
#define __BCMTF__H

/*!
 * \class BCMTF
 * \brief A class for ...
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * \version 1.0
 * \date 04.2012
 * \detail
 *
 *
 *
 *
 */

/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "../../BAT/BCModel.h"

#include <TH1D.h>

class BCMTFChannel;
class BCMTFProcess;
class BCMTFSystematic;
class TF1;

// ---------------------------------------------------------
class BCMTF : public BCModel
{

   public:

      // Constructors and destructor
      BCMTF();
      ~BCMTF();

      // setters

      // getters

      // return number of channels
      int GetNChannels()
         { return fNChannels; };

      // return number of processes
      int GetNProcesses()
         { return fNProcesses; };

      // return number of systematics
      int GetNSystematics()
         { return fNSystematics; };

      // return channel index
      int GetChannelIndex(const char * name);

      // return process index
      int GetProcessIndex(const char * name);

      // return systematic index
      int GetSystematicIndex(const char * name);

      // return parameter index of process
      int GetParIndexProcess(int index)
         { return fProcessParIndexContainer.at(index); };

      // return parameter index of systematic
      int GetParIndexSystematic(int index)
         { return fSystematicParIndexContainer.at(index); };

      // return a channel
      BCMTFChannel * GetChannel(int index)
         { return fChannelContainer.at(index); };

      // return a process
      BCMTFProcess * GetProcess(int index)
         { return fProcessContainer.at(index); };

      // return a systematic
      BCMTFSystematic * GetSystematic(int index)
         { return fSystematicContainer.at(index); };

      // misc

      // add a channel
      int AddChannel(const char * name);

      // add a process
      int AddProcess(const char * name, double nmin = 0., double nmax = 1.);

      // add a systematic uncertainty
      int AddSystematic(const char * name, double min = -5., double max = 5.);

      // set template
      int SetTemplate(const char * channelname, const char * processname, TH1D hist, double efficiency = 1.);

      // set template
      int SetTemplate(const char * channelname, const char * processname, std::vector<TF1 *> * funccont, int nbins, double efficiency = 1.);

      // set systematic variation
      int SetSystematicVariation(const char * channelname, const char * processname,  const char * systematicname, double variation_up, double variation_down);

      // set systematic variation
      int SetSystematicVariation(const char * channelname, const char * processname,  const char * systematicname, TH1D hist_up, TH1D hist_down);

      // set systematic variation
      int SetSystematicVariation(const char * channelname, const char * processname,  const char * systematicname, TH1D hist, TH1D hist_up, TH1D hist_down);

      // set data
      int SetData(const char * channelname, TH1D hist);

      // set flag (efficiency constraint on (true) or off (false)
      void SetFlagEfficiencyConstraint(bool flag)
         { fFlagEfficiencyConstraint = flag; };

      // set a function for the expectation
      void SetExpectationFunction(int parindex, TF1 * func)
         { fExpectationFunctionContainer[parindex] = func; };

      // print summary
      int PrintSummary(const char * filename = "summary.txt");

      // return expectation value for a channel and bin
      double Expectation(int channelindex, int binindex, const std::vector<double> & parameters);

      // return a function for the expecation
      double ExpectationFunction(int parindex, int channelindex, int processindex, const std::vector<double> & parameters);

      // return efficiency for a channel, process and bin
      double Efficiency(int channelindex, int processindex, int binindex, const std::vector<double> & parameters);

      // return probability for a channel, process and bin
      double Probability(int channelindex, int processindex, int binindex, const std::vector<double> & parameters);

      // print stack
      int PrintStack(int channelindex, const std::vector<double> & parameters, const char * filename = "stack.eps", const char * options = "");

      // print stack
      int PrintStack(const char * channelname, const std::vector<double> & parameters, const char * filename = "stack.eps", const char * options = "");

      // calculate chi2 for single channel
      double CalculateChi2(int channelindex, const std::vector<double> & parameters);

      // calculate chi2 for all channels
      double CalculateChi2(const std::vector<double> & parameters);

      // calculate Cash statistics for single channel
      double CalculateCash(int channelindex, const std::vector<double> & parameters);

      // calculate Cash statistics for all channels
      double CalculateCash(const std::vector<double> & parameters);

      // BAT

      // the log likelihood
      double LogLikelihood(const std::vector<double> & parameters);

    // void MCMCIterationInterface();

 private:

      // a container of channels
      std::vector<BCMTFChannel *> fChannelContainer;

      // a container of processes
      std::vector<BCMTFProcess *> fProcessContainer;

      // a container of systematic sources
      std::vector<BCMTFSystematic *> fSystematicContainer;

      // number of channels
      int fNChannels;

      // number of processes
      int fNProcesses;

      // number of systematics
      int fNSystematics;

      // a container of parameter indeces for the process normalization
      std::vector<int> fProcessParIndexContainer;

      // a container of parameter indeces for the systematics
      std::vector<int> fSystematicParIndexContainer;

      // flag: efficiency within 0 and 1 (true) or not (false)
      bool fFlagEfficiencyConstraint;

      // a container of functions for the expectation
      std::vector<TF1 *> fExpectationFunctionContainer;

};
// ---------------------------------------------------------

#endif

