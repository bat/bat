/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <THStack.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphAsymmErrors.h>

#include "../../BAT/BCMath.h"
#include "../../BAT/BCLog.h"

#include "BCMTFChannel.h"
#include "BCMTFProcess.h"
#include "BCMTFTemplate.h"
#include "BCMTFSystematic.h"
#include "BCMTFSystematicVariation.h"

#include "BCMTF.h"

// ---------------------------------------------------------
BCMTF::BCMTF()
   : BCModel("Multi-template Fitter")
   , fNChannels(0)
   , fNProcesses(0)
   , fNSystematics(0)
   , fFlagEfficiencyConstraint(false)
{}

// ---------------------------------------------------------
BCMTF::BCMTF(const char * name)
   : BCModel(name)
   , fNChannels(0)
   , fNProcesses(0)
   , fNSystematics(0)
   , fFlagEfficiencyConstraint(false)
{}

// ---------------------------------------------------------
BCMTF::~BCMTF()
// default destructor
{
   for (int i = 0; i < fNChannels; ++i)
      delete fChannelContainer.at(i);
}

// ---------------------------------------------------------
int BCMTF::GetChannelIndex(const char * name)
{
   // loop over all channels and compare names
   for (int i = 0; i < fNChannels; ++i) {
      // get channel
      BCMTFChannel * channel = GetChannel(i);

      // compare names
      if (!channel->GetName().compare(name))
         return i;
   }

   // if channel does not exist, return -1
   return -1;
}

// ---------------------------------------------------------
int BCMTF::GetProcessIndex(const char * name)
{
   // loop over all processs and compare names
   for (int i = 0; i < fNProcesses; ++i) {
      // get process
      BCMTFProcess * process = GetProcess(i);

      // compare names
      if (!process->GetName().compare(name))
         return i;
   }

   // if process does not exist, return -1
   return -1;
}

// ---------------------------------------------------------
int BCMTF::GetSystematicIndex(const char * name)
{
   // loop over all systematics and compare names
   for (int i = 0; i < fNSystematics; ++i) {
      // get systematic
      BCMTFSystematic * systematic = GetSystematic(i);

      // compare names
      if (!systematic->GetName().compare(name))
         return i;
   }

   // if process does not exist, return -1
   return -1;
}

// ---------------------------------------------------------
int BCMTF::SetTemplate(const char * channelname, const char * processname, TH1D hist, double efficiency, double norm)
{
   // get channel index
   int channelindex = GetChannelIndex(channelname);

   // check if channel exists
   if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
      return -1;
   }

   // get process index
   int processindex = GetProcessIndex(processname);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Process does not exist.");
      return -1;
   }

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // get template
   BCMTFTemplate * bctemplate = channel->GetTemplate(processindex);

   // remove statistics box
   hist.SetStats(kFALSE);

   int color = GetProcess(processindex)->GetHistogramColor();
   if (color < 0)
     color = 2 + processindex;
   int fillstyle = GetProcess(processindex)->GetHistogramFillStyle();
   if (fillstyle < 0)
     fillstyle = 1001;
   int linestyle = GetProcess(processindex)->GetHistogramLineStyle();
   if (linestyle < 0)
     linestyle = 1;

   // set color and fill style
   hist.SetFillColor(color);
   hist.SetFillStyle(fillstyle);
   hist.SetLineStyle(linestyle);

   // create new histogram
   TH1D * temphist = new TH1D(hist);

   // set histogram
   bctemplate->SetHistogram(temphist, norm);

   // set efficiency
   bctemplate->SetEfficiency(efficiency);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::SetTemplate(const char * channelname, const char * processname, std::vector<TF1 *> * funccont, int nbins, double efficiency)
{
   // get channel index
   int channelindex = GetChannelIndex(channelname);

   // check if channel exists
   if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
      return -1;
   }

   // get process index
   int processindex = GetProcessIndex(processname);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Process does not exist.");
      return -1;
   }

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // get template
   BCMTFTemplate * bctemplate = channel->GetTemplate(processindex);

   // set histogram
   bctemplate->SetFunctionContainer(funccont, nbins);

   // set efficiency
   bctemplate->SetEfficiency(efficiency);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::SetData(const char * channelname, TH1D hist, double minimum, double maximum)
{
   int channelindex = GetChannelIndex(channelname);

   // check if channel exists
   if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
      return -1;
   }

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // get template
   BCMTFTemplate * data = channel->GetData();

   // remove statistics box
   hist.SetStats(kFALSE);

   // set marker
   hist.SetMarkerStyle(20);
   hist.SetMarkerSize(1.1);

   // set divisions
   hist.SetNdivisions(509);

   // remove old data set if it exists
   if (data->GetHistogram()) {
      delete data->GetHistogram();
      data->SetHistogram(0);
   }

   // remove old uncertainty histograms if they exist
   if (channel->GetHistUncertaintyBandExpectation()) {
      delete channel->GetHistUncertaintyBandExpectation();
      channel->SetHistUncertaintyBandExpectation(0);
   }
   if (channel->GetHistUncertaintyBandPoisson()) {
      delete channel->GetHistUncertaintyBandPoisson();
      channel->SetHistUncertaintyBandPoisson(0);
   }

   // create new histograms for uncertainty bands
   //	 double minimum = floor(TMath::Max(0., hist.GetMinimum() - 7.*sqrt(hist.GetMinimum())));
   if (minimum==-1)
      minimum = 0;
   if (maximum==-1)
      maximum = ceil(hist.GetMaximum() + 5.*sqrt(hist.GetMaximum()));

   std::vector<double> a(hist.GetNbinsX()+1);
   for (int i = 0; i < hist.GetNbinsX()+1; ++i) {
      a[i] = hist.GetXaxis()->GetBinLowEdge(i+1);
   }

   TH2D* hist_uncbandexp = new TH2D(TString::Format("UncertaintyBandExpectation_%i", BCLog::GetHIndex()), "",
                                    hist.GetNbinsX(), &a[0], 1000, minimum, maximum);
   hist_uncbandexp->SetStats(kFALSE);

   TH2D* hist_uncbandpoisson = new TH2D(TString::Format("UncertaintyBandPoisson_%i", BCLog::GetHIndex()), "",
                                        hist.GetNbinsX(), &a[0], int(maximum-minimum), minimum, maximum);
   hist_uncbandpoisson->SetStats(kFALSE);

   // set histograms
   data->SetHistogram(new TH1D(hist), hist.Integral());
   channel->SetHistUncertaintyBandExpectation(hist_uncbandexp);
   channel->SetHistUncertaintyBandPoisson(hist_uncbandpoisson);

   // set y-range for printing
   channel->SetRangeY(minimum, maximum);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::AddChannel(const char * name)
{
   // check if channel exists
   for (int i = 0; i < fNChannels; ++i) {
      // compare names
      if (GetChannelIndex(name) >= 0) {
         BCLog::OutWarning("BCMultitemplateFitter::AddChannel() : Channel with this name exists already.");
         return -1;
      }
   }

   // create new channel
   BCMTFChannel * channel = new BCMTFChannel(name);

   // create new data
   BCMTFTemplate * bctemplate = new BCMTFTemplate(channel->GetName().c_str(), "data");

   // add data
   channel->SetData(bctemplate);

   // add process templates
   for (int i = 0; i < fNProcesses; ++i) {
      // get process
      BCMTFProcess * process = GetProcess(i);

      // create new template
      BCMTFTemplate * bctemplate = new BCMTFTemplate(name, process->GetName().c_str());

      // add template
      channel->AddTemplate(bctemplate);
   }

   // loop over all systematics
   for (int i = 0; i < fNSystematics; ++i) {
      // get systematic
      BCMTFSystematic * systematic = GetSystematic(i);

      // create new systematic variation
      BCMTFSystematicVariation * variation = new BCMTFSystematicVariation(name, systematic->GetName().c_str(), fNProcesses);

      // add systematic variation
      channel->AddSystematicVariation(variation);
   }

   // add channel
   fChannelContainer.push_back(channel);

   // increase number of channels
   fNChannels++;

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::AddProcess(const char * name, double nmin, double nmax, int color, int fillstyle, int linestyle)
{
   // check if process exists
   for (int i = 0; i < fNProcesses; ++i) {
      // compare names
      if (GetProcessIndex(name) >= 0) {
         BCLog::OutWarning("BCMultitemplateFitter::AddProcess() : Process with this name exists already.");
         return -1;
      }
   }

   // create new process
   BCMTFProcess * process = new BCMTFProcess(name);
   process->SetHistogramColor(color);
   process->SetHistogramFillStyle(fillstyle);
   process->SetHistogramLineStyle(linestyle);

   // add process
   fProcessContainer.push_back(process);

   // add process templates
   for (int i = 0; i < fNChannels; ++i) {
      // get channel
      BCMTFChannel * channel = GetChannel(i);

      // create new template
      BCMTFTemplate * bctemplate = new BCMTFTemplate(channel->GetName().c_str(), name);

      // add template
      channel->AddTemplate(bctemplate);

      // loop over all systematic
      for (int j = 0; j < fNSystematics; ++j) {
         // get systematic variation
         BCMTFSystematicVariation * variation = channel->GetSystematicVariation(j);

         // add histogram
         variation->AddHistograms(0, 0);
      }
   }

   // increase number of processes
   fNProcesses++;

   // add parameter index to container
   fProcessParIndexContainer.push_back(GetNParameters());

   // add parameter
   AddParameter(name, nmin, nmax);

   // add a functional form for the expectation
   fExpectationFunctionContainer.push_back(0);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::AddSystematic(const char * name, double min, double max)
{
   // check if systematic exists
   for (int i = 0; i < fNSystematics; ++i) {
      // compare names
      if (GetSystematicIndex(name) >= 0) {
         BCLog::OutWarning("BCMultitemplateFitter::AddSystematic() : Systematic with this name exists already.");
         return -1;
      }
   }

   // create new systematic
   BCMTFSystematic * systematic = new BCMTFSystematic(name);

   // add systematic
   fSystematicContainer.push_back(systematic);

   // add systematic variations
   for (int i = 0; i < fNChannels; ++i) {
      // get channel
      BCMTFChannel * channel = GetChannel(i);

      // create new systematic variation
      BCMTFSystematicVariation * variation = new BCMTFSystematicVariation(channel->GetName().c_str(), name, fNProcesses);

      // add systematic variation
      channel->AddSystematicVariation(variation);
   }
   // ...

   // increase number of systematices
   fNSystematics++;

   // add parameter index to container
   fSystematicParIndexContainer.push_back(GetNParameters());

   // add parameter
   AddParameter(name, min, max);

   // add a functional form for the expectation
   fExpectationFunctionContainer.push_back(0);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::SetSystematicVariation(const char * channelname, const char * processname,  const char * systematicname, double variation_up, double variation_down)
{

   // get channel index
   int channelindex = GetChannelIndex(channelname);

   // check if channel exists
   if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
      return -1;
   }

   // get process index
   int processindex = GetProcessIndex(processname);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Process does not exist.");
      return -1;
   }

   // get systematic index
   int systematicindex = GetSystematicIndex(systematicname);

   // check if systematic exists
   if (systematicindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Systematic does not exist.");
      return -1;
   }

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   BCMTFTemplate * bctemplate = channel->GetTemplate(processindex);

   TH1D * hist_template = bctemplate->GetHistogram();

   TH1D hist_up = TH1D(*hist_template);
   TH1D hist_down = TH1D(*hist_template);

   int nbins = hist_up.GetNbinsX();

   // loop over all bins
   for (int ibin = 1; ibin <= nbins; ++ibin) {
      hist_up.SetBinContent(ibin, variation_up);
      hist_down.SetBinContent(ibin, variation_down);
   }

   // get systematic variation
   BCMTFSystematicVariation * variation = channel->GetSystematicVariation(systematicindex);

   // set histogram
   variation->SetHistograms(processindex, new TH1D(hist_up), new TH1D(hist_down));

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::SetSystematicVariation(const char * channelname, const char * processname,  const char * systematicname, TH1D hist_up, TH1D hist_down)
{
   // get channel index
   int channelindex = GetChannelIndex(channelname);

   // check if channel exists
   if (channelindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Channel does not exist.");
      return -1;
   }

   // get process index
   int processindex = GetProcessIndex(processname);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Process does not exist.");
      return -1;
   }

   // get systematic index
   int systematicindex = GetSystematicIndex(systematicname);

   // check if systematic exists
   if (systematicindex < 0) {
      BCLog::OutWarning("BCMultitemplateFitter::SetTemplate() : Systematic does not exist.");
      return -1;
   }

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // get systematic variation
   BCMTFSystematicVariation * variation = channel->GetSystematicVariation(systematicindex);

   // set histogram
   variation->SetHistograms(processindex, new TH1D(hist_up), new TH1D(hist_down));

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTF::SetSystematicVariation(const char * channelname, const char * processname,  const char * systematicname, TH1D hist, TH1D hist_up, TH1D hist_down)
{
   // get number of bins
   int nbins = hist.GetNbinsX();

   TH1D * hist_up_rel = new TH1D(hist);
   TH1D * hist_down_rel = new TH1D(hist);

   // loop over all bins
   for (int ibin = 1; ibin <= nbins; ++ibin) {
      hist_up_rel->SetBinContent(ibin, (hist_up.GetBinContent(ibin) - hist.GetBinContent(ibin)) / hist.GetBinContent(ibin));
      hist_down_rel->SetBinContent(ibin, (hist.GetBinContent(ibin) - hist_down.GetBinContent(ibin)) / hist.GetBinContent(ibin));
   }

   // set the systematic variation
   return SetSystematicVariation(channelname, processname, systematicname, *hist_up_rel, *hist_down_rel);
}

// ---------------------------------------------------------
int BCMTF::PrintSummary(const char * filename)
{
   // open file
   std::ofstream ofi(filename);
   ofi.precision(3);

   // check if file is open
   if(!ofi.is_open()) {
      BCLog::OutWarning(Form("BCMultitemplateFitter::PrintSummary() : Could not open file %s", filename));
      return 0;
   }

   ofi
      << " Multi template fitter summary " << std::endl
      << " ----------------------------- " << std::endl
      << std::endl
      << " Number of channels      : " << fNChannels << std::endl
      << " Number of processes     : " << fNProcesses << std::endl
      << " Number of systematics   : " << fNSystematics << std::endl
      << std::endl;

   ofi
      << " Channels :" << std::endl;
   for (int i = 0; i < GetNChannels(); ++i) {
      ofi
         << " " << i
         << " : \"" << GetChannel(i)->GetName().c_str() << "\""
         << std::endl;
   }
   ofi
      << std::endl;

   ofi
      << " Processes :" << std::endl;
   for (int i = 0; i < GetNProcesses(); ++i) {
      ofi
         << " " << i
         << " : \"" << GetProcess(i)->GetName().c_str()  << "\""
         << " (par index " << GetParIndexProcess(i) << ")"
         << std::endl;
   }
   ofi
      << std::endl;

   ofi
      << " Systematics :" << std::endl;
   for (int i = 0; i < GetNSystematics(); ++i) {
      ofi
         << " " << i
         << " : \"" << GetSystematic(i)->GetName().c_str()  << "\""
         << " (par index " << GetParIndexSystematic(i) << ")"
         << std::endl;
   }
   ofi
      << std::endl;
   if (GetNSystematics() == 0)
      ofi
         << " - none - " << std::endl;

   ofi
      << " Goodness-of-fit: " << std::endl;
   for (int i = 0; i < GetNChannels(); ++i) {
      ofi
         << " i : \"" << GetChannel(i)->GetName().c_str() << "\" : chi2 = "
         << CalculateChi2( i, GetBestFitParameters() )
         << std::endl;
   }
   ofi
      << std::endl;


   // close file
   ofi.close();

   // no error
   return 1;
}

// ---------------------------------------------------------
double BCMTF::Expectation(int channelindex, int binindex, const std::vector<double> & parameters)
{
   double expectation = 0.;

   // loop over all processes
   for (int i = 0; i < fNProcesses; ++i) {
      // get efficiency
      double efficiency = Efficiency(channelindex, i, binindex, parameters);

      // get probability
      double probability = Probability(channelindex, i, binindex, parameters);

      // get parameter index
      int parindex = fProcessParIndexContainer[i];

      // add to expectation
      expectation += ExpectationFunction(parindex, channelindex, i, parameters)
         * efficiency
         * probability;
   }

   // check if expectation is positive
   if (expectation < 0)
      expectation = 0.;

   return expectation;
}

// ---------------------------------------------------------
double BCMTF::ExpectationFunction(int parindex, int channelindex, int processindex, const std::vector<double> & parameters)
{
   // get function container
   std::vector<TF1 *> * funccont = fChannelContainer[channelindex]->GetTemplate(processindex)->GetFunctionContainer();

   if (funccont->size()>0)
      return 1.;

   else if (!fExpectationFunctionContainer[parindex])
      return parameters[parindex];

   else {
      TF1 * func = fExpectationFunctionContainer[parindex];
      return func->Eval(parameters[parindex]);
   }
}

// ---------------------------------------------------------
double BCMTF::Efficiency(int channelindex, int processindex, int binindex, const std::vector<double> & parameters)
{
   // get channel
   BCMTFChannel * channel = fChannelContainer[channelindex];

   double efficiency = channel->GetTemplate(processindex)->GetEfficiency();

   double defficiency = 1.;

   // loop over all systematics
   for (int i = 0; i < fNSystematics; ++i) {
      if (!(fSystematicContainer[i]->GetFlagSystematicActive()))
         continue;

      // get parameter index
      int parindex = fSystematicParIndexContainer[i];

      // get histogram
      TH1D * hist = 0;

      if (parameters[parindex] > 0)
         hist = channel->GetSystematicVariation(i)->GetHistogramUp(processindex);
      else
         hist = channel->GetSystematicVariation(i)->GetHistogramDown(processindex);

      // check if histogram exists
      if (!hist)
         continue;

      // multiply efficiency
      defficiency += parameters[parindex] * hist->GetBinContent(binindex);
   }

   // calculate efficiency
   efficiency *= defficiency;

   // make sure efficiency is between 0 and 1
   if (fFlagEfficiencyConstraint) {
      if (efficiency < 0.)
         efficiency = 0.;
      if (efficiency > 1.)
         efficiency = 1.;
   }

   return efficiency;
}

// ---------------------------------------------------------
double BCMTF::Probability(int channelindex, int processindex, int binindex, const std::vector<double> & parameters)
{
   // get histogram
   TH1D * hist = fChannelContainer[channelindex]->GetTemplate(processindex)->GetHistogram();

   // get function container
   std::vector<TF1 *> * funccont = fChannelContainer[channelindex]->GetTemplate(processindex)->GetFunctionContainer();

   // this needs to be fast
   if (!hist && !(funccont->size()>0))
      return 0.;

   if (hist)
      return hist->GetBinContent(binindex);
   else {
      int parindex = fProcessParIndexContainer[processindex];
      return funccont->at(binindex-1)->Eval(parameters[parindex]);
   }
}

// ---------------------------------------------------------
int BCMTF::PrintStack(const char * channelname, const std::vector<double> & parameters, const char * filename, const char * options)
{
   int index = GetChannelIndex(channelname);

   return PrintStack(index, parameters, filename, options);
}

// ---------------------------------------------------------
int BCMTF::PrintStack(int channelindex, const std::vector<double> & parameters, const char * filename, const char * options)
{
   // todo:
   // - add difference/ratio/significance plot below
   // - check for b0/1 if the mcmc was run

   // check if parameters are filled
   if (!parameters.size())
      return -1;

   // check options
   bool flag_logx   = false; // plot x-axis in log-scale
   bool flag_logy   = false; // plot y-axis in log-scale
   bool flag_bw     = false; // plot in black and white

   bool flag_sum    = false; // plot sum of all templates
   bool flag_stack  = false; // plot stack of templates

   bool flag_e0     = false; // do not draw error bars on data
   bool flag_e1     = false; // draw sqrt(N) error bars on data

   bool flag_b0     = false; // draw an error band on the expectation
   bool flag_b1     = false; // draw an error band on the number of events

   if (std::string(options).find("logx") < std::string(options).size())
      flag_logx = true;

   if (std::string(options).find("logy") < std::string(options).size())
      flag_logy = true;

   if (std::string(options).find("bw") < std::string(options).size())
      flag_bw = true;

   if (std::string(options).find("sum") < std::string(options).size())
      flag_sum = true;

   if (std::string(options).find("stack") < std::string(options).size())
      flag_stack = true;

   if (std::string(options).find("e0") < std::string(options).size())
      flag_e0 = true;

   if (std::string(options).find("e1") < std::string(options).size())
      flag_e1 = true;

   if (std::string(options).find("b0") < std::string(options).size())
      flag_b0 = true;

   if (std::string(options).find("b1") < std::string(options).size())
      flag_b1 = true;

   if (!flag_e0)
      flag_e1=true;

   // check if MCMC ran
   if (!(GetMarginalizationMethod() == BCIntegrate::kMargMetropolis)) {
     flag_b0 = false;
     flag_b1 = false;
     BCLog::OutWarning("BCMTF::PrintStack : Did not run MCMC. Error bands are not available.");
   }

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // create canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // set log or linear scale
   if (flag_logx)
      c1->SetLogx();

   if (flag_logy)
      c1->SetLogy();

   // get data histogram
   TH1D* hist_data = channel->GetData()->GetHistogram();

   // get number of bins
   int nbins = hist_data->GetNbinsX();

   // define sum of templates
   TH1D* hist_sum = new TH1D(*hist_data);
   hist_sum->SetLineColor(kBlack);
   for (int i = 1; i <= nbins; ++i)
      hist_sum->SetBinContent(i, 0);

   // define error band
   TH1D* hist_error_band = new TH1D(*hist_data);
   hist_error_band->SetFillColor(kBlack);
   hist_error_band->SetFillStyle(3005);
   hist_error_band->SetLineWidth(1);
   hist_error_band->SetStats(kFALSE);
   hist_error_band->SetMarkerSize(0);

   TGraphAsymmErrors * graph_error_exp = new TGraphAsymmErrors(nbins);
   //   graph_error_exp->SetLineWidth(2);
   graph_error_exp->SetMarkerStyle(0);
   graph_error_exp->SetFillColor(kBlack);
   graph_error_exp->SetFillStyle(3005);

   // get histogram for uncertainty band
   TH2D* hist_uncbandexp = channel->GetHistUncertaintyBandExpectation();

   // fill error band
   if (flag_b0) {
      for (int i = 1; i <= nbins; ++i) {
         TH1D * proj = hist_uncbandexp->ProjectionY("_py", i, i);
         if (proj->Integral() > 0)
            proj->Scale(1.0 / proj->Integral());
         double quantiles[3];
         double sums[3] = {0.16, 0.5, 0.84};
         proj->GetQuantiles(3, quantiles, sums);
         graph_error_exp->SetPoint(i-1, hist_data->GetBinCenter(i), quantiles[1]);
         graph_error_exp->SetPointError(i-1, 0.0, 0.0, quantiles[1] - quantiles[0], quantiles[2]-quantiles[1]);
         hist_error_band->SetBinContent(i, 0.5*(quantiles[2]+quantiles[0]));
         hist_error_band->SetBinError(i, 0, 0.5*(quantiles[2]-quantiles[0]));
         delete proj;
      }
   }

   // create stack
   THStack * stack = new THStack("", "");

   // create a container of temporary histograms
   std::vector<TH1D *> histcontainer;

   // get number of templates
   unsigned int ntemplates = GetNProcesses();

   // loop over all templates
   for (unsigned int i = 0; i < ntemplates; ++i) {

      // get histogram
      TH1D * temphist = channel->GetTemplate(i)->GetHistogram();

      // get function container
      std::vector<TF1 *> * funccont = channel->GetTemplate(i)->GetFunctionContainer();

      // create new histogram
      TH1D * hist(0);

      if (temphist)
         hist = new TH1D( *(temphist) );
      else if (funccont)
         hist = new TH1D( *(channel->GetData()->GetHistogram()));
      else
         continue;

      // set histogram style
      int color = GetProcess(i)->GetHistogramColor();
      if (color < 0)
        color = 2 + i;
      int fillstyle = GetProcess(i)->GetHistogramFillStyle();
      if (fillstyle < 0)
        fillstyle = 1001;
      int linestyle = GetProcess(i)->GetHistogramLineStyle();
      if (linestyle < 0)
        linestyle = 1;

      // set color and fill style
      hist->SetFillColor(color);
      hist->SetFillStyle(fillstyle);
      hist->SetLineStyle(linestyle);

      if (flag_bw) {
        hist->SetFillColor(0);
      }

      // scale histogram
      for (int ibin = 1; ibin <= nbins; ++ibin) {

         // get efficiency
         double efficiency = Efficiency(channelindex, i, ibin, parameters);

         // get probability
         double probability = Probability(channelindex, i, ibin, parameters);

         // get parameter index
         int parindex = GetParIndexProcess(i);

         // add to expectation
         double expectation = parameters[parindex] * efficiency * probability;

         // set bin content
         hist->SetBinContent(ibin, expectation);

         // add bin content
         hist_sum->SetBinContent(ibin, hist_sum->GetBinContent(ibin) + expectation);
      }

      // add histogram to container (for memory management)
      histcontainer.push_back(hist);

      // add histogram to stack
      stack->Add(hist);
   }

   //draw data
   hist_data->Draw("P0");

   // define variable for maximum in y-direction
   double ymax = 0;;

   if (flag_e1)
      ymax = hist_data->GetMaximum() + sqrt(hist_data->GetMaximum());
   else
      ymax = hist_data->GetMaximum();

   // set range user
   hist_data->GetYaxis()->SetRangeUser(channel->GetRangeYMin(), channel->GetRangeYMax());

   // draw stack
   if (flag_stack) {
      stack->Draw("SAMEHIST");
      if (stack->GetMaximum() > ymax)
         ymax = stack->GetMaximum();
   }

   // draw error band on number of observed events
   if (flag_b1) {
      channel->CalculateHistUncertaintyBandPoisson();
      TH1D* hist_temp = channel->CalculateUncertaintyBandPoisson(0.001, 0.999, kRed);
      hist_temp->Draw("SAMEE2");
      channel->CalculateUncertaintyBandPoisson(0.023, 0.977, kOrange)->Draw("SAMEE2");
      channel->CalculateUncertaintyBandPoisson(0.159, 0.841, kGreen)->Draw("SAMEE2");

      // get bin with maximum
      int ymaxbin = hist_temp->GetMaximumBin();

      if (hist_temp->GetBinContent(ymaxbin)+hist_temp->GetBinError(ymaxbin)> ymax)
         ymax = hist_temp->GetBinContent(ymaxbin)+hist_temp->GetBinError(ymaxbin);
   }

   // draw error band on expectation
   if (flag_b0) {
      hist_error_band->Draw("SAMEE2");
   }

   if (flag_sum)
      hist_sum->Draw("SAME");

   //draw data again
   if (flag_e0)
      hist_data->Draw("SAMEP0");

   if (flag_e1)
      hist_data->Draw("SAMEP0E");

   hist_data->GetYaxis()->SetRangeUser(0., 1.1* ymax);

   // redraw the axes
   gPad->RedrawAxis();

   // print
   c1->Print(filename);

   // free memory
   for (unsigned int i = 0; i < histcontainer.size(); ++i) {
      TH1D * hist = histcontainer.at(i);
      delete hist;
   }
   delete stack;
   delete c1;
   delete graph_error_exp;
   delete hist_error_band;
   delete hist_sum;

   // no error
   return 1;
}

// ---------------------------------------------------------
double BCMTF::CalculateChi2(int channelindex, const std::vector<double> & parameters)
{
   if (parameters.size() == 0)
      return -1;

   double chi2 = 0;

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // get data
   BCMTFTemplate * data = channel->GetData();

   // get histogram
   TH1D * hist = data->GetHistogram();

   // check if histogram exists
   if (hist) {
      // get number of bins in data
      int nbins = hist->GetNbinsX();

      // loop over all bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         // get expectation value
         double expectation = Expectation(channelindex, ibin, parameters);

         // get observation
         double observation = hist->GetBinContent(ibin);

         // add Poisson term
         chi2 += (expectation - observation) * (expectation - observation) / expectation;
      }
   }

   // return chi2;
   return chi2;
}

// ---------------------------------------------------------
double BCMTF::CalculateChi2(const std::vector<double> & parameters)
{
   if (parameters.size() == 0)
      return -1;

   double chi2 = 0;

   // get number of channels
   int nchannels = GetNChannels();

   // loop over all channels
   for (int i = 0; i < nchannels; ++i) {
      chi2 += CalculateChi2(i, parameters);
   }

   // return chi2
   return chi2;
}

// ---------------------------------------------------------
double BCMTF::CalculateCash(int channelindex, const std::vector<double> & parameters)
{
   if (parameters.size() == 0)
      return -1;

   double cash = 0;

   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // get data
   BCMTFTemplate * data = channel->GetData();

   // get histogram
   TH1D * hist = data->GetHistogram();

   // check if histogram exists
   if (hist) {
      // get number of bins in data
      int nbins = hist->GetNbinsX();

      // loop over all bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         // get expectation value
         double expectation = Expectation(channelindex, ibin, parameters);

         // get observation
         double observation = hist->GetBinContent(ibin);

         // calculate Cash statistic
         cash += 2. * (expectation - observation);

         // check negative log
         if (observation > 0)
            cash += 2. * observation * log (observation/expectation);
      }
   }

   // return cash;
   return cash;

}

// ---------------------------------------------------------
double BCMTF::CalculateCash(const std::vector<double> & parameters)
{
   if (parameters.size() == 0)
      return -1;

   double cash = 0;

   // get number of channels
   int nchannels = GetNChannels();

   // loop over all channels
   for (int i = 0; i < nchannels; ++i) {
      cash += CalculateCash(i, parameters);
   }

   // return cash
   return cash;
}

// ---------------------------------------------------------
double BCMTF::CalculatePValue(int channelindex, const std::vector<double> & parameters)
{
   // get channel
   BCMTFChannel * channel = GetChannel(channelindex);

   // get data histogram
   TH1D * hist = channel->GetData()->GetHistogram();

   // check if histogram exists
   if (!hist) {
      return -1;
   }

   // get number of bins in data
   int nbins = hist->GetNbinsX();

   // copy observed and expected values
   std::vector<unsigned> observation(nbins);
   std::vector<double> expectation(nbins);

   // loop over all bins
   for (int ibin = 0; ibin < nbins; ++ibin) {
      // get expectation value
      expectation[ibin] = Expectation(channelindex, ibin + 1, parameters);

      // get observation
      observation[ibin]= unsigned(hist->GetBinContent(ibin + 1));
   }

   // create pseudo experiments
   static const unsigned nIterations = unsigned (1e5);
   return BCMath::FastPValue(observation, expectation, nIterations, fRandom->GetSeed());
}

// ---------------------------------------------------------
double BCMTF::CalculatePValue(const std::vector<double> & parameters)
{

   // copy observed and expected values
   std::vector<unsigned> observation;
   std::vector<double> expectation;

   // loop over all channels
   for (int ichannel = 0; ichannel < fNChannels; ++ichannel) {

      // get channel
      BCMTFChannel * channel = fChannelContainer[ichannel];

      // check if channel is active
      if (!(channel->GetFlagChannelActive()))
         continue;

      // get data histogram
      TH1D * hist = channel->GetData()->GetHistogram();

      // check if histogram exists
      if (!hist) {
         return -1;
      }

      // get number of bins in data
      int nbins = hist->GetNbinsX();

      // loop over all bins
      for (int ibin = 0; ibin < nbins; ++ibin) {
         // get expectation value
         expectation.push_back(Expectation(ichannel, ibin + 1, parameters));

         // get observation
         observation.push_back(unsigned(hist->GetBinContent(ibin + 1)));
      }
   }

   // create pseudo experiments
   static const unsigned nIterations = unsigned(1e5);
   fPValue = BCMath::FastPValue(observation, expectation, nIterations, fRandom->GetSeed());
   fPValueNDoF = BCMath::CorrectPValue(fPValue, parameters.size(), observation.size());

   return fPValue;
}

// ---------------------------------------------------------
double BCMTF::LogLikelihood(const std::vector<double> & parameters)
{
   double logprob = 0.;

   // loop over all channels
   for (int ichannel = 0; ichannel < fNChannels; ++ichannel) {

      // get channel
      BCMTFChannel * channel = fChannelContainer[ichannel];

      // check if channel is active
      if (!(channel->GetFlagChannelActive()))
         continue;

      // get data
      BCMTFTemplate * data = channel->GetData();

      // get histogram
      TH1D * hist = data->GetHistogram();

      // check if histogram exists
      if (!hist)
         continue;

      // get number of bins in data
      int nbins = data->GetNBins();

      // loop over all bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {

         // get expectation value
         double expectation = Expectation(ichannel, ibin, parameters);

         // get observation
         double observation = hist->GetBinContent(ibin);

         // add Poisson term
         logprob += BCMath::LogPoisson(observation, expectation);
      }
   }

   return logprob;
}

// ---------------------------------------------------------
void BCMTF::MCMCUserIterationInterface()
{
   // loop over all channels
   for (int ichannel = 0; ichannel < fNChannels; ++ichannel) {

      // get channel
      BCMTFChannel * channel = fChannelContainer[ichannel];

      // check if channel is active
      if (!(channel->GetFlagChannelActive()))
         continue;

      // get data
      BCMTFTemplate * data = channel->GetData();

      // get histogram
      TH1D * hist_data = data->GetHistogram();

      // check if histogram exists
      if (!hist_data)
         continue;

      // get histogram for uncertainty band
      TH2D* hist_uncbandexp = channel->GetHistUncertaintyBandExpectation();

      // check if histogram exists
      if (!hist_uncbandexp)
         continue;

      // get number of bins in data
      int nbins = hist_data->GetNbinsX();

      // loop over all bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {

         // get expectation value
         double expectation = Expectation(ichannel, ibin, fMCMCx);

         // fill uncertainty band on expectation
         hist_uncbandexp->Fill(hist_data->GetBinCenter(ibin), expectation);
      }
   }

}

// ---------------------------------------------------------
