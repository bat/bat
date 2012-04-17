/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <sys/stat.h>
#include <iostream>

#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TFile.h>

#include "../../BAT/BCLog.h"
#include "../../BAT/BCH1D.h"

#include "BCMTF.h"
#include "BCMTFChannel.h"
#include "BCMTFTemplate.h"
#include "BCMTFSystematic.h"
#include "BCMTFComparisonTool.h"

#include "BCMTFAnalysisFacility.h"

// ---------------------------------------------------------
BCMTFAnalysisFacility::BCMTFAnalysisFacility(BCMTF * mtf)
 : fRandom(new TRandom3(0))
 , fFlagMCMC(false)
{
   fMTF = mtf;
};

// ---------------------------------------------------------
BCMTFAnalysisFacility::~BCMTFAnalysisFacility()
{
};

// ---------------------------------------------------------
std::vector<TH1D> BCMTFAnalysisFacility::BuildEnsemble(std::vector<double> parameters)
{
   // create vector of histograms
   std::vector<TH1D> histograms;

   // loop over channels
   for (int ichannel = 0; ichannel < fMTF->GetNChannels(); ++ichannel) {

      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // create new histogram
      TH1D hist( *(channel->GetData()->GetHistogram()) );

      // get number of bins
      int nbins = hist.GetNbinsX();

      // loop over all bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         double expectation = fMTF->Expectation(ichannel, ibin, parameters);

         double observation = gRandom->Poisson(expectation);

         hist.SetBinContent(ibin, observation);
      }

      // add histogram
      histograms.push_back(hist);
   }

   // return histograms
   return histograms;
}

// ---------------------------------------------------------
TTree * BCMTFAnalysisFacility::BuildEnsembles(TTree * tree, int nensembles)
{
   // get number of parameters
   int nparameters = fMTF->GetNParameters();

   // create tree variables for input tree
   std::vector<double> parameters(nparameters);

   // set branch addresses
   for (int i = 0; i < nparameters; ++i) {
      tree->SetBranchAddress(Form("Parameter%i", i), &(parameters[i]));
   }

   // create tree
   TTree * tree_out = new TTree("ensembles", "ensembles");

   // get number of channels
   int nchannels = fMTF->GetNChannels();

   // create matrix of number of bins
   std::vector< std::vector<double> > nbins_matrix;

   // prepare the tree variables
   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // get number of bins
      int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

      // create new matrix row
      std::vector<double> nbins_column(nbins);

      // add matrix
      nbins_matrix.push_back(nbins_column);
   }

   std::vector<double> in_parameters(nparameters);

   // create branches
   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // get number of bins
      int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

      // loop over bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         // create branches
         tree_out->Branch(Form("channel_%i_bin_%i", ichannel, ibin),
                                  &(nbins_matrix[ichannel])[ibin-1],
                                  "n/D");
      }
   }

   for (int i = 0; i < nparameters; ++i) {
      tree_out->Branch(Form("parameter_%i", i), &in_parameters[i], Form("parameter_%i/D", i));
   }

   // create vector of histograms
   std::vector<TH1D> histograms;

   // loop over ensembles
   for (int iensemble = 0; iensemble < nensembles; ++iensemble) {
      // get random event from tree
      int index = (int) fRandom->Uniform(tree->GetEntries());
      tree->GetEntry(index);

      // create ensembles
      histograms = BuildEnsemble(parameters);

      // copy information from histograms into tree variables
      // loop over channels
      for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
         // get channel
         BCMTFChannel * channel = fMTF->GetChannel(ichannel);

         // get number of bins
         int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

         // loop over bins
         for (int ibin = 1; ibin <= nbins; ++ibin) {
            // fill tree variable
            (nbins_matrix[ichannel])[ibin-1] = histograms.at(ichannel).GetBinContent(ibin);
         }
      }

      // copy parameter information
      for (int i = 0; i < nparameters; ++i) {
         in_parameters[i] = parameters.at(i);
      }

      // fill tree
      tree_out->Fill();
   }

   // return tree
   return tree_out;
}

// ---------------------------------------------------------
TTree * BCMTFAnalysisFacility::BuildEnsembles(std::vector<double> parameters, int nensembles)
{
   // create tree
   TTree * tree = new TTree("ensembles", "ensembles");

   // get number of channels
   int nchannels = fMTF->GetNChannels();

   // create matrix of number of bins
   std::vector< std::vector<double> > nbins_matrix;

   // prepare the tree variables
   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // get number of bins
      int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

      // create new matrix row
      std::vector<double> nbins_column(nbins);

      // add matrix
      nbins_matrix.push_back(nbins_column);
   }

   // get number of parameters
   int nparameters = fMTF->GetNParameters();

   std::vector<double> in_parameters(nparameters);

   // create branches
   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // get number of bins
      int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

      // loop over bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         // create branches
         tree->Branch(Form("channel_%i_bin_%i", ichannel, ibin),
                            &(nbins_matrix[ichannel])[ibin-1],
                            "n/D");
      }
   }

   for (int i = 0; i < nparameters; ++i) {
      tree->Branch(Form("parameter_%i", i), &in_parameters[i], Form("parameter_%i/D", i));
   }

   // create vector of histograms
   std::vector<TH1D> histograms;

   // loop over ensembles
   for (int iensemble = 0; iensemble < nensembles; ++iensemble) {
      // create ensembles
      histograms = BuildEnsemble(parameters);

      // copy information from histograms into tree variables
      // loop over channels
      for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
         // get channel
         BCMTFChannel * channel = fMTF->GetChannel(ichannel);

         // get number of bins
         int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

         // loop over bins
         for (int ibin = 1; ibin <= nbins; ++ibin) {
            // fill tree variable
            (nbins_matrix[ichannel])[ibin-1] = histograms.at(ichannel).GetBinContent(ibin);
         }
      }

      // copy parameter information
      for (int i = 0; i < nparameters; ++i) {
         in_parameters[i] = parameters.at(i);
      }

      // fill tree
      tree->Fill();
   }

   // return tree
   return tree;
}

// ---------------------------------------------------------
TTree * BCMTFAnalysisFacility::PerformEnsembleTest(std::vector<double> parameters, int nensembles)
{
   // create new tree
   TTree * tree = 0;

   // create ensembles
   tree = BuildEnsembles(parameters, nensembles);

   // perform ensemble test
   return PerformEnsembleTest(tree, nensembles);
}

// ---------------------------------------------------------
TTree * BCMTFAnalysisFacility::PerformEnsembleTest(TTree * tree, int nensembles, int start)
{
   // set log level to nothing
   BCLog::LogLevel ll = BCLog::GetLogLevelScreen();
   BCLog::SetLogLevel(BCLog::nothing);

   // get number of channels
   int nchannels = fMTF->GetNChannels();

   // define set of histograms for the original data sets
   std::vector<TH1D *> histograms_data(nchannels);

   // create matrix of number of bins
   std::vector< std::vector<double> > nbins_matrix;

   // prepare the tree
   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // get number of bins
      int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

      // create new matrix row
      std::vector<double> nbins_column(nbins);

      // add matrix
      nbins_matrix.push_back(nbins_column);
   }

   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // get number of bins
      int nbins = channel->GetData()->GetHistogram()->GetNbinsX();

      // loop over bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         // create branches
         tree->SetBranchAddress(Form("channel_%i_bin_%i", ichannel, ibin),
                                           &(nbins_matrix[ichannel])[ibin-1]);
      }
   }

   // get number of parameters
   int nparameters = fMTF->GetNParameters();

   // define tree variables
   std::vector<double> out_parameters(nparameters);

   for (int i = 0; i < nparameters; ++i) {
      tree->SetBranchAddress(Form("parameter_%i", i), &out_parameters[i]);
   }

   // copy the original data sets
   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // set data pointer
      histograms_data[ichannel] = channel->GetData()->GetHistogram();
   }

   // create output tree
   TTree * tree_out = new TTree("ensemble_test", "ensemble test");

   // define tree variables
   std::vector<double> out_mode_global(nparameters);
   std::vector<double> out_std_global(nparameters);
   std::vector<double> out_mode_marginalized(nparameters);
   std::vector<double> out_mean_marginalized(nparameters);
   std::vector<double> out_median_marginalized(nparameters);
   std::vector<double> out_5quantile_marginalized(nparameters);
   std::vector<double> out_10quantile_marginalized(nparameters);
   std::vector<double> out_16quantile_marginalized(nparameters);
   std::vector<double> out_84quantile_marginalized(nparameters);
   std::vector<double> out_90quantile_marginalized(nparameters);
   std::vector<double> out_95quantile_marginalized(nparameters);
   std::vector<double> out_std_marginalized(nparameters);
   std::vector<double> out_chi2_generated(nchannels);
   std::vector<double> out_chi2_mode(nchannels);
   std::vector<double> out_cash_generated(nchannels);
   std::vector<double> out_cash_mode(nchannels);
   std::vector<int> out_nevents(nchannels);
   double out_chi2_generated_total;
   double out_chi2_mode_total;
   double out_cash_generated_total;
   double out_cash_mode_total;
   int out_nevents_total;

   // create branches
   for (int i = 0; i < nparameters; ++i) {
      tree_out->Branch(Form("parameter_%i", i), &out_parameters[i], Form("parameter %i/D", i));
      tree_out->Branch(Form("mode_global_%i", i), &out_mode_global[i], Form("global mode of par. %i/D", i));
      tree_out->Branch(Form("std_global_%i", i), &out_std_global[i], Form("global std of par. %i/D", i));
      if (fFlagMCMC) {
         tree_out->Branch(Form("mode_marginalized_%i", i), &out_mode_marginalized[i], Form("marginalized mode of par. %i/D", i));
         tree_out->Branch(Form("mean_marginalized_%i", i), &out_mean_marginalized[i], Form("marginalized mean of par. %i/D", i));
         tree_out->Branch(Form("median_marginalized_%i", i), &out_median_marginalized[i], Form("median of par. %i/D", i));
         tree_out->Branch(Form("5quantile_marginalized_%i", i), &out_5quantile_marginalized[i], Form("marginalized 5 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("10quantile_marginalized_%i", i), &out_10quantile_marginalized[i], Form("marginalized 10 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("16quantile_marginalized_%i", i), &out_16quantile_marginalized[i], Form("marginalized 16 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("84quantile_marginalized_%i", i), &out_84quantile_marginalized[i], Form("marginalized 84 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("90quantile_marginalized_%i", i), &out_90quantile_marginalized[i], Form("marginalized 90 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("95quantile_marginalized_%i", i), &out_95quantile_marginalized[i], Form("marginalized 95 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("std_marginalized_%i", i), &out_std_marginalized[i], Form("marginalized std of par. %i/D", i));
      }
   }
   for (int i = 0; i < nchannels; ++i) {
      tree_out->Branch(Form("chi2_generated_%i", i), &out_chi2_generated[i], Form("chi2 (generated par.) in channel %i/D", i));
      tree_out->Branch(Form("chi2_mode_%i", i), &out_chi2_mode[i], Form("chi2 (mode of par.) in channel %i/D", i));
      tree_out->Branch(Form("cash_generated_%i", i), &out_cash_generated[i], Form("cash statistic (generated par.) in channel %i/D", i));
      tree_out->Branch(Form("cash_mode_%i", i), &out_cash_mode[i], Form("cash statistic (mode of par.) in channel %i/D", i));
      tree_out->Branch(Form("nevents_%i", i), &out_nevents[i], Form("nevents in channel %i/I", i));
   }
   tree_out->Branch("chi2_generated_total", &out_chi2_generated_total, "chi2 (generated par.) in all channels/D");
   tree_out->Branch("chi2_mode_total", &out_chi2_mode_total, "chi2 (mode of par.) in all channels/D");
   tree_out->Branch("cash_generated_total", &out_cash_generated_total, "cash statistic (generated par.) in all channels/D");
   tree_out->Branch("cash_mode_total", &out_cash_mode_total, "cash statistic (mode of par.) in all channels/D");
   tree_out->Branch("nevents_total", &out_nevents_total, "total number of events/I");

   // loop over ensembles
   for (int iensemble = 0; iensemble < nensembles; ++iensemble) {
      // print status
      if ((iensemble+1) % 100 == 0 && iensemble > 0)
         std::cout << "Fraction of ensembles analyzed: " << double(iensemble+1) / double(nensembles) * 100 << "%" << std::endl;

      // get next (commented out: random) event from tree
      //      int index = (int) fRandom->Uniform(tree->GetEntries());
      int index = iensemble + start;
      tree->GetEntry(index);

      // define histograms
      std::vector<TH1D> histograms;

      // transform matrix into histograms
      histograms = MatrixToHistograms(nbins_matrix);

      // loop over channels
      for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
         // get channel
         BCMTFChannel * channel = fMTF->GetChannel(ichannel);

         // overwrite data
         if (iensemble == 0)
            channel->GetData()->SetHistogram(0);

         // set data histogram
         fMTF->SetData(channel->GetName().c_str(), histograms.at(ichannel));
      }

      // check if MCMC should be run and perform analysis
      if (fFlagMCMC) {
         // work-around: force initialization
         fMTF->MCMCResetResults();

         // run mcmc
         fMTF->MarginalizeAll();

         // find mode
         fMTF->FindMode( fMTF->GetBestFitParameters() );
      }
      else {
         // find mode
         fMTF->FindMode();
      }

      // fill tree variables
      out_mode_global = fMTF->GetBestFitParameters();
      out_std_global = fMTF->GetBestFitParameterErrors();
      out_chi2_generated_total = 0;
      out_chi2_mode_total = 0;
      out_cash_generated_total = 0;
      out_cash_mode_total = 0;
      out_nevents_total = 0;

      for (int i = 0; i < nparameters; ++i) {
         if (fFlagMCMC) {
            BCH1D * hist = fMTF->GetMarginalized( fMTF->GetParameter(i) );
            out_mode_marginalized[i] = hist->GetMode();
            out_mean_marginalized[i] = hist->GetMean();
            out_median_marginalized[i] = hist->GetMedian();
            out_5quantile_marginalized[i] = hist->GetQuantile(0.05);
            out_10quantile_marginalized[i] = hist->GetQuantile(0.10);
            out_16quantile_marginalized[i] = hist->GetQuantile(0.16);
            out_84quantile_marginalized[i] = hist->GetQuantile(0.84);
            out_90quantile_marginalized[i] = hist->GetQuantile(0.90);
            out_95quantile_marginalized[i] = hist->GetQuantile(0.95);
            out_std_marginalized[i]=hist->GetRMS();
         }
      }

      for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
         // get channel
         BCMTFChannel * channel = fMTF->GetChannel(ichannel);

         // get number of events
         out_nevents[ichannel] = (int) channel->GetData()->GetHistogram()->Integral();

         // calculate chi2
         out_chi2_generated[ichannel] = fMTF->CalculateChi2( ichannel, out_parameters );
         out_chi2_mode[ichannel] = fMTF->CalculateChi2( ichannel, fMTF->GetBestFitParameters() );

         // calculate cash statistic
         out_cash_generated[ichannel] = fMTF->CalculateCash( ichannel, out_parameters );
         out_cash_mode[ichannel] = fMTF->CalculateCash( ichannel, fMTF->GetBestFitParameters() );

         // increase the total number of events
         out_nevents_total += out_nevents[ichannel];

         // increase the total chi2
         out_chi2_generated_total += out_chi2_generated[ichannel];
         out_chi2_mode_total += out_chi2_mode[ichannel];

         // increase the total cash statistic
         out_cash_generated_total += out_cash_generated[ichannel];
         out_cash_mode_total += out_cash_mode[ichannel];
      }

      // fill tree
      tree_out->Fill();
   }

   // put the original data back in place
   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // set data pointer
      channel->GetData()->SetHistogram(histograms_data.at(ichannel));
   }

   // reset log level
   BCLog::SetLogLevel(ll);

   // work-around: force initialization
   fMTF->MCMCResetResults();

   // return output tree
   return tree_out;
}

// ---------------------------------------------------------
std::vector<TH1D> BCMTFAnalysisFacility::MatrixToHistograms(std::vector< std::vector<double> > matrix)
{
   // create vector of histograms
   std::vector<TH1D> histograms;

   // get number of channels
   int nchannels = matrix.size();;

   // loop over channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // get column
      std::vector<double> nbins_column = matrix[ichannel];

      // create new histogram
      TH1D hist( *(channel->GetData()->GetHistogram()) );

      // get number of bins
      int nbins = hist.GetNbinsX();

      // fill bin content
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         hist.SetBinContent(ibin, nbins_column.at(ibin-1));
      }

      // add histogram to container
      histograms.push_back(hist);
   }

   // return histograms
   return histograms;

}

// ---------------------------------------------------------
int BCMTFAnalysisFacility::PerformSingleChannelAnalyses(const char * dirname, const char * options)
{
   // ---- create new directory ---- //

   mkdir(dirname, 0777);
   chdir(dirname);

   // ---- check options ---- //

   bool flag_syst = true;
   bool flag_mcmc = true;

   if (std::string(options).find("nosyst") < std::string(options).size())
      flag_syst = false;

   if (std::string(options).find("mcmc") < std::string(options).size())
      flag_mcmc = true;

   // get number of channels
   int nchannels = fMTF->GetNChannels();

   // get number of systematics
   int nsystematics = fMTF->GetNSystematics();

   // container of flags for channels
   std::vector<bool> flag_channel(nchannels);

   // container of flags for systematics
   std::vector<bool> flag_systematic(nsystematics);

   // create new container of comparison tools
   std::vector<BCMTFComparisonTool *> ctc;
   std::vector<BCMTFComparisonTool *> ctc_mcmc;

   // get number of parameters
   int nparameters = fMTF->GetNParameters();

   // ---- add one comparison tool for each parameter ---- //

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // create new comparison tool
      BCMTFComparisonTool * ct = new BCMTFComparisonTool(fMTF->GetParameter(i)->GetName().c_str());
      BCMTFComparisonTool * ct_mcmc = new BCMTFComparisonTool(fMTF->GetParameter(i)->GetName().c_str());

      // add comparison tool to container
      ctc.push_back(ct);
      ctc_mcmc.push_back(ct_mcmc);
   }

   // ---- switch on/off all systematics ---- //

   // loop over all systematics
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      BCMTFSystematic * systematic = fMTF->GetSystematic(isystematic);

      // remember old setting
      flag_systematic[isystematic] = systematic->GetFlagSystematicActive();

      // switch off systematic
      if (flag_systematic[isystematic])
         systematic->SetFlagSystematicActive(flag_syst);
   }

   // ---- switch on all channels ---- //

   // loop over all channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // remember old setting
      flag_channel[ichannel] = channel->GetFlagChannelActive();

      // switch channel on/off
      channel->SetFlagChannelActive(true);
   }

   // ---- perform analysis for combination ---- //
   if (flag_mcmc) {
      // work-around: force initialization
      fMTF->MCMCResetResults();

      // run mcmc
      fMTF->MarginalizeAll();

      // find mode
      fMTF->FindMode( fMTF->GetBestFitParameters() );
   }
   else {
      // find mode
      fMTF->FindMode();
   }

   // print results
   if (flag_mcmc)
      fMTF->PrintAllMarginalized("marginalized_combined.ps");
   fMTF->PrintResults("results_combined.txt");

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // get comparison tool
      BCMTFComparisonTool * ct = ctc.at(i);
      BCMTFComparisonTool * ct_mcmc = ctc_mcmc.at(i);

      ct->AddContribution("all channels",
                                    fMTF->GetBestFitParameters().at(i),
                                    fMTF->GetBestFitParameterErrors().at(i));
      if (flag_mcmc) {
         BCH1D * hist = fMTF->GetMarginalized( fMTF->GetParameter(i) );

         ct_mcmc->AddContribution("all channels",
                                              hist->GetMean(),
                                              hist->GetRMS());
      }
   }

   // ---- switch off all channels ----//

   // loop over all channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // switch off channel
      channel->SetFlagChannelActive(false);
   }

   // ---- perform analysis on all channels separately ---- //

   // loop over all channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // switch on one channel
      channel->SetFlagChannelActive(true);

      // perform analysis

      if (flag_mcmc) {
         // work-around: force initialization
         fMTF->MCMCResetResults();

         // run mcmc
         fMTF->MarginalizeAll();

         // find mode
         fMTF->FindMode( fMTF->GetBestFitParameters() );
      }
      else {
         // find mode
         fMTF->FindMode();
      }

      // print results
      if (flag_mcmc)
         fMTF->PrintAllMarginalized(Form("marginalized_channel_%i.ps", ichannel));
      fMTF->PrintResults(Form("results_channel_%i.txt", ichannel));

      // ---- update comparison tools ---- //

      // loop over all parameters
      for (int i = 0; i < nparameters; ++ i) {
         // get comparison tool
         BCMTFComparisonTool * ct = ctc.at(i);
         BCMTFComparisonTool * ct_mcmc = ctc_mcmc.at(i);

         ct->AddContribution(channel->GetName().c_str(),
                                       fMTF->GetBestFitParameters().at(i),
                                       fMTF->GetBestFitParameterErrors().at(i));
         if (flag_mcmc) {
            BCH1D * hist = fMTF->GetMarginalized( fMTF->GetParameter(i) );

            ct_mcmc->AddContribution(channel->GetName().c_str(),
                                                 hist->GetMean(),
                                                 hist->GetRMS());
         }

         // switch off channel
         channel->SetFlagChannelActive(false);
      }
   }

   // ---- reset all systematics ---- //
   // loop over all systematics
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      BCMTFSystematic * systematic = fMTF->GetSystematic(isystematic);

      // switch off systematic
      if (flag_systematic[isystematic])
         systematic->SetFlagSystematicActive(flag_systematic[isystematic]);
   }

   // ---- reset all channels ---- //

   // loop over all channels
   for (int ichannel = 0; ichannel < nchannels; ++ichannel) {
      // get channel
      BCMTFChannel * channel = fMTF->GetChannel(ichannel);

      // switch channel on/off
      channel->SetFlagChannelActive(flag_channel[ichannel]);
   }

   // ---- workaround: reset MCMC ---- //
   fMTF->MCMCResetResults();

   // ---- print everything ---- //
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // draw first one
   BCMTFComparisonTool * ct =  ctc.at(0);
   ct->DrawOverview();
   //   c1->Print((std::string(filename)+std::string("(")).c_str());
   c1->Print((std::string("overview.ps")+std::string("(")).c_str());

   // loop over all parameters
   for (int i = 1; i < nparameters-1; ++ i) {
      // create new comparison tool
      BCMTFComparisonTool * ct = ctc.at(i);

      ct->DrawOverview();
      c1->Print("overview.ps");
   }

   // draw last one
   ct =  ctc.at(nparameters-1);
   ct->DrawOverview();
   c1->Print((std::string("overview.ps")+std::string(")")).c_str());

   // ---- print everything (mcmc) ---- //
   if (flag_mcmc) {
      TCanvas * c2 = new TCanvas();
      c2->cd();

      // draw first one
      BCMTFComparisonTool * ct_mcmc =  ctc_mcmc.at(0);
      ct_mcmc->DrawOverview();
      //   c2->Print((std::string(filename)+std::string("(")).c_str());
      c2->Print((std::string("overview_mcmc.ps")+std::string("(")).c_str());

      // loop over all parameters
      for (int i = 1; i < nparameters-1; ++ i) {
         // create new comparison tool
         BCMTFComparisonTool * ct_mcmc = ctc_mcmc.at(i);

         ct_mcmc->DrawOverview();
         c2->Print("overview_mcmc.ps");
      }

      // draw last one
      ct_mcmc =  ctc_mcmc.at(nparameters-1);
      ct_mcmc->DrawOverview();
      c2->Print((std::string("overview_mcmc.ps")+std::string(")")).c_str());
      delete c2;
   }

   // ---- free memory ---- //
   for (int i = 0; i < nparameters; ++i) {
      BCMTFComparisonTool * ct = ctc[i];
      BCMTFComparisonTool * ct_mcmc = ctc_mcmc[i];
      delete ct;
      delete ct_mcmc;
   }
   ctc.clear();
   ctc_mcmc.clear();

   delete c1;

   // ---- change directory ---- //

   chdir("../");

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTFAnalysisFacility::PerformSingleSystematicAnalyses(const char * dirname, const char * options)
{
   // ---- create new directory ---- //

   mkdir(dirname, 0777);
   chdir(dirname);

   // ---- check options ---- //

   bool flag_mcmc = true;

   if (std::string(options).find("mcmc") < std::string(options).size())
      flag_mcmc = true;

   // get number of channels
   int nchannels = fMTF->GetNChannels();

   // get number of systematics
   int nsystematics = fMTF->GetNSystematics();

   // container of flags for channels
   std::vector<bool> flag_channel(nchannels);

   // container of flags for systematics
   std::vector<bool> flag_systematic(nsystematics);

   // create new container of comparison tools
   std::vector<BCMTFComparisonTool *> ctc;

   // get number of parameters
   int nparameters = fMTF->GetNParameters();

   // ---- add one comparison tool for each systematic ---- //

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // create new comparison tool
      BCMTFComparisonTool * ct = new BCMTFComparisonTool(fMTF->GetParameter(i)->GetName().c_str());

      // add comparison tool to container
      ctc.push_back(ct);
   }

   // ---- switch on all systematics ---- //

   for (int isystematic = 0; isystematic < nsystematics; ++ isystematic) {
      // get systematic
      BCMTFSystematic * systematic = fMTF->GetSystematic(isystematic);

      // remember old setting
      flag_systematic[isystematic] = systematic->GetFlagSystematicActive();

      // switch on
      systematic->SetFlagSystematicActive(true);
   }

   if (flag_mcmc) {
      // work-around: force initialization
      fMTF->MCMCResetResults();

      // run mcmc
      fMTF->MarginalizeAll();

      // find mode
      fMTF->FindMode( fMTF->GetBestFitParameters() );
   }
   else {
      // find mode
      fMTF->FindMode();
   }

   // print results
   if (flag_mcmc)
      fMTF->PrintAllMarginalized("marginalized_all.ps");
   fMTF->PrintResults("results_all.txt");

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // get comparison tool
      BCMTFComparisonTool * ct = ctc.at(i);

      ct->AddContribution("all systematics",
                                    fMTF->GetBestFitParameters().at(i),
                                    fMTF->GetBestFitParameterErrors().at(i));
   }

   // ---- switch off all systematics ---- //

   // loop over all systematics
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      BCMTFSystematic * systematic = fMTF->GetSystematic(isystematic);

      // switch off
      systematic->SetFlagSystematicActive(false);
   }

   // ---- perform analysis with all systematics separately ---- //

   // loop over all channels
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      BCMTFSystematic * systematic = fMTF->GetSystematic(isystematic);

      // switch on systematic
      systematic->SetFlagSystematicActive(true);

      // perform analysis
      if (flag_mcmc) {
         // work-around: force initialization
         fMTF->MCMCResetResults();

         // run mcmc
         fMTF->MarginalizeAll();

         // find mode
         fMTF->FindMode( fMTF->GetBestFitParameters() );
      }
      else {
         // find mode
         fMTF->FindMode();
      }

      // print results
      if (flag_mcmc)
         fMTF->PrintAllMarginalized(Form("marginalized_systematic_%i.ps", isystematic));
      fMTF->PrintResults(Form("results_systematic_%i.txt", isystematic));

      // ---- update comparison tools ---- //

      // loop over all parameters
      for (int i = 0; i < nparameters; ++ i) {
         // get comparison tool
         BCMTFComparisonTool * ct = ctc.at(i);

         ct->AddContribution(systematic->GetName().c_str(),
                                       fMTF->GetBestFitParameters().at(i),
                                       fMTF->GetBestFitParameterErrors().at(i));
      }

      // switch off systematic
      systematic->SetFlagSystematicActive(false);
   }

   // ---- analysis without any systematic uncertainty ---- //

   if (flag_mcmc) {
      // work-around: force initialization
      fMTF->MCMCResetResults();

      // run mcmc
      fMTF->MarginalizeAll();

      // find mode
      fMTF->FindMode( fMTF->GetBestFitParameters() );
   }
   else {
      // find mode
      fMTF->FindMode();
   }

   // print results
   if (flag_mcmc)
      fMTF->PrintAllMarginalized("marginalized_none.ps");
   fMTF->PrintResults("results_none.txt");

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // get comparison tool
      BCMTFComparisonTool * ct = ctc.at(i);

      ct->AddContribution("no systematics",
                                    fMTF->GetBestFitParameters().at(i),
                                    fMTF->GetBestFitParameterErrors().at(i));
   }



   // ---- reset all systematics ---- //

   // loop over all systematics
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      BCMTFSystematic * systematic = fMTF->GetSystematic(isystematic);

      // switch off systematic
      if (flag_systematic[isystematic])
         systematic->SetFlagSystematicActive(flag_systematic[isystematic]);
   }

   // ---- workaround: reset MCMC ---- //
   fMTF->MCMCResetResults();

   // ---- print everything ---- //
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // draw first one
   BCMTFComparisonTool * ct =  ctc.at(0);
   ct->DrawOverview();
   //   c1->Print((std::string(filename)+std::string("(")).c_str());
   c1->Print((std::string("overview.ps")+std::string("(")).c_str());

   // loop over all parameters
   for (int i = 1; i < nparameters-1; ++i) {
      // create new comparison tool
      BCMTFComparisonTool * ct = ctc.at(i);

      ct->DrawOverview();
      c1->Print("overview.ps");
   }

   // draw last one
   ct =  ctc.at(nparameters-1);
   ct->DrawOverview();
   c1->Print((std::string("overview.ps")+std::string(")")).c_str());

   // ---- free memory ---- //
   for (int i = 0; i < nparameters; ++i) {
      BCMTFComparisonTool * ct = ctc[i];
      delete ct;
   }
   ctc.clear();

   delete c1;

   // ---- change directory ---- //

   chdir("../");

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCMTFAnalysisFacility::PerformCalibrationAnalysis(const char * dirname, std::vector<double> default_parameters, int index, std::vector<double> parametervalues, int nensembles)
{

   // ---- create new directory ---- //

   mkdir(dirname, 0777);
   chdir(dirname);

   // ---- loop over parameter values and perform analysis  ---- //

   int nvalues = int(parametervalues.size());
   for (int ivalue = 0; ivalue < nvalues; ++ivalue) {

      // open file
      TFile * file = new TFile(Form("ensemble_%i.root", ivalue), "RECREATE");
      file->cd();

      // set parameters
      std::vector<double> parameters = default_parameters;
      parameters[index] = parametervalues.at(ivalue);

      // create ensemble
      TTree * tree = PerformEnsembleTest(parameters, nensembles);

      // write tree
      tree->Write();

      // close file
      file->Close();

      // free memory
      delete file;
   }

   // ---- change directory ---- //

   chdir("../");

   // no error
   return 1;
}

// ---------------------------------------------------------
