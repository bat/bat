/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPPDiagnostics.h"

#include <cmath>

#include <TROOT.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>

#include <BAT/BCMath.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
BCPPDiagnostics::BCPPDiagnostics() : BCPostProcessor()
                                   , fBatchLength(1000)
{
}

// ---------------------------------------------------------
BCPPDiagnostics::~BCPPDiagnostics()
{
}

// ---------------------------------------------------------
int BCPPDiagnostics::GetNBatches()
{
  return int(fNSamplesMainRun / fBatchLength);
}

// ---------------------------------------------------------
void BCPPDiagnostics::PrintLogProbability(std::string filename, std::string options)
{
  // check if file extension does not exist or is not pdf or ps
  if ( (filename.find_last_of(".") == std::string::npos) or
       ((filename.substr(filename.find_last_of(".")+1) != "pdf") and
        (filename.substr(filename.find_last_of(".")+1) != "ps"))) {
    // make it a PDF file
    filename += ".pdf";
  }

  // create canvas
  TCanvas* c1 = new TCanvas("c1");
  c1->cd();

  // option flags
  bool flag_all = false;
  bool flag_sum = false;

  // check if options is all
  if (options.find("all") < options.size()) {
    flag_all = true;
  }

  // check if options is all
  if (options.find("sum") < options.size()) {
    flag_sum = true;
  }

  // define histograms
  TH1D* hist_sum = new TH1D();
  std::vector<TH1D*> hist_all;

  // loop over all chains and fill histograms
  for (int i = 0; i < fNTrees; ++i) {
    TTree* tree = fTrees.at(i);

    TH1D* temphist = new TH1D(Form("hist_logprob%i", i), "", 100, fLogProbabilityMin - 0.1*(fLogProbabilityMax - fLogProbabilityMin), fLogProbabilityMax + 0.1*(fLogProbabilityMax - fLogProbabilityMin));
    temphist->GetXaxis()->SetTitle("log(probability)");
    temphist->GetYaxis()->SetTitle("dN/dlog(probability)");
    temphist->SetStats(kFALSE);
    temphist->SetLineColor(1+i);

    tree->Draw(Form("LogProbability>>hist_logprob%i", i), "Phase==2");

    hist_all.push_back(temphist);

    if (i == 0)
      hist_sum = (TH1D*) temphist->Clone();
    else
      hist_sum->Add(temphist);
  }

  // draw histograms
  if (flag_all) {
    if (flag_sum)
      hist_sum->Draw();
    for (int i = 0; i < fNTrees; ++i) {
      if (i > 0)
        hist_all.at(i)->Draw("SAME");
      else {
        if (!flag_sum)
          hist_all.at(i)->Draw();
      }
    }
  }
  else if (flag_sum) {
    hist_sum->Draw();
  }

  // print
  c1->Print(filename.c_str());

  // free memory
  delete c1;
  delete hist_sum;
  for (int i = 0; i < fNTrees; ++i)
    delete hist_all[i];
  hist_all.clear();
}

// ---------------------------------------------------------
void BCPPDiagnostics::CalculateBatchQuantities()
{
  // clear old entries
  fParameterBatchMean.clear();
  fParameterBatchVariance.clear();
  fParameterBatchStdDev.clear();

  // loop over parameters
  for (int ipar = 0; ipar < GetNParameters(); ++ipar) {

    // loop over chains
    for (int ichain = 0; ichain < GetNChains(); ++ichain) {

      // loop over batches
      for (int ibatch = 0; ibatch < GetNBatches(); ++ibatch) {

        // start and stop numbers
        int start = ibatch * fBatchLength;
        int stop  = (ibatch + 1) * fBatchLength;

        // helper variables
        double sum = 0;
        double sum2 = 0;

        // loop over samples and calculate helper quantities
        for (int i = start; i < stop; ++i) {
          double temp = GetParameterValue(ipar, ichain, i, false);

          sum += temp;
          sum2 += temp*temp;
        }

        // calculate batch quantities
        double mean     = sum / double(stop-start+1);
        double variance = sum2 / double(stop-start+1) - mean*mean;
        double stddev    = sqrt(variance);

        // write quantities
        fParameterBatchMean.push_back(mean);
        fParameterBatchVariance.push_back(variance);
        fParameterBatchStdDev.push_back(stddev);
      }
    }
  }
}

// ---------------------------------------------------------
void BCPPDiagnostics::PrintBatchQuantities(std::string filename)
{
  // check if file extension does not exist or is not pdf or ps
  if ( (filename.find_last_of(".") == std::string::npos) or
       ((filename.substr(filename.find_last_of(".")+1) != "pdf") and
        (filename.substr(filename.find_last_of(".")+1) != "ps"))) {
    // make it a PDF file
    filename += ".pdf";
  }

  // helper
  int npar = GetNParameters();
  int nbatches = GetNBatches();
  int nchains = GetNChains();

  // create canvas
  TCanvas* c1 = new TCanvas("c1", "", 900., 300.);
  c1->SetLeftMargin(0.06);
  c1->SetRightMargin(0.02);
  c1->cd();

  // --------------------------------------
  // Draw parameter mean values
  // --------------------------------------

  // loop over parameters
  for (int ipar = 0; ipar < npar; ++ipar) {

    // define the graphs
    std::vector<TGraph*> graphs(fNTrees);

    // minimum and maximum
    double ymin = 0;
    double ymax = 0;

    // loop over chains
    for (int ichain = 0; ichain < nchains; ++ichain) {
      graphs[ichain] = new TGraph(nbatches);
      graphs[ichain]->SetLineColor(1+ichain);

      // loop over batches
      for (int ibatch = 0; ibatch < nbatches; ++ibatch) {
        double value = fParameterBatchMean[ibatch + ichain*nbatches + ipar*nchains*nbatches];
        graphs[ichain]->SetPoint(ibatch, ibatch, value);

        // find minimum and maximum
        if (ichain == 0 && ibatch == 0) {
          ymin = value;
          ymax = value;
        }
        if (value < ymin)
          ymin = value;
        if (value > ymax)
          ymax = value;
      }
    }

    // draw axes
    TH2D* hist_axes = new TH2D(Form("hist_axes_%i", BCLog::GetHIndex()), Form(";Batch number; Mean value of parameter %i", ipar), nbatches, -5.5, nbatches+4.5, 1, ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin));
    hist_axes->SetTitleOffset(0.7, "Y");
    hist_axes->SetStats(kFALSE);
    hist_axes->Draw();

    // draw graphs
    for (int ichain = 0; ichain < nchains; ++ichain) {
      graphs[ichain]->Draw("L");
    }
    // print
    if (ipar == 0)
      c1->Print(Form("%s(", filename.c_str()));
    else
      c1->Print(filename.c_str());

    // free memory
    delete hist_axes;
    for (int ichain = 0; ichain < nchains; ++ichain)
      delete graphs[ichain];
    graphs.clear();
  }

  // --------------------------------------
  // Draw parameter standard deviation
  // --------------------------------------

  // loop over parameters
  for (int ipar = 0; ipar < npar; ++ipar) {

    // define the graphs
    std::vector<TGraph*> graphs(fNTrees);

    // minimum and maximum
    double ymin = 0;
    double ymax = 0;

    // loop over chains
    for (int ichain = 0; ichain < nchains; ++ichain) {
      graphs[ichain] = new TGraph(nbatches);
      graphs[ichain]->SetLineColor(1+ichain);

      // loop over batches
      for (int ibatch = 0; ibatch < nbatches; ++ibatch) {
        double value = fParameterBatchStdDev[ibatch + ichain*nbatches + ipar*nchains*nbatches];
        graphs[ichain]->SetPoint(ibatch, ibatch, value);

        // find minimum and maximum
        if (ichain == 0 && ibatch == 0) {
          ymin = value;
          ymax = value;
        }
        if (value < ymin)
          ymin = value;
        if (value > ymax)
          ymax = value;
      }
    }

    // draw axes
    TH2D* hist_axes = new TH2D(Form("hist_axes_%i", BCLog::GetHIndex()), Form(";Batch number; Standard deviation of parameter %i", ipar), nbatches, -5.5, nbatches+4.5, 1, ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin));
    hist_axes->SetTitleOffset(0.7, "Y");
    hist_axes->SetStats(kFALSE);
    hist_axes->Draw();

    // draw graphs
    for (int ichain = 0; ichain < nchains; ++ichain) {
      graphs[ichain]->Draw("L");
    }
    // print
    c1->Print(filename.c_str());

    // free memory
    delete hist_axes;
    for (int ichain = 0; ichain < nchains; ++ichain)
      delete graphs[ichain];
    graphs.clear();
  }

  // --------------------------------------
  // Draw parameter R values
  // --------------------------------------

  // define the graphs
  std::vector<TGraph*> graphs_r(npar);

  // minimum and maximum
  double ymin = 0;
  double ymax = 0;

  // loop over parameters
  for (int ipar = 0; ipar < npar; ++ipar) {

    // define graphs
    graphs_r[ipar] = new TGraph(nbatches);
    graphs_r[ipar]->SetLineColor(1+ipar);

    // vector or mean values and variances
    std::vector<double> meanvalues;
    std::vector<double> variances;

    // loop over batches
    for (int ibatch = 0; ibatch < nbatches; ++ibatch) {

      // loop over chains
      for (int ichain = 0; ichain < nchains; ++ichain) {
        meanvalues.push_back(fParameterBatchMean[ibatch + ichain*nbatches + ipar*nchains*nbatches]);
        variances.push_back(fParameterBatchVariance[ibatch + ichain*nbatches + ipar*nchains*nbatches]);
      }

      // calculate R value
      double value = BCMath::Rvalue(meanvalues, variances, fBatchLength, true);

      graphs_r[ipar]->SetPoint(ibatch, ibatch, value);

      // find minimum and maximum
      if (ipar == 0 && ibatch == 0) {
        ymin = value;
        ymax = value;
      }
      if (value < ymin)
        ymin = value;
      if (value > ymax)
        ymax = value;
    }
  }

  // draw axes
  TH2D* hist_axes_r = new TH2D(Form("hist_axes_%i", BCLog::GetHIndex()), ";Batch number; R values", nbatches, -5.5, nbatches+4.5, 1, ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin));
  hist_axes_r->SetTitleOffset(0.7, "Y");
  hist_axes_r->SetStats(kFALSE);
  hist_axes_r->Draw();

  // draw graphs
  for (int ipar = 0; ipar < npar; ++ipar) {
    graphs_r[ipar]->Draw("L");
  }

  // print
  c1->Print(Form("%s)", filename.c_str()));

  // free memory
  delete hist_axes_r;
  for (int ipar = 0; ipar < npar; ++ipar) {
    delete graphs_r[ipar];
  }
  graphs_r.clear();


  // free memory
  delete c1;
}

// ---------------------------------------------------------
