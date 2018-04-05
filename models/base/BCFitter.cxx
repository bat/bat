/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include <config.h>

#include "BCFitter.h"

#include <BAT/BCDataPoint.h>
#include <BAT/BCH1D.h>
#include <BAT/BCLog.h>
#include <BAT/BCModel.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>

#include <stdexcept>

// ---------------------------------------------------------
BCFitter::BCFitter(const TF1& f, const std::string& name)
    : BCModel(name),
      fFitFunction(1, f),
      fFlagFillErrorBand(true),
      fFitFunctionIndexX(-1),
      fFitFunctionIndexY(-1),
      fErrorBandContinuous(true),
      fErrorBandNbinsX(100),
      fErrorBandNbinsY(500),
      fErrorBandExtensionLowEdgeX(0),
      fErrorBandExtensionUpEdgeX(0),
      fErrorBandExtensionLowEdgeY(0),
      fErrorBandExtensionUpEdgeY(0)
{
    if (f.GetNdim() != 1)
        throw std::invalid_argument("BCFitter and descendants only support 1D problems");

    // set name if name empty
    if (name.empty())
        SetName(std::string("fitter_model_") + fFitFunction.front().GetName());

    // fill parameter set
    fParameters = BCParameterSet();
    for (int i = 0; i < fFitFunction.front().GetNpar(); ++i) {
        double xmin;
        double xmax;

        fFitFunction.front().GetParLimits(i, xmin, xmax);

        AddParameter(fFitFunction.front().GetParName(i), xmin, xmax);
    }

    // create a data set. this is necessary in order to calculate the
    // error band.
    SetDataSet(&fFitterDataSet);

    // by default set all priors constant
    fParameters.SetPriorConstantAll();
}

// ---------------------------------------------------------
BCFitter::~BCFitter()
{
}

// ---------------------------------------------------------
void BCFitter::MCMCUserInitialize()
{
    // add or remove copies
    fFitFunction.resize(fMCMCNChains, fFitFunction.front());
}

// ---------------------------------------------------------
void BCFitter::MCMCUserIterationInterface()
{
    // fill error band
    if (fFlagFillErrorBand)
        FillErrorBand();
}

// ---------------------------------------------------------
void BCFitter::MarginalizePreprocess()
{
    // prepare function fitting
    double dx = 0.;
    double dy = 0.;

    if (GetDataSet() && fFitFunctionIndexX >= 0 && fFitFunctionIndexY >= 0) {

        dx = (GetDataSet()->GetRangeWidth(fFitFunctionIndexX)
              + fErrorBandExtensionLowEdgeX
              + fErrorBandExtensionUpEdgeX)
             / fErrorBandNbinsX;


        dy = (GetDataSet()->GetRangeWidth(fFitFunctionIndexY)
              + fErrorBandExtensionLowEdgeY
              + fErrorBandExtensionUpEdgeY)  / fErrorBandNbinsY;

        double xRangeLow = GetDataSet()->GetLowerBound(fFitFunctionIndexX) - fErrorBandExtensionLowEdgeX - dx / 2;
        double xRangeHigh = GetDataSet()->GetUpperBound(fFitFunctionIndexX) + fErrorBandExtensionUpEdgeX + dx / 2;

        double yRangeLow = GetDataSet()->GetLowerBound(fFitFunctionIndexY) - fErrorBandExtensionLowEdgeY - dy / 2;
        double yRangeHigh = GetDataSet()->GetUpperBound(fFitFunctionIndexY) + fErrorBandExtensionUpEdgeY + dy / 2;

        fErrorBandXY = TH2D(TString::Format("errorbandxy_%s", GetSafeName().data()), "",
                            fErrorBandNbinsX,
                            xRangeLow,
                            xRangeHigh,
                            fErrorBandNbinsY,
                            yRangeLow,
                            yRangeHigh);
        fErrorBandXY.SetStats(kFALSE);

        // why are we doing this?
        for (unsigned ix = 1; ix <= fErrorBandNbinsX; ++ix)
            for (unsigned iy = 1; iy <= fErrorBandNbinsX; ++iy)
                fErrorBandXY.SetBinContent(ix, iy, 0.);
    }

}

// ---------------------------------------------------------
void BCFitter::FillErrorBand()
{
    // function fitting
    if (fFitFunctionIndexX < 0)
        return;

    // loop over all possible x values ...
    if (fErrorBandContinuous) {
        double x = 0;
        for (unsigned ix = 1; ix <= fErrorBandNbinsX; ++ix) {
            // calculate x
            x = fErrorBandXY.GetXaxis()->GetBinCenter(ix);

            // calculate y
            std::vector<double> xvec;
            xvec.push_back(x);

            // loop over all chains
            for (unsigned ichain = 0; ichain < GetNChains(); ++ichain) {
                // calculate y
                double y = FitFunction(xvec, Getx(ichain));

                // fill histogram
                fErrorBandXY.Fill(x, y);
            }

            xvec.clear();
        }
    }
    // ... or evaluate at the data point x-values
    else {
        unsigned ndatapoints = fErrorBandX.size();
        double x = 0;

        for (unsigned ix = 0; ix < ndatapoints; ++ix) {
            // calculate x
            x = fErrorBandX.at(ix);

            // calculate y
            std::vector<double> xvec;
            xvec.push_back(x);

            // loop over all chains
            for (unsigned ichain = 0; ichain < GetNChains(); ++ichain) {
                // calculate y
                double y = FitFunction(xvec, Getx(ichain));

                // fill histogram
                fErrorBandXY.Fill(x, y);
            }

            xvec.clear();
        }
    }
}

// ---------------------------------------------------------
double BCFitter::FitFunction(const std::vector<double>& x, const std::vector<double>& params)
{
    // update parameters in right TF1 and evaluate
    TF1& f = GetFitFunction();
    f.SetParameters(&params[0]);
    return f.Eval(x[0]);
}

double BCFitter::Integral(const std::vector<double>& params, const double xmin, const double xmax)
{
    TF1& f = GetFitFunction();

    f.SetParameters(&params[0]);

    // use ROOT's TH1::Integral method
    if (fFlagIntegration) {
        return f.Integral(xmin, xmax);
    } // use linear interpolation
    else {
        return (f.Eval(xmin) + f.Eval(xmax)) * (xmax - xmin) / 2.;
    }
}

// ---------------------------------------------------------
void BCFitter::PrintShortFitSummary()
{
    BCModel::PrintShortFitSummary();
    if (GetPValue() >= 0) {
        BCLog::OutSummary("   Goodness-of-fit test:");
        BCLog::OutSummary(Form("      p-value = %.3g", GetPValue()));
        BCLog::OutSummary("---------------------------------------------------");
    }
}

// ---------------------------------------------------------
std::vector<double> BCFitter::GetErrorBand(double level) const
{
    std::vector<double> errorband;

    int nx = fErrorBandXY.GetNbinsX();
    errorband.assign(nx, 0.);

    // loop over x and y bins
    for (int ix = 1; ix <= nx; ix++) {
        TH1D* temphist = fErrorBandXY.ProjectionY("temphist", ix, ix);

        int nprobSum = 1;
        double q[1];
        double probSum[1];
        probSum[0] = level;

        temphist->GetQuantiles(nprobSum, q, probSum);

        errorband[ix - 1] = q[0];
        delete temphist;
    }

    return errorband;
}

// ---------------------------------------------------------
TGraph* BCFitter::GetErrorBandGraph(double level1, double level2) const
{
    // define new graph
    int nx = fErrorBandXY.GetNbinsX();

    TGraph* graph = new TGraph(2 * nx);
    graph->SetFillStyle(1001);
    graph->SetFillColor(kYellow);
    graph->SetLineWidth(0);

    // get error bands
    std::vector<double> ymin = GetErrorBand(level1);
    std::vector<double> ymax = GetErrorBand(level2);

    for (int i = 0; i < nx; i++) {
        graph->SetPoint(i, fErrorBandXY.GetXaxis()->GetBinCenter(i + 1),
                        ymin[i] * GraphCorrection(i + 1));
        graph->SetPoint(nx + i, fErrorBandXY.GetXaxis()->GetBinCenter(nx - i),
                        ymax[nx - i - 1] * GraphCorrection(i + 1));
    }

    return graph;
}

// ---------------------------------------------------------
TH2* BCFitter::GetGraphicalErrorBandXY(double level, int nsmooth, bool overcoverage) const
{
    int nx = fErrorBandXY.GetNbinsX();
    int ny = fErrorBandXY.GetNbinsY();

    // copy existing histogram
    TH2* hist_tempxy = (TH2*) fErrorBandXY.Clone(TString::Format("%s_sub_%f.2", fErrorBandXY.GetName(), level));
    hist_tempxy->Reset();
    hist_tempxy->SetFillColor(kYellow);

    // loop over x bins
    for (int ix = 1; ix < nx; ix++) {
        BCH1D hist_temp(fErrorBandXY.ProjectionY("temphist", ix, ix));
        if (nsmooth > 0)
            hist_temp.Smooth(nsmooth);
        std::vector<std::pair<double, double> > bound = hist_temp.GetSmallestIntervalBounds(std::vector<double>(1, level), overcoverage);
        for (int iy = 1; iy <= ny; ++iy)
            if (hist_temp.GetHistogram()->GetBinContent(iy) >= bound.front().first)
                hist_tempxy->SetBinContent(ix, iy, 1);
    }

    return hist_tempxy;
}


// ---------------------------------------------------------
TGraph* BCFitter::GetFitFunctionGraph(const std::vector<double>& parameters)
{
    // define new graph
    const int nx = fErrorBandXY.GetNbinsX();
    TGraph* graph = new TGraph(nx);

    std::vector<double> xvec(1);

    // loop over x values
    for (int i = 0; i < nx; i++) {
        xvec[0] = fErrorBandXY.GetXaxis()->GetBinCenter(i + 1);
        const double y = FitFunction(xvec, parameters) * GraphCorrection(i + 1);

        graph->SetPoint(i, xvec[0], y);
    }

    return graph;
}

// ---------------------------------------------------------
TGraph* BCFitter::GetFitFunctionGraph(const std::vector<double>& parameters, double xmin, double xmax, int n)
{
    // define new graph
    TGraph* graph = new TGraph(n + 1);

    double dx = (xmax - xmin) / (double) n;
    std::vector<double> xvec(1);

    // loop over x values
    for (int i = 0; i <= n; i++) {
        xvec[0] = (double) i * dx + xmin;
        const double y = FitFunction(xvec, parameters) * GraphCorrection(i + 1);

        graph->SetPoint(i, xvec[0], y);
    }

    return graph;
}

// ---------------------------------------------------------
void BCFitter::FixDataAxis(unsigned int index, bool fixed)
{
    fFitterDataSet.Fix(index, fixed);
}

// ---------------------------------------------------------
bool BCFitter::GetFixedDataAxis(unsigned int index) const
{
    return fFitterDataSet.IsFixed(index);
}

// ---------------------------------------------------------
void BCFitter::SetErrorBandContinuous(bool flag)
{
    fErrorBandContinuous = flag;

    if (flag)
        return;

    // copy data x-values
    fErrorBandX = fDataSet->GetDataComponents(fFitFunctionIndexX);
}
