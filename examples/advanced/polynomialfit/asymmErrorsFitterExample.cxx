#include "PolAsymm.h"

#include <BAT/BCAux.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCLog.h>
#include <BAT/BCModelManager.h>

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TLegend.h>

// ---------------------------------------------------------

int main()
{
    // set nice style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file with default level of logging
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // 1st order polynomial
    PolAsymm m1("Pol1");
    m1.AddParameter("p0", -0.2, 1.2);
    m1.AddParameter("p1", 0.015, 0.045);
    m1.SetPriorConstantAll();
    m1.PrintSummary();

    // 2nd order polynomial
    PolAsymm m2("Pol2");
    m2.AddParameter("p0", 0., 2.2);
    m2.AddParameter("p1", -0.1, 0.05);
    m2.AddParameter("p2", 0.0, 0.001);
    m2.SetPriorConstantAll();
    m2.PrintSummary();

    // create model manager
    BCModelManager mgr;

    // add models to it, each with 50% of the a priori probability
    mgr.AddModel(&m1, 0.5);
    mgr.AddModel(&m2, 0.5);

    // create a new data set and fill it with data from a textfile
    // with four columns: x, y, y_error_low, y_error_high
    BCDataSet data;
    if (!data.ReadDataFromFileTxt("data.txt", 4)) {
        BCLog::OutError("You should run CreateData.C macro in ROOT to create the data file.");
        return -1;
    }

    data.Dump();

    // assigns the data set to the model manager.
    mgr.SetDataSet(&data);

    // normalize all models and calculate the Bayes factors
    //	mgr.SetIntegrationMethod(BCIntegrate::kIntMonteCarlo);
    mgr.Integrate();

    // one can manually set the data bounds with the following three lines
    // but BAT will automatically know the data bounds
    // mgr.GetDataSet()->SetBounds(0, 0.0, 100.0); // possible x-values
    // mgr.GetDataSet()->SetBounds(1, 0.0, 5.0);   // possible y-values
    // mgr.GetDataSet()->SetBounds(2, 0.2, 5.2);   // possible sigmas

    // in this example the data points are not chosen arbitrarily but
    // they have fixed values on the x-axis, i.e, 5, 15, ... . also, the
    // resolution is fixed. in general this can be solved by fixing the
    // data axes.
    mgr.GetDataSet()->Fix(0); // fixes x-values
    mgr.GetDataSet()->Fix(2); // fixes resolution
    mgr.GetDataSet()->Fix(3); // fixes resolution

    // during the marginalization, the error propagation can be done
    // thus, the x- and y-indices of the data values have to be set.
    m1.SetFitFunctionIndices(0, 1);
    m2.SetFitFunctionIndices(0, 1);

    // turns on the calculation of the error band
    m1.SetFillErrorBand();
    m2.SetFillErrorBand();

    // set Metropolis algorithm for both models
    mgr.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    mgr.MCMCSetPrecision(BCEngineMCMC::kMedium);

    // marginalizes the probability density with respect to all
    // parameters subsequently for all models
    mgr.MarginalizeAll();

    // the MarginalizeAll method will also perform a mode finding
    // but we run minuit on top of that starting in the mode
    // found by MCMC to get more precise result
    for (unsigned int i = 0; i < mgr.GetNModels(); i++)
        mgr.GetModel(i)->FindMode(mgr.GetModel(i)->GetGlobalMode());

    // write marginalized distributions to a root file
    m1.WriteMarginalizedDistributions("Pol1.root", "RECREATE");
    m2.WriteMarginalizedDistributions("Pol2.root", "RECREATE");

    // print all marginalized distributions to a postscript file
    m1.PrintAllMarginalized("Pol1-marg.pdf");
    m2.PrintAllMarginalized("Pol2-marg.pdf");

    // write the summary
    mgr.PrintResults();
    mgr.PrintModelComparisonSummary("model-comparison.log");

    // print data with best fit result from both models
    TCanvas* c1 = new TCanvas();

    // defines a data graph with asymmetric errors
    TGraphAsymmErrors* gdata = new TGraphAsymmErrors();
    gdata->SetMarkerStyle(20);
    gdata->SetMarkerColor(kBlack);

    // sets the points of the graph to the data points read in
    // previously. loops over all entries.
    for (unsigned i = 0; i < mgr.GetDataSet()->GetNDataPoints(); i++) {
        gdata->SetPoint(i, mgr.GetDataSet()->GetDataPoint(i)->GetValue(0), mgr.GetDataSet()->GetDataPoint(i)->GetValue(1));
        gdata->SetPointEYlow(i, mgr.GetDataSet()->GetDataPoint(i)->GetValue(2));
        gdata->SetPointEYhigh(i, mgr.GetDataSet()->GetDataPoint(i)->GetValue(3));
    }

    mgr.PrintSummary("summary.txt");

    // defines a histogram for the axes and draws it.
    TH2D* haxes = new TH2D("haxes", ";x;y", 1, 0., 100., 1, 0., 5.);
    haxes->SetStats(false);
    haxes->Draw();

    // draw data
    gdata->Draw("p same");

    // gets the best fit function graphs and draws them
    TGraph* gm1_best = m1.GetFitFunctionGraph();
    gm1_best->SetLineColor(kBlue);

    TGraph* gm2_best = m2.GetFitFunctionGraph();
    gm2_best->SetLineColor(kRed);

    gm1_best->Draw("c same");
    gm2_best->Draw("c same");

    // defines and draw a legend.
    TLegend* leg = new TLegend(0.2, 0.60, 0.6, 0.8);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    leg->AddEntry(gdata, "Data", "p");
    leg->AddEntry(gm1_best, "1st order polynomial", "l");
    leg->AddEntry(gm2_best, "2nd order polynomial", "l");
    leg->Draw("same");

    // print the canvas to a .pdf file
    c1->Print("data_allmodels.pdf");

    // defines a new canvas
    TCanvas* c2 = new TCanvas("c2", "", 800, 500);
    c2->Divide(2, 1);

    // draw best fit for model 1 with errorband
    c2->cd(1);
    haxes->Draw();
    // draw central 68% errorband
    m1.GetErrorBandGraph(.16, .84)->Draw("f");
    gdata->Draw("p same");
    gm1_best->SetLineColor(kRed);
    gm1_best->Draw("c same");
    gPad->RedrawAxis();

    // draw best fit for model 2 with errorband
    c2->cd(2);
    haxes->Draw();
    // draw central 68% errorband
    m2.GetErrorBandGraph(.16, .84)->Draw("f");
    gdata->Draw("p same");
    gm2_best->SetLineColor(kRed);
    gm2_best->Draw("c same");
    gPad->RedrawAxis();

    // print the canvas to a .pdf file
    c2->Print("data_errorbands.pdf");

    // close log file
    BCLog::CloseLog();

    // all ok
    // exit with success
    return 0;
}

// ---------------------------------------------------------
