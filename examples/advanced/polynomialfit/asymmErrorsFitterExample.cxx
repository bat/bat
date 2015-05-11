#include "PolAsymm.h"

#include <BAT/BCAux.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCLog.h>
#include <BAT/BCModelManager.h>

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TH2.h>
#include <TLegend.h>

// ---------------------------------------------------------

int main()
{
    // set nice style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file with default level of logging
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::summary);

    // 1st order polynomial
    PolAsymm m1("1st-Order Polynomial");
    m1.AddParameter("p0", -0.1, 2.5);
    m1.AddParameter("p1", -0.1, 0.1);
    m1.SetPriorConstantAll();
    m1.PrintSummary();

    // 2nd order polynomial
    PolAsymm m2("2nd-Order Polynomial");
    m2.AddParameter("p0", -0.1, 2.5);
    m2.AddParameter("p1", -0.1, 0.1);
    m2.AddParameter("p2", -1e-4, 15e-4);
    m2.SetPriorConstantAll();
    m2.PrintSummary();

    // create model manager
    BCModelManager mgr;

    // add models to it, each with 50% of the a priori probability
    mgr.AddModel(&m1, 1. / 2);
    mgr.AddModel(&m2, 1. / 2);

    // create a new data set and fill it with data from a textfile
    // with four columns: x, y, y_error_low, y_error_high
    BCDataSet data;
    if (!data.ReadDataFromFileTxt("data.txt", 4)) {
        BCLog::OutError("You should run CreateData.C macro in ROOT to create the data file.");
        return -1;
    }

    // in this example the data points are not chosen arbitrarily but
    // they have fixed values on the x-axis, i.e, 5, 15, ... . also, the
    // resolution is fixed. So these data axes are fixed:
    data.Fix(0); // fixes x-values
    data.Fix(2); // fixes resolution
    data.Fix(3); // fixes resolution

    // adjust y (1) for 2.*sigma uncertainties
    // from y_error_low (2) and y_error_high (3)
    data.AdjustBoundForUncertainties(1, 2., 2, 3);

    data.PrintSummary();

    // assigns the data set to the model manager.
    mgr.SetDataSet(&data);

    // normalize all models
    mgr.Integrate();

    // for each model (as a BCFitter)
    for (unsigned i = 0; i < mgr.GetNModels(); ++i) {
        BCFitter* m = dynamic_cast<BCFitter*>(mgr.GetModel(i));
        // tell model which data axes correspond to x and y axes
        m->SetFitFunctionIndices(0, 1);
        // turn on calculaton of error band
        m->SetFillErrorBand();
    }

    // set Metropolis algorithm for both models
    mgr.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    mgr.MCMCSetPrecision(BCEngineMCMC::kMedium);

    // the manager will marginalize all models
    mgr.MarginalizeAll();

    // For each model
    for (unsigned i = 0; i < mgr.GetNModels(); ++i) {
        BCModel* m = mgr.GetModel(i);
        // find mode, starting from global mode found by marginalizing
        m->FindMode(m->GetGlobalMode());
        // write distrubutions to .root file
        m->WriteMarginalizedDistributions(m->GetSafeName() + "_plots.root", "RECREATE");
        // write distributions to .pdf file
        m->PrintAllMarginalized(m->GetSafeName() + "_plots.pdf");
    }

    // print the model comparison summary to the log
    mgr.PrintModelComparisonSummary();

    // print data with best fit result from both models
    TCanvas* c1 = new TCanvas();

    // Get data set as TGraphAsymmErrors, with x-axis errors set 0
    TGraphAsymmErrors* gdata = data.GetGraph(0, 1, -1, -1, 2, 3);
    gdata->SetMarkerStyle(20);
    gdata->SetMarkerColor(kBlack);

    // defines a histogram for the axes and draws it.
    TH2* haxes = data.CreateH2("haxes", ";x;y", 0, 1);
    haxes->SetStats(false);
    haxes->Draw();

    // define the legend
    TLegend* leg = new TLegend(0.2, 0.60, 0.6, 0.8);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);

    // add data to legend
    leg->AddEntry(gdata, "Data", "PE");

    int Colors[3] = {kBlue, kRed, kGreen};

    // draw best fit functions
    // for each model (as BCFitter)
    for (unsigned i = 0; i < mgr.GetNModels(); ++i) {
        BCFitter* m = dynamic_cast<BCFitter*>(mgr.GetModel(i));

        TGraph* g_err = m->GetErrorBandGraph(0.16, 0.84);
        g_err->SetFillColor(Colors[i] - 7);
        g_err->SetFillStyle(3002);
        g_err->Draw("f");

        TGraph* g_best = m->GetFitFunctionGraph();
        g_best->SetLineColor(Colors[i]);
        g_best->SetFillColor(Colors[i] - 7);
        g_best->SetFillStyle(3002);

        g_best->Draw("c same");
        leg->AddEntry(g_best, mgr.GetModel(i)->GetName().data(), "LF");
    }

    gdata->Draw("p same");

    // draw the legend
    leg->Draw();

    // print the canvas to a .pdf file
    c1->Print("data_allmodels.pdf");

    // close log file
    BCLog::CloseLog();

    // all ok
    // exit with success
    return 0;
}

// ---------------------------------------------------------
