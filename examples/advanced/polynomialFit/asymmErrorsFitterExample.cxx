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
    // open log file with default level of logging
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::summary);

    // empty parameter set
    BCParameterSet pars;

    // 1st order polynomial
    pars.Add("p0", -0.1, 2.5);
    pars.Add("p1", -0.1, 0.1);

    PolAsymm m1("1st-Order Polynomial", pars);
    m1.SetPriorConstantAll();
    m1.PrintSummary();

    // 2nd order polynomial
    pars.Add("p2", -1e-4, 15e-4);
    PolAsymm m2("2nd-Order Polynomial", pars);
    m2.SetPriorConstantAll();
    m2.PrintSummary();

    // create model manager
    BCModelManager mgr;

    // add models to it, each with 50% prior probability
    mgr.AddModel(&m1, 1. / 2);
    mgr.AddModel(&m2, 1. / 2);

    // keep track of BCFitter objects
    BCFitter* models[] = { &m1, &m2 };

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
    for (unsigned i = 0; i < 2; ++i) {
        // tell model which data axes correspond to x and y axes
        models[i]->SetFitFunctionIndices(0, 1);
        // turn on calculaton of error band
        models[i]->SetFillErrorBand();
    }

    // set Metropolis algorithm for both models
    mgr.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    mgr.SetPrecision(BCEngineMCMC::kQuick);

    // the manager will marginalize all models
    mgr.MarginalizeAll();

    // each model is available as a BCModel* from the manager
    for (unsigned i = 0; i < mgr.GetNModels(); ++i) {
        BCModel* m = mgr.GetModel(i);
        // find mode, starting from global mode found by marginalizing
        m->FindMode();
        // write distrubutions to .root file
        m->WriteMarginalizedDistributions(m->GetSafeName() + "_plots.root", "RECREATE");
        // write distributions to .pdf file
        m->PrintAllMarginalized(m->GetSafeName() + "_plots.pdf");
    }

    // print the model comparison summary to the log
    mgr.PrintModelComparisonSummary();

    // print data with best fit result from both models
    TCanvas c1;

    // keep heap-allocated objects alive until plot is saved
    // => put them in the trash that is emptied at the end of the main function
    BCAux::BCTrash<TObject> trash;

    // Get data set as TGraphAsymmErrors, with x-axis errors set 0
    TGraphAsymmErrors* gdata = data.GetGraph(0, 1, -1, -1, 2, 3);
    trash.Put(gdata);
    gdata->SetMarkerStyle(20);
    gdata->SetMarkerColor(kBlack);

    // defines a histogram for the axes and draws it.
    TH2* haxes = data.CreateH2("haxes", ";x;y", 0, 1);
    trash.Put(haxes);
    haxes->SetStats(false);
    haxes->Draw();

    // define the legend
    TLegend leg(0.2, 0.60, 0.6, 0.8);
    leg.SetFillColor(kWhite);
    leg.SetBorderSize(0);

    // add data to legend
    leg.AddEntry(gdata, "Data", "PE");

    int Colors[3] = {kBlue, kRed, kGreen};

    // draw best fit functions
    // for each model (as BCFitter)
    for (unsigned i = 0; i < 2; ++i) {
        TGraph* g_err = models[i]->GetErrorBandGraph(0.16, 0.84);
        trash.Put(g_err);
        g_err->SetFillColor(Colors[i] - 7);
        g_err->SetFillStyle(3002);
        g_err->Draw("f");


        TGraph* g_best = models[i]->GetFitFunctionGraph();
        trash.Put(g_best);
        g_best->SetLineColor(Colors[i]);
        g_best->SetFillColor(Colors[i] - 7);
        g_best->SetFillStyle(3002);
        g_best->Draw("c same");
        leg.AddEntry(g_best, mgr.GetModel(i)->GetName().data(), "LF");
    }

    gdata->Draw("p same");

    // draw the legend
    leg.Draw();

    // print the canvas to a pdf file
    c1.Print("data_allmodels.pdf");

    // close log file
    BCLog::CloseLog();

    // all ok
    // exit with success
    return 0;
}

// ---------------------------------------------------------
