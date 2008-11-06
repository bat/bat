#include <BCModelPol0.h>
#include <BCModelPol1.h>
#include <BCModelPol2.h>
#include <BCModelManager.h>
#include <BCLog.h>
#include <BCModelOutput.h>

#include "style.c"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLegend.h>

// ---------------------------------------------------------

int main()
{

	// ---------------------------------------------------------
	// set style
	// ----------------------------------------------------------

	// calls a function which defines a nicer style than the ROOT
	// default (included in file "style.c").
	SetStyle();

	// ---------------------------------------------------------
	// open log file
	// ---------------------------------------------------------
	BCLog::OpenLog("log.txt");

	// ---------------------------------------------------------
	// model definition
	// ---------------------------------------------------------

	// create three models
	BCModelPol0 * fModelPol0 = new BCModelPol0("ModelPol0");
	BCModelPol1 * fModelPol1 = new BCModelPol1("ModelPol1");
	BCModelPol2 * fModelPol2 = new BCModelPol2("ModelPol2");

	// add parameters to the models.
	fModelPol0 -> AddParameter("constant", 0.0, 4.0); // index 0

	// add parameters to the models.
	fModelPol1 -> AddParameter("constant", 0.0, 4.0); // index 0
	fModelPol1 -> AddParameter("slope", 0., 0.04); // index 1

	// add parameters to the models.
	fModelPol2 -> AddParameter("constant",  0.0, 4.0); // index 0
	fModelPol2 -> AddParameter("slope", -0.05, 0.05); // index 1
	fModelPol2 -> AddParameter("quadratic", 0.0, 0.008); // index 1

	// ---------------------------------------------------------
	// define model manager
	// ---------------------------------------------------------

	// defines a model manager. this manager will be used to make sure
	// all models use the same data set.
	BCModelManager* fModelManager = new BCModelManager();

	// adds all three models to the manager and passes the a priori
	// probability to the manager.
	fModelManager -> AddModel(fModelPol0, 1./3.);
	fModelManager -> AddModel(fModelPol1, 1./3.);
	fModelManager -> AddModel(fModelPol2, 1./3.);

	// ---------------------------------------------------------
	// read data from file
	// ---------------------------------------------------------

	// creates a new data set.
	BCDataSet * fDataSet = new BCDataSet();

	// data is read in from a text file. three values per data point are
	// read in: x, y, and the uncertainty on y. if the file is not found
	// or corrupt, the program returns -1.
	if (fDataSet -> ReadDataFromFileTxt("./data/data.txt", 3) != 0)
	{
		BCLog::Out(BCLog::error,BCLog::error,"You should run CreateData.c macro in ROOT to create the data file.");
		return -1;
	}

	// assigns the data set to the model manager.
	fModelManager -> SetDataSet(fDataSet);

	// the allowed range of data values has to be defined for error
	// propagation and fitting.
	fModelManager -> SetDataBoundaries(0, 0.0, 100.0); // possible x-values
	fModelManager -> SetDataBoundaries(1, 0.0, 5.0); // possible y-values
	fModelManager -> SetDataBoundaries(2, 0.2, 5.2); // possible sigmas

	// in this example the data points are not chosen arbitrarily but
	// they have fixed values on the x-axis, i.e, 5, 15, ... . also, the
	// resolution is fixed. in general this can be solved by fixing the
	// "axes".
	fModelManager -> FixDataAxis(0, true); // fixes x-values
	fModelManager -> FixDataAxis(2, true); // fixes resolution

	// ---------------------------------------------------------
	// optimization
	// ---------------------------------------------------------

	// set the optimization algorithm to the ROOT-version of Minuit.
	fModelManager -> SetOptimizationMethod(BCIntegrate::kOptMinuit);

	// perform optimization. the best fit parameter values will be
	// stored and can later be retrieved.
	fModelManager -> FindMode();

	// ---------------------------------------------------------
	// marginalize
	// ---------------------------------------------------------

	// during the marginalization, the error propagation can be done.
	// thus, the x- and y-indices of the data values have to be set.
	fModelManager -> SetFitFunctionIndices(0, 1);

	// turn on the calculation of the error band
	fModelManager -> SetFillErrorBand();

	// marginalizes the probability density with respect to all
	// parameters, i.e. constant and slope, and with respect to all
	// combinations of two parameters, in this case constant-slope.
	fModelManager -> MarginalizeAll();

	// ---------------------------------------------------------
	// Do goodness-of-fit test
	// ---------------------------------------------------------

	// once the most probable model was found the goodness-of-fit can be
	// tested. the question to be answered is: given the best fit
	// parameters what is the probability to obtain a better agreement
	// (conditional probability)? the answer can be given by integrating
	// over data.

	// calculate the p-value for the set of best-fit parameters.
	fModelManager -> CalculatePValue();

	// ---------------------------------------------------------
	// summarize
	// ---------------------------------------------------------

	// prints the results of the analysis into a file
	fModelManager -> PrintResults();

	// ---------------------------------------------------------
	// Print data with best fit result
	// ---------------------------------------------------------

	TCanvas * canvas_summary = new TCanvas();
	canvas_summary -> cd();

	// defines a graph with errors
	TGraphErrors * graph_data = new TGraphErrors();
	graph_data -> SetMarkerStyle(20);
	graph_data -> SetMarkerColor(kBlack);

	// sets the points of the graph to the data points read in
	// previously. loops over all entries.
	for (int i = 0; i < fModelManager -> GetNDataPoints(); i++)
	{
		graph_data -> SetPoint(i,
				fModelManager -> GetDataPoint(i) -> GetValue(0),
				fModelManager -> GetDataPoint(i) -> GetValue(1));
		graph_data -> SetPointError(i,
				0.0,
				fModelManager -> GetDataPoint(i) -> GetValue(2));
	}

	// defines a histogram for the axes and draws it.
	TH2D * hist_axes = new TH2D("hist_axes", ";x;y", 1, 0.0, 100.0, 1, 0.0, 6.0);
	hist_axes -> SetStats(false);
	hist_axes -> Draw();

	// draw data
	graph_data -> Draw("p same");

	// define best fit function graphs
	TGraph * graph_bestfit_pol0 = fModelPol0 -> GetFitFunctionGraph();
	graph_bestfit_pol0 -> SetLineColor(kRed);

	TGraph * graph_bestfit_pol1 = fModelPol1 -> GetFitFunctionGraph();
	graph_bestfit_pol1 -> SetLineColor(kBlack);

	TGraph * graph_bestfit_pol2 = fModelPol2 -> GetFitFunctionGraph();
	graph_bestfit_pol2 -> SetLineColor(kGreen);

	graph_bestfit_pol0 -> Draw("c same");
	graph_bestfit_pol1 -> Draw("c same");
	graph_bestfit_pol2 -> Draw("c same");

	// define and draw a legend.

	TLegend * legend = new TLegend(0.65, 0.70, 0.95, 0.95);
	legend -> SetFillColor(kWhite);
	legend -> SetBorderSize(0);
	legend -> AddEntry(graph_data, "Data", "P");
	legend -> AddEntry(graph_bestfit_pol0, "ModelPol0", "L");
	legend -> AddEntry(graph_bestfit_pol1, "ModelPol1", "L");
	legend -> AddEntry(graph_bestfit_pol2, "ModelPol2", "L");
	legend -> Draw("same");

	// print the canvas to a .ps file
	canvas_summary -> Print("data_allmodels.ps");

	// defines a new canvas
	TCanvas * canvas_bestfit = new TCanvas();
	canvas_bestfit -> cd();

	if (fModelPol0 -> GetPValue() > fModelPol1 -> GetPValue() &&
			fModelPol0 -> GetPValue() > fModelPol2 -> GetPValue())
	{
		hist_axes -> Draw();
		fModelPol0 -> GetErrorBandGraph(0.16, 0.84) -> Draw("f");
		graph_bestfit_pol0 -> SetLineColor(kBlack);
		graph_bestfit_pol0 -> Draw("c same");
	}

	if (fModelPol1 -> GetPValue() > fModelPol0 -> GetPValue() &&
			fModelPol1 -> GetPValue() > fModelPol2 -> GetPValue())
	{
		hist_axes -> Draw();
		fModelPol1 -> GetErrorBandGraph(0.16, 0.84) -> Draw("f");
		graph_bestfit_pol1 -> SetLineColor(kBlack);
		graph_bestfit_pol1 -> Draw("c same");
	}

	if (fModelPol2 -> GetPValue() > fModelPol0 -> GetPValue() &&
			fModelPol2 -> GetPValue() > fModelPol1 -> GetPValue())
	{
		hist_axes -> Draw();
		fModelPol2 -> GetErrorBandGraph(0.16, 0.84) -> Draw("f");
		graph_bestfit_pol2 -> SetLineColor(kBlack);
		graph_bestfit_pol2 -> Draw("c same");
	}

	// draw data points
	graph_data -> Draw("p same");

	// print the canvas to a .ps file
	canvas_bestfit -> Print("data_bestfit.ps");

	// ---------------------------------------------------------
	// close log file
	// ---------------------------------------------------------
	BCLog::CloseLog();

	return 1;

}

// ---------------------------------------------------------
