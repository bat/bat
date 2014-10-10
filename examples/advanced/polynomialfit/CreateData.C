/*
 * This macro generates a dataset consisting of 4 values for
 * each data point: x, y, error on y low, error on y high
 * The data are generated according to 2nd order polynomial with
 * the asymmetric spread corresponding to two half-gaussians with
 * different widths. Below you can adjust the parameters of the
 * function (par0, par1, par2), spread of the down side gausian
 * (sigmay1), spread of the up side gausian (sigmay2), the
 * range of x-values of the dataset (xmin, xmax) as well as
 * the number of points in the dataset (npoints)
 *
 * The random seed for the generation (RandomSeed) of the data
 * has to be set to 0 if you want to generate different datasets.
 * Otherwise always the same dataset will be generated.
 *
 * You can invoke this macro from the command line by calling
 *    $ root -l CreateData.c
 * and from  within ROOT by calling
 *    root[1] .x CreateData.c
 * or
 *    root[1] .L CreateData.c
 *    root[2] CreateData()
 *
 * To store the data to a different file then the default,
 * you can call
 *    $ root -l 'CreateData.c("data01.txt")'
 * or
 *    root[1] .x CreateData("data01.txt")
 * or
 *    root[1] .L CreateData.c
 *    root[2] CreateData("data01.txt"")
 */

#include <fstream>

void CreateData(const char * fname = "data.txt")
{
	// parameters
	int npoints = 10; // number of data points
	float par0 = 1.; // offset
	float par1 = -0.02; // slope
	float par2 = 0.0005; // quad
	float sigmay1 = .15; // down uncertainty on each measurement
	float sigmay2 = .4; // up uncertainty on each measurement

	// range of x-values
	float xmin = 0.;
	float xmax = 100.;

	// initial seed for the random number generator
	// set to 0 if you want to generate different dataset in every call
	int RandomSeed = 1000;

	cout <<
		"Generating dataset using function   f(x) = ax^2 + bx + c" << endl <<
		" - Range of x-values:      " << xmin << " < x < " << xmax << endl <<
		" - Number of data points:  " << npoints << endl <<
		" - Parameter a:            " << par2 << endl <<
		" - Paramerer b:            " << par1 << endl <<
		" - Paramerer c:            " << par0 << endl <<
		" - Uncertainty down:       " << sigmay1 << endl <<
		" - Uncertainty up:         " << sigmay2 << endl;

	// initialize random number generator
	TRandom3 * fRandom = new TRandom3(RandomSeed);

	// open file
	fstream file_data;
	char filename[200];
	file_data.open(fname, std::fstream::out);

	// create graph to draw the generated data
	TGraphAsymmErrors * g = new TGraphAsymmErrors(npoints);

	// loop over points
	for (int i = 0; i < npoints; i++)
	{
		// get x value
		float x = (xmax - xmin) / float(npoints) * (float(i) + .5);

		// get y value
		float mean = par0 + par1 * x + par2 * x * x;
		float y = fRandom -> Gaus(mean, sigmay1);
		float y2 = fRandom -> Gaus(mean, sigmay2);
		if (y>mean)
		{
			if(y2<mean)
				y2 = mean + mean - y2;
			y = y2;
		}

		// write to file
		file_data << x << " " << y << " " << sigmay1 << " " << sigmay2;

		// add point to graph
		g -> SetPoint(i, x, y);
		g -> SetPointEYlow(i, sigmay1);
		g -> SetPointEYhigh(i, sigmay2);

		file_data << endl;
	}

	// close file
	file_data.close();

	cout << "Data have been recorded to file  " << fname << endl;

	// draw true function and the genetated data
	TCanvas * c = new TCanvas();
	g -> SetMarkerStyle(20);
	g -> SetMarkerSize(1.);
	g -> Draw("a p");
	g -> GetXaxis() -> SetTitle("x");
	g -> GetYaxis() -> SetTitle("y");
	g -> SetTitle("");

	TF1 * f2d = new TF1("f2d","[0] + [1]*x + [2]*x*x", xmin, xmax);
	f2d -> SetParameters(par0, par1, par2);
	f2d -> SetLineWidth(2);
	f2d -> SetLineStyle(2);
	f2d -> SetLineColor(kBlue);
	f2d -> Draw("c same");

	TLegend * leg = new TLegend(0.15, 0.65, 0.5, 0.85);
	leg -> SetFillStyle(0);
	leg -> AddEntry(f2d, "True function", "l");
	leg -> AddEntry(g, "Generated data", "p");
	leg -> Draw();
}

