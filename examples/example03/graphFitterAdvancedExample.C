//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// a graph using a function defined by the user. TGraphErrors
// with defined y uncertainties has to be defined. In the fit the
// uncertainties are considered to be Gaussian.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x graphFitterAdvancedExample.C
//
// or
//
//    root[1] .L graphFitterAdvancedExample.C
//    root[2] graphFitterAdvancedExample()
//
// or from the command line
//
//    $ root graphFitterAdvancedExample.C
//
// By default the data in the file data/datax.txt are fitted. These
// are the same data that are used for the example in the BAT paper
// (see BAT webpage). To generate new data and fit them use
// CreateDataGraph() function.
//

TGraphErrors * CreateDataGraph(int n=100, double xmin=0.1, double xmax=19.9);

//
// The data are generated according to second order polynomial plus
// a gaussian peak. The parameters and the smearing can be set below

	const double p0 =  0.;
	const double p1 =  0.5;
	const double p2 =  0.02;
	const double a  = 15.;
	const double m  =  5.;
	const double s  =  0.5;

	const double sigmay = 4.;

//
// The macro performs fits with four different functions (models):
//   2nd order polynomial
//   gaussian peak + constant
//   gaussian peak + straight line
//   gaussian peak + 2nd order polynomial

#include <TGraphErrors.h>
#include <TF1.h>
#include <TRandom.h>
#include <TCanvas.h>

#include <fstream>
#include <vector>

// ----------------------------------------------------------------------
void graphFitterAdvancedExample()
{
	// set up nice style for plots
	BCAux::SetStyle();

	// uncomment the this command if you want more printouts during the run
	BCLog::SetLogLevel(BCLog::detail);

	// uncomment the next command if you want the data
	// to be generated randomly
//	TGraphErrors * gr = CreateDataGraph();

	// uncomment the next command if you want the data
	// to be read data from file
	// the file supplied with BAT distribution
	// is data/datax.txt and contains the data which were
	// used in the BAT paper
	TGraphErrors * gr = ReadDataGraph();

	gr -> SetMarkerStyle(20);
	gr -> SetMarkerSize(.5);

	// prepare fitting functions
	// 2nd order polynomial
	TF1 * f1 = new TF1("f1", "[0]+[1]*x+[2]*x*x", 0., 100.);
	f1 -> SetParLimits(0,   0.,    5.);  // offset
	f1 -> SetParLimits(1,   0.,    1.2); // slope
	f1 -> SetParLimits(2,  -0.1,   0.1); // quad

	// constant + gaussian
	TF1 * f2 = new TF1("f2", "[0]+[1]/(sqrt(2.*3.141592)*[3]) * exp(-(x-[2])*(x-[2])/(2.*[3]*[3]))", 0., 100.);
	f2 -> SetParLimits(0,   0.,   10.); // offset
	f2 -> SetParLimits(1,   0.,  200.); // A_gauss
	f2 -> SetParLimits(2,   2.,   18.); // mean
	f2 -> SetParLimits(3,   .2,    4.); // sigma

	// straight line + gaussian
	TF1 * f3 = new TF1("f3", "[0]+[1]*x+[2]/(sqrt(2.*3.141592)*[4]) * exp(-(x-[3])*(x-[3])/(2.*[4]*[4]))", 0., 100.);
	f3 -> SetParLimits(0,   0.,   10.); // offset
	f3 -> SetParLimits(1,   0.,    2.); // slope
	f3 -> SetParLimits(2,   0.,  200.); // A_gauss
	f3 -> SetParLimits(3,   2.,   18.); // mean
	f3 -> SetParLimits(4,   .2,    4.); // sigma

	// define fit function
	TF1 * f4 = new TF1("f4", "[0]+[1]*x+[2]*x*x+[3]/(sqrt(2.*3.141592)*[5]) * exp(-(x-[4])*(x-[4])/(2.*[5]*[5]))", 0., 100.);
	f4 -> SetParLimits(0,   0.,   10.); // offset
	f4 -> SetParLimits(1,   0.,    2.); // slope
	f4 -> SetParLimits(2,   0.,    .5); // quad
	f4 -> SetParLimits(3,   0.,  200.); // A_gauss
	f4 -> SetParLimits(4,   2.,   18.); // mean
	f4 -> SetParLimits(5,   .2,    4.); // sigma

	// initialize graph fitter and perform the fit
	BCGraphFitter * gf1 = new BCGraphFitter(gr,f1);
	gf1 -> MCMCSetNIterationsRun(20000);
	gf1 -> Fit();

	BCGraphFitter * gf2 = new BCGraphFitter(gr,f2);
	gf2 -> MCMCSetNIterationsRun(20000);
	gf2 -> Fit();

	BCGraphFitter * gf3 = new BCGraphFitter(gr,f3);
	gf3 -> MCMCSetNIterationsRun(20000);
	gf3 -> Fit();

	BCGraphFitter * gf4 = new BCGraphFitter(gr,f4);
	gf4 -> MCMCSetNIterationsRun(20000);
	gf4 -> Fit();

	// draw fit including the error band
	TCanvas * c = new TCanvas();
	c -> Divide(2,2);
	c -> cd(1);
	gf1 -> DrawFit();
	c -> cd(2);
	gf2 -> DrawFit();
	c -> cd(3);
	gf3 -> DrawFit();
	c -> cd(4);
	gf4 -> DrawFit();
	c -> Print("data-all-band.eps");
	delete c;

	// draw all fits in the same plot (w/o error bands)
	c = new TCanvas();
	gr -> Draw("ap");
	f1 -> SetLineColor(1);
	f1 -> SetLineWidth(2);
	f1 -> Draw("l same");
	f2 -> SetLineColor(2);
	f2 -> SetLineWidth(2);
	f2 -> Draw("l same");
	f3 -> SetLineColor(3);
	f3 -> SetLineWidth(2);
	f3 -> Draw("l same");
	f4 -> SetLineColor(4);
	f4 -> SetLineWidth(2);
	f4 -> Draw("l same");
	c -> Print("data-all.eps");

	// draw marginalized posterior distributions for all parameters
	// and print them into a postscript file
	gf1 -> PrintAllMarginalized("plots-1.ps");
	gf2 -> PrintAllMarginalized("plots-2.ps");
	gf3 -> PrintAllMarginalized("plots-3.ps");
	gf4 -> PrintAllMarginalized("plots-4.ps");
}


// ----------------------------------------------------------------------
TGraphErrors * CreateDataGraph(int n, double xmin, double xmax)
{
	const double pi =  3.141592;

	// initialize random number generator
	TRandom3 * ran = new TRandom3(0);

	double * xx = new double[n];
	double * yy = new double[n];
	double * err = new double[n];

	double dx = (xmax-xmin)/(double)(n-1);

	// loop over points
	for (int i=0;i<n;i++)
	{
		// get x value
		double x = xmin + (double)i*dx;

		// get y value
		double yexp = p0 + p1*x + p2*x*x + a/(sqrt(2.*pi)*s)*exp(-.5*(x-m)*(x-m)/(s*s));
		double y = ran->Gaus(yexp*100., sigmay*100.) / 100.;

		xx[i]=x;
		yy[i]=y;
		err[i]=sigmay;
	}

	TGraphErrors * g = new TGraphErrors(n,xx,yy,0,err);

	delete [] xx;
	delete [] yy;
	delete [] err;

	return g;
}

// ----------------------------------------------------------------------
TGraphErrors * ReadDataGraph(const char * file = "data/datax.txt")
{
	ifstream ifi(file);
	if (!ifi.is_open())
	{
		std::cerr<<"Couldn't open file "<<file<<std::endl;
		return 0;
	}

	std::vector<double> x(0);
	std::vector<double> y(0);
	std::vector<double> ey(0);
	while (!ifi.eof())
	{
		double tx,ty,tey;
		ifi>>tx>>ty>>tey;
		x.push_back(tx);
		y.push_back(ty);
		ey.push_back(tey);
	}
	ifi.close();

	int n = x.size();
	double * xx = new double[n];
	double * yy = new double[n];
	double * err = new double[n];

	// loop over points
	for (int i=0;i<n;i++)
	{
		xx[i]=x[i];
		yy[i]=y[i];
		err[i]=ey[i];
	}

	x.clear();
	y.clear();
	ey.clear();

	TGraphErrors * g = new TGraphErrors(n,xx,yy,0,err);

	delete [] xx;
	delete [] yy;
	delete [] err;

	return g;
}

// ----------------------------------------------------------------------
