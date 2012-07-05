#include "DataGen.cxx"

class datagen;

int test()
{
	// initialize random number generator
	gRandom = new TRandom3(1001);
	
	// create data generator
	DataGen* datagen = new DataGen("test", 10, 0., 10., 40, 0., 10.);

	// define truth spectrum
	TF1* f_truth = new TF1("f_truth", "0.5*(gaus(0) + gaus(3))", 0., 10.);
	f_truth->SetParameter(0, 1.);
	f_truth->SetParameter(1, 2.5);
	f_truth->SetParameter(2, 0.5);
	f_truth->SetParameter(3, 1.);
	f_truth->SetParameter(4, 7.5);
	f_truth->SetParameter(5, 0.5);

	// define background spectrum
	TF1* f_background = new TF1("f_background", "0.1+0.01*x", 0., 10.);

	// define resolution
	TF1* f_resolution = new TF1("f_resolution", "gaus(0)", -3., 3.);
	f_resolution->SetParameter(0, 1.);
	f_resolution->SetParameter(1, 0.);
	f_resolution->SetParameter(2, 1.0);

	// set functions
	datagen->SetFuncTruthSignal(f_truth);
	datagen->SetFuncBackground(f_background);
	datagen->SetFuncResolution(f_resolution);

	// fill histograms
	datagen->FillHistograms(100000, 1000, 1000);

	// write histograms to file
	datagen->Write("histograms.root");

	// clean up memory
	delete datagen;

	// no error
	return 0;
}
