#include "../datagenerator/DataGen.h"
#include "../datagenerator/DataGen.cxx"

class datagen;

int CreateHistograms()
{
	// initialize random number generator
	gRandom = new TRandom3(1001);
	
	// create data generator
	DataGen* datagen = new DataGen("test", 10, 0., 10., 40, 0., 10.);

	// define truth spectrum
	TF1* f_truth = new TF1("f_truth", "gaus(0)", 0., 10.);
	f_truth->SetParameter(0, 1.);
	f_truth->SetParameter(1, 5.0);
	f_truth->SetParameter(2, 1.0);

	// define resolution
	TF1* f_resolution = new TF1("f_resolution", "gaus(0)", -3., 3.);
	f_resolution->SetParameter(0, 1.);
	f_resolution->SetParameter(1, 0.);
	f_resolution->SetParameter(2, 1.5);

	// set functions
	datagen->SetFuncTruthSignal(f_truth);
	datagen->SetFuncResolution(f_resolution);

	// fill histograms
	datagen->FillHistograms(100000, 1000, 0);

	// write histograms to file
	datagen->Write("histograms.root");

	// clean up memory
	delete datagen;

	// no error
	return 0;
}
