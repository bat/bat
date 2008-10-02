#include <fstream.h>

/*
 * This macro generates a dataset consisting of 3 values for
 * each data point: x, y, error on y
 * The data are generated according to 2nd order polynomial with
 * gausian spread. Below you can adjust the parameters of the
 * function (par0, par1, par2), spread of the gausian (sigmay), the
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

void CreateData(const char * fname = "./data/data.txt")
{
	// parameters
	int npoints = 10; // number of data points
	float par0 = 1.; // offset
	float par1 = 0.; // slope
	float par2 = 0.0003; // quad
	float sigmay = .2; // uncertainty on each measurement

	// range of x-values
	float xmin = 0.;
	float xmax = 100.;

	// random seed
	int RandomSeed = 1000;

	cout <<
		"Generating dataset function   f(x) = ax^2 + bx + c" << endl <<
		" - Range of x-values:      " << xmin << " < x < " << xmax << endl <<
		" - Number of data points:  " << npoints << endl <<
		" - Parameter a:            " << par2 << endl <<
		" - Paramerer b:            " << par1 << endl <<
		" - Paramerer c:            " << par0 << endl <<
		" - Uncertainty:            " << sigmay << endl;


	// initialize random number generator
	TRandom3 * fRandom = new TRandom3(RandomSeed);

	// open file
	fstream file_data;
	char filename[200];
	file_data.open(fname, std::fstream::out);

	// loop over points
	for (int i = 0; i < npoints; i++)
	{
		// get x value
		float x = (xmax - xmin) / float(npoints) * (float(i) + .5);

		// get y value
		float mean = par0 + par1 * x + par2 * x * x;
		float y = fRandom -> Gaus(mean, sigmay);

		// write to file
		file_data << x << " " << y << " " << sigmay;

		if (i < npoints - 1)
			file_data << endl;
	}

	// close file
	file_data.close();

	cout << "Data have been recorded to file  " << fname << endl;

}
