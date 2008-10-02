#include <fstream.h>

/*
 * This macro generates a dataset consisting of 3 values for
 * each data point: x, y, error on y
 * The data are generated according to a linear function with
 * gausian spread. Below you can adjust the parameters of the
 * function (par0, par1), spread of the gausian (sigmay), the
 * range of x-values of the dataset (xmin, xmax) as well as
 * the number of points in the dataset (npoints)
 */

void CreateData()
{
	// parameters
	int npoints = 10; // number of data points
	float par0 = 2.; // offset
	float par1 = .01; // slope
	float sigmay = .2; // uncertainty on each measurement

	// range of x-values
	float xmin = 0.;
	float xmax = 100.;

	// initialize random number generator
	TRandom3 * fRandom = new TRandom3(1000);

	// open file
	std::fstream file_data;
	char filename[200];
	file_data.open("./data/data.txt", std::fstream::out);

	// loop over points
	for (int i = 0; i < npoints; i++)
	{
		// get x value
		float x = (xmax - xmin) / float(npoints) * (float(i) + 0.5);

		// get y value
		float mean = par0 + par1 * x;
		float y = fRandom -> Gaus(mean, sigmay);

		// write to file
		file_data << x << " " << y << " " << sigmay;

		if (i < npoints - 1)
			file_data << std::endl;
	}

	// close file
	file_data.close();
}
