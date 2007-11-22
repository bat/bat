#include "BCModelPol1.h" 

// --------------------------------------------------------- 

BCModelPol1::BCModelPol1() : BCModel()
{

	// define parameters 

	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

BCModelPol1::BCModelPol1(const char* name) : BCModel(name)
{

	// define parameters 

	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

void BCModelPol1::DefineParameters()
{

	// adds the parameters which define the model and their limits. 

	this -> AddParameter("constant", 1.0, 3.0);   // index 0 
	this -> AddParameter("slope",    -0.03, 0.03); // index 1 

}

// --------------------------------------------------------- 

double BCModelPol1::LogAPrioriProbability(std::vector <double> parameters)
{

	double logprob = 0.; 

	// in this case we define flat a priopi probability across the whole parameter space

	double offset_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double offset_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	double slope_lower = this -> GetParameter(1) -> GetLowerLimit(); 
	double slope_upper = this -> GetParameter(1) -> GetUpperLimit(); 

	logprob -= log(offset_upper - offset_lower);
	logprob -= log(slope_upper - slope_lower);

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelPol1::LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{

	// define data values 

	//	double x       = datapoint -> GetValue(0); 
	double y       = datapoint -> GetValue(1); 
	double sigma_y = datapoint -> GetValue(2); 

	// calculate expectation value 

	double yex = this -> FitFunction(datapoint -> GetValues(), parameters); 

	// calculate probability for a single measurement 

	return BCMath::LogGaus(y, yex, sigma_y, true); 

}

// --------------------------------------------------------- 

double BCModelPol1::LogPoissonProbability(int nentries, std::vector <double> parameters)
{

	// Poisson term is 1. => Log of it is 0.

	return 0.;

}

// --------------------------------------------------------- 

void BCModelPol1::GetRandomVectorMetro(std::vector <double> &x)
{
  
  x[0] = fRandom -> Gaus(0.0, 0.05); 
  
  while (x[0] < -1.0 || x[0] > 1.0)
    x[0] = fRandom -> Gaus(0.0, 0.05); 
  
  x[1] = fRandom -> Gaus(0.0, 0.05); 
  
  while (x[1] < -1.0 || x[1] > 1.0)
    x[1] = fRandom -> Gaus(0.0, 0.05); 

  //  double * randx = new double[fNvar];

  //  fRandom -> RndmArray(fNvar, randx);

  //  for(int i=0;i<fNvar;i++)
  //    x[i] = randx[i];

  //  delete[] randx;
  //  randx = 0;
  
  //    cout << x[0] << " " << x[1] << endl; 
  
}

// --------------------------------------------------------- 

double BCModelPol1::FitFunction(std::vector <double> x, std::vector <double> parameters)
{
	
	// get parameter values

	double offset = parameters.at(0);
	double slope  = parameters.at(1);

	return offset + x.at(0) * slope; 

}

// --------------------------------------------------------- 

