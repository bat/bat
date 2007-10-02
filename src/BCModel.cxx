#include "BCModel.h" 
#include "BCLog.h"
#include "BCErrorCodes.h"
#include "BCMath.h"

#include <fstream.h>

// --------------------------------------------------------- 

BCModel::BCModel(const char * name) : BCIntegrate()
{

	fNormalization = -1.0; 

	fDataSet = 0; 

	fParameterSet = new BCParameterSet; 

	fIndex = -1; 
	fPValue = -1; 

	fName = (char *) name; 

	flag_ConditionalProbabilityEntry = true; 

	fErrorBandXY = 0; 

}

// --------------------------------------------------------- 

BCModel::BCModel() : BCIntegrate()
{

	fNormalization = -1.0; 

	fDataSet = 0; 

	fParameterSet = new BCParameterSet(); 

	fIndex = -1; 
	fPValue = -1; 

	fName = "model"; 

	flag_ConditionalProbabilityEntry = true; 

}

// --------------------------------------------------------- 

BCModel::~BCModel()
{

	delete fParameterSet; 

	if (fDataSet) 
		delete fDataSet; 

	if (fDataPointLowerBoundaries)
		delete fDataPointLowerBoundaries; 

	if (fDataPointUpperBoundaries)
		delete fDataPointUpperBoundaries; 

}

// --------------------------------------------------------- 

int BCModel::GetNDataPoints() 
{

	int npoints = 0; 

	if (fDataSet) 
		npoints = fDataSet -> GetNDataPoints(); 
	else
	  {
	    BCLog::Out(BCLog::warning, BCLog::warning,"BCModel::GetNDataPoints(). No data set defined."); 	 
	    return ERROR_NOEVENTS; 
	  }

	return npoints; 

}

// --------------------------------------------------------- 

BCDataPoint * BCModel::GetDataPoint(int index) 
{

	if (fDataSet) 
		return fDataSet -> GetDataPoint(index); 

	BCLog::Out(BCLog::warning, BCLog::warning,"BCModel::GetDataPoint. No data set defined."); 
	return 0; 

}

// --------------------------------------------------------- 

BCParameter * BCModel::GetParameter(int index)
{
	if (!fParameterSet) 
		return 0; 

	if (index < 0 || index >= this -> GetNParameters())
	{
		BCLog::Out(BCLog::warning, BCLog::warning, 
				Form("BCModel::GetParameter. Parameter index %d not within range.", index)); 
		return 0; 
	}

	return fParameterSet -> at(index); 
}

// --------------------------------------------------------- 

BCParameter * BCModel::GetParameter(const char * name)
{
	if (!fParameterSet) 
		return 0; 

	int index = -1; 

	for (int i = 0; i < this->GetNParameters(); i++)
		if (this->CompareStrings(name, this->GetParameter(i)->GetName()) == 0) 
			index = i; 

	if (index<0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
				Form("BCModel::GetParameter. Model %s has no parameter named %s.", this->GetName(), name));
		return 0;
	}

	return this->GetParameter(index);
}

// --------------------------------------------------------- 

std::vector <double> BCModel::GetErrorBand(double level) 
{

	std::vector <double> errorband; 

	if (!fErrorBandXY) 
		return errorband; 

	int nx = fErrorBandXY -> GetNbinsX(); 
	int ny = fErrorBandXY -> GetNbinsY(); 

	errorband.assign(nx, 0.0); 
	
	// loop over x and y bins 

	for (int ix = 1; ix <= nx; ix++) 
		{
			double sum = 0.0; 

			for (int iy = 1; iy <= ny; iy++)
				{
					double sumplusone = sum + fErrorBandXY -> GetBinContent(ix, iy); 

					if (sumplusone > level && sum < level)
						{
							errorband[ix-1] = fErrorBandXY -> GetYaxis() -> GetBinLowEdge(iy) + 
								(level - sumplusone) / fErrorBandXY -> GetBinContent(ix, iy) * fErrorBandXY -> GetYaxis() -> GetBinWidth(iy); 
						}

					sum += fErrorBandXY -> GetBinContent(ix, iy); 
				}
		}
	
	return errorband; 

}

// --------------------------------------------------------- 

TGraph * BCModel::GetErrorBandGraph(double level1, double level2) 
{
	
	if (!fErrorBandXY)
		return 0; 

	// define new graph 

	int nx = fErrorBandXY -> GetNbinsX(); 

	TGraph * graph = new TGraph(2 * nx); 
	graph -> SetFillStyle(1001);
  graph -> SetFillColor(kYellow);

	// get error bands 

	std::vector <double> ymin = this -> GetErrorBand(level1); 
	std::vector <double> ymax = this -> GetErrorBand(level2); 

	for (int i = 0; i < nx; i++)
		{
			graph -> SetPoint(i, fErrorBandXY -> GetXaxis() -> GetBinCenter(i + 1), ymin.at(i)); 
			graph -> SetPoint(nx + i, fErrorBandXY -> GetXaxis() -> GetBinCenter(nx - i), ymax.at(nx - i - 1)); 
		}

	return graph; 

}

// --------------------------------------------------------- 

TGraph * BCModel::GetFitFunctionGraph(std::vector <double> parameters)
{
	
	if (!fErrorBandXY)
		return 0; 

	// define new graph 

	int nx = fErrorBandXY -> GetNbinsX(); 

	TGraph * graph = new TGraph(nx); 

	// loop over x values 

	for (int i = 0; i < nx; i++)
		{
			double x = fErrorBandXY -> GetXaxis() -> GetBinCenter(i + 1);
			
			std::vector <double> xvec; 
			xvec.push_back(x); 

			double y = this -> FitFunction(xvec, parameters); 

			xvec.clear(); 

			graph -> SetPoint(i, x, y); 
		}

	return graph; 

}

// --------------------------------------------------------- 

void BCModel::SetDataBoundaries(int index, double lowerboundary, double upperboundary)
{

	// check if data set exists 

	if (!fDataSet)
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 "BCModel::SetDataBoundaries. Need to define data set first."); 

			return; 
		}

	// check if index is within range 

	if (index < 0 || index > fDataSet -> GetDataPoint(0) -> GetNValues())
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 "BCModel::SetDataBoundaries. Index out of range."); 

			return; 
		}

	// check if boundary data points exist 

	if (!fDataPointLowerBoundaries)
		fDataPointLowerBoundaries = new BCDataPoint(fDataSet -> GetDataPoint(0) -> GetNValues()); 

	if (!fDataPointUpperBoundaries)
		fDataPointUpperBoundaries = new BCDataPoint(fDataSet -> GetDataPoint(0) -> GetNValues()); 

	// set boundaries 

	fDataPointLowerBoundaries -> SetValue(index, lowerboundary); 
	fDataPointUpperBoundaries -> SetValue(index, upperboundary); 

}

// --------------------------------------------------------- 

int BCModel::AddParameter(const char * name, double lowerlimit, double upperlimit) 
{

	// create new parameter 

	BCParameter * parameter = new BCParameter(name, lowerlimit, upperlimit);   

	int flag_ok = 0; 

	flag_ok = this -> AddParameter(parameter); 

	if (flag_ok != 0) 
		delete parameter; 

	return flag_ok; 

}

// --------------------------------------------------------- 

int BCModel::AddParameter(BCParameter * parameter) 
{ 
  
	// check if parameter set exists 

	if (!fParameterSet) 
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 "BCModel::AddParameter. Parameter set does not exist"); 

			return ERROR_PARAMETERSETDOESNOTEXIST; 
		}

	// check if parameter with same name exists 

	int flag_exists = 0; 

	for (int i = 0; i < this -> GetNParameters(); i++)
		if (this -> CompareStrings(parameter -> GetName(), this -> GetParameter(i) -> GetName()) == 0)
			flag_exists = -1; 

	if (flag_exists < 0)
		{
			BCLog::Out(BCLog::warning, BCLog::warning, 
								 Form("BCModel::AddParameter. Parameter with name %s exists already. ", parameter -> GetName())); 

			return ERROR_PARAMETEREXISTSALREADY; 
		}

	// define index of new parameter 

	int index = int(fParameterSet -> size()); 

	parameter -> SetIndex(index); 

	// add parameter to parameter container 

	fParameterSet -> push_back(parameter); 

	// add parameters to integation methods 

	this -> SetParameters(fParameterSet); 

	return 0; 

}; 

// --------------------------------------------------------- 

double BCModel::LogProbabilityNN(std::vector <double> parameters) 
{

	// add log of conditional probability

	double logprob = this -> LogLikelihood(parameters);

	// add log of prior probability

	logprob += this -> LogAPrioriProbability(parameters);

	return logprob;

}

// --------------------------------------------------------- 

double BCModel::LogProbability(std::vector <double> parameters) 
{

	// check if normalized

	if (fNormalization<0. || fNormalization==0.)
		{
			BCLog::Out(BCLog::warning, BCLog::warning,
								 "BCModel::LogProbability. Normalization not available or zero.");
			return 0.;
		}

	return this -> LogProbabilityNN(parameters) - log(fNormalization);

}
  
// --------------------------------------------------------- 

double BCModel::LogLikelihood(std::vector <double> parameters)
{

	int ndatapoints = fDataSet -> GetNDataPoints();

	// get the log of poisson term first

	double logprob = this -> LogPoissonProbability(ndatapoints, parameters);

	// add log of conditional probabilities event-by-event

	for (int i=0;i<ndatapoints;i++)
		{
			BCDataPoint * datapoint = this -> GetDataPoint(i);
			logprob += this -> LogConditionalProbabilityEntry(datapoint, parameters);
		}

	return logprob;

}

// --------------------------------------------------------- 

double BCModel::LogEval(std::vector <double> parameters)
{

	return this -> LogProbabilityNN(parameters); 

}

// --------------------------------------------------------- 

double BCModel::EvalSampling(std::vector <double> parameters)
{

	return this -> SamplingFunction(parameters);

}

// --------------------------------------------------------- 

void BCModel::CalculateErrorBandXY(int nx, double xmin, double xmax, int ny, double ymin, double ymax, int niter)
{

	// calculate number of elements in the Markov chain if set negative or zero 
	
	if (niter <= 0) 
		niter = fNbins * fNbins * fNSamplesPer2DBin * fNvar;

	// define 2d histogram 

	if (!fErrorBandXY)
		delete fErrorBandXY; 

	double dx = (xmax - xmin) / double(nx); 
	double dy = (ymax - ymin) / double(ny); 

	fErrorBandXY = new TH2D("errorbandxy", "", 
													nx + 1, 
													xmin - 0.5 * dx, 
													xmax + 0.5 * dx, 
													ny + 1, 
													ymin - 0.5 * dy, 
													ymax + 0.5 * dy); 
	fErrorBandXY -> SetStats(kFALSE); 

	// initialize Markov chain 

	std::vector <double> randparameters;
	randparameters.assign(fNvar, 0.0);

	this -> InitMetro(); 

	// run Markov chain 

	for(int i = 0; i <= niter; i++)
		{
			// get point in Markov chain 

			GetRandomPointMetro(randparameters);

			// loop over x values 

			double x = 0; 

			for (int ix = 0; ix < nx; ix++)
				{
						// calculate x 

						x = fErrorBandXY -> GetXaxis() -> GetBinCenter(ix + 1); 

						// calculate y 

						std::vector <double> xvec; 
						xvec.push_back(x); 

						double y = this -> FitFunction(xvec, randparameters); 

						xvec.clear(); 

						// fill histogram 

						fErrorBandXY -> Fill(x, y); 
				}
		}
	
	fErrorBandXY -> Scale(1.0/fErrorBandXY -> Integral() * fErrorBandXY -> GetNbinsX()); 

}

// --------------------------------------------------------- 

double BCModel::SamplingFunction(std::vector <double> parameters)
{

	double probability = 1.0; 

	for (std::vector<BCParameter*>::const_iterator it = fParameterSet -> begin(); it != fParameterSet -> end(); ++it)
		{
			probability *= 1.0 / ((*it) -> GetUpperLimit() - 
														(*it) -> GetLowerLimit()); 
  }

	return probability; 

} 

// --------------------------------------------------------- 

std::vector<double> BCModel::FindMode()
{

	BCLog::Out(BCLog::summary, BCLog::summary, Form("Model \'%s\': Finding mode using Simulated Annealing algorithm.",this->GetName()));

	this -> FindModeSA();



/*
  BCIntegrate::BCIntegrationType integrationtypetemp = this -> GetIntegrationMethod(); 

  double normalizationtemp = this -> GetNormalization(); 

  this -> SetIntegrationMethod(BCIntegrate::kIMonteCarlo); 

  this -> Integrate(); 

  this -> SetNormalization(normalizationtemp); 

  this -> SetIntegrationMethod(integrationtypetemp); 
*/

	return this -> GetBestFitParameters(); 

}

// --------------------------------------------------------- 

double BCModel::Normalize()
{
	BCLog::Out(BCLog::summary, BCLog::summary, Form("Model \'%s\': Normalizing probability",this->GetName()));

	int n = this->GetNvar();
	
	// initialize BCIntegrate if not done already

	if (n == 0)
		{
			this->SetParameters(fParameterSet);
			n = this->GetNvar();
		}

	// integrate and get best fit parameters
	// maybe we have to remove the mode finding from here in the future

	fNormalization = this -> Integrate();

	return fNormalization;

}

// --------------------------------------------------------- 

int BCModel::CheckParameters(std::vector <double> parameters)
{

	// check if vectors are of equal size 

	if (!parameters.size() == fParameterSet -> size())
		return  ERROR_INVALIDNUMBEROFPARAMETERS; 

	// check if parameters are within limits 

	for (int i = 0; i < int(fParameterSet -> size()); i++)
		{
			BCParameter * modelparameter = fParameterSet -> at(i); 
  
			if (modelparameter -> GetLowerLimit() > parameters.at(i) ||
					modelparameter -> GetUpperLimit() < parameters.at(i)) 
				{
					BCLog::Out(BCLog::warning, BCLog::warning, 
										 Form("BCModel::CheckParameters. Parameter %s not within limits.", fParameterSet -> at(i) -> GetName())); 

					return ERROR_PARAMETERNOTWITHINRANGE; 
				}
		}

	return 0; 
  
}

// --------------------------------------------------------- 

BCH1D * BCModel::MarginalizeProbability(BCParameter * parameter)
{
	// print log

	BCLog::Out(BCLog::summary, BCLog::summary, Form("Marginalize probability with respect to %s", parameter -> GetName()));

	BCH1D * hprobability = new BCH1D();

	// get histogram

	TH1D * hist = this -> Marginalize(parameter);

	// set axis labels

	hist -> SetName(Form("hist_%s_%s", this -> GetName(),parameter -> GetName()));
	hist -> SetXTitle(parameter -> GetName());
	hist -> SetYTitle(Form("p(%s|data)", parameter -> GetName()));
	hist -> SetStats(kFALSE);

	// set histogram

	hprobability -> SetHistogram(hist);

	// set best fit parameter

	double bestfit = hprobability -> GetMode();

	int index = parameter -> GetIndex();

	if (fBestFitParametersMarginalized.size() == 0)
		for (int i = 0; i < this -> GetNParameters(); i++)
			fBestFitParametersMarginalized.push_back(0.0);

	fBestFitParametersMarginalized[index] = bestfit;

	return hprobability;

}

// --------------------------------------------------------- 

BCH2D * BCModel::MarginalizeProbability(BCParameter * parameter1, BCParameter * parameter2)
{
	// print log

	BCLog::Out(BCLog::summary, BCLog::summary,
		Form("Marginalize probability with respect to %s and %s",
			parameter1 -> GetName(),
			parameter2 -> GetName()));

	BCH2D * hprobability = new BCH2D();

	// get histogram

	TH2D * hist = this -> Marginalize(parameter1, parameter2);
	hist -> SetXTitle(Form("%s", parameter1 -> GetName()));
	hist -> SetYTitle(Form("%s", parameter2 -> GetName()));
	hist -> SetStats(kFALSE);

	// set histogram

	hprobability -> SetHistogram(hist);

	return hprobability;

}

// --------------------------------------------------------- 

int BCModel::MarginalizeAll()
{

	BCLog::Out(BCLog::summary, BCLog::summary,
		Form("Model \'%s\': Marginalizing all probability distributions.",this->GetName()));

	return this -> MarginalizeAllByMetro(this->GetName());

}


// --------------------------------------------------------- 

BCH1D * BCModel::GetMarginalized(BCParameter * parameter)
{

	if(fHProb1D.size()==0)
		{
			BCLog::Out(BCLog::warning, BCLog::warning,
								 "BCModel::GetMarginalized. MarginalizeAll() has to be run prior to this.");
			return 0;
		}

	int index = parameter -> GetIndex();

	// get histogram

	TH1D * hist = this -> GetH1D(index);

	if(!hist)
		return 0;

	BCH1D * hprob = new BCH1D();

	// set axis labels

	hist -> SetName(Form("hist_%s_%s", this -> GetName(), parameter -> GetName()));
	hist -> SetXTitle(parameter -> GetName());
	hist -> SetYTitle(Form("p(%s|data)", parameter -> GetName()));
	hist -> SetStats(kFALSE);

	// set histogram

	hprob -> SetHistogram(hist);

	// set best fit parameter

	double bestfit = hprob -> GetMode();

	if (fBestFitParametersMarginalized.size() == 0)
		for (int i = 0; i < this -> GetNParameters(); i++)
			fBestFitParametersMarginalized.push_back(0.0);

	fBestFitParametersMarginalized[index] = bestfit;

	return hprob;

}

// --------------------------------------------------------- 

BCH2D * BCModel::GetMarginalized(BCParameter * parameter1, BCParameter * parameter2)
{

	if(fHProb2D.size()==0)
		{
			BCLog::Out(BCLog::warning, BCLog::warning,
								 "BCModel::GetMarginalized. MarginalizeAll() has to be run prior to this.");
			return 0;
		}

	int index1 = parameter1 -> GetIndex();
	int index2 = parameter2 -> GetIndex();

	// get histogram

	TH2D * hist = this -> GetH2D(index1,index2);

	if(hist==0)
		return 0;

	BCH2D * hprob = new BCH2D();

	// set axis labels

	hist -> SetName(Form("hist_%s_%s_%s", this -> GetName(), parameter1 -> GetName(), parameter2 -> GetName()));
	hist -> SetXTitle(Form("%s", parameter1 -> GetName()));
	hist -> SetYTitle(Form("%s", parameter2 -> GetName()));
	hist -> SetStats(kFALSE);

	// set histogram

	hprob -> SetHistogram(hist);

	return hprob;

}

// --------------------------------------------------------- 

void BCModel::CreateData(int ndatasets, std::vector <double> parameters)
{

	// print log 

	BCLog::Out(BCLog::summary, BCLog::summary, "CreateData"); 

	std::vector <bool> grid; 
	std::vector <double> limits; 

	for (int i = 0; i < this -> GetDataSet() -> GetDataPoint(0) -> GetNValues(); i++)
		grid.push_back(false); 

	return this -> CreateDataGrid(ndatasets, parameters, grid, limits); 

}

// --------------------------------------------------------- 

void BCModel::CreateDataGrid(int ndatasets, std::vector <double> parameters, std::vector <bool> grid, std::vector <double> limits) 
{

	// define data stream 

	std::fstream stream_data; 

	// define stream for list of data files 

	std::fstream stream_list; 

	char listname[200]; 

	sprintf(listname, "./data/list_%s.txt", this -> GetName()); 

	stream_list.open(listname, std::fstream::out); 

	if (!stream_list.is_open())
		{
			BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::CreateDataGrid. Couldn't open filelist: %s", listname)); 
			return; 
		}

	// get number of data values 

	int nvalues = this -> GetDataSet() -> GetDataPoint(0) -> GetNValues(); 

	// calculate number of grid axes 

	int ngridaxes = 0; 

	for (int i = 0; i < nvalues; i++)
		if (grid.at(i) == true)
			ngridaxes++; 

	// initialize data creation 

	double pmaximum = 0.0; 
	double pmaximumsingle = 0.0; 

	int ninit = 1000; 

	for (int idataset = 0; idataset < ninit; idataset++)
		{
			std::vector <double> x; 

			// loop over values 

			for (int ivalue = 0; ivalue < nvalues; ivalue++)
				{
					// randomly choose value ... 
  
					if (this -> GetDataPointLowerBoundary(ivalue) != this -> GetDataPointUpperBoundary(ivalue))
						x.push_back(fRandom -> Uniform(this -> GetDataPointLowerBoundary(ivalue), 
																					 this -> GetDataPointUpperBoundary(ivalue))); 
					else 
						x.push_back(this -> GetDataPointUpperBoundary(ivalue)); 
				}

			// correlate data point values 

			this -> CorrelateDataPointValues(x);

			// define data object 

			BCDataPoint * datapoint = new BCDataPoint(x); 
      
			// calculate probability 
    
			double p = this -> ConditionalProbabilityEntry(datapoint, parameters); 

			// check if probability larger 

			if (p > pmaximumsingle) 
				pmaximumsingle = p; 

			// delete data object 
    
			delete datapoint; 

			// clear vector 

			x.clear(); 
		}

	// calculate maximum probability 

	pmaximum = BCMath::Min(1.0, pmaximumsingle * this -> GetNDataPointsMaximum()); 

	// loop over data sets 

	for (int idataset = 0; idataset < ndatasets; idataset++) 
		{
			// open new stream 

			char filename[200]; 

			sprintf(filename, "./data/data_%s_%d.txt", this -> GetName(), idataset); 

			stream_data.open(filename, std::fstream::out); 

			if (!stream_data.is_open())
				{
					BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::CreateDataGrid. Couldn't open file: %s", filename)); 
					return; 
				}
    
			// calculate number of entries 

			int nentries = 0; 

			// random number of entries ... 

			if (ngridaxes == 0) 
				{
					double pnentries = 0.0; 
					double temp_rand = 1.0; 

					while (temp_rand > pnentries)
						{	
							if (this -> GetNDataPointsMinimum() != this -> GetNDataPointsMaximum())
								nentries = BCMath::Nint(fRandom -> Uniform(this -> GetNDataPointsMinimum(), 
																													 this -> GetNDataPointsMaximum())); 
							else
								nentries = this -> GetNDataPointsMaximum(); 

							pnentries = this -> PoissonProbability(nentries, parameters); 
							temp_rand = fRandom -> Uniform(0, 1); 
						}
				}

			// ... or fixed number of entries on a grid 

			else 
				{
					nentries = 1; 

					for (int i = 0; i < ngridaxes; i++)
						nentries *= BCMath::Nint(limits.at(2 + i * 3)); 
				}
    
			// loop over entries 
    
			for (int ientry = 0; ientry < nentries; ientry++)
				{
					double p = 0.0; 
					double temp_rand = 1.0; 

					// calculate random data object 

					std::vector <double> x; 

					while (temp_rand > p)
						{
							// clear vector 

							x.clear(); 

							// loop over values 

							for (int ivalue = 0; ivalue < nvalues; ivalue++)
								{
									// randomly choose value ... 

									if (grid.at(ivalue) == false) 
										{
											if (this -> GetDataPointLowerBoundary(ivalue) != this -> GetDataPointUpperBoundary(ivalue))
												x.push_back(fRandom -> Uniform(this -> GetDataPointLowerBoundary(ivalue), 

																											 this -> GetDataPointUpperBoundary(ivalue))); 
											else
												x.push_back(this -> GetDataPointUpperBoundary(ivalue)); 
										}
	  
									// ... or calculate on a grid 

									else
										{
											double xgrid = 0; 

											if (ngridaxes == 1) 
												xgrid = limits.at(0) + double(ientry) * limits.at(1); 

											if (ngridaxes == 2) 
												{
													int ngrid = 0; 
		  
													for (int j = 0; j < ivalue; j++)
														if (grid.at(j) == true)
															ngrid++; 

													int ix = ientry % BCMath::Nint(limits.at(1 + ngrid * 3)); 
													xgrid = limits.at(ngrid * 3) + double(ix) * limits.at(1 + ngrid * 3); 
												}
											x.push_back(xgrid); 
										}
								}
      
							// correlate data point values 
      
							this -> CorrelateDataPointValues(x);

							// define data object 

							BCDataPoint * datapoint = new BCDataPoint(x); 
      
							// calculate probability 

							p = this -> ConditionalProbabilityEntry(datapoint, parameters); 

							// check if limit is set correctly 

							if (p > pmaximum) 
								{ 
									BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::CreateDataGrid. Probability larger than expected. Set limit to 1.0.")); 
									pmaximum = 1.0; 
								} 

							// delete data object 

							delete datapoint; 

							// calculate random number 

							temp_rand = pmaximum * fRandom -> Uniform(); 
						}

					// write data to file 

					for (int ivalue = 0; ivalue < nvalues; ivalue++)
						{
							stream_data << x.at(ivalue); 
							if (ivalue < nvalues-1)
								stream_data << " "; 
						}

					if (ientry < nentries - 1) 
						stream_data << std::endl; 

					// clear x 

					x.clear(); 
				}

			// write filename into filelist 

			stream_list << filename; 

			if (idataset < ndatasets - 1) 
				stream_list << std::endl; 

			// close data stream 

			stream_data.close(); 
		}

	// close list stream 

	stream_list.close(); 

}

// --------------------------------------------------------- 

BCH1D * BCModel::GoodnessOfFitTest(const char * filename, std::vector <double> parameters)
{

	// create marginalized probability 

	BCH1D * probability = new BCH1D(); 

	// vector containing the conditional probabilities 

	std::vector<double> likelihoodcontainer; 

	// open list file 

	std::fstream file;

	file.open(filename, std::fstream::in); 

	// check if file is open 

	if (file.is_open() == false)
		{
			BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::GoodnessOfFitTest. Couldn't open filelist: %s", filename)); 
			return probability; 
		} 

	// temporarily store data object container 

	BCDataSet * fDataSetTemp = fDataSet; 

	double NormTemp = fNormalization; 

	// loop over filenames 

	BCLog::Out(BCLog::detail, BCLog::detail, "BCModel::GoodnessOfFitTest. Loop over files"); 

	while (!file.eof())
		{
			char filename_temp[100]; 
    
			// read filename 
    
			file >> filename_temp; 

			// create new data set 

			BCDataSet * dataset = new BCDataSet(); 

			// open file 

			dataset -> ReadDataFromFileTxt(filename_temp, fDataSetTemp -> GetDataPoint(0) -> GetNValues()); 

			// set new data set 

			this -> SetDataSet(dataset); 

			// define event weight 

			double weight = 1.0; 

			// calculate weight 

			weight = this -> Likelihood(parameters); 

			// fill container  

			likelihoodcontainer.push_back(log10(weight)); 
		}

	// restore original data object container 

	this -> SetDataSet(fDataSetTemp); 

	fNormalization = NormTemp; 

	// find minimum and maximum 

	double minimum = 0.0; 
	double maximum = 0.0; 

	for (int i = 0; i < int(likelihoodcontainer.size()); i++)
		{
			if (likelihoodcontainer.at(i) < minimum)
				minimum = likelihoodcontainer.at(i); 

			if (likelihoodcontainer.at(i) > maximum || i == 0)
				maximum = likelihoodcontainer.at(i); 
		}

	// create histogram 

	TH1D * hist = new TH1D(Form("GOF_%s", this -> GetName()), "", 100, minimum - 0.1 * fabs(minimum), BCMath::Min(0.0, maximum + 0.1 * fabs(minimum))); 
	hist -> SetXTitle("log_{10}y=log_{10}p(data|#lambda^{*})"); 
	hist -> SetYTitle("1/N dN/dlog_{10}y"); 
	hist -> SetStats(kFALSE); 

	// fill histogram 

	for (int i = 0; i < int(likelihoodcontainer.size()); i++)
		hist -> Fill(likelihoodcontainer.at(i)); 

	// normalize to unity 

	hist -> Scale(1.0 / hist -> Integral()); 

	// set histogram 

	probability -> SetHistogram(hist); 

	// print summary to screen 

	double likelihood = this -> Likelihood(parameters); 
	// double fPValue =  probability -> GetIntegral(log10(likelihood), 0.0); 
	fPValue =  probability -> GetPValue(log10(likelihood)); 

	std::cout << std::endl; 
	std::cout << " Goodness-of-fit : " << std::endl; 
	std::cout << std::endl; 
	std::cout << " Model : " << this -> GetName() << std::endl; 
	std::cout << std::endl; 
	std::cout << " Parameters : " << std::endl; 

	for (int i = 0; i < int(parameters.size()); i++)
		std::cout << " Parameter : " 
							<< fParameterSet -> at(i) -> GetName() 
							<< " = " << parameters.at(i) << std::endl; 
	std::cout << std::endl; 
	std::cout << " Conditional probability p(data|lambda*) = " << likelihood << std::endl; 
	std::cout << " p-value = " << fPValue << std::endl; 
	std::cout << std::endl; 

	// clear container 

	likelihoodcontainer.clear(); 

return probability; 

}

// --------------------------------------------------------- 

BCH1D * BCModel::DoGoodnessOfFitTest(int ndatasets, std::vector<double> parameters, std::vector <bool> grid, std::vector <double> limits)
{

	// check if conditional probability has been defined on entry basis 

	if (flag_ConditionalProbabilityEntry == false) 
		{
			BCLog::Out(BCLog::warning, BCLog::warning, "BCModel::DoGoodnessOfFitTest. The method ConditionalProbabilityEntry has not been overloaded"); 
			return 0; 
		} 

	// print log 

	BCLog::Out(BCLog::summary, BCLog::summary, "Do goodness-of-fit-test"); 

	// create data sets 

	this -> CreateDataGrid(ndatasets, parameters, grid, limits); 

	// do goodness-of-fit test 

	BCH1D * gof = this -> GoodnessOfFitTest(Form("./data/list_%s.txt", this -> GetName()), parameters); 

	return gof; 

}

// --------------------------------------------------------- 

BCH1D * BCModel::DoGoodnessOfFitTest(int ndatasets, std::vector<double> parameters)
{

	// check if conditional probability has been defined on entry basis 

	if (flag_ConditionalProbabilityEntry == false) 
		{
			BCLog::Out(BCLog::warning, BCLog::warning, "BCModel::DoGoodnessOfFitTest. The method ConditionalProbabilityEntry has not been overloaded"); 
			return 0; 
		} 

	// print log 

	BCLog::Out(BCLog::summary, BCLog::summary, "Do goodness-of-fit-test"); 

	// create data sets 

	this -> CreateData(ndatasets, parameters); 

	// do goodness-of-fit test 

	BCH1D * gof = this -> GoodnessOfFitTest(Form("./data/list_%s.txt", this -> GetName()), parameters); 

	return gof; 

}

// --------------------------------------------------------- 

BCH1D * BCModel::DoGoodnessOfFitTest(int ndatasets)
{

	// check if conditional probability has been defined on entry basis 

	if (flag_ConditionalProbabilityEntry == false) 
		{
			BCLog::Out(BCLog::warning, BCLog::warning, "BCModel::DoGoodnessOfFitTest. The method ConditionalProbabilityEntry has not been overloaded"); 
			return 0; 
		} 

	std::vector<bool> grid; 
	std::vector<double> limits; 

	return this -> DoGoodnessOfFitTest(ndatasets, this -> GetBestFitParameters(), grid, limits); 

}

// --------------------------------------------------------- 

BCH1D * BCModel::DoGoodnessOfFitTest(const char * filename, std::vector<double> parameters)
{

	// print log 

	BCLog::Out(BCLog::summary, BCLog::summary, "Do goodness-of-fit-test"); 

	// do goodness-of-fit test 

	BCH1D * gof = this -> GoodnessOfFitTest(filename, parameters); 

	return gof; 

}

// --------------------------------------------------------- 

BCH1D * BCModel::DoGoodnessOfFitTest(const char * filename)
{

	return this -> DoGoodnessOfFitTest(filename, this -> GetBestFitParameters());  

}

// --------------------------------------------------------- 

void BCModel::CorrelateDataPointValues(vector<double> &x) 
{
  
}

// --------------------------------------------------------- 

double BCModel::HessianMatrixElement(BCParameter * parameter1, BCParameter * parameter2, std::vector<double> point)
{

	// check number of entries in vector 

	if (int(point.size()) != this -> GetNParameters())
		{
			BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::HessianMatrixElement. Invalid number of entries in the vector.")); 
			return -1; 
		}

	// define steps 

	double nsteps = 1e5; 

	double dx1 = (parameter1 -> GetUpperLimit() - parameter1 -> GetLowerLimit()) / nsteps; 
	double dx2 = (parameter2 -> GetUpperLimit() - parameter2 -> GetLowerLimit()) / nsteps ; 

	// define points at which to evaluate 

	std::vector<double> xpp; 
	std::vector<double> xpm; 
	std::vector<double> xmp; 
	std::vector<double> xmm; 

	xpp = point; 
	xpm = point; 
	xmp = point; 
	xmm = point; 

	xpp.at(parameter1 -> GetIndex()) = xpp.at(parameter1 -> GetIndex()) + dx1; 
	xpp.at(parameter2 -> GetIndex()) = xpp.at(parameter2 -> GetIndex()) + dx2; 

	xpm.at(parameter1 -> GetIndex()) = xpm.at(parameter1 -> GetIndex()) + dx1; 
	xpm.at(parameter2 -> GetIndex()) = xpm.at(parameter2 -> GetIndex()) - dx2; 

	xmp.at(parameter1 -> GetIndex()) = xmp.at(parameter1 -> GetIndex()) - dx1; 
	xmp.at(parameter2 -> GetIndex()) = xmp.at(parameter2 -> GetIndex()) + dx2; 

	xmm.at(parameter1 -> GetIndex()) = xmm.at(parameter1 -> GetIndex()) - dx1; 
	xmm.at(parameter2 -> GetIndex()) = xmm.at(parameter2 -> GetIndex()) - dx2; 

	// calculate probability at these points 

	double ppp = this -> Likelihood(xpp); 
	double ppm = this -> Likelihood(xpm); 
	double pmp = this -> Likelihood(xmp); 
	double pmm = this -> Likelihood(xmm); 

	// calculate derivative 

	double derivative = (ppp + pmm - ppm - pmp) / (4.0 * dx1 * dx2); 

	return derivative; 

}

// --------------------------------------------------------- 

void BCModel::PrintSummary()
{
	int nparameters = this -> GetNParameters(); 

	// model summary 
	cout
		<<endl
		<<"   ---------------------------------"<<endl
		<<"    Model : " << fName <<endl
		<<"   ---------------------------------"<<endl
		<<"     Index                : "<< fIndex <<endl
		<<"     Number of parameters : "<< nparameters <<endl
		<<endl
		<<"     - Parameters : " <<endl
		<<endl;

	// parameter summary
	for (int i=0; i<nparameters; i++)
		fParameterSet -> at(i) -> PrintSummary();

	// best fit parameters 
	if (this -> GetBestFitParameters().size() > 0) 
	{
		cout
			<<endl
			<<"     - Best fit parameters :"<<endl
			<<endl;

		for (int i=0; i<nparameters; i++)
		{
			cout
				<<"       "<< fParameterSet -> at(i) -> GetName()
				<<" = "<< this -> GetBestFitParameter(i)
				<<" (overall)"<<endl;
			if ((int)fBestFitParametersMarginalized.size() == nparameters) 
				cout
					<<"       "<< fParameterSet -> at(i) -> GetName()
					<<" = "<< this -> GetBestFitParameterMarginalized(i)
					<<" (marginalized)"<<endl;
		}
	}

	cout<<endl;

	// model testing 
	if (fPValue >= 0)
	{
		double likelihood = this -> Likelihood(this -> GetBestFitParameters()); 

		cout
			<<"   - Model testing:"<<endl
			<<endl
			<<"       p(data|lambda*) = "<< likelihood <<endl
			<<"       p-value         = "<< fPValue <<endl
			<<endl;
	}
}

// --------------------------------------------------------- 

void BCModel::PrintHessianMatrix(std::vector<double> parameters)
{

	// check number of entries in vector 

	if (int(parameters.size()) != this -> GetNParameters())
		{
			BCLog::Out(BCLog::warning, BCLog::warning, Form("BCModel::PrintHessianMatrix. Invalid number of entries in the vector.")); 
			return; 
		}

	// print to screen 

	std::cout << std::endl; 
	std::cout << " Hessian matrix elements: " << std::endl; 
	std::cout << " Point: "; 

	for (int i = 0; i < int(parameters.size()); i++)
		std::cout << parameters.at(i) << " "; 
	std::cout << endl; 

	// loop over all parameter pairs 

	for (int i = 0; i < this -> GetNParameters(); i++)
		for (int j = 0; j < i; j++) 
			{
				// calculate Hessian matrix element 

				double hessianmatrixelement = this -> HessianMatrixElement(fParameterSet -> at(i), 
																																	 fParameterSet -> at(j), 
																																	 parameters); 

				// print to screen 

				std::cout << " " << i << " " << j << " : " << hessianmatrixelement << std::endl; 
      }

}

// --------------------------------------------------------- 

BCDataPoint * BCModel::VectorToDataPoint(std::vector<double> data) 
{

	int sizeofvector = int(data.size()); 

	BCDataPoint * datapoint = new BCDataPoint(sizeofvector); 

	datapoint -> SetValues(data); 

	return datapoint; 

}

// --------------------------------------------------------- 

int BCModel::CompareStrings(const char * string1, const char * string2) 
{

	int flag_same = 0; 

	if (strlen(string1) != strlen(string2))
		return -1; 

	for (int i = 0; i < int(strlen(string1)); i++)
		if (string1[i] != string2[i]) 
			flag_same = -1; 

	return flag_same; 

}

// --------------------------------------------------------- 

