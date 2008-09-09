#include "BCModelTest.h" 

// --------------------------------------------------------- 

BCModelTest::BCModelTest(const char* name) : BCModel(name)
{

	// set original data set to zero 
	fTemporaryDataSet = 0; 

	// set test mode to zero 
	fTestModel = 0; 

	// reset pvalue and counter 
	fPValue = 0; 
	fPValueAbove = 0; 
	fPValueBelow = 0; 

	// reset loglikelihood and range
	fLogLikelihood = 0; 
	fLogLikelihoodMin = 1e99; 
	fLogLikelihoodMax = -1e99; 

	// define new histogram 
	fHistogramLogProb = 0; 

}

// --------------------------------------------------------- 

BCModelTest::~BCModelTest()
{

	// restore original data set 
	*(fTestModel -> GetDataSet()) = *fTemporaryDataSet; 

	// restore data point limits 
	for (int i = 0; i < this -> GetNParameters(); ++i)
		fTestModel -> SetDataBoundaries(fMapDataValue[i], 
																		this -> GetParameter(i) -> GetLowerLimit(), 
																		this -> GetParameter(i) -> GetUpperLimit()); 
	
	// delete temporary data set 
	//	delete fTemporaryDataSet; 

}

// --------------------------------------------------------- 

double BCModelTest::LogLikelihood(std::vector <double> parameters)
{

	// set the original data set to the new parameters 
	for (int i = 0; i < int(parameters.size()); ++i)
		{
			fTestModel -> GetDataSet() -> GetDataPoint(fMapDataPoint[i]) -> SetValue(fMapDataValue[i], parameters.at(i)); 
		}

	// calculate likelihood at the point of the original parameters 
	double loglikelihood = fTestModel -> LogLikelihood(fDataSet -> GetDataPoint(0) -> GetValues()); 
	
	// return likelihood 
	return loglikelihood; 

}

// --------------------------------------------------------- 

void BCModelTest::MCMCUserInterface()
{
	
	int nchains = this -> MCMCGetNChains(); 

	for (int i = 0; i < nchains; ++i)
		{
			// get likelihood at the point of the original parameters 
			double loglikelihood = this -> MCMCGetLogProbx(i); 
	
			// calculate pvalue 
			if (loglikelihood < fLogLikelihood)
				fPValueBelow++; 
			else
				fPValueAbove++; 

			// if histogram exists already, then fill it ... 
			if (fHistogramLogProb)
				fHistogramLogProb -> Fill(loglikelihood);  

			// ...otherwise find range 
			else
				{
					if (loglikelihood > fLogLikelihoodMax)
						fLogLikelihoodMax = loglikelihood;
					else if (loglikelihood < fLogLikelihoodMin)
						fLogLikelihoodMin = loglikelihood; 
				}
		}

}

// --------------------------------------------------------- 

int BCModelTest::SetTestPoint(std::vector<double> parameters)
{

	// check if the boundaries of the original data set exist. 
	if (!fTestModel -> GetFlagBoundaries())
		{
			BCLog::Out(BCLog::warning, BCLog::warning,"BCModelTest::SetTestDataPoint(). Boundaries of the original data set are not defined.");    

			return 0; 
		}

	// reset histogram 
	if (fHistogramLogProb)
		{
			delete fHistogramLogProb; 
			fHistogramLogProb = 0;
		}

	// reset variables 
	fPValue = 0; 
	fPValueAbove = 0; 
	fPValueBelow = 0; 

	// create temporary data set ... 
	fTemporaryDataSet = new BCDataSet(); 

	// ... and fill with the origianl one 
	*fTemporaryDataSet = *(fTestModel -> GetDataSet()); 

	// get number of data points and values 
	int ndatapoints = fTemporaryDataSet -> GetNDataPoints(); 
	int ndatavalues = fTemporaryDataSet -> GetDataPoint(0) -> GetNValues(); 
	
	// clear maps 
	fMapDataPoint.clear(); 
	fMapDataValue.clear(); 

	int counter = 0; 
	
	// remove parameters 
	fParameterSet -> clear(); 
	delete fParameterSet; 
	fParameterSet = new BCParameterSet; 

	// loop through data points and values 
	for (int i = 0; i < ndatapoints; ++i)
		for (int j = 0; j < ndatavalues; ++j)
			{
				if (fTestModel -> GetFixedDataAxis(j))
					continue; 

				// add parameter to this model 
				this -> AddParameter(Form("parameter_%i", counter), 
														 fTestModel -> GetDataPointLowerBoundary(j), 
														 fTestModel -> GetDataPointUpperBoundary(j)); 
				
				// add another element to the maps 
				fMapDataPoint.push_back(i); 
				fMapDataValue.push_back(j); 

				// increase counter 
				counter ++; 
			}

	// check if there are any non-fixed data values left 
	if (counter == 0)
		{
			BCLog::Out(BCLog::warning, BCLog::warning,"BCModelTest::SetTestDataPoint(). No non-fixed data values left.");    

			return 0; 
		}

	// create a new data set containing the vector of parameters which
	// are to be tested
	BCDataPoint * datapoint = new BCDataPoint(parameters); 
	BCDataSet * dataset = new BCDataSet(); 
	dataset -> AddDataPoint(datapoint); 

	// calculate likelihood of the original data set 
	fLogLikelihood = fTestModel -> LogLikelihood(parameters); 
	
	// if data set has been set before, delete 
	if (fDataSet)
		delete fDataSet; 

	// set data set of this model 
	fDataSet = dataset; 

	// put proper range to new data set 
	for (int i = 0; i < int(parameters.size()); ++i)
		this -> SetDataBoundaries(i, 
															fTestModel -> GetParameter(i) -> GetLowerLimit(), 
															fTestModel -> GetParameter(i) -> GetUpperLimit()); 
	
	return 1; 

}

// --------------------------------------------------------- 

double BCModelTest::GetCalculatedPValue(bool flag_histogram)
{

	// set histogram point to null 
	fHistogramLogProb = 0; 

	if (flag_histogram)
		{
			// modify MCMC for first run 
			this -> MCMCSetNIterationsMax(100000); 
			this -> MCMCSetNIterationsRun(10000); 
			
			// perform first run to obtain limits for the log(likelihood)
			this -> MarginalizeAll(); 

			// modify MCMC for first run 
			this -> MCMCSetNIterationsMax(100000); 
			this -> MCMCSetNIterationsRun(100000); 

			// create histogram 
			double D = fLogLikelihoodMax - fLogLikelihoodMin; 
			fHistogramLogProb = new TH1D(Form("hist_%s_logprob", this -> GetName().data()), ";N;log(prob)", 100, fLogLikelihoodMin - 0.1*D, fLogLikelihoodMax + 0.1*D); 
			fHistogramLogProb -> SetStats(kFALSE); 
		}
	else
		{
			// modify MCMC
			this -> MCMCSetNIterationsMax(100000); 
			this -> MCMCSetNIterationsRun(100000); 
		}

	// run MCMC 
	this -> MarginalizeAll(); 

	// check for convergence 
	if (this -> MCMCGetNIterationsConvergenceGlobal() < 0.)
		{
			BCLog::Out(BCLog::detail, BCLog::detail, 
								 " --> MCMC did not converge in evaluation of the p-value."); 

			return -1; 
		}

	// calculate the p-value 
	fPValue = double(fPValueBelow) / double(fPValueBelow + fPValueAbove); 

	return fPValue; 

}


// --------------------------------------------------------- 
