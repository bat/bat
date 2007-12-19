#include "BCEngineMCMC.h" 

// --------------------------------------------------------- 

BCEngineMCMC::BCEngineMCMC(int n) 
{

 	fMCMCNParameters         = 0; 
	fMCMCNChains             = n; 
	fMCMCNIterationsMax      = 100000; 
	fMCMCNIterationsBurnIn   = 100000; 
	fMCMCNIterationsPCA      = 100000; 
	fMCMCFlagIterationsAuto  = true; 
	fMCMCTrialFunctionScale  = 1.0; 
	fMCMCFlagInitialPosition = 0; 
	fMCMCRValueCriterion     = 0.25; 
	fMCMCNIterationsConvergenceGlobal = -1; 
	fMCMCRValue              = 100; 

	fMCMCRandom              = new TRandom3(0); 

}

// --------------------------------------------------------- 

BCEngineMCMC::~BCEngineMCMC() 
{

	if (fMCMCRandom)
		delete fMCMCRandom; 

	if (fMCMCPCA) 
		delete fMCMCPCA; 

}

// --------------------------------------------------------- 

BCEngineMCMC::BCEngineMCMC(const BCEngineMCMC & enginemcmc)
{

	enginemcmc.Copy(*this); 

}

// --------------------------------------------------------- 

BCEngineMCMC & BCEngineMCMC::operator = (const BCEngineMCMC & enginemcmc)
{

	if (this != &enginemcmc) 
		enginemcmc.Copy(* this); 

	return * this; 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCSetInitialPosition(int chain, std::vector<double> x0) 
{

	// check consistency 

	if (chain >= fMCMCNChains)
		{
			fMCMCInitialPosition.clear(); 

			return; 
		}

	if (int(fMCMCInitialPosition.size()) != (fMCMCNChains * fMCMCNParameters) && 
			int(fMCMCInitialPosition.size()) == fMCMCNParameters)
		{
			fMCMCInitialPosition.clear(); 
			for (int j = 0; j < fMCMCNChains; ++j)
				for (int i = 0; i < fMCMCNParameters; ++i)	
					fMCMCInitialPosition.push_back(x0.at(i)); 
		}

	bool flagin = true; 
	
	for (int j = 0; j < fMCMCNChains; ++j)
		for (int i = 0; i < fMCMCNParameters; ++i)
			if (fMCMCInitialPosition[j * fMCMCNParameters + i] < fMCMCBoundaryMin[i] || 
					fMCMCInitialPosition[j * fMCMCNParameters + i] > fMCMCBoundaryMax[i])
			flagin = false; 

	if (flagin == false)
		{
			fMCMCInitialPosition.clear(); 

			return; 			
		}

	// copy this intial position into the particular Markov chain

	else
		{
			for (int i = 0; i < fMCMCNParameters; ++i)
				fMCMCInitialPosition[chain * fMCMCNParameters + i] = x0.at(i); 
		}

	// use this intial position for the Markov chain 

	this -> MCMCSetFlagInitialPosition(2); 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCSetInitialPosition(std::vector<double> x0) 
{

	for (int i = 0; i < fMCMCNChains; ++i)
		this -> MCMCSetInitialPosition(i, x0); 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCSetInitialPositions(std::vector<double> x0s)
{

	if (int(fMCMCInitialPosition.size()) != (fMCMCNChains * fMCMCNParameters))
		return; 

	if (int(x0s.size()) != (fMCMCNChains * fMCMCNParameters))
		return; 
	
// copy these intial positions into the Markov chains 

	for (int i = 0; i < (fMCMCNChains * fMCMCNParameters); ++i)
		fMCMCInitialPosition[i] = x0s.at(i); 

	// use these intial positions for the Markov chain 

	this -> MCMCSetFlagInitialPosition(2); 

}	

// --------------------------------------------------------

void BCEngineMCMC::MCMCUserInterface()
{

}

// --------------------------------------------------------

void BCEngineMCMC::Copy(BCEngineMCMC & enginemcmc) const 
{

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCTrialFunction(std::vector <double> &x)
{

	// use uniform distribution for now 

	for (int i = 0; i < fMCMCNParameters; ++i) 
	  x[i] = 2.0 * (0.5 - fMCMCRandom -> Rndm()); 

}


// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetProposalPoint(int chain, std::vector <double> &x, bool pca)
{

	// get unscaled random point. this point might not be in the correct
	// volume.

	this -> MCMCTrialFunction(x); 

	// shift the point to the old point (x0) and scale it. 

	if (pca == false)
		for (int i = 0; i < fMCMCNParameters; ++i) 
			x[i] = fMCMCx[chain * fMCMCNParameters + i] + fMCMCTrialFunctionScale * x[i] * (fMCMCBoundaryMax.at(i) - fMCMCBoundaryMin.at(i)); 
	
	else 
		{
			double * data = new double[fMCMCNParameters]; 

			//			for (int i = 0; i < fMCMCNParameters; i++)
			//				data[i] = 

			delete [] data; 
		}

	// check if the point is in the correct volume. 
	
	for (int i = 0; i < fMCMCNParameters; ++i) 	
		if ((x[i] < fMCMCBoundaryMin[i]) || (x[i] > fMCMCBoundaryMax[i]))
			return false; 

	return true; 

}

// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetNewPointMetropolis(int chain, bool pca)
{
  
  // calculate index 

	int index = chain * fMCMCNParameters; 

 	// get proposal point 

	int counter = 0; 
	
 	while (!this -> MCMCGetProposalPoint(chain, fMCMCxLocal, pca) && counter < 1000)
		counter++; 

	// calculate probabilities of the old and new points 

	double p0 = fMCMCLogProbx[chain]; 
	double p1 = this -> LogEval(fMCMCxLocal); 
     
	// flag for accept 

 	bool accept = false;

	// if the new point is more probable, keep it ... 

	if (p1 >= p0)
	  accept = true;

	// ... or else throw dice. 

	else
 		{
 			double r = log(fMCMCRandom -> Rndm());

 			if(r < p1 - p0)
			  accept = true;
 		}

 	// fill the new point 

 	if(accept)
 	  {	    
			// increase counter 

			fMCMCNTrialsTrue[chain]++; 

	    // copy the point 	   

 	    for(int i = 0; i < fMCMCNParameters; ++i)
 	      {
					// save the point 
					
					fMCMCx[index + i] = fMCMCxLocal[i]; 
					
					// save the probability of the point 
					
					fMCMCLogProbx[chain] = p1; 
	      }
	  }

	else
		{
			// increase counter 

			fMCMCNTrialsFalse[chain]++; 
		}

	// increase counter 

	fMCMCNIterations[chain]++; 
	
	return accept; 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCUpdateConvergence()
{

	double sum = 0; 
	double sum2 = 0; 
	double sumv = 0; 

	// loop over chains 

	for (int i = 0; i < fMCMCNChains; ++i)
		{
			// calculate mean value of each chain 
			
			fMCMCMean[i]     += (fMCMCLogProbx[i] - fMCMCMean[i]) / double(fMCMCNIterations[i]); 

			// calculate variance of each chain 

			if (fMCMCNIterations[i] > 1) 
				fMCMCVariance[i] = (1.0 - 1/double(fMCMCNIterations[i])) * fMCMCVariance[i] 
					+ (fMCMCLogProbx[i] - fMCMCMean[i]) * (fMCMCLogProbx[i] - fMCMCMean[i]) / double(fMCMCNIterations[i] - 1); 

			sum  += fMCMCMean[i]; 
			sum2 += fMCMCMean[i] * fMCMCMean[i]; ; 
			sumv += fMCMCVariance[i]; 
		}

	// caluclate r-value 

	double mean            = sum / double(fMCMCNChains); 
	double varianceofmeans = sum2 / double(fMCMCNChains) - mean * mean; 
	double meanofvariance  = sumv / double(fMCMCNChains); 

	if (meanofvariance > 0)
		fMCMCRValue = varianceofmeans / meanofvariance; 

	if (fMCMCNIterationsConvergenceGlobal == -1 && fMCMCRValue < fMCMCRValueCriterion)
		fMCMCNIterationsConvergenceGlobal = fMCMCNIterations[0]; 

}

// --------------------------------------------------------

double BCEngineMCMC::LogEval(std::vector <double> parameters)
{

	// test function for now 
	// this will be overloaded by the user 

	double logprob = 0.0; 

	double x = parameters.at(0); 
	double y = parameters.at(1); 

	double cphi = cos(20.0/180*3.1416); 
	double sphi = sin(20.0/180*3.1416); 

	double xnew = - x * sphi + y * cphi; 
	double ynew =   x * cphi + y * sphi; 

	logprob += - (xnew - 70.0) * (xnew - 70.0)/ (2.0 * 15.0 * 15.0); 
	logprob += - (ynew - 70.0) * (ynew - 70.0)/ (2.0 * 7.0 * 7.0); 
	//	logprob += - (parameters.at(1) * parameters.at(1) - 70.0) * (parameters.at(1) - 70.0)/ 100.0; 

	return logprob; 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCPCARun()
{

	// reset run statistics 

	this -> MCMCResetRunStatistics(); 

	// create new TPrincipal 

	fMCMCPCA = new TPrincipal(fMCMCNParameters, "D"); 

	double * dataall = new double[fMCMCNParameters * fMCMCNIterationsPCA]; 

	double * sum = new double[fMCMCNParameters]; 
	double * sum2 = new double[fMCMCNParameters]; 

	for (int i = 0; i < fMCMCNIterationsPCA; ++i)
		{
			this -> MCMCGetNewPointMetropolis(0, false); 

			double * data = new double[fMCMCNParameters]; 

			for (int j = 0; j < fMCMCNParameters; ++j)
				{
					data[j]                           = fMCMCx[j]; 
					dataall[i * fMCMCNParameters + j] = fMCMCx[j]; 
				}

			fMCMCPCA -> AddRow(data); 

			delete [] data; 
		}

	// perform PCA 

	fMCMCPCA -> MakePrincipals();

	// re-run over data points to gain a measure for the spread of the variables 

	for (int i = 0; i < fMCMCNIterationsPCA; ++i)
		{
			double * data = new double[fMCMCNParameters]; 
			double * p    = new double[fMCMCNParameters]; 

			for (int j = 0; j < fMCMCNParameters; ++j)
				data[j] = dataall[i * fMCMCNParameters + j]; 

			fMCMCPCA -> X2P(data, p); 
			
			for (int j = 0; j < fMCMCNParameters; ++j)
				{
					sum[j]  += p[j]; 
					sum2[j] += p[j] * p[j]; 
				}

			delete [] data; 
			delete [] p; 
		}

	delete [] dataall; 

	// debug
	double * data = new double[fMCMCNParameters]; 
	double * p    = new double[fMCMCNParameters]; 
	data[0] = 40.0; 
	data[1] = 80.0; 
	fMCMCPCA -> X2P(data, p); 
	cout << p[0] << " " << p[1] << endl; 

	p[0] = 0.0; 
	p[1] = 0.0; 
	fMCMCPCA -> P2X(p, data, 2); 
	cout << data[0] << " " << data[1] << endl; 

	fMCMCPCAMean.clear(); 
	fMCMCPCAVariance.clear(); 
	
	for (int j = 0; j < fMCMCNParameters; ++j)
		{
			fMCMCPCAMean.push_back(sum[j] / double(fMCMCNIterationsPCA)); 
			fMCMCPCAVariance.push_back(sum2[j] / double(fMCMCNIterationsPCA) - fMCMCPCAMean[j] * fMCMCPCAMean[j]); 

			// debug
			cout << fMCMCPCAMean[j] << " " << sqrt(fMCMCPCAVariance[j]) << endl; 
		}

	// debug 
	fMCMCPCA -> Print(); 

	// reset run statistics 

	this -> MCMCResetRunStatistics(); 

}

// --------------------------------------------------------

int BCEngineMCMC::MCMCMetropolis()
{

	// initialize Markov chain 

	this -> MCMCInitialize(); 

	// perform burn-in run 

	for (int i = 0; i < fMCMCNIterationsBurnIn; ++i)
		for (int j = 0; j < fMCMCNChains; ++j)
			this -> MCMCGetNewPointMetropolis(j, false);

	// reset run statistics 

	this -> MCMCResetRunStatistics(); 

	// perform PCA run 

	this -> MCMCPCARun(); 

	// perform run 

	for (int i = 0; i < fMCMCNIterationsMax; ++i)
		{
			for (int j = 0; j < fMCMCNChains; ++j)
				{
					// get new point and increase counters 

					this -> MCMCGetNewPointMetropolis(j, false); 
				}

			// update the convergence calculation 

			this -> MCMCUpdateConvergence(); 

			// call user interface 

			this -> MCMCUserInterface(); 

		}

	return 0; 
	
}

// --------------------------------------------------------

void BCEngineMCMC::MCMCResetRunStatistics()
{

	for (int j = 0; j < fMCMCNChains; ++j)
		{
			fMCMCNIterations[j]  = 0; 
			fMCMCNTrialsTrue[j]  = 0; 
			fMCMCNTrialsFalse[j] = 0; 
			fMCMCMean[j]         = 0;
			fMCMCVariance[j]     = 0;
		}

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCInitialRun() 
{

	// burn-in run 

	for (int i = 0; i < fMCMCNIterationsBurnIn; ++i)
		cout << "here" << endl; 

}

// --------------------------------------------------------

int BCEngineMCMC::MCMCAddParameter(double min, double max) 
{

	// add the boundaries to the corresponding vectors 

	fMCMCBoundaryMin.push_back(min); 
	fMCMCBoundaryMax.push_back(max); 

	// increase the number of parameters by one 

	fMCMCNParameters++; 

	// return the number of parameters 

	return fMCMCNParameters; 

}

// --------------------------------------------------------

int BCEngineMCMC::MCMCInitialize()
{

	// reset variables 

	fMCMCNIterations.clear(); 
	fMCMCNIterationsConvergenceLocal.clear(); 

	fMCMCNTrialsTrue.clear(); 
	fMCMCNTrialsFalse.clear(); 

	fMCMCMean.clear(); 
	fMCMCVariance.clear(); 

	fMCMCSum.clear(); 
	fMCMCSum2.clear(); 

	fMCMCx.clear(); 

	fMCMCLogProbx.clear(); 

	fMCMCMinimumPoints.clear(); 

	fMCMCMinimumLogProb.clear(); 

	fMCMCNIterationsConvergenceGlobal = -1; 

	fMCMCTrees.clear(); 

	for (int i = 0; i < fMCMCNChains; ++i)
		{
			fMCMCNIterations.push_back(0); 
			fMCMCNIterationsConvergenceLocal.push_back(-1); 
			fMCMCNTrialsTrue.push_back(0); 
			fMCMCNTrialsFalse.push_back(0); 
			fMCMCMean.push_back(0); 
			fMCMCVariance.push_back(0); 
			fMCMCLogProbx.push_back(-1.0); 
			fMCMCMinimumLogProb.push_back(-1.0); 

			TTree * tree = new TTree(Form("MarkovChainTree_%i", i), Form("MarkovChainTree_%i", i)); 
 			fMCMCTrees.push_back(tree); 

			for (int j = 0; j < fMCMCNParameters; ++j)
				{
					fMCMCMinimumPoints.push_back(0.0); 
				}
		}

	// set initial position 

	for (int j = 0; j < fMCMCNChains; ++j)
		for (int i = 0; i < fMCMCNParameters; ++i)
			{
				switch(fMCMCFlagInitialPosition)
					{

						// use the center of the region 

					case 0 :
					  
					  fMCMCx.push_back(fMCMCBoundaryMin[i] + 0.5 * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i])); 

					  break; 
					  
					  // use a random variable in the valid region 
						
					case 1 : 
						
					  fMCMCx.push_back(fMCMCBoundaryMin[i] + fMCMCRandom -> Rndm() * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i])); 

					  break; 
					  
					  // use user-defined value 
					  
					case 2 : 
					  
						if (int(fMCMCInitialPosition.size()) == fMCMCNParameters * fMCMCNChains)
							fMCMCx.push_back(fMCMCInitialPosition[j * fMCMCNParameters + i]); 
						else
							fMCMCx.push_back(fMCMCBoundaryMin[i] + 0.5 * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i])); 

					  break; 
					}
			}
	
	std::vector <double> x0; 
	
	for (int j = 0; j < fMCMCNChains; ++j)
		{
			x0.clear(); 

			for (int i = 0; i < fMCMCNParameters; ++i)
				{
					x0.push_back(fMCMCx[j * fMCMCNParameters + i]); 
				}

			fMCMCLogProbx[j] = this -> LogEval(x0); 
		}

	x0.clear(); 

	fMCMCxLocal.clear(); 

	for (int i = 0; i < fMCMCNParameters; ++i)
		{
			fMCMCxLocal.push_back(fMCMCx[i]); 
		}

	return 0; 

}

// --------------------------------------------------------- 
