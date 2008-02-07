#include "BCEngineMCMC.h" 
#include "BCLog.h" 

// debug
#include <TCanvas.h> 

#define DEBUG 0

// --------------------------------------------------------- 

BCEngineMCMC::BCEngineMCMC() 
{

  // default settings 

 	fMCMCNParameters          = 0; 
	fMCMCNChains              = 5; 
	fMCMCNIterationsMax       = 1000000; 
	fMCMCNIterationsRun       = 100000; 
	fMCMCNIterationsBurnIn    = 0; 
	fMCMCNIterationsPCA       = 10000; 
	fMCMCFlagIterationsAuto   = false; 
	fMCMCTrialFunctionScale   = 1.0; 
	fMCMCFlagInitialPosition  = 1; 
	fMCMCRValueCriterion      = 0.01; 
	fMCMCNIterationsConvergenceGlobal = -1; 
	fMCMCRValue               = 100; 
	fMCMCFlagPCA              = false; 
	fMCMCSimulatedAnnealingT0 = 100; 
	fMCMCH1NBins              = 100; 
	fMCMCH2NBins              = 100; 
	fMCMCFlagPCATruncate      = false; 
	fMCMCPCAMinimumRatio      = 1e-7; 
	fMCMCNIterationsUpdate    = 1000; 
	fMCMCFlagWriteChainToFile = false; 
	fMCMCFlagPreRun           = false; 
	
	// set pointer to control histograms to NULL 

	for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
	  fMCMCH1Marginalized[i] = 0; 

	for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
	  fMCMCH2Marginalized[i] = 0; 

	fMCMCH1RValue = 0; 
	fMCMCH1Efficiency = 0; 

	// initialize random number generator 

	fMCMCRandom = new TRandom3(0); 

	// initialize 

	this -> MCMCInitialize(); 

}

// --------------------------------------------------------- 

BCEngineMCMC::BCEngineMCMC(int n)
{

  // set number of chains to n 

  fMCMCNChains = n; 

  // call default constructor 

  BCEngineMCMC(); 

}

// --------------------------------------------------------- 

BCEngineMCMC::~BCEngineMCMC() 
{

  // delete random number generator 

	if (fMCMCRandom)
		delete fMCMCRandom; 

	// delete PCA object

	if (fMCMCPCA) 
		delete fMCMCPCA; 

	// delete constrol histograms and plots 

	for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
		if (fMCMCH1Marginalized[i])
			delete fMCMCH1Marginalized[i]; 

	fMCMCH1Marginalized.clear(); 

	for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
		if (fMCMCH2Marginalized[i])
			delete fMCMCH2Marginalized[i]; 

	fMCMCH2Marginalized.clear(); 

	if (fMCMCH1RValue)
	  delete fMCMCH1RValue; 

	if (fMCMCH1Efficiency)
	  delete fMCMCH1Efficiency; 

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

TH2D * BCEngineMCMC::MCMCGetH2Marginalized(int index1, int index2) 
{

	int counter = 0; 
	int index = 0;

	// search for correct combination 

	for(int i = 0; i < fMCMCNParameters; i++)
		for (int j = 0; j < i; ++j)
			{
				if(j == index1 && i == index2)
					index = counter; 
				counter++;
		}	

	return fMCMCH2Marginalized[index]; 

}

// --------------------------------------------------------

std::vector <double> BCEngineMCMC::MCMCGetMaximumPoint(int i) 
{

  // create a new vector with the lenght of fMCMCNParameters 

  std::vector <double> x; 

  // check if i is in range 

  if (i < 0 || i >= fMCMCNChains)
    return x; 

  // copy the point in the ith chain into the temporary vector 

  for (int j = 0; j < fMCMCNParameters; ++j)
    x.push_back(fMCMCMaximumPoints.at(i * fMCMCNParameters + j)); 

  return x; 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCSetNChains(int n)
{

  fMCMCNChains = n; 

  // re-initialize 

  this -> MCMCInitialize(); 

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

	if (int(x0.size()) == fMCMCNParameters)
		{
			fMCMCInitialPosition.clear(); 
			for (int j = 0; j < fMCMCNChains; ++j)
				for (int i = 0; i < fMCMCNParameters; ++i)	
					fMCMCInitialPosition.push_back(x0.at(i)); 
		}

	else 
		return; 

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

void BCEngineMCMC::MCMCSetMarkovChainTrees(std::vector <TTree *> trees)
{

  // clear vector 

	fMCMCTrees.clear(); 

	// copy tree 

	for (int i = 0; i < int(trees.size()); ++i)
		fMCMCTrees.push_back(trees[i]); 

} 

// --------------------------------------------------------

void BCEngineMCMC::MCMCUserInterface()
{

  // this function has to be overloaded by the user. 

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
	  x[i] = fMCMCTrialFunctionScale * 2.0 * (0.5 - fMCMCRandom -> Rndm()); 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCTrialFunctionSingle(int ichain, int iparameter, std::vector <double> &x)
{

  // no check of range for performance reasons   
  
  // use uniform distribution 
  
  x[iparameter] = fMCMCTrialFunctionScaleFactor[ichain * fMCMCNParameters + iparameter] * 2.0 * (0.5 - fMCMCRandom -> Rndm()); 
  
}

// --------------------------------------------------------

double BCEngineMCMC::MCMCTrialFunctionIndependent(std::vector <double> &xnew, std::vector <double> &xold, bool newpoint)
{

  // use uniform distribution for now 

  if (newpoint)
    for (int i = 0; i < fMCMCNParameters; ++i) 
      xnew[i] = fMCMCRandom -> Rndm()* (fMCMCBoundaryMax.at(i) - fMCMCBoundaryMin.at(i)); 

  double prob = 1.0; 
  for (int i = 0; i < fMCMCNParameters; ++i) 
    prob *= 1.0 / (fMCMCBoundaryMax.at(i) - fMCMCBoundaryMin.at(i)); 
  
  return prob; 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCTrialFunctionAuto(std::vector <double> &x)
{

  // use uniform distribution for now 

  for (int i = 0; i < fMCMCNParameters; ++i) 
    // debug: is that correct? 
		x[i] = fMCMCRandom -> Gaus(fMCMCPCAMean[i], sqrt(fMCMCPCAVariance[i])); 
		//    x[i] = fMCMCRandom -> Gaus(0.0 , sqrt(fMCMCPCAVariance[i])); 
}

// --------------------------------------------------------

std::vector <double> BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(int i)
{
  
  // create a new vector with the length of fMCMCNParameters 

  std::vector <double> x; 

  // check if i is in range 

  if (i < 0 || i >= fMCMCNChains)
    return x; 

  // copy the scale factors into the temporary vector 

  for (int j = 0; j < fMCMCNParameters; ++j)
    x.push_back(fMCMCTrialFunctionScaleFactor.at(i * fMCMCNParameters + j)); 

  return x; 

}

// --------------------------------------------------------

double BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(int i, int j)
{

  // check if i is in range 

  if (i < 0 || i >= fMCMCNChains)
    return 0; 

  // check if j is in range 

  if (j < 0 || j >= fMCMCNParameters)
    return 0; 

  // return component of jth point in the ith chain 

  return fMCMCTrialFunctionScaleFactor.at(i *  fMCMCNChains +j);

}

// --------------------------------------------------------

std::vector <double> BCEngineMCMC::MCMCGetx(int i)
{

  // create a new vector with the length of fMCMCNParameters 

  std::vector <double> x; 

  // check if i is in range 

  if (i < 0 || i >= fMCMCNChains)
    return x; 

  // copy the point in the ith chain into the temporary vector 

  for (int j = 0; j < fMCMCNParameters; ++j)
    x.push_back(fMCMCx.at(i * fMCMCNParameters + j)); 

  return x; 
 
}

// --------------------------------------------------------

double BCEngineMCMC::MCMCGetx(int i, int j)
{

  // check if i is in range 

  if (i < 0 || i >= fMCMCNChains)
    return 0; 

  // check if j is in range 

  if (j < 0 || j >= fMCMCNParameters)
    return 0; 

  // return component of jth point in the ith chain 

  return fMCMCx.at(i *  fMCMCNParameters + j);
 
}

// --------------------------------------------------------

double BCEngineMCMC::MCMCGetLogProbx(int i)
{

  // check if i is in range 

  if (i < 0 || i >= fMCMCNChains)
    return -1; 

  // return log of the probability at the current point in the ith chain 

  return fMCMCLogProbx.at(i); 

}

// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetProposalPointMetropolis(int chain, std::vector <double> &x, bool pca)
{

	// get unscaled random point. this point might not be in the correct
	// volume.

	this -> MCMCTrialFunction(x); 

	// shift the point to the old point (x0) and scale it. 

	if (pca == false)
		{
			// get a proposal point from the trial function and scale it 

			for (int i = 0; i < fMCMCNParameters; ++i) 
				x[i] = fMCMCx[chain * fMCMCNParameters + i] + x[i] * (fMCMCBoundaryMax.at(i) - fMCMCBoundaryMin.at(i)); 
		}

	else 
		{
		  // create temporary points in x and p space 

			double * newp = new double[fMCMCNParameters]; 
			double * newx = new double[fMCMCNParameters]; 

			for (int i = 0; i < fMCMCNParameters; i++)
				{
					newp[i] = 0.0; 
					newx[i] = 0.0; 
				}

			// get a new trial point 

			this -> MCMCTrialFunctionAuto(x); 

			// get the old point in x space 

			for (int i = 0; i < fMCMCNParameters; i++)
				newx[i] = fMCMCx[chain * fMCMCNParameters + i]; 

			// transform the old point into p space 

			fMCMCPCA -> X2P(newx, newp); 

			// add new trial point to old point 

			for (int i = 0; i < fMCMCNParameters; i++)
				newp[i] += fMCMCTrialFunctionScale * x[i]; 

			// transform new point back to x space 

			//			fMCMCPCA -> P2X(newp, newx, fMCMCNParameters); 
			fMCMCPCA -> P2X(newp, newx, fMCMCNDimensionsPCA); 
			
			// copy point into vector 

			for (int i = 0; i < fMCMCNParameters; ++i) 
				x[i] = newx[i]; 

			delete [] newp; 
			delete [] newx; 
		}
	
	// check if the point is in the correct volume. 
	
	for (int i = 0; i < fMCMCNParameters; ++i) 	
		if ((x[i] < fMCMCBoundaryMin[i]) || (x[i] > fMCMCBoundaryMax[i]))
			return false; 

	return true; 

}
// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetProposalPointMetropolis(int chain, int parameter, std::vector <double> &x, bool pca)
{

  // get unscaled random point in the dimension of the chosen
  // parameter. this point might not be in the correct volume.
  
  this -> MCMCTrialFunctionSingle(chain, parameter, x); 
  
  // shift the point to the old point (x0) and scale it. 

  if (pca == false)
    {
      // copy the old point into the new 
      
      for (int i = 0; i < fMCMCNParameters; ++i)
	if (i != parameter)
	  x[i] = fMCMCx[chain * fMCMCNParameters + i];

      // modify the parameter under study 

      x[parameter] = fMCMCx[chain * fMCMCNParameters + parameter] + x[parameter] * (fMCMCBoundaryMax.at(parameter) - fMCMCBoundaryMin.at(parameter)); 
    }
  
	else 
		{
		  // create temporary points in x and p space 

			double * newp = new double[fMCMCNParameters]; 
			double * newx = new double[fMCMCNParameters]; 

			for (int i = 0; i < fMCMCNParameters; i++)
				{
					newp[i] = 0.0; 
					newx[i] = 0.0; 
				}

			// get a new trial point 

			this -> MCMCTrialFunctionSingle(chain, parameter, x); 

			// get the old point in x space 

			for (int i = 0; i < fMCMCNParameters; i++)
				newx[i] = fMCMCx[chain * fMCMCNParameters + i]; 

			// transform the old point into p space 

			fMCMCPCA -> X2P(newx, newp); 

			// add new trial point to old point 

			newp[parameter] += x[parameter] * sqrt(fMCMCPCAVariance[parameter]); 

			// transform new point back to x space 

			//			fMCMCPCA -> P2X(newp, newx, fMCMCNParameters); 
			fMCMCPCA -> P2X(newp, newx, fMCMCNDimensionsPCA); 
			
			// copy point into vector 

			for (int i = 0; i < fMCMCNParameters; ++i) 
				x[i] = newx[i]; 

			delete [] newp; 
			delete [] newx; 
		}
	
	// check if the point is in the correct volume. 
	
	for (int i = 0; i < fMCMCNParameters; ++i) 	
		if ((x[i] < fMCMCBoundaryMin[i]) || (x[i] > fMCMCBoundaryMax[i]))
			return false; 

	return true; 

}

// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetProposalPointMetropolisHastings(int chain, std::vector <double> &xnew, std::vector <double> &xold)
{

  // get a scaled random point.

  this -> MCMCTrialFunctionIndependent(xnew, xold, true); 

  // check if the point is in the correct volume. 
	
  for (int i = 0; i < fMCMCNParameters; ++i) 	
    if ((xnew[i] < fMCMCBoundaryMin[i]) || (xnew[i] > fMCMCBoundaryMax[i]))
      return false; 
  
  return true; 

}

// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetNewPointMetropolis(int chain, int parameter, bool pca)
{

  // calculate index 

	int index = chain * fMCMCNParameters; 

 	// get proposal point 

	int counter = 0; 
	
 	while (!this -> MCMCGetProposalPointMetropolis(chain, parameter, fMCMCxLocal, pca) && counter < 1000)
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
	    
	    fMCMCNTrialsTrue[chain * fMCMCNParameters + parameter]++; 

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

		  fMCMCNTrialsFalse[chain * fMCMCNParameters + parameter]++; 
		}

	// increase counter 

	fMCMCNIterations[chain]++; 
	
	return accept; 

}

// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetNewPointMetropolis(int chain, bool pca)
{

  // calculate index 

	int index = chain * fMCMCNParameters; 

 	// get proposal point 

	int counter = 0; 
	
 	while (!this -> MCMCGetProposalPointMetropolis(chain, fMCMCxLocal, pca) && counter < 1000)
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

	    for (int i = 0; i < fMCMCNParameters; ++i)
	      fMCMCNTrialsTrue[chain * fMCMCNParameters + i]++; 

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

		  for (int i = 0; i < fMCMCNParameters; ++i)
		    fMCMCNTrialsFalse[chain * fMCMCNParameters + i]++; 

		}

	// increase counter 

	fMCMCNIterations[chain]++; 
	
	return accept; 

}

// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetNewPointMetropolisHastings(int chain)
{
  
  // calculate index 

  int index = chain * fMCMCNParameters; 

  // save old point

  std::vector <double> xold; 

  for (int i = 0; i < fMCMCNParameters; ++i)
    xold.push_back(fMCMCxLocal.at(i)); 

  // get proposal point 

  int counter = 0; 
	
  while (!this -> MCMCGetProposalPointMetropolisHastings(chain, fMCMCxLocal, xold) && counter < 1000)
    counter++; 

  // calculate probabilities of the old and new points 

  double p0 = fMCMCLogProbx[chain] + log(this -> MCMCTrialFunctionIndependent(xold, fMCMCxLocal, false)); 
  double p1 = this -> LogEval(fMCMCxLocal) + log(this -> MCMCTrialFunctionIndependent(xold, fMCMCxLocal, false)); 
     
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

		  for (int i = 0; i < fMCMCNParameters; ++i)
		    fMCMCNTrialsTrue[chain * fMCMCNParameters + i]++; 

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

		  for (int i = 0; i < fMCMCNParameters; ++i)
		    fMCMCNTrialsFalse[chain * fMCMCNParameters + i]++; 

		}

	// increase counter 

	fMCMCNIterations[chain]++; 
	
	return accept; 

}

// --------------------------------------------------------

bool BCEngineMCMC::MCMCGetNewPointSimulatedAnnealing(int chain, bool pca)
{
  
  // calculate index 

	int index = chain * fMCMCNParameters; 

 	// get proposal point 

	int counter = 0; 
	
 	while (!this -> MCMCGetProposalPointMetropolis(chain, fMCMCxLocal, pca) && counter < 1000)
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

 			if(r < (p1 - p0) / this -> MCMCAnnealingSchedule(chain))
			  accept = true;
 		}

 	// fill the new point 

 	if(accept)
 	  {	    
			// increase counter 

		  for (int i = 0; i < fMCMCNParameters; ++i)
		    fMCMCNTrialsTrue[chain * fMCMCNParameters + i]++; 


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

		  for (int i = 0; i < fMCMCNParameters; ++i)
		    fMCMCNTrialsFalse[chain * fMCMCNParameters + i]++; 
		}

	// increase counter 

	fMCMCNIterations[chain]++; 
	
	return accept; 

}

// --------------------------------------------------------

double BCEngineMCMC::MCMCAnnealingSchedule(int chain)
{

  // this function can be overloaded by the user

  return fMCMCSimulatedAnnealingT0 * pow(1.0/fMCMCSimulatedAnnealingT0, double(fMCMCNIterations[chain])/fMCMCNIterationsMax);    
  
}

// --------------------------------------------------------

void BCEngineMCMC::MCMCUpdateStatisticsCheckMaximum()
{

	// loop over all chains 
  
  for (int i = 0; i < fMCMCNChains; ++i)
    {
      if (fMCMCLogProbx[i] > fMCMCMaximumLogProb[i] || fMCMCNIterations[i] == 1)
				{
					fMCMCMaximumLogProb[i] = fMCMCLogProbx[i]; 
					
					for (int j = 0; j < fMCMCNParameters; ++j)
						fMCMCMaximumPoints[i * fMCMCNParameters + j] = fMCMCx[i * fMCMCNParameters + j]; 
				}
    }	
	
}

// --------------------------------------------------------

void BCEngineMCMC::MCMCUpdateStatisticsFillHistograms()
{

	// loop over chains 

  for (int i = 0; i < fMCMCNChains; ++i)
    {
      // fill each 1-dimensional histogram 
      
      for (int j = 0; j < fMCMCNParameters; ++j)
				fMCMCH1Marginalized[i * fMCMCNParameters + j] -> Fill(fMCMCx[i * fMCMCNParameters + j]); 
			
      // fill each 2-dimensional histogram 
      
      int counter = 0; 
      
      for (int j = 0; j < fMCMCNParameters; ++j)
				for (int k = 0; k < j; ++k)
					{
						fMCMCH2Marginalized[counter] -> Fill(fMCMCx[i * fMCMCNParameters + k], fMCMCx[i * fMCMCNParameters + j]); 
						
						counter ++; 
					}
    }
	
}

// --------------------------------------------------------

void BCEngineMCMC::MCMCUpdateStatisticsTestConvergenceAllChains()
{

  if (fMCMCNChains > 1)
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
      
      // calculate r-value 
      
      double mean            = sum / double(fMCMCNChains); 
      double varianceofmeans = sum2 / double(fMCMCNChains) - mean * mean; 
      double meanofvariance  = sumv / double(fMCMCNChains); 
      
      if (meanofvariance > 0)
				fMCMCRValue = varianceofmeans / meanofvariance; 
      
      if (fMCMCNIterationsConvergenceGlobal == -1 && fMCMCRValue < fMCMCRValueCriterion)
				fMCMCNIterationsConvergenceGlobal = fMCMCNIterations[0] / fMCMCNParameters; 
      
      int dn = int(double(fMCMCNIterationsMax) / 100.0); 
      
      if (fMCMCNIterations[0] % dn == 0)
				fMCMCH1RValue -> SetBinContent(fMCMCNIterations[0] / dn + 1, log(fMCMCRValue)); 
    }
	
}

// --------------------------------------------------------

void BCEngineMCMC::MCMCUpdateStatisticsWriteChains()
{

	// loop over all chains

	for (int i = 0; i < fMCMCNChains; ++i)
		fMCMCTrees[i] -> Fill(); 
	
}

// --------------------------------------------------------

void BCEngineMCMC::MCMCUpdateStatistics()
{

  // check if maximum is reached

	this -> MCMCUpdateStatisticsCheckMaximum(); 
  
  // fill histograms of marginalized distributions 

	this -> MCMCUpdateStatisticsFillHistograms(); 
  
  // test convergence of single Markov chains
  //  --> not implemented yet 
  
  // test convergence of all Markov chains 
  
	this -> MCMCUpdateStatisticsTestConvergenceAllChains(); 

  // fill Markov chains
  
  if (fMCMCFlagWriteChainToFile)
		this -> MCMCUpdateStatisticsWriteChains(); 

}

// --------------------------------------------------------

double BCEngineMCMC::LogEval(std::vector <double> parameters)
{

	// test function for now 
	// this will be overloaded by the user 

  return 0.0; 

}

// --------------------------------------------------------

void BCEngineMCMC::MCMCPCARun()
{

	// reset run statistics 

	this -> MCMCResetRunStatistics(); 

	// create new TPrincipal 

	fMCMCPCA = new TPrincipal(fMCMCNParameters, "D"); 

	// create buffer of sampled points 

	double * dataall = new double[fMCMCNParameters * fMCMCNIterationsPCA]; 

	// create buffer for the sum and the sum of the squares 

	double * sum = new double[fMCMCNParameters]; 
	double * sum2 = new double[fMCMCNParameters]; 

	// reset the buffers

	for (int i = 0; i < fMCMCNParameters; ++i)
		{
			sum[i]  = 0.0; 
			sum2[i] = 0.0; 
		}

	// loop over iterations 

	for (int i = 0; i < fMCMCNIterationsPCA; ++i)
		{
		  // get new point 
		  
			this -> MCMCGetNewPointMetropolis(0, false); 

			// create buffer for the new point 

			double * data = new double[fMCMCNParameters]; 

			// copy the current point 

			for (int j = 0; j < fMCMCNParameters; ++j)
				{
					data[j]                           = fMCMCx[j]; 
					dataall[i * fMCMCNParameters + j] = fMCMCx[j]; 
				}

			// add point to PCA object 

			fMCMCPCA -> AddRow(data); 

			// delete buffer 

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

	fMCMCPCAMean.clear(); 
	fMCMCPCAVariance.clear(); 
	
	for (int j = 0; j < fMCMCNParameters; ++j)
		{
			fMCMCPCAMean.push_back(sum[j] / double(fMCMCNIterationsPCA)); 
			fMCMCPCAVariance.push_back(sum2[j] / double(fMCMCNIterationsPCA) - fMCMCPCAMean[j] * fMCMCPCAMean[j]); 
		}

	// debug 

	if (DEBUG)
		{
			fMCMCPCA -> Print("MSEV"); 
	
			for (int j = 0; j < fMCMCNParameters; ++j)
				cout << fMCMCPCAMean.at(j) << " " << sqrt(fMCMCPCAVariance.at(j)) << endl; 
		}

	// reset run statistics 

	this -> MCMCResetRunStatistics(); 

}

// --------------------------------------------------------

int BCEngineMCMC::MCMCMetropolisPreRun()
{

	// print on screen 

	BCLog::Out(BCLog::summary, BCLog::summary, "Pre-run Metropolis MCMC."); 

	// initialize Markov chain 
	
	this -> MCMCInitialize(); 

	this -> MCMCInitializeMarkovChains(); 
	
	// perform burn-in run 

	if (fMCMCNIterationsBurnIn > 0)
		BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Start burn-in run with %i iterations.", fMCMCNIterationsBurnIn)); 
	
	for (int i = 0; i < fMCMCNIterationsBurnIn; ++i)
		for (int j = 0; j < fMCMCNChains; ++j)
      for (int k = 0; k < fMCMCNParameters; ++k)
	this -> MCMCGetNewPointMetropolis(j, k, false);
  
  // reset run statistics 
  
  this -> MCMCResetRunStatistics(); 
	
  // perform PCA run 

  if (fMCMCFlagPCA) 
    {
      BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Start PCA run with %i iterations.", fMCMCNIterationsPCA)); 
      
      this -> MCMCPCARun(); 

      if (fMCMCFlagPCATruncate == true) 
				{
					fMCMCNDimensionsPCA = 0; 
					
					const double * ma = fMCMCPCA -> GetEigenValues() -> GetMatrixArray(); 
					
					for (int i = 0; i < fMCMCNParameters; ++i)
					  if (ma[i]/ma[0] > fMCMCPCAMinimumRatio)
					    fMCMCNDimensionsPCA++; 
					
					BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Use %i out of %i dimensions.", fMCMCNDimensionsPCA, fMCMCNParameters)); 
					
				}
      
      else
				fMCMCNDimensionsPCA = fMCMCNParameters; 
      
    }
  
  else
    BCLog::Out(BCLog::detail, BCLog::detail, " --> No PCA run."); 
  
  // reset run statistics 
  
  this -> MCMCResetRunStatistics(); 
  
  // perform run 
  
  BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Perform MCMC run with %i chains each with %i iterations maximum.", fMCMCNChains, fMCMCNIterationsMax)); 

  bool tempflag_writetofile = fMCMCFlagWriteChainToFile; 

  // don't write to file during pre run 

  fMCMCFlagWriteChainToFile = false; 
  
  int counter = 0;
  int counterupdate = 0; 
  bool convergence = false;
  bool flagefficiency = false; 
	
	std::vector <double> efficiency; 

  for (int i = 0; i < fMCMCNParameters; ++i)
		for (int j = 0; j < fMCMCNChains; ++j)
			efficiency.push_back(0.0); 

  while (counter < 100 || (counter < fMCMCNIterationsMax && !(convergence && flagefficiency)))
    {
      // loop over parameters 
			
      for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters)
				{
					// loop over chains 
					
					for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
						{
							this -> MCMCGetNewPointMetropolis(ichains, iparameters, fMCMCFlagPCA); 
						}
					
					// search for global maximum
					
					this -> MCMCUpdateStatisticsCheckMaximum(); 
					
					// check convergence status
					
					this -> MCMCUpdateStatisticsTestConvergenceAllChains(); 
				}
			
      // update scale factors 
			
      if (counterupdate % fMCMCNIterationsUpdate == 0 && counterupdate > 0)
				{
					// set flag

					flagefficiency = true; 

					// loop over chains
					
					for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
						{
							// loop over parameters 
							
							for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter)
								{
									// calculate efficiency 

									efficiency[ichains * fMCMCNParameters + iparameter] = double(fMCMCNTrialsTrue[ichains * fMCMCNParameters + iparameter]) / double(fMCMCNTrialsTrue[ichains * fMCMCNParameters + iparameter] + fMCMCNTrialsFalse[ichains * fMCMCNParameters + iparameter]); 
									
									// adjust scale factors if efficiency is too low 

									if (efficiency[ichains * fMCMCNParameters + iparameter] < 0.1)
										{
											fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] = fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] / 2.0; 
											BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Efficiency of parameter %i dropped below 10%% (eps = %.2lf%%) in chain %i. Set scale to %.2lf%%. ", iparameter, 100.0 * efficiency[ichains * fMCMCNParameters + iparameter], ichains, 100.0 * fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter])); 
											
										}
									
									// adjust scale factors if efficiency is too high 

									if (efficiency[ichains * fMCMCNParameters + iparameter] > 0.5 && (fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] < 1.0 || fMCMCFlagPCA))
										{
											fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] = fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] * 2.0; 
											BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Efficiency of parameter %i above 50%% (eps = %.2lf%%) in chain %i. Set scale to %.2lf%%. ", iparameter, 100.0 * efficiency[ichains * fMCMCNParameters + iparameter], ichains, 100.0 * fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter])); 
										}
									
									// reset counters 
									
									counterupdate = 0; 
									fMCMCNTrialsTrue[ichains * fMCMCNParameters + iparameter] = 0;
									fMCMCNTrialsFalse[ichains * fMCMCNParameters + iparameter] = 0;

									// check flag 

									if (efficiency[ichains * fMCMCNParameters + iparameter] < 0.1 || efficiency[ichains * fMCMCNParameters + iparameter] > 0.5) 
										flagefficiency = false; 
								}
						}
				}
			
      // increase counter 
			
      counter++;      
      counterupdate++;       

			// check convergence 

      if (fMCMCNIterationsConvergenceGlobal > 0)
				convergence = true; 
			
    }

  // print convergence status 

	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Markov chains ran for %i iterations. ", counter)); 

  if (fMCMCNIterationsConvergenceGlobal > 0) 
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains converged within %i iterations. ", fMCMCNChains, fMCMCNIterationsConvergenceGlobal)); 
  
  else
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains did not converge within %i iterations. ", fMCMCNChains, fMCMCNIterationsMax)); 
  // print efficiencies 

  std::vector <double> efficiencies; 

  for (int i = 0; i < fMCMCNParameters; ++i)
    efficiencies.push_back(0.0); 
	
  for (int i = 0; i < fMCMCNParameters; ++i)
    {
      for (int j = 0; j < fMCMCNChains; ++j)
				{
					efficiencies[i] += efficiency[j * fMCMCNParameters + i] / double(fMCMCNChains); 
				}
			
      BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Average efficiency for parameter %i: %.2lf%%. ", i, 100.0 * efficiencies[i])); 
    }
	
  // print scale factors 
	
  std::vector <double> scalefactors; 
	
  for (int i = 0; i < fMCMCNParameters; ++i)
    scalefactors.push_back(0.0); 
	
  for (int i = 0; i < fMCMCNParameters; ++i)
    {
      for (int j = 0; j < fMCMCNChains; ++j)
				{
					scalefactors[i] += fMCMCTrialFunctionScaleFactor[j * fMCMCNParameters + i] / double(fMCMCNChains); 
				}
			
      BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Average scale factor for parameter %i: %.2lf%%. ", i, 100 * scalefactors[i])); 
    }
  
	// fill efficiency plot 
  
  for (int i = 0; i < fMCMCNParameters; ++i)
    fMCMCH1Efficiency -> SetBinContent(i + 1, efficiencies[i]); 

	// reset flag 

  fMCMCFlagWriteChainToFile = tempflag_writetofile; 

	// set flag 

	fMCMCFlagPreRun = true; 

  return 0; 
  
}

// --------------------------------------------------------

int BCEngineMCMC::MCMCMetropolis()
{

	// check if prerun has been performed 

	if (!fMCMCFlagPreRun)
		this -> MCMCMetropolisPreRun(); 

	// print to screen 	

  BCLog::Out(BCLog::summary, BCLog::summary, "Run Metropolis MCMC."); 
    
  // reset run statistics 
  
  this -> MCMCResetRunStatistics(); 
	
  // perform run 
  
  BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Perform MCMC run with %i chains each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun)); 

  int counterupdate = 0; 
  bool convergence = false;
	bool flagefficiency = false; 
	
	std::vector <double> efficiency; 

  for (int i = 0; i < fMCMCNParameters; ++i)
		for (int j = 0; j < fMCMCNChains; ++j)
			efficiency.push_back(0.0); 

  for (int iiterations = 0; iiterations < fMCMCNIterationsRun; ++iiterations)
    {
      // loop over parameters 
			
      for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters)
				{
					// loop over chains 
					
					for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
						{
							this -> MCMCGetNewPointMetropolis(ichains, iparameters, fMCMCFlagPCA); 
						}
					
					// update statistics
      
					this -> MCMCUpdateStatistics(); 
      
					// update function fitting 

					this -> MCMCUpdateFunctionFitting(); 

					// call user interface 
      
					this -> MCMCUserInterface(); 
				}

      // update scale factors 
			
      if (counterupdate % fMCMCNIterationsUpdate == 0 && counterupdate > 0)
				{
					// set flag

					flagefficiency = true; 

					// loop over chains
					
					for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
						{
							// loop over parameters 
							
							for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter)
								{
									// calculate efficiency 

									efficiency[ichains * fMCMCNParameters + iparameter] = double(fMCMCNTrialsTrue[ichains * fMCMCNParameters + iparameter]) / double(fMCMCNTrialsTrue[ichains * fMCMCNParameters + iparameter] + fMCMCNTrialsFalse[ichains * fMCMCNParameters + iparameter]); 
									
									// adjust scale factors if efficiency is too low 

									if (efficiency[ichains * fMCMCNParameters + iparameter] < 0.1)
										{
											fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] = fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] / 2.0; 
											BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Efficiency of parameter %i dropped below 10%% (eps = %.2lf%%) in chain %i. Set scale to %.2lf%%. ", iparameter, 100.0 * efficiency[ichains * fMCMCNParameters + iparameter], ichains, 100.0 * fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter])); 
											
										}
									
									// adjust scale factors if efficiency is too high 

									if (efficiency[ichains * fMCMCNParameters + iparameter] > 0.5 && (fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] < 1.0 || fMCMCFlagPCA))
										{
											fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] = fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter] * 2.0; 
											BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Efficiency of parameter %i above 50%% (eps = %.2lf%%) in chain %i. Set scale to %.2lf%%. ", iparameter, 100.0 * efficiency[ichains * fMCMCNParameters + iparameter], ichains, 100.0 * fMCMCTrialFunctionScaleFactor[ichains * fMCMCNParameters + iparameter])); 
										}
									
									// reset counters 
									
									counterupdate = 0; 
									fMCMCNTrialsTrue[ichains * fMCMCNParameters + iparameter] = 0;
									fMCMCNTrialsFalse[ichains * fMCMCNParameters + iparameter] = 0;

									// check flag 

									if (efficiency[ichains * fMCMCNParameters + iparameter] < 0.1 || efficiency[ichains * fMCMCNParameters + iparameter] > 0.5) 
										flagefficiency = false; 
								}
						}
				}
			
      // increase counter 
			
      counterupdate++;       

			// check convergence 

      if (fMCMCNIterationsConvergenceGlobal > 0)
				convergence = true; 
			
    }

  // print convergence status 

	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Markov chains ran for %i iterations. ", fMCMCNIterationsRun)); 

  if (fMCMCNIterationsConvergenceGlobal > 0) 
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains converged within %i iterations. ", fMCMCNChains, fMCMCNIterationsConvergenceGlobal)); 
  
  else
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains did not converge within %i iterations. ", fMCMCNChains, fMCMCNIterationsRun)); 

  // print efficiencies 

  std::vector <double> efficiencies; 

  for (int i = 0; i < fMCMCNParameters; ++i)
    efficiencies.push_back(0.0); 
	
  for (int i = 0; i < fMCMCNParameters; ++i)
    {
      for (int j = 0; j < fMCMCNChains; ++j)
				{
					double efficiency = double(fMCMCNTrialsTrue[j * fMCMCNParameters + i]) / double(fMCMCNTrialsTrue[j * fMCMCNParameters + i] + fMCMCNTrialsFalse[j * fMCMCNParameters + i]);
					efficiencies[i] += efficiency / double(fMCMCNChains); 
				}
			
      BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Average efficiency for parameter %i: %.2lf%%. ", i, 100.0 * efficiencies[i])); 
    }
	
  // print scale factors 
	
  std::vector <double> scalefactors; 
	
  for (int i = 0; i < fMCMCNParameters; ++i)
    scalefactors.push_back(0.0); 
	
  for (int i = 0; i < fMCMCNParameters; ++i)
    {
      for (int j = 0; j < fMCMCNChains; ++j)
				{
					scalefactors[i] += fMCMCTrialFunctionScaleFactor[j * fMCMCNParameters + i] / double(fMCMCNChains); 
				}
			
      BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Average scale factor for parameter %i: %.2lf%%. ", i, 100.0 * scalefactors[i])); 
    }
  
	// fill efficiency plot 
  
  for (int i = 0; i < fMCMCNParameters; ++i)
    fMCMCH1Efficiency -> SetBinContent(i + 1, efficiencies[i]); 
 
	// set flag 

	fMCMCFlagPreRun = false; 
 
  return 0; 
  
}

// --------------------------------------------------------

int BCEngineMCMC::MCMCMetropolisHastings()
{

  BCLog::Out(BCLog::summary, BCLog::summary, "Run Metropolis-Hastings MCMC."); 
  
  // initialize Markov chain 
  
  this -> MCMCInitializeMarkovChains(); 
  
  // perform burn-in run 
  
  BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Start burn-in run with %i iterations.", fMCMCNIterationsBurnIn)); 
  
  for (int i = 0; i < fMCMCNIterationsBurnIn; ++i)
    for (int j = 0; j < fMCMCNChains; ++j)
      this -> MCMCGetNewPointMetropolisHastings(j);
  
  // reset run statistics 
  
  this -> MCMCResetRunStatistics(); 
	
  // perform run 
  
  BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Perform MCMC run with %i iterations.", fMCMCNIterationsMax)); 
  
  for (int i = 0; i < fMCMCNIterationsMax; ++i)
    {
      for (int j = 0; j < fMCMCNChains; ++j)
				{
					// get new point and increase counters 
	  
					this -> MCMCGetNewPointMetropolisHastings(j); 
				}
      
      // update statistics
      
      this -> MCMCUpdateStatistics(); 
      
      // call user interface 
      
      this -> MCMCUserInterface(); 
    }
  
  if (fMCMCNIterationsConvergenceGlobal > 0) 
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains converged within %i iterations. ", fMCMCNChains, fMCMCNIterationsConvergenceGlobal)); 
  
  else
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains did not converge within %i iterations. ", fMCMCNChains, fMCMCNIterationsMax)); 

  // debug

  if (DEBUG)
    {
      for (int i = 0; i < fMCMCNChains; ++i)
				{
					cout << i << " "  
							 << fMCMCMaximumLogProb[i] << endl;
	  
					for (int j = 0; j < fMCMCNParameters; ++j)
						cout << fMCMCMaximumPoints[i * fMCMCNParameters + j] << " "; 
					cout << endl; 
				}
    }

  return 0; 
  
}

// --------------------------------------------------------

int BCEngineMCMC::MCMCSimulatedAnnealing()
{

  BCLog::Out(BCLog::summary, BCLog::summary, "Run Simulated Annealing MCMC."); 
  
  // initialize Markov chain 
  
  this -> MCMCInitializeMarkovChains(); 
  
  // perform burn-in run 
  
  BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Start burn-in run with %i iterations.", fMCMCNIterationsBurnIn)); 
  
  for (int i = 0; i < fMCMCNIterationsBurnIn; ++i)
    for (int j = 0; j < fMCMCNChains; ++j)
      this -> MCMCGetNewPointSimulatedAnnealing(j, false);
  
  // reset run statistics 
  
  this -> MCMCResetRunStatistics(); 
	
  // perform PCA run 

  if (fMCMCFlagPCA) 
    {
      BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Start PCA run with %i iterations.", fMCMCNIterationsPCA)); 
      
      this -> MCMCPCARun(); 
    }
  
  else
    BCLog::Out(BCLog::detail, BCLog::detail, " --> No PCA run."); 
  
  // reset run statistics 
  
  this -> MCMCResetRunStatistics(); 
  
  // perform run 
  
  BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Perform MCMC run with %i chains each with %i iterations.", fMCMCNChains, fMCMCNIterationsMax)); 
  
  for (int i = 0; i < fMCMCNIterationsMax; ++i)
    {
      for (int j = 0; j < fMCMCNChains; ++j)
				{
					// get new point and increase counters 
	  
					this -> MCMCGetNewPointSimulatedAnnealing(j, fMCMCFlagPCA); 
				}
      
      // update statistics
      
      this -> MCMCUpdateStatistics(); 
      
      // call user interface 
      
      this -> MCMCUserInterface(); 
    }
  
  if (fMCMCNIterationsConvergenceGlobal > 0) 
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains converged within %i iterations. ", fMCMCNChains, fMCMCNIterationsConvergenceGlobal)); 
  
  else
    BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Set of %i Markov chains did not converge within %i iterations. ", fMCMCNChains, fMCMCNIterationsMax)); 
  
  // debug

  if (DEBUG)
    {
      for (int i = 0; i < fMCMCNChains; ++i)
				{
					cout << i << " "  
							 << fMCMCMaximumLogProb[i] << endl;
	  
					for (int j = 0; j < fMCMCNParameters; ++j)
						cout << fMCMCMaximumPoints[i * fMCMCNParameters + j] << " "; 
					cout << endl; 
				}
    }

  return 0; 
  
}

// --------------------------------------------------------

void BCEngineMCMC::MCMCResetRunStatistics()
{

	fMCMCNIterationsConvergenceGlobal = -1; 

	for (int j = 0; j < fMCMCNChains; ++j)
		{
			fMCMCNIterations[j]  = 0; 
			fMCMCNTrialsTrue[j]  = 0; 
			fMCMCNTrialsFalse[j] = 0; 
			fMCMCMean[j]         = 0;
			fMCMCVariance[j]     = 0;

			for (int k = 0; k < fMCMCNParameters; ++k)
			  {
			    fMCMCNTrialsTrue[j * fMCMCNParameters + k]  = 0; 
			    fMCMCNTrialsFalse[j * fMCMCNParameters + k] = 0; 
			  }
		}

	// reset marginalized distributions 

	for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
	  fMCMCH1Marginalized[i] -> Reset(); 

	for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
		fMCMCH2Marginalized[i] -> Reset(); 

	if (fMCMCH1RValue)
	  fMCMCH1RValue -> Reset(); 

	if (fMCMCH1Efficiency)
	  fMCMCH1Efficiency -> Reset(); 

	fMCMCRValue = 100; 

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

void BCEngineMCMC::MCMCInitializeMarkovChains()
{
  
	// evaluate function at the starting point 

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

}

// --------------------------------------------------------

int BCEngineMCMC::MCMCInitialize()
{

	// reset variables 

	fMCMCNIterations.clear(); 
	fMCMCNIterationsConvergenceLocal.clear(); 

	fMCMCNTrialsTrue.clear(); 
	fMCMCNTrialsFalse.clear(); 

	fMCMCTrialFunctionScaleFactor.clear();

	fMCMCMean.clear(); 
	fMCMCVariance.clear(); 

	fMCMCSum.clear(); 
	fMCMCSum2.clear(); 

	fMCMCx.clear(); 

	fMCMCLogProbx.clear(); 

	fMCMCMaximumPoints.clear(); 

	fMCMCMaximumLogProb.clear(); 

	fMCMCNIterationsConvergenceGlobal = -1; 

	for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)          
	  if (fMCMCH1Marginalized[i])
			delete fMCMCH1Marginalized[i]; 

	for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
	  if (fMCMCH2Marginalized[i]) 
	    delete fMCMCH2Marginalized[i]; 

	// clear plots 

	fMCMCH1Marginalized.clear(); 
	fMCMCH2Marginalized.clear(); 

	//	if (fMCMCH1RValue)
	//	  delete fMCMCH1RValue; 

	//	if (fMCMCH1Efficiency) 
	//	  delete fMCMCH1Efficiency; 

	// free memory for vectors 

	for (int i = 0; i < fMCMCNChains; ++i)
		{
			fMCMCNIterations.push_back(0); 
			fMCMCNIterationsConvergenceLocal.push_back(-1); 
			fMCMCMean.push_back(0); 
			fMCMCVariance.push_back(0); 
			fMCMCLogProbx.push_back(-1.0); 
			fMCMCMaximumLogProb.push_back(-1.0); 

			for (int j = 0; j < fMCMCNParameters; ++j)
				{
					fMCMCMaximumPoints.push_back(0.0); 
					fMCMCTrialFunctionScaleFactor.push_back(1.0);
					fMCMCNTrialsTrue.push_back(0); 
					fMCMCNTrialsFalse.push_back(0); 
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
				{
					fMCMCx.push_back(fMCMCInitialPosition[j * fMCMCNParameters + i]); 
				}
		  else
				{
					fMCMCx.push_back(fMCMCBoundaryMin[i] + 0.5 * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i])); 
				}

			break; 
		  
		  // use the center of the region as default 

		default: 
		  
		  fMCMCx.push_back(fMCMCBoundaryMin[i] + 0.5 * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i])); 
		  break; 
		}
	    }

	// copy the point of the first chain 
		
	fMCMCxLocal.clear(); 
	
	for (int i = 0; i < fMCMCNParameters; ++i)
	  {
	    fMCMCxLocal.push_back(fMCMCx[i]); 
	  }
	
	// define 1-dimensional histograms for marginalization 
	
	for (int j = 0; j < fMCMCNChains; ++j)
	  for(int i = 0; i < fMCMCNParameters; ++i)
	    {
	      double hmin1 = fMCMCBoundaryMin.at(i); 
	      double hmax1 = fMCMCBoundaryMax.at(i); 
	      
	      TH1D * h1 = new TH1D(Form("h1_chain_%i_parameter_%i", j, i), "",
				   fMCMCH1NBins, hmin1, hmax1);
	      
	      fMCMCH1Marginalized.push_back(h1); 
	    }
	
	for (int j = 0; j < fMCMCNChains; ++j)
	  for(int i = 0; i < fMCMCNParameters; ++i)
	    for (int k = 0; k < i; ++k)
	      {
					double hmin1 = fMCMCBoundaryMin.at(k); 
					double hmax1 = fMCMCBoundaryMax.at(k); 
					
					double hmin2 = fMCMCBoundaryMin.at(i); 
					double hmax2 = fMCMCBoundaryMax.at(i); 
		
					TH2D * h2 = new TH2D(Form("h2_chain_%i_parameters_%i_vs_%i", j, i, k), "",
															 fMCMCH2NBins, hmin1, hmax1, 
															 fMCMCH2NBins, hmin2, hmax2);
					
					fMCMCH2Marginalized.push_back(h2); 
	      }
	
	// define plot for R-value 
	
	if (fMCMCNChains > 1)
	  {
	    fMCMCH1RValue = new TH1D("h1_rvalue", ";Iteration;R-value", 
				     100, 0, double(fMCMCNIterationsMax)); 
	    fMCMCH1RValue -> SetStats(false); 
	  }
	
	fMCMCH1Efficiency = new TH1D("h1_efficiency", ";Chain;Efficiency",
				     fMCMCNParameters, -0.5, fMCMCNParameters - 0.5); 
	fMCMCH1Efficiency -> SetStats(false); 
	
	return 0; 

}

// --------------------------------------------------------- 
