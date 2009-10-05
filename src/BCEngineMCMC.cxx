/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BAT/BCEngineMCMC.h"

#include "BAT/BCParameter.h"
#include "BAT/BCDataPoint.h"
#include "BAT/BCLog.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TRandom3.h>

#include <iostream>
#include <math.h>

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC()
{
	// set default parameters for the mcmc
	this -> MCMCSetValuesDefault();

	// initialize random number generator
	fMCMCRandom = new TRandom3(0);
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
void BCEngineMCMC::MCMCSetValuesDefault()
{
	fMCMCNParameters          = 0;
	fMCMCFlagWriteChainToFile = false;
	fMCMCFlagPreRun           = false;
	fMCMCFlagRun              = false;
	fMCMCFlagFillHistograms   = true;
	fMCMCEfficiencyMin        = 0.15;
	fMCMCEfficiencyMax        = 0.50;
	fMCMCFlagInitialPosition  = 1;
	fMCMCNLag                 = 1;
	fMCMCNIterationsUpdateMax = 10000;

	this -> MCMCSetValuesDetail();
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetValuesQuick()
{
	fMCMCNChains              = 1;
	fMCMCNIterationsMax       = 1000;
	fMCMCNIterationsRun       = 10000;
	fMCMCNIterationsBurnIn    = 0;
	fMCMCNIterationsPreRunMin = 0;
	fMCMCFlagInitialPosition  = 1;
	fMCMCRValueCriterion      = 0.1;
	fMCMCRValueParametersCriterion = 0.1;
	fMCMCNIterationsConvergenceGlobal = -1;
	fMCMCFlagConvergenceGlobal = false;
	fMCMCRValue               = 100;
	fMCMCNIterationsUpdate    = 1000;
	fMCMCFlagOrderParameters  = true;
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetValuesDetail()
{
	fMCMCNChains              = 5;
	fMCMCNIterationsMax       = 1000000;
	fMCMCNIterationsRun       = 100000;
	fMCMCNIterationsBurnIn    = 0;
	fMCMCNIterationsPreRunMin = 500;
	fMCMCFlagInitialPosition  = 1;
	fMCMCRValueCriterion      = 0.1;
	fMCMCRValueParametersCriterion = 0.1;
	fMCMCNIterationsConvergenceGlobal = -1;
	fMCMCFlagConvergenceGlobal = false;
	fMCMCRValue               = 100;
	fMCMCNIterationsUpdate    = 1000;
	fMCMCFlagOrderParameters  = true;
}

// ---------------------------------------------------------
BCEngineMCMC::~BCEngineMCMC()
{
	// delete random number generator
	if (fMCMCRandom)
		delete fMCMCRandom;

	// delete 1-d marginalized distributions
	for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
		if (fMCMCH1Marginalized[i])
			delete fMCMCH1Marginalized[i];
	fMCMCH1Marginalized.clear();

	// delete 2-d marginalized distributions
	for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
		if (fMCMCH2Marginalized[i])
			delete fMCMCH2Marginalized[i];
	fMCMCH2Marginalized.clear();
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
TH1D * BCEngineMCMC::MCMCGetH1Marginalized(int index)
{
		// check index
	if (index < 0 || index >= int(fMCMCH1Marginalized.size()))
	{
		BCLog::OutWarning("BCEngineMCMC::MCMCGetH1Marginalized. Index out of range.");
		return 0;
	}

	return fMCMCH1Marginalized[index];
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

	// check index
	if (index < 0 || index >= int(fMCMCH2Marginalized.size()))
	{
		BCLog::OutWarning("BCEngineMCMC::MCMCGetH2Marginalized. Index out of range.");
		return 0;
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
void BCEngineMCMC::MCMCSetInitialPositions(std::vector<double> x0s)
{
	// clear the existing initial position vector
	fMCMCInitialPosition.clear();

	// copy the initial positions
	int n = int(x0s.size());

	for (int i = 0; i < n; ++i)
		fMCMCInitialPosition.push_back(x0s.at(i));

	// use these intial positions for the Markov chain
	this -> MCMCSetFlagInitialPosition(2);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetInitialPositions(std::vector< std::vector<double> > x0s)
{
	// create new vector
	std::vector <double> y0s;

	// loop over vector elements
	for (int i = 0; i < int(x0s.size()); ++i)
		for (int j = 0; j < int((x0s.at(i)).size()); ++j)
			y0s.push_back((x0s.at(i)).at(j));

	this -> MCMCSetInitialPositions(y0s);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetFlagFillHistograms(bool flag)
{
	fMCMCFlagFillHistograms = flag;

	for (int i = 0; i < fMCMCNParameters; ++i)
		fMCMCFlagsFillHistograms[i] = flag;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetFlagFillHistograms(int index, bool flag)
{
	// check if index is within range
	if (index < 0 || index > fMCMCNParameters)
	{
		BCLog::OutWarning("BCEngineMCMC :MCMCSetFlagFillHistograms. Index out of range.");
		return;
	}

	// set flag
	fMCMCFlagsFillHistograms[index] = flag;
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
void BCEngineMCMC::Copy(BCEngineMCMC & enginemcmc) const
{}

// --------------------------------------------------------
void BCEngineMCMC::MCMCTrialFunction(int chain, std::vector <double> &x)
{
	// use uniform distribution for now
//	for (int i = 0; i < fMCMCNParameters; ++i)
//		x[i] = fMCMCTrialFunctionScaleFactor[i] * 2.0 * (0.5 - fMCMCRandom -> Rndm());

	// Breit-Wigner with adjustable width
	for (int i = 0; i < fMCMCNParameters; ++i)
		x[i] = fMCMCRandom -> BreitWigner(0.0, fMCMCTrialFunctionScaleFactor[chain * fMCMCNParameters + i]);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCTrialFunctionSingle(int ichain, int iparameter, std::vector <double> &x)
{
	// no check of range for performance reasons

	// use uniform distribution
//	x[iparameter] = fMCMCTrialFunctionScaleFactor[ichain * fMCMCNParameters + iparameter] * 2.0 * (0.5 - fMCMCRandom -> Rndm());

	// Breit-Wigner width adjustable width
	x[iparameter] = fMCMCRandom -> BreitWigner(0.0, fMCMCTrialFunctionScaleFactor[ichain * fMCMCNParameters + iparameter]);
}

// --------------------------------------------------------
std::vector <double> BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(int ichain)
{
	// create a new vector with the length of fMCMCNParameters
	std::vector <double> x;

	// check if ichain is in range
	if (ichain < 0 || ichain >= fMCMCNChains)
		return x;

	// copy the scale factors into the temporary vector
	for (int j = 0; j < fMCMCNParameters; ++j)
		x.push_back(fMCMCTrialFunctionScaleFactor.at(ichain * fMCMCNParameters + j));

	return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(int ichain, int ipar)
{
	// check if ichain is in range
	if (ichain < 0 || ichain >= fMCMCNChains)
		return 0;

	// check if ipar is in range
	if (ipar < 0 || ipar >= fMCMCNParameters)
		return 0;

	// return component of ipar point in the ichain chain
	return fMCMCTrialFunctionScaleFactor.at(ichain *  fMCMCNChains + ipar);
}

// --------------------------------------------------------
std::vector <double> BCEngineMCMC::MCMCGetx(int ichain)
{
	// create a new vector with the length of fMCMCNParameters
	std::vector <double> x;

	// check if ichain is in range
	if (ichain < 0 || ichain >= fMCMCNChains)
		return x;

	// copy the point in the ichain chain into the temporary vector
	for (int j = 0; j < fMCMCNParameters; ++j)
		x.push_back(fMCMCx.at(ichain * fMCMCNParameters + j));

	return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetx(int ichain, int ipar)
{
	// check if ichain is in range
	if (ichain < 0 || ichain >= fMCMCNChains)
		return 0;

	// check if ipar is in range
	if (ipar < 0 || ipar >= fMCMCNParameters)
		return 0;

	// return component of jth point in the ith chain
	return fMCMCx.at(ichain *  fMCMCNParameters + ipar);
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetLogProbx(int ichain)
{
	// check if ichain is in range
	if (ichain < 0 || ichain >= fMCMCNChains)
		return -1;

	// return log of the probability at the current point in the ith chain
	return fMCMCLogProbx.at(ichain);
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(int chain, std::vector <double> &x)
{
	// get unscaled random point. this point might not be in the correct volume.
	this -> MCMCTrialFunction(chain, x);

	// get a proposal point from the trial function and scale it
	for (int i = 0; i < fMCMCNParameters; ++i)
		x[i] = fMCMCx[chain * fMCMCNParameters + i] + x[i] * (fMCMCBoundaryMax.at(i) - fMCMCBoundaryMin.at(i));

	// check if the point is in the correct volume.
	for (int i = 0; i < fMCMCNParameters; ++i)
		if ((x[i] < fMCMCBoundaryMin[i]) || (x[i] > fMCMCBoundaryMax[i]))
			return false;

	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(int chain, int parameter, std::vector <double> &x)
{
	// get unscaled random point in the dimension of the chosen
	// parameter. this point might not be in the correct volume.
	this -> MCMCTrialFunctionSingle(chain, parameter, x);

	// copy the old point into the new
	for (int i = 0; i < fMCMCNParameters; ++i)
		if (i != parameter)
			x[i] = fMCMCx[chain * fMCMCNParameters + i];

	// modify the parameter under study
	x[parameter] = fMCMCx[chain * fMCMCNParameters + parameter] + x[parameter] * (fMCMCBoundaryMax.at(parameter) - fMCMCBoundaryMin.at(parameter));

	// check if the point is in the correct volume.
	for (int i = 0; i < fMCMCNParameters; ++i)
		if ((x[i] < fMCMCBoundaryMin[i]) || (x[i] > fMCMCBoundaryMax[i]))
			return false;

	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(int chain, int parameter)
{
	// calculate index
	int index = chain * fMCMCNParameters;

	// increase counter
	fMCMCNIterations[chain]++;

	// get proposal point
	if (!this -> MCMCGetProposalPointMetropolis(chain, parameter, fMCMCxLocal))
	{
		// increase counter
		fMCMCNTrialsFalse[chain * fMCMCNParameters + parameter]++;

		return false;
	}

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

	return accept;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(int chain)
{
	// calculate index
	int index = chain * fMCMCNParameters;

	// increase counter
	fMCMCNIterations[chain]++;

	// get proposal point
	if (!this -> MCMCGetProposalPointMetropolis(chain, fMCMCxLocal))
	{
		// increase counter
		for (int i = 0; i < fMCMCNParameters; ++i)
			fMCMCNTrialsFalse[chain * fMCMCNParameters + i]++;

		return false;
	}

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

	return accept;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCUpdateStatisticsCheckMaximum()
{
	// loop over all chains
	for (int i = 0; i < fMCMCNChains; ++i)
	{
		// check if new maximum is found or chain is at the beginning
		if (fMCMCLogProbx[i] > fMCMCMaximumLogProb[i] || fMCMCNIterations[i] == 1)
		{
			// copy maximum value
			fMCMCMaximumLogProb[i] = fMCMCLogProbx[i];

			// copy mode of chain
			for (int j = 0; j < fMCMCNParameters; ++j)
				fMCMCMaximumPoints[i * fMCMCNParameters + j] = fMCMCx[i * fMCMCNParameters + j];
		}
	}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCUpdateStatisticsFillHistograms()
{
	// check if histograms are supposed to be filled
	if (!fMCMCFlagFillHistograms)
		return;

	// loop over chains
	for (int i = 0; i < fMCMCNChains; ++i)
	{
		// fill each 1-dimensional histogram (if supposed to be filled)
		for (int j = 0; j < fMCMCNParameters; ++j)
			if (fMCMCFlagsFillHistograms.at(j))
				fMCMCH1Marginalized[j] -> Fill(fMCMCx[i * fMCMCNParameters + j]);

		// fill each 2-dimensional histogram (if supposed to be filled)
		int counter = 0;

		for (int j = 0; j < fMCMCNParameters; ++j)
			for (int k = 0; k < j; ++k)
			{
				if (fMCMCFlagsFillHistograms.at(j) && fMCMCFlagsFillHistograms.at(k))
					fMCMCH2Marginalized[counter]->Fill(fMCMCx[i*fMCMCNParameters+k],fMCMCx[i* fMCMCNParameters+j]);
				counter ++;
			}
	}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCUpdateStatisticsTestConvergenceAllChains()
{
	// calculate number of entries in this part of the chain
	int npoints = fMCMCNTrialsTrue[0] + fMCMCNTrialsFalse[0];

	if (fMCMCNChains > 1 && npoints > 1)
	{
		// define flag for convergence
		bool flag_convergence = true;

		// loop over parameters
		for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters)
		{
			double sum = 0;
			double sum2 = 0;
			double sumv = 0;

			// loop over chains
			for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
			{
				// get parameter index
				int index = ichains * fMCMCNParameters + iparameters;

				// calculate mean value of each parameter in the chain for this part
				fMCMCMeanx[index] += (fMCMCx[index] - fMCMCMeanx[index]) / double(npoints);

				// calculate variance of each chain for this part
				fMCMCVariancex[index] = (1.0 - 1/double(npoints)) * fMCMCVariancex[index]
						+ (fMCMCx[index] - fMCMCMeanx[index]) * (fMCMCx[index] - fMCMCMeanx[index]) / double(npoints - 1);

				sum  += fMCMCMeanx[index];
				sum2 += fMCMCMeanx[index] * fMCMCMeanx[index];
				sumv += fMCMCVariancex[index];
			}

			// calculate r-value for each parameter
			double mean = sum / double(fMCMCNChains);
			double B = (sum2 / double(fMCMCNChains) - mean * mean) * double(fMCMCNChains) / double(fMCMCNChains-1) * double(npoints);
			double W = sumv * double(npoints) / double(npoints - 1) / double(fMCMCNChains);

			double r = 100.0;

			if (W > 0)
			{
				r = sqrt( ( (1-1/double(npoints)) * W + 1/double(npoints) * B ) / W);
				fMCMCRValueParameters[iparameters] = r;
			}
			// set flag to false if convergence criterion is not fulfilled for the parameter
			if (! ((fMCMCRValueParameters[iparameters]-1.0) < fMCMCRValueParametersCriterion))
				flag_convergence = false;
		}

		// convergence criterion applied on the log-likelihood
		double sum = 0;
		double sum2 = 0;
		double sumv = 0;

		// loop over chains
		for (int i = 0; i < fMCMCNChains; ++i)
		{
			// calculate mean value of each chain for this part
			fMCMCMean[i] += (fMCMCLogProbx[i] - fMCMCMean[i]) / double(npoints);

			// calculate variance of each chain for this part
			if (npoints > 1)
				fMCMCVariance[i] = (1.0 - 1/double(npoints)) * fMCMCVariance[i]
						+ (fMCMCLogProbx[i] - fMCMCMean[i]) * (fMCMCLogProbx[i] - fMCMCMean[i]) / double(npoints - 1);

			sum  += fMCMCMean[i];
			sum2 += fMCMCMean[i] * fMCMCMean[i]; ;
			sumv += fMCMCVariance[i];
		}

		// calculate r-value for log-likelihood
		double mean = sum / double(fMCMCNChains);
		double B = (sum2 / double(fMCMCNChains) - mean * mean) * double(fMCMCNChains) / double(fMCMCNChains-1) * double(npoints);
		double W = sumv * double(npoints) / double(npoints - 1) / double(fMCMCNChains);
		double r = 100.0;

		if (W > 0)
		{
			r = sqrt( ( (1-1/double(npoints)) * W + 1/double(npoints) * B ) / W);
			fMCMCRValue = r;
		}

		// set flag to false if convergence criterion is not fulfilled for the log-likelihood
		if (! ((fMCMCRValue-1.0) < fMCMCRValueCriterion))
			flag_convergence = false;

		// remember number of iterations needed to converge
		if (fMCMCNIterationsConvergenceGlobal == -1 && flag_convergence == true)
			fMCMCNIterationsConvergenceGlobal = fMCMCNIterations[0] / fMCMCNParameters;
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
double BCEngineMCMC::LogEval(std::vector <double> parameters)
{
	// test function for now
	// this will be overloaded by the user
	return 0.0;
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCMetropolisPreRun()
{
	// print on screen
	BCLog::OutSummary("Pre-run Metropolis MCMC...");

	// initialize Markov chain
	this -> MCMCInitialize();
	this -> MCMCInitializeMarkovChains();

	// helper variable containing number of digits in the number of parameters
	int ndigits = (int)log10(fMCMCNParameters) +1;
	if(ndigits<4)
		ndigits=4;

	// perform burn-in run
	if (fMCMCNIterationsBurnIn > 0)
		BCLog::OutDetail(Form(" --> Start burn-in run with %i iterations", fMCMCNIterationsBurnIn));

	// if the flag is set then run over the parameters one after the other.
	if (fMCMCFlagOrderParameters)
	{
		for (int i = 0; i < fMCMCNIterationsBurnIn; ++i)
			for (int j = 0; j < fMCMCNChains; ++j)
				for (int k = 0; k < fMCMCNParameters; ++k)
					this -> MCMCGetNewPointMetropolis(j, k);
	}
	// if the flag is not set then run over the parameters at the same time.
	else
	{
		for (int i = 0; i < fMCMCNIterationsBurnIn; ++i)
			// loop over chains
			for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
				// get new point
				this -> MCMCGetNewPointMetropolis(ichains);
	}

	// reset run statistics
	this -> MCMCResetRunStatistics();
	fMCMCNIterationsConvergenceGlobal = -1;

	// perform run
	BCLog::OutSummary(Form(" --> Perform MCMC pre-run with %i chains, each with maximum %i iterations", fMCMCNChains, fMCMCNIterationsMax));


	// don't write to file during pre run
	bool tempflag_writetofile = fMCMCFlagWriteChainToFile;
	fMCMCFlagWriteChainToFile = false;

	// initialize counter variables and flags
	int counter = 0;
	int counterupdate = 0;
	bool convergence = false;     // convergence reached?
	bool flagefficiency = false;  // efficiency reached?

	// array of efficiencies
	std::vector <double> efficiency;
	efficiency.assign(fMCMCNParameters * fMCMCNChains, 0.0);

	// how often to check convergence and efficiencies?
	int updateLimit = ( fMCMCNIterationsUpdateMax<fMCMCNIterationsUpdate*(fMCMCNParameters+1)  && fMCMCNIterationsUpdateMax>0 ) ?
			fMCMCNIterationsUpdateMax : fMCMCNIterationsUpdate*(fMCMCNParameters+1);

	// run chain ...
	// (a) for at least a minimum number of iterations,
	// (b) until a maximum number of iterations is reached,
	// (c) or until convergence is reached and the efficiency is in the
	//     specified region
	while (counter < fMCMCNIterationsPreRunMin || (counter < fMCMCNIterationsMax && !(convergence && flagefficiency)))
	{
		// set convergence to false by default
		convergence = false;

		// set number of iterations needed to converge to negative
		fMCMCNIterationsConvergenceGlobal = -1;

		// if the flag is set then run over the parameters one after the other.
		if (fMCMCFlagOrderParameters)
		{
			// loop over parameters
			for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters)
			{
				// loop over chains
				for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
					this -> MCMCGetNewPointMetropolis(ichains, iparameters);

				// search for global maximum
				this -> MCMCUpdateStatisticsCheckMaximum();
			} // end loop over parameters
		} // end if fMCMCFlagOrderParameters

		// if the flag is not set then run over the parameters at the same time.
		else
		{
			// loop over chains
			for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
				// get new point
				this -> MCMCGetNewPointMetropolis(ichains);

			// search for global maximum
			this -> MCMCUpdateStatisticsCheckMaximum();
		}

		// progress printout
		if ( counter > 0 && counter % fMCMCNIterationsUpdate == 0 )
			BCLog::OutDetail(Form(" --> Iteration %i", fMCMCNIterations[0]/fMCMCNParameters - 1));

		//-------------------------------------------
		// update scale factors and check convergence
		//-------------------------------------------
		if (counterupdate % updateLimit == 0 && counterupdate > 0 && counter >= fMCMCNIterationsPreRunMin)
		{
			bool rvalues_ok = true;
			static bool has_converged = false;

			// -----------------------------
			// check convergence status
			// -----------------------------

			// reset the number of iterations needed for convergence to
			// negative
			fMCMCNIterationsConvergenceGlobal = -1;

			this -> MCMCUpdateStatisticsTestConvergenceAllChains();

			// check convergence
			if (fMCMCNIterationsConvergenceGlobal > 0)
				convergence = true;

			// print convergence status:
			if (convergence && fMCMCNChains > 1)
				BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains converged within %i iterations.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));
			else if (!convergence && fMCMCNChains > 1)
			{
				BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains did not converge after %i iterations.", fMCMCNChains, counter));

				BCLog::OutDetail("       - R-Values:");
				for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter)
				{
					if(fabs(fMCMCRValueParameters[iparameter]-1.) < fMCMCRValueParametersCriterion)
						BCLog::OutDetail(Form( TString::Format("         parameter %%%di :  %%.06f",ndigits), iparameter, fMCMCRValueParameters.at(iparameter)));
					else
					{
						BCLog::OutDetail(Form( TString::Format("         parameter %%%di :  %%.06f <--",ndigits), iparameter, fMCMCRValueParameters.at(iparameter)));
						rvalues_ok = false;
					}
				}
				if(fabs(fMCMCRValue-1.) < fMCMCRValueCriterion)
					BCLog::OutDetail(Form("         log-likelihood :  %.06f", fMCMCRValue));
				else
				{
					BCLog::OutDetail(Form("         log-likelihood :  %.06f <--", fMCMCRValue));
					rvalues_ok = false;
				}
			}

			if(!has_converged)
				if(rvalues_ok)
					has_converged = true;

			// -----------------------------
			// check efficiency status
			// -----------------------------

			// set flag
			flagefficiency = true;

			bool flagprintefficiency = true;

			// loop over chains
			for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
			{
				// loop over parameters
				for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter)
				{
					// global index of the parameter (throughout all the chains)
					int ipc = ichains * fMCMCNParameters + iparameter;

					// calculate efficiency
					efficiency[ipc] = double(fMCMCNTrialsTrue[ipc]) / double(fMCMCNTrialsTrue[ipc] + fMCMCNTrialsFalse[ipc]);

					// adjust scale factors if efficiency is too low
					if (efficiency[ipc] < fMCMCEfficiencyMin && fMCMCTrialFunctionScaleFactor[ipc] > .01)
					{
						if (flagprintefficiency)
						{
							BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range."));
							BCLog::OutDetail(Form("       - Efficiencies:"));
							flagprintefficiency = false;
						}

						double fscale=2.;
						if(has_converged && fMCMCEfficiencyMin/efficiency[ipc] > 2.)
							fscale = 4.;
						fMCMCTrialFunctionScaleFactor[ipc] /= fscale;

						BCLog::OutDetail(Form("         Efficiency of parameter %i dropped below %.2lf%% (eps = %.2lf%%) in chain %i. Set scale to %.4g",
								iparameter, 100. * fMCMCEfficiencyMin, 100. * efficiency[ipc], ichains, fMCMCTrialFunctionScaleFactor[ipc]));
					}

					// adjust scale factors if efficiency is too high
					else if (efficiency[ipc] > fMCMCEfficiencyMax && fMCMCTrialFunctionScaleFactor[ipc] < 1.0)
					{
						if (flagprintefficiency)
						{
							BCLog::OutDetail(Form("   * Efficiency status: Efficiencies not within pre-defined ranges."));
							BCLog::OutDetail(Form("     - Efficiencies:"));
							flagprintefficiency = false;
						}

						fMCMCTrialFunctionScaleFactor[ipc] *= 2.;

						BCLog::OutDetail(Form("         Efficiency of parameter %i above %.2lf%% (eps = %.2lf%%) in chain %i. Set scale to %.4g",
								iparameter, 100.0 * fMCMCEfficiencyMax, 100.0 * efficiency[ipc], ichains, fMCMCTrialFunctionScaleFactor[ipc]));
					}

					// check flag
					if ((efficiency[ipc] < fMCMCEfficiencyMin && fMCMCTrialFunctionScaleFactor[ipc] > .01)
							|| (efficiency[ipc] > fMCMCEfficiencyMax && fMCMCTrialFunctionScaleFactor[ipc] < 1.))
						flagefficiency = false;

					// reset counters
					counterupdate = 0;
					fMCMCNTrialsTrue[ipc] = 0;
					fMCMCNTrialsFalse[ipc] = 0;
				} // end of running over all parameters

				// reset counters
				fMCMCMean[ichains] = 0;
				fMCMCVariance[ichains] = 0;
			} // end of running over all chains

			// print to scrren
			if (flagefficiency)
				BCLog::OutDetail(Form("     * Efficiency status: Efficiencies within pre-defined ranges."));

		} // end if update scale factors and check convergence

		// increase counter
		counter++;
		counterupdate++;
	} // end of running

	// did we check convergence at least once ?
	if (counter-1<updateLimit)
	{
		BCLog::OutWarning(" Convergence never checked !");
		BCLog::OutWarning("   Increase maximum number of iterations in the pre-run /MCMCSetNIterationsMax()/");
		BCLog::OutWarning("   or decrease maximum number of iterations for update  /MCMCSetNIterationsUpdateMax()/");
	}

	// ---------------
	// after chain run
	// ---------------

	// define convergence status
	if (fMCMCNIterationsConvergenceGlobal > 0)
		fMCMCFlagConvergenceGlobal = true;
	else
		fMCMCFlagConvergenceGlobal = false;

	// print convergence status
	if (fMCMCFlagConvergenceGlobal && fMCMCNChains > 1 && !flagefficiency)
		BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations but could not adjust scales.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));

	else if (fMCMCFlagConvergenceGlobal && fMCMCNChains > 1 && flagefficiency)
		BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations and could adjust scales.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));

	else if (!fMCMCFlagConvergenceGlobal && (fMCMCNChains > 1) && flagefficiency)
		BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations.", fMCMCNChains, fMCMCNIterationsMax));

	else if (!fMCMCFlagConvergenceGlobal && (fMCMCNChains > 1) && !flagefficiency)
		BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations and could not adjust scales.", fMCMCNChains, fMCMCNIterationsMax));

	else if(fMCMCNChains == 1)
		BCLog::OutSummary(" --> No convergence criterion for a single chain defined.");

	else
		BCLog::OutSummary(" --> Only one Markov chain. No global convergence criterion defined.");

	BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", counter));


	// print efficiencies
	std::vector <double> efficiencies;

	for (int i = 0; i < fMCMCNParameters; ++i)
		efficiencies.push_back(0.);

	BCLog::OutDetail(" --> Average efficiencies:");
	for (int i = 0; i < fMCMCNParameters; ++i)
	{
		for (int j = 0; j < fMCMCNChains; ++j)
			efficiencies[i] += efficiency[j * fMCMCNParameters + i] / double(fMCMCNChains);

		BCLog::OutDetail(Form(TString::Format(" -->      parameter %%%di :  %%.02f%%%%",ndigits), i, 100. * efficiencies[i]));
	}


	// print scale factors
	std::vector <double> scalefactors;

	for (int i = 0; i < fMCMCNParameters; ++i)
		scalefactors.push_back(0.0);

	BCLog::OutDetail(" --> Average scale factors:");
	for (int i = 0; i < fMCMCNParameters; ++i)
	{
		for (int j = 0; j < fMCMCNChains; ++j)
			scalefactors[i] += fMCMCTrialFunctionScaleFactor[j * fMCMCNParameters + i] / double(fMCMCNChains);

		BCLog::OutDetail(Form( TString::Format(" -->      parameter %%%di :  %%.02f%%%%",ndigits), i, 100. * scalefactors[i]));
	}

	// reset flag
	fMCMCFlagWriteChainToFile = tempflag_writetofile;

	// set pre-run flag
	fMCMCFlagPreRun = true;

	// no error
	return 1;
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCMetropolis()
{
	// check if prerun has been performed
	if (!fMCMCFlagPreRun)
		this -> MCMCMetropolisPreRun();

	// print to screen
	BCLog::OutSummary( "Run Metropolis MCMC...");

	// reset run statistics
	this -> MCMCResetRunStatistics();

	// perform run
	BCLog::OutSummary(Form(" --> Perform MCMC run with %i chains, each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun));

//	int counterupdate = 0;
//	bool convergence = false;
//	bool flagefficiency = false;

//	std::vector <double> efficiency;

//	for (int i = 0; i < fMCMCNParameters; ++i)
//		for (int j = 0; j < fMCMCNChains; ++j)
//			efficiency.push_back(0.0);

	int nwrite = fMCMCNIterationsRun/10;
	if(nwrite < 100)
		nwrite=100;
	else if(nwrite < 500)
		nwrite=1000;
	else if(nwrite < 10000)
		nwrite=1000;
	else
		nwrite=10000;

	// start the run
	for (int iiterations = 0; iiterations < fMCMCNIterationsRun; ++iiterations)
	{
		if ( (iiterations+1)%nwrite == 0 )
			BCLog::OutDetail(Form(" --> iteration number %i (%.2f%%)", iiterations+1, (double)(iiterations+1)/(double)fMCMCNIterationsRun*100.));

		// if the flag is set then run over the parameters one after the other.
		if (fMCMCFlagOrderParameters)
		{
			// loop over parameters
			for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters)
			{
				// loop over chains
				for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
					this -> MCMCGetNewPointMetropolis(ichains, iparameters);

				// update search for maximum
				this -> MCMCUpdateStatisticsCheckMaximum();

				// check if the current iteration is consistent with the lag
				if ( (fMCMCNParameters * iiterations + iparameters) % (fMCMCNLag * fMCMCNParameters) == 0)
				{
					// fill histograms
					this -> MCMCUpdateStatisticsFillHistograms();

					// write chain to file
					if (fMCMCFlagWriteChainToFile)
						this -> MCMCUpdateStatisticsWriteChains();

					// do anything interface
					this -> MCMCIterationInterface();
				}

			} // end loop over all parameters
		}
		// if the flag is not set then run over the parameters at the same time.
		else
		{
			// loop over chains
			for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
				// get new point
				this -> MCMCGetNewPointMetropolis(ichains);

			// update search for maximum
			this -> MCMCUpdateStatisticsCheckMaximum();

			// check if the current iteration is consistent with the lag
			if (iiterations % fMCMCNLag == 0)
			{
				// fill histograms
				this -> MCMCUpdateStatisticsFillHistograms();

				// write chain to file
				if (fMCMCFlagWriteChainToFile)
					this -> MCMCUpdateStatisticsWriteChains();

				// do anything interface
				this -> MCMCIterationInterface();
			}
		}

	} // end run

	// print convergence status
	BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCNIterationsRun));

	// print modes

	// find global maximum
	double probmax = fMCMCMaximumLogProb.at(0);
	int probmaxindex = 0;

	// loop over all chains and find the maximum point
	for (int i = 1; i < fMCMCNChains; ++i)
		if (fMCMCMaximumLogProb.at(i) > probmax)
		{
			probmax = fMCMCMaximumLogProb.at(i);
			probmaxindex = i;
		}

	BCLog::OutDetail(" --> Global mode from MCMC:");
	int ndigits = (int) log10(fMCMCNParameters);
	for (int i = 0; i < fMCMCNParameters; ++i)
		BCLog::OutDetail(Form( TString::Format(" -->      parameter %%%di:   %%.4g", ndigits+1),
				i, fMCMCMaximumPoints[probmaxindex * fMCMCNParameters + i]));

	// set flags
	fMCMCFlagPreRun = false;
	fMCMCFlagRun = true;

	return 1;
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

		for (int k = 0; k < fMCMCNParameters; ++k)
		{
			fMCMCNTrialsTrue[j * fMCMCNParameters + k]  = 0;
			fMCMCNTrialsFalse[j * fMCMCNParameters + k] = 0;
		}
	}

	// reset marginalized distributions
	for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
		if (fMCMCH1Marginalized[i])
			fMCMCH1Marginalized[i] -> Reset();

	for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
		if (fMCMCH2Marginalized[i])
			fMCMCH2Marginalized[i] -> Reset();

	fMCMCRValue = 100;
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCAddParameter(double min, double max)
{
	// add the boundaries to the corresponding vectors
	fMCMCBoundaryMin.push_back(min);
	fMCMCBoundaryMax.push_back(max);

	// set flag for individual parameters
	fMCMCFlagsFillHistograms.push_back(true);

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
			x0.push_back(fMCMCx[j * fMCMCNParameters + i]);
		fMCMCLogProbx[j] = this -> LogEval(x0);
	}

	x0.clear();
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCInitialize()
{
	// reset variables
	fMCMCNIterations.clear();
	fMCMCNTrialsTrue.clear();
	fMCMCNTrialsFalse.clear();
	fMCMCTrialFunctionScaleFactor.clear();
	fMCMCMean.clear();
	fMCMCVariance.clear();
	fMCMCMeanx.clear();
	fMCMCVariancex.clear();
	fMCMCx.clear();
	fMCMCLogProbx.clear();
	fMCMCMaximumPoints.clear();
	fMCMCMaximumLogProb.clear();
	fMCMCNIterationsConvergenceGlobal = -1;
	fMCMCFlagConvergenceGlobal = false;
	fMCMCRValueParameters.clear();

	for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
		if (fMCMCH1Marginalized[i])
			delete fMCMCH1Marginalized[i];

	for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
		if (fMCMCH2Marginalized[i])
			delete fMCMCH2Marginalized[i];

	// clear plots
	fMCMCH1Marginalized.clear();
	fMCMCH2Marginalized.clear();

// free memory for vectors
	fMCMCNIterations.assign(fMCMCNChains, 0);
	fMCMCMean.assign(fMCMCNChains, 0);
	fMCMCVariance.assign(fMCMCNChains, 0);
	fMCMCLogProbx.assign(fMCMCNChains, -1.0);
	fMCMCMaximumLogProb.assign(fMCMCNChains, -1.0);

	fMCMCNTrialsTrue.assign(fMCMCNChains * fMCMCNParameters, 0);
	fMCMCNTrialsFalse.assign(fMCMCNChains * fMCMCNParameters, 0);
	fMCMCMaximumPoints.assign(fMCMCNChains * fMCMCNParameters, 0.);
	fMCMCMeanx.assign(fMCMCNChains * fMCMCNParameters, 0);
	fMCMCVariancex.assign(fMCMCNChains * fMCMCNParameters, 0);

	fMCMCRValueParameters.assign(fMCMCNParameters, 0.);

	if (fMCMCTrialFunctionScaleFactorStart.size() == 0)
		fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains * fMCMCNParameters, 1.0);
	else
		for (int i = 0; i < fMCMCNChains; ++i)
			for (int j = 0; j < fMCMCNParameters; ++j)
				fMCMCTrialFunctionScaleFactor.push_back(fMCMCTrialFunctionScaleFactorStart.at(j));

	// set initial position
	if (fMCMCFlagInitialPosition == 2) // user defined points
	{
		// define flag
		bool flag = true;

		// check the length of the array of initial positions
		if (int(fMCMCInitialPosition.size()) != (fMCMCNChains * fMCMCNParameters))
		{
			BCLog::OutError("BCEngine::MCMCInitialize : Length of vector containing initial positions does not have required length.");
			flag = false;
		}

		// check the boundaries
		if (flag)
		{
			for (int j = 0; j < fMCMCNChains; ++j)
				for (int i = 0; i < fMCMCNParameters; ++i)
					if (fMCMCInitialPosition[j * fMCMCNParameters + i] < fMCMCBoundaryMin[i] ||
							fMCMCInitialPosition[j * fMCMCNParameters + i] > fMCMCBoundaryMax[i])
					{
						BCLog::OutError("BCEngine::MCMCInitialize : Initial position out of boundaries.");
						flag = false;
					}
		}

		// check flag
		if (!flag)
			fMCMCFlagInitialPosition = 1;
	}

	if (fMCMCFlagInitialPosition == 0) // center of the interval
		for (int j = 0; j < fMCMCNChains; ++j)
			for (int i = 0; i < fMCMCNParameters; ++i)
				fMCMCx.push_back(fMCMCBoundaryMin[i] + .5 * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i]));

	else if (fMCMCFlagInitialPosition == 2) // user defined
	{
		for (int j = 0; j < fMCMCNChains; ++j)
			for (int i = 0; i < fMCMCNParameters; ++i)
				fMCMCx.push_back(fMCMCInitialPosition.at(j * fMCMCNParameters + i));
	}

	else
		for (int j = 0; j < fMCMCNChains; ++j) // random number (default)
			for (int i = 0; i < fMCMCNParameters; ++i)
				fMCMCx.push_back(fMCMCBoundaryMin[i] + fMCMCRandom->Rndm() * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i]));

	// copy the point of the first chain
	fMCMCxLocal.clear();
	for (int i = 0; i < fMCMCNParameters; ++i)
		fMCMCxLocal.push_back(fMCMCx[i]);

	// define 1-dimensional histograms for marginalization
	for(int i = 0; i < fMCMCNParameters; ++i)
	{
		TH1D * h1 = 0;
		if (fMCMCFlagsFillHistograms.at(i))
			h1 = new TH1D(TString::Format("h1_%d_parameter_%i", BCLog::GetHIndex() ,i), "",
					fMCMCH1NBins[i], fMCMCBoundaryMin[i], fMCMCBoundaryMax[i]);
		fMCMCH1Marginalized.push_back(h1);
	}

	for(int i = 0; i < fMCMCNParameters; ++i)
		for (int k = 0; k < i; ++k)
		{
			TH2D * h2 = 0;
			if (fMCMCFlagsFillHistograms.at(i) && fMCMCFlagsFillHistograms.at(k))
				h2 = new TH2D(Form("h2_%d_parameters_%i_vs_%i", BCLog::GetHIndex(), i, k), "",
						fMCMCH1NBins[k], fMCMCBoundaryMin[k], fMCMCBoundaryMax[k],
						fMCMCH1NBins[i], fMCMCBoundaryMin[i], fMCMCBoundaryMax[i] );
			fMCMCH2Marginalized.push_back(h2);
		}

	fMCMCFlagPreRun = false;
	fMCMCFlagRun = false;

	return 1;
}

// ---------------------------------------------------------
int BCEngineMCMC::SetMarginalized(int index, TH1D * h)
{
	if((int)fMCMCH1Marginalized.size()<=index)
		return 0;

	if(h==0)
		return 0;

	if((int)fMCMCH1Marginalized.size()==index)
		fMCMCH1Marginalized.push_back(h);
	else
		fMCMCH1Marginalized[index]=h;

	return index;
}

// ---------------------------------------------------------
int BCEngineMCMC::SetMarginalized(int index1, int index2, TH2D * h)
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

	if((int)fMCMCH2Marginalized.size()<=index)
		return 0;

	if(h==0)
		return 0;

	if((int)fMCMCH2Marginalized.size()==index)
		fMCMCH2Marginalized.push_back(h);
	else
		fMCMCH2Marginalized[index]=h;

	return index;
}

// ---------------------------------------------------------

