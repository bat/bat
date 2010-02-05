#ifndef __ENSEMBLETESTTOOL__H
#define __ENSEMBLETESTTOOL__H

/*!
 * \class EnsembleTestTool
 * This class can be used for ensemble tests in template method
 * analyses.  For the template fit, the "StackTool" is used. The
 * fitting can be done with Minuit or with Markov Chains. These
 * Methods are provided by the BAT package.
 *
 * \brief A class for doing ensemble tests.
 * \author Andrea Knue, Kevin Kroeninger
 * \version 0.1
 * \date 25.12.2009
 *
 */

// --------------------------------------------------------------------------------------------

#include <BAT/BCModel.h>
#include <StackModel.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h> 

// --------------------------------------------------------------------------------------------

class EnsembleTestTool 
{
 public: 

  /**
   * The constructor.
   */ 
  EnsembleTestTool();

  /**
   * The destructor.
   */
  ~EnsembleTestTool();

  /**
   * Set the template used to generate the ensembles. 
   */
  int SetEnsembleTemplate(TH1D hist); 

	/**
	 * Set the StackModel used to analyze the ensembles.
	 */ 
	void SetStackModel(StackModel* model)
	{ fStackModel = model; }; 

  /**
   * A function to perform an ensemble test for each data set in the container. 
   */
  int PerformEnsembleTest();

  /**
   * A function to define the number of events per ensemble.
   */
  void SetEnsembleExpectation(double expectation)
	{ fEnsembleExpectation = expectation; };

  /**
   * A function to define the number of ensembles per data set.
   */
  void SetNEnsembles(int n)
	{ fNEnsembles = n; }; 

  /**
   * A function to set the MCMC flag.
   */
  void SetFlagMCMC(bool flag)
	{ fFlagMCMC = flag; }; 

	/**
	 * Write tree to file 
	 */ 
	int Write(const char* filename); 

	/**
	 * Prepare tree.
	 */ 
	int PrepareTree(); 

 private:
	
  /**
   * Create a new ensemble.
	 * @return A histogram with the new ensemble.
   */
  TH1D BuildEnsemble();

 protected:

  /**
   * The template used for the generation of ensembles.
   */
	TH1D fEnsembleTemplate; 

	/** 
	 * The stack model used to analyze the ensembles.
	 */ 
	StackModel* fStackModel; 

	/**
	 * Output file.
	 */ 
	TFile* fFile; 

	/**
	 * Output tree.
	 */ 
	TTree* fTree; 

	/**
	 * A counter for the number of ensembles. 
	 */ 
	int fEnsembleCounter; 

  /**
   * Exepectation value 
   */
  double fEnsembleExpectation;

  /**
   * Number of ensembles per data set.
   */
  int fNEnsembles;

  /**
   * A flag to turn the Markov Chains on.
   */
  bool fFlagMCMC;

	/**
	 * The random number generator. 
	 */ 
	TRandom3* fRandom; 

	/**
	 * Tree variable: the global mode
	 */ 
	std::vector<double> fOutModeGlobal; 

	/**
	 * Tree variable: positive uncertainty on global mode 
	 */ 
	std::vector<double> fOutErrorUpGlobal; 

	/**
	 * Tree variable: negative uncertainty on global mode
	 */ 
	std::vector<double> fOutErrorDownGlobal; 

	/**
	 * Tree variable: marginalized mode
	 */ 
	std::vector<double> fOutModeMarg; 

	/**
	 * Tree variable: median
	 */ 
	std::vector<double> fOutMedianMarg;

	/**
	 * Tree variable: mean
	 */ 
	std::vector<double> fOutMeanMarg;

	/**
	 * Tree variable: rms
	 */ 
	std::vector<double> fOutRMSMarg;

	/**
	 * Tree variable: 16% quantile
	 */ 
	std::vector<double> fOutErrorUpMarg;

	/**
	 * Tree variable: 84% quantile
	 */ 
	std::vector<double> fOutErrorDownMarg;

	/**
	 * Tree variable: 5% quantile
	 */ 
	std::vector<double> fOutQuantile5Marg;

	/**
	 * Tree variable: 10% quantile
	 */ 
	std::vector<double> fOutQuantile10Marg;

	/**
	 * Tree variable: 90% quantile
	 */ 
	std::vector<double> fOutQuantile90Marg;

	/**
	 * Tree variable: 95% quantile
	 */ 
	std::vector<double> fOutQuantile95Marg;

	/**
	 * Tree variable: chi2
	 */ 
	double fOutChi2;

	/**
	 * Tree variable: ndf
	 */ 
	int fOutNDF; 

	/**
	 * Tree variable: chi2-probability
	 */ 
	double fOutChi2Prob;

	/**
	 * Tree variable: KL probability
	 */ 
	double fOutKSProb;

	/**
	 * Tree variable: p-value
	 */ 
	double fOutPValue;

	/**
	 * Tree variable: number of events in the data
	 */ 
	int fOutNEvents;
};

#endif  
