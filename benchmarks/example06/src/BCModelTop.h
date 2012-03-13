// --------------------------------------------------------- 
//
// This class derives from BCModel. It describes a linear correlation
// relation between measured points. Two parameters are defined within
// the model, an offset, b, and a slope, m. The data are points (x,y)
// with an uncertainty on y, s_y. The uncertainty is assumed to be
// Gaussian. 
// 
// --------------------------------------------------------- 

#ifndef __BCMODELTOP__H
#define __BCMODELTOP__H

#include "BCModel.h" 
#include "BCDataPoint.h" 
#include "BCDataSet.h" 

#include <vector> 

#include <TLorentzVector.h> 

// --------------------------------------------------------- 

class BCModelTop : public BCModel 
{

 public: 

	// constructor 

	BCModelTop(); 

	BCModelTop(const char* name); 

	// destructor 

	~BCModelTop()
		{ ;};  

	// getters 

	TLorentzVector * GetLorentzVector(int index) 
		{ return fLorentzVectorContainer.at(index); }; 

	// methods 

	void DefineParameters(); 

	double LogAPrioriProbability(std::vector<double> parameters); 

	virtual double LogLikelihood(std::vector<double> parameters);

	double EnergyResolutionBJets(double Emeasured, double Eestimated); 
	double EnergyResolutionLightJets(double Emeasured, double Eestimated);
	double EnergyResolutionElectrons(double Emeasured, double Eestimated); 

	void MCMCUserInterface();

	void InitializeEvent(BCDataSet * dataset, int index); 

	void SetPermutation(int index); 
	void CalculateLorentzVectors(std::vector<double> parameters); 

	int fPermutationTable[12][4]; 

	TLorentzVector fLV_bhad; 
	TLorentzVector fLV_blep; 
	TLorentzVector fLV_qup; 
	TLorentzVector fLV_qdown; 
	TLorentzVector fLV_electron;
	TLorentzVector fLV_neutrino; 
	TLorentzVector fLV_Whad; 
	TLorentzVector fLV_Wlep; 
	TLorentzVector fLV_Tophad; 
	TLorentzVector fLV_Toplep; 

	std::vector<TLorentzVector *> fLorentzVectorContainer; 

	double fMb; 
	double fMW; 
	double fGammaW; 

	BCDataPoint * fDataPoint; 

	BCDataPoint * fDataPointFull; 

}; 

// --------------------------------------------------------- 

#endif 
