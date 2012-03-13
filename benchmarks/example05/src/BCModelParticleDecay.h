// --------------------------------------------------------- 
//
// This class derives from BCModel. It describes a linear correlation
// relation between measured points. Two parameters are defined within
// the model, an offset, b, and a slope, m. The data are points (x,y)
// with an uncertainty on y, s_y. The uncertainty is assumed to be
// Gaussian. 
// 
// --------------------------------------------------------- 

#ifndef __BCMODELPARTICLEDECAY__H
#define __BCMODELPARTICLEDECAY__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelParticleDecay : public BCModel 
{

 public: 

	// constructor 

	BCModelParticleDecay(); 

	BCModelParticleDecay(const char* name); 

	// destructor 

	~BCModelParticleDecay()
		{ ;};  

	// methods 

	virtual void DefineParameters(); 

	virtual double LogAPrioriProbability(std::vector<double> parameters); 

	virtual double LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector<double> parameters); 

	virtual double LogPoissonProbability(int nentries, std::vector<double> parameters); 

 private: 

	double logEnergyResolution(double Emeasured, double Etrue); 

	double logMassDistribution(double mass, double polemass, double width); 

}; 

// --------------------------------------------------------- 

#endif 
