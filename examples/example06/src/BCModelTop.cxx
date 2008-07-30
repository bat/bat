#include "BCModelTop.h" 

// --------------------------------------------------------- 

BCModelTop::BCModelTop() : BCModel()
{

	// define parameters 
	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

BCModelTop::BCModelTop(const char* name) : BCModel(name)
{

	// define parameters 
	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

void BCModelTop::DefineParameters()
{

	fLorentzVectorContainer.push_back(&fLV_bhad); 
	fLorentzVectorContainer.push_back(&fLV_blep); 
	fLorentzVectorContainer.push_back(&fLV_qup); 
	fLorentzVectorContainer.push_back(&fLV_qdown); 
	fLorentzVectorContainer.push_back(&fLV_electron); 
	fLorentzVectorContainer.push_back(&fLV_neutrino); 
	fLorentzVectorContainer.push_back(&fLV_Whad); 
	fLorentzVectorContainer.push_back(&fLV_Wlep); 
	fLorentzVectorContainer.push_back(&fLV_Tophad); 
	fLorentzVectorContainer.push_back(&fLV_Toplep); 

	// create permutation table 
	int table_temp[12][4] = { 
		{0,1,2,3}, // all correct 
		{2,1,3,0}, // leptonic correct
		{3,1,2,0}, // leptonic correct 
		{1,0,2,3}, // hadronic W correct
		{0,3,1,2}, // all wrong 
		{0,2,3,1}, // all wrong 
		{1,3,0,2}, // all wrong 
		{1,2,3,0}, // all wrong 
		{2,0,1,3}, // all wrong
		{2,3,0,1}, // all wrong 
		{3,0,1,2}, // all wrong 
		{3,2,0,1}, // all wrong 
	};

	// fill permutation table 
	for (int i = 0; i < 12; ++i)
		for (int j = 0; j < 4; ++j)
			fPermutationTable[i][j] = table_temp[i][j]; 

	// set constants 
	fMb = 4.7; 
	fMW = 80.4; 
	fGammaW = 2.1; 

	// initialize variables 
	fDataPoint = 0; 

}

// --------------------------------------------------------- 

double BCModelTop::LogAPrioriProbability(std::vector <double> parameters)
{

	double logprob = 0.; 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelTop::LogLikelihood(std::vector <double> parameters)
{

	this -> CalculateLorentzVectors(parameters); 

	double logprob = 0.0; 
	
	logprob += this -> EnergyResolutionBJets(fDataPoint -> GetValue(0), 
																					 parameters.at(0)); 
 	logprob += this -> EnergyResolutionBJets(fDataPoint -> GetValue(4), 
 																					 parameters.at(1)); 
	logprob += this -> EnergyResolutionLightJets(fDataPoint -> GetValue(8) / parameters.at(6), 
																							 parameters.at(2)); 
	logprob += this -> EnergyResolutionLightJets(fDataPoint -> GetValue(12) / parameters.at(6), 
																							 parameters.at(3)); 
	logprob += this -> EnergyResolutionElectrons(fDataPoint -> GetValue(16), 
																							 parameters.at(4)); 

 	logprob += BCMath::LogBreitWignerNonRel(fLV_Whad.M(), fMW, fGammaW, true); 
 	logprob += BCMath::LogBreitWignerNonRel(fLV_Wlep.M(), fMW, fGammaW, true); 
	logprob += BCMath::LogGaus(fLV_Tophad.M() - fLV_Toplep.M(), 0.0, 1.5, true); 

	return logprob; 

}

// --------------------------------------------------------- 

void BCModelTop::InitializeEvent(BCDataSet * dataset, int index)
{

	// set data point 
	this -> SetSingleDataPoint(dataset, index); 

	// create data point 
	if (fDataPoint != 0)
		delete fDataPoint; 
	fDataPoint = new BCDataPoint(dataset -> GetDataPoint(index) -> GetNValues()); 

	// get point of current data point 
	*fDataPoint = *dataset -> GetDataPoint(index); 

	// chose measured values as starting values 
	std::vector <double> xstart; 
	xstart.push_back(fDataPoint -> GetValue(0)); 
	xstart.push_back(fDataPoint -> GetValue(4)); 
	xstart.push_back(fDataPoint -> GetValue(8)); 
	xstart.push_back(fDataPoint -> GetValue(12)); 
	xstart.push_back(fDataPoint -> GetValue(16)); 
	xstart.push_back(0.0); 
	xstart.push_back(1.0); 

	this -> MCMCSetInitialPositions(xstart); 

}

// --------------------------------------------------------- 

void BCModelTop::SetPermutation(int index)
{

	// copy temporary vectors from permutated table 
	double temp_bhad_E  = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][0]*4); 
	double temp_bhad_px = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][0]*4+1); 
	double temp_bhad_py = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][0]*4+2); 
	double temp_bhad_pz = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][0]*4+3); 

	double temp_blep_E  = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][1]*4); 
	double temp_blep_px = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][1]*4+1); 
	double temp_blep_py = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][1]*4+2); 
	double temp_blep_pz = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][1]*4+3); 

	double temp_qup_E  = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][2]*4); 
	double temp_qup_px = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][2]*4+1); 
	double temp_qup_py = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][2]*4+2); 
	double temp_qup_pz = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][2]*4+3); 

	double temp_qdown_E  = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][3]*4); 
	double temp_qdown_px = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][3]*4+1); 
	double temp_qdown_py = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][3]*4+2); 
	double temp_qdown_pz = this -> GetDataPoint(0) -> GetValue(fPermutationTable[index][3]*4+3); 

	// fill current data point 
	fDataPoint -> SetValue(0, temp_bhad_E); 
	fDataPoint -> SetValue(1, temp_bhad_px); 
	fDataPoint -> SetValue(2, temp_bhad_py); 
	fDataPoint -> SetValue(3, temp_bhad_pz); 

	fDataPoint -> SetValue(4, temp_blep_E); 
	fDataPoint -> SetValue(5, temp_blep_px); 
	fDataPoint -> SetValue(6, temp_blep_py); 
	fDataPoint -> SetValue(7, temp_blep_pz); 

	fDataPoint -> SetValue(8, temp_qup_E); 
	fDataPoint -> SetValue(9, temp_qup_px); 
	fDataPoint -> SetValue(10, temp_qup_py); 
	fDataPoint -> SetValue(11, temp_qup_pz); 

	fDataPoint -> SetValue(12, temp_qdown_E); 
	fDataPoint -> SetValue(13, temp_qdown_px); 
	fDataPoint -> SetValue(14, temp_qdown_py); 
	fDataPoint -> SetValue(15, temp_qdown_pz); 

}

// --------------------------------------------------------- 

void BCModelTop::CalculateLorentzVectors(std::vector <double> parameters) 
{

	// hadronic b quark 
	double E = parameters.at(0); 
	double px = fDataPoint -> GetValue(1); 
	double py = fDataPoint -> GetValue(2); 
	double pz = fDataPoint -> GetValue(3); 
	double pold = sqrt(px*px + py*py + pz*pz); 
	double pnew = sqrt(E*E - fMb*fMb); 

	fLV_bhad.SetPxPyPzE(px/pold*pnew, py/pold*pnew, pz/pold*pnew, E); 

	// leptonic b quark 
	E = parameters.at(1); 
	px = fDataPoint -> GetValue(5); 
	py = fDataPoint -> GetValue(6); 
	pz = fDataPoint -> GetValue(7); 
	pold = sqrt(px*px + py*py + pz*pz); 
	pnew = sqrt(E*E - fMb*fMb); 
	
	fLV_blep.SetPxPyPzE(px/pold*pnew, py/pold*pnew, pz/pold*pnew, E); 

	// light up quark 
	E = parameters.at(2); 
	px = fDataPoint -> GetValue(9); 
	py = fDataPoint -> GetValue(10); 
	pz = fDataPoint -> GetValue(11); 
	pold = sqrt(px*px + py*py + pz*pz); 
	pnew = E; 
	
	fLV_qup.SetPxPyPzE(px/pold*pnew, py/pold*pnew, pz/pold*pnew, E); 

	// light down quark 
	E = parameters.at(3); 
	px = fDataPoint -> GetValue(13); 
	py = fDataPoint -> GetValue(14); 
	pz = fDataPoint -> GetValue(15); 
	pold = sqrt(px*px + py*py + pz*pz); 
	pnew = E; 
	
	fLV_qdown.SetPxPyPzE(px/pold*pnew, py/pold*pnew, pz/pold*pnew, E); 

	// electron
	E = parameters.at(4); 
	px = fDataPoint -> GetValue(17); 
	py = fDataPoint -> GetValue(18); 
	pz = fDataPoint -> GetValue(19); 
	pold = sqrt(px*px + py*py + pz*pz); 
	pnew = E; 
	
	fLV_electron.SetPxPyPzE(px/pold*pnew, py/pold*pnew, pz/pold*pnew, E); 

	// neutrino
	// calculate from pT balance 
	px = -(fLV_bhad.Px() + 
				 fLV_blep.Px() + 
				 fLV_qup.Px() + 
				 fLV_qdown.Px() + 
				 fLV_electron.Px());

	py = -(fLV_bhad.Py() + 
				 fLV_blep.Py() + 
				 fLV_qup.Py() + 
				 fLV_qdown.Py() + 
				 fLV_electron.Py());
				 
	pz = parameters.at(5); 
	E =  sqrt(px*px + py*py + pz*pz); 

	fLV_neutrino.SetPxPyPzE(px, py, pz, E); 

	// calculate composite particles 
	fLV_Whad = fLV_qup + fLV_qdown; 
	fLV_Wlep = fLV_electron + fLV_neutrino; 
	fLV_Tophad = fLV_Whad + fLV_bhad; 
	fLV_Toplep = fLV_Wlep + fLV_blep; 

}

// --------------------------------------------------------- 

double BCModelTop::EnergyResolutionBJets(double Emeasured, double Eestimated) 
{

	double logprob = 0.0; 

	double sigma = 0.5 * sqrt(Eestimated); 

	logprob = BCMath::LogGaus(Eestimated, Emeasured, sigma, true); 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelTop::EnergyResolutionLightJets(double Emeasured, double Eestimated) 
{

	double logprob = 0.0; 

	double sigma = 0.5 * sqrt(Eestimated); 

	logprob = BCMath::LogGaus(Eestimated, Emeasured, sigma, true); 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelTop::EnergyResolutionElectrons(double Emeasured, double Eestimated) 
{

	double logprob = 0.0; 

	double sigma = 0.1 * sqrt(Eestimated); 

	logprob = BCMath::LogGaus(Eestimated, Emeasured, sigma, true); 

	return logprob; 

}

// --------------------------------------------------------- 

void BCModelTop::MCMCUserInterface() 
{
	
}

// --------------------------------------------------------- 

