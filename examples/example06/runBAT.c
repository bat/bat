#include <BCModelTop.h>
#include <BCLog.h>
#include <BCModelOutput.h> 

#include "style.c" 

#include <TCanvas.h> 
#include <TGraphErrors.h>
#include <TH2D.h> 
#include <TF1.h> 
#include "BCParameter.h" 

// ---------------------------------------------------------
  
const double sigmas = 2.0; 
const double cut_whad = 5.0; 
const double cut_wlep = 5.0; 
const double cut_deltat = 5.0; 

// ---------------------------------------------------------

int main()
{

	// ---------------------------------------------------------
	// set style  
	// ----------------------------------------------------------

	// calls a function which defines a nice style. 
	SetStyle(); 

	// ---------------------------------------------------------
	// histograms definitions 
	// ---------------------------------------------------------

	// initialize histograms 
	TH1D * hist_mass_Whad = new TH1D("hist_mass_Whad", ";m_{W} (hadr.) [GeV/c^{2};N", 100, 60.0, 100.0); 
	hist_mass_Whad -> SetStats(kFALSE); 

	TH1D * hist_mass_Wlep = new TH1D("hist_mass_Wlep", ";m_{W} (lept.) [GeV/c^{2};N", 100, 60.0, 100.0); 
	hist_mass_Wlep -> SetStats(kFALSE); 

	TH1D * hist_mass_Tophad = new TH1D("hist_mass_Tophad", ";m_{top} (hadr.) [GeV/c^{2};N", 100, 140.0, 200.0); 
	hist_mass_Tophad -> SetStats(kFALSE); 

	TH1D * hist_mass_Toplep = new TH1D("hist_mass_Toplep", ";m_{top} (lept.) [GeV/c^{2};N", 100, 140.0, 200.0); 
	hist_mass_Toplep -> SetStats(kFALSE); 

	TH1D * hist_JES_all = new TH1D("hist_JES_all", ";JES (light);N", 101, 0.495, 1.505); 
	hist_JES_all -> SetStats(kFALSE); 

	TH2D * hist_JES_all_mass_Tophad = new TH2D("hist_JES_all_mass_Tophad", ";JES (light);m_{top} (hadr.) [GeV/c^{2}]", 101, 0.495, 1.505, 61, 139.5, 200.5); 
	hist_JES_all_mass_Tophad -> SetStats(kFALSE); 

	TH1D * hist_combinations = new TH1D("hist_combinations", ";combination;N", 12, -0.5, 11.5); 
	hist_combinations -> SetStats(kFALSE); 

	TH2D * hist_sign_comb = new TH2D("hist_sign_comb", ";combination;significance", 12, -0.5, 11.5, 101, -0.005, 1.005); 
	hist_sign_comb -> SetStats(kFALSE); 

	TH2D * hist_best_comb = new TH2D("hist_best_comb", ";combination;ln p(x)_{max}", 12, -0.5, 11.5, 100, -20., 0.0); 
	hist_best_comb -> SetStats(kFALSE); 

	// ---------------------------------------------------------
	// open log file 
	// ---------------------------------------------------------

	// opens the log file. 
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail); 

	/*
	Need: 
	- energy resolution of all jets and the electron
	- mass resolution and offset 
	- incorporate JES 
	- b-tagging 
	- anything else? 
	*/

	// ---------------------------------------------------------
	// read data from file 
	// ---------------------------------------------------------

	// creates a new data set. 
	BCDataSet* fDataSet = new BCDataSet(); 

	// read in data from file 
	if (fDataSet -> ReadDataFromFileTree("input.root", "fTree", 
																			 //																			 "bhad_E,bhad_px,bhad_py,bhad_pz,blep_E,blep_px,blep_py,blep_pz,qup_E,qup_px,qup_py,qup_pz,qdown_E,qdown_px,qdown_py,qdown_pz,lcharged_E,lcharged_px,lcharged_py,lcharged_pz,lneutral_E,lneutral_px,lneutral_py,lneutral_pz,true_bhad_E,true_blep_E,true_qup_E,true_qdown_E,true_lcharged_E,true_lneutral_E") != 0)
																			 "bhad_E,bhad_px,bhad_py,bhad_pz,blep_E,blep_px,blep_py,blep_pz,qup_E,qup_px,qup_py,qup_pz,qdown_E,qdown_px,qdown_py,qdown_pz,lcharged_E,lcharged_px,lcharged_py,lcharged_pz,lneutral_E,lneutral_px,lneutral_py,lneutral_pz") != 0)
		//																			 "true_bhad_E,true_bhad_px,true_bhad_py,true_bhad_pz,true_blep_E,true_blep_px,true_blep_py,true_blep_pz,true_qup_E,true_qup_px,true_qup_py,true_qup_pz,true_qdown_E,true_qdown_px,true_qdown_py,true_qdown_pz,true_lcharged_E,true_lcharged_px,true_lcharged_py,true_lcharged_pz,true_lneutral_E,true_lneutral_px,true_lneutral_py,true_lneutral_pz,") != 0)
		return -1; 

	// define parameters 
	BCParameter * par_bhad_E = new BCParameter("bhad_E", 0.0, 700.0); 
	BCParameter * par_blep_E = new BCParameter("blep_E", 0.0, 700.0); 
	BCParameter * par_qup_E = new BCParameter("qup_E", 0.0, 700.0); 
	BCParameter * par_qdown_E = new BCParameter("qdown_E",0.0, 700.0); 
	BCParameter * par_electron_E = new BCParameter("electron_E", 0.0, 700.0); 
	// debug
	//	BCParameter * par_neutrino_pz = new BCParameter("neutrino_pz", -700.0, 700.0); 
	BCParameter * par_neutrino_pz = new BCParameter("neutrino_pz", -50.0, 200.0); 
	BCParameter * par_JES_all = new BCParameter("JES_all",  0.5, 1.5); 

	// reset counter 
	int nevents_processed = 0; 
	int nevents_minuiterror = 0; 
	int nevents_accepted = 0; 
	int nevents_correctcombination = 0; 

	// get number of events 
	int nevents = fDataSet -> GetNDataPoints();

	// debug
	nevents = 10; 

	// loop over events 
	for (int ievent = 0; ievent < nevents; ++ievent)
		{
			// debug
			std::cout << "pvalue - event number " << ievent << std::endl; 

			// increase counter of prcoessed events 
			nevents_processed++; 

			// define model 
			BCModelTop* fModelTop = new BCModelTop("myModelTop"); 

			// get current event
			fModelTop -> InitializeEvent(fDataSet, ievent); 

			// set boundaries on data 
			fModelTop -> SetDataBoundaries( 0, 0.0, 1000.0, false); 
			fModelTop -> SetDataBoundaries( 1, fDataSet -> GetDataPoint(ievent) -> GetValue(1), fDataSet -> GetDataPoint(ievent) -> GetValue(1), true); 
			fModelTop -> SetDataBoundaries( 2, fDataSet -> GetDataPoint(ievent) -> GetValue(2), fDataSet -> GetDataPoint(ievent) -> GetValue(2), true); 
			fModelTop -> SetDataBoundaries( 3, fDataSet -> GetDataPoint(ievent) -> GetValue(3), fDataSet -> GetDataPoint(ievent) -> GetValue(3), true); 

			fModelTop -> SetDataBoundaries( 4, 0.0, 1000.0, false); 
			fModelTop -> SetDataBoundaries( 5, fDataSet -> GetDataPoint(ievent) -> GetValue(5), fDataSet -> GetDataPoint(ievent) -> GetValue(5), true); 
			fModelTop -> SetDataBoundaries( 6, fDataSet -> GetDataPoint(ievent) -> GetValue(6), fDataSet -> GetDataPoint(ievent) -> GetValue(6), true); 
			fModelTop -> SetDataBoundaries( 7, fDataSet -> GetDataPoint(ievent) -> GetValue(7), fDataSet -> GetDataPoint(ievent) -> GetValue(7), true); 

			fModelTop -> SetDataBoundaries( 8, 0.0, 1000.0, false); 	
			fModelTop -> SetDataBoundaries( 9, fDataSet -> GetDataPoint(ievent) -> GetValue(9), fDataSet -> GetDataPoint(ievent) -> GetValue(9), true); 
			fModelTop -> SetDataBoundaries(10, fDataSet -> GetDataPoint(ievent) -> GetValue(10), fDataSet -> GetDataPoint(ievent) -> GetValue(10), true); 
			fModelTop -> SetDataBoundaries(11, fDataSet -> GetDataPoint(ievent) -> GetValue(11), fDataSet -> GetDataPoint(ievent) -> GetValue(11), true); 

			fModelTop -> SetDataBoundaries(12, 0.0, 1000.0, false); 
			fModelTop -> SetDataBoundaries(13, fDataSet -> GetDataPoint(ievent) -> GetValue(13), fDataSet -> GetDataPoint(ievent) -> GetValue(13), true); 
			fModelTop -> SetDataBoundaries(14, fDataSet -> GetDataPoint(ievent) -> GetValue(14), fDataSet -> GetDataPoint(ievent) -> GetValue(14), true); 
			fModelTop -> SetDataBoundaries(15, fDataSet -> GetDataPoint(ievent) -> GetValue(15), fDataSet -> GetDataPoint(ievent) -> GetValue(15), true); 

			fModelTop -> SetDataBoundaries(16, 0.0, 1000.0, false); 
			fModelTop -> SetDataBoundaries(17, fDataSet -> GetDataPoint(ievent) -> GetValue(17), fDataSet -> GetDataPoint(ievent) -> GetValue(17), true); 
			fModelTop -> SetDataBoundaries(18, fDataSet -> GetDataPoint(ievent) -> GetValue(18), fDataSet -> GetDataPoint(ievent) -> GetValue(8), true); 
			fModelTop -> SetDataBoundaries(19, fDataSet -> GetDataPoint(ievent) -> GetValue(19), fDataSet -> GetDataPoint(ievent) -> GetValue(19), true); 

			fModelTop -> SetDataBoundaries(20, 0.0, 1000.0, false); 
			fModelTop -> SetDataBoundaries(21, fDataSet -> GetDataPoint(ievent) -> GetValue(21), fDataSet -> GetDataPoint(ievent) -> GetValue(21), true); 
			fModelTop -> SetDataBoundaries(22, fDataSet -> GetDataPoint(ievent) -> GetValue(22), fDataSet -> GetDataPoint(ievent) -> GetValue(22), true); 
			fModelTop -> SetDataBoundaries(23, fDataSet -> GetDataPoint(ievent) -> GetValue(23), fDataSet -> GetDataPoint(ievent) -> GetValue(23), true); 

			// adjust settings 
			fModelTop -> MCMCSetNIterationsMax(100000); 
			fModelTop -> MCMCSetNIterationsRun(100000); 

			// adjust limits 
			double x    = fDataSet -> GetDataPoint(ievent) -> GetValue(0); 
			double xlow = BCMath::Max(0.0, x - sigmas * sqrt(x)); 
			double xhi  = x + sigmas * sqrt(x); 
			par_bhad_E -> SetLowerLimit(xlow); 
			par_bhad_E -> SetUpperLimit(xhi); 

			x    = fDataSet -> GetDataPoint(ievent) -> GetValue(4); 
			xlow = BCMath::Max(0.0, x - sigmas * sqrt(x)); 
			xhi  = x + sigmas * sqrt(x); 
			par_blep_E -> SetLowerLimit(xlow); 
			par_blep_E -> SetUpperLimit(xhi); 

			x    = fDataSet -> GetDataPoint(ievent) -> GetValue(8); 
			xlow = BCMath::Max(0.0, x - sigmas * sqrt(x)); 
			xhi  = x + sigmas * sqrt(x); 
			par_qup_E -> SetLowerLimit(xlow); 
			par_qup_E -> SetUpperLimit(xhi); 

			x    = fDataSet -> GetDataPoint(ievent) -> GetValue(12); 
			xlow = BCMath::Max(0.0, x - sigmas * sqrt(x)); 
			xhi  = x + sigmas * sqrt(x); 
			par_qdown_E -> SetLowerLimit(xlow); 
			par_qdown_E -> SetUpperLimit(xhi); 

			x    = fDataSet -> GetDataPoint(ievent) -> GetValue(16); 
			xlow = BCMath::Max(0.0, x - sigmas * sqrt(x)); 
			xhi  = x + sigmas * sqrt(x); 
			par_electron_E -> SetLowerLimit(xlow); 
			par_electron_E -> SetUpperLimit(xhi); 

			// add parameters 
			fModelTop -> AddParameter(par_bhad_E); 
			fModelTop -> AddParameter(par_blep_E); 
			fModelTop -> AddParameter(par_qup_E); 
			fModelTop -> AddParameter(par_qdown_E); 
			fModelTop -> AddParameter(par_electron_E); 
			fModelTop -> AddParameter(par_neutrino_pz); 
			fModelTop -> AddParameter(par_JES_all); 

			// reset max likelihood
			double maxlike = -1; 
			double maxloglike = -1;
			int maxlikeindex = -1; 
			double sumlike = 0.0; 

			std::vector <double> tempvector; 
			std::vector <double> bestfitvector; 

			// loop over permutations 
			// debug
			for (int iperm = 0; iperm < 12; ++iperm)
			//			for (int iperm = 0; iperm < 1; ++iperm)
				{
					// change to current permutation 
					fModelTop -> SetPermutation(iperm); 
					

					// perform minimization, etc. 
					//					fModelTop -> GetMinuit() -> SetPrintLevel(-1);
					// debug
					//					fModelTop -> FindMode();

					/*
					if (fModelTop -> GetMinuitErrorFlag() != 0)
						{
							nevents_minuiterror++; 
							std::cout << " Minuit error: " << fModelTop -> GetMinuitErrorFlag() << std::endl; 
						}
					*/

					fModelTop -> MarginalizeAll(); 
						

					fModelTop -> PrintAllMarginalized(Form("output/plots_%i_%i.ps", ievent, iperm)); 

					// calculate best fitting lorentz vectors 
					tempvector = fModelTop -> GetBestFitParameters(); 
					
					fModelTop -> CalculateLorentzVectors(tempvector); 

					// remove combinations with unrealistic masses 
					/*
					if (fabs(fModelTop -> fLV_Whad.M() - 80.4) > cut_whad)
						{ 
							//							std::cout << " hadronic W-mass not within range : " << fModelTop -> fLV_Whad.M() - 80.4 << " GeV/c2" << std::endl; 
							continue; 
						}

					if (fabs(fModelTop -> fLV_Wlep.M() - 80.4) > cut_wlep)
						{ 
							//							std::cout << " leptonic W-mass not within range : " << fModelTop -> fLV_Wlep.M() - 80.4 << " GeV/c2" << std::endl; 
							continue; 
						}

					if (fabs(fModelTop -> fLV_Tophad.M() - fModelTop -> fLV_Toplep.M()) > cut_deltat)
						{
							//							std::cout << " top mass difference too large : " << fModelTop -> fLV_Tophad.M() - fModelTop -> fLV_Toplep.M() << " GeV/c2" << std::endl; 
							continue; 
						}
					*/
					// calculate best likelihood 
					double like = fModelTop -> Likelihood(tempvector); 
					double loglike = fModelTop -> LogLikelihood(tempvector); 

					// sum over best likelihoods 
					sumlike += like; 

					// check if current permutation is more likely 
					if (loglike > maxloglike)
						{
							maxlike = like; 
							maxloglike = loglike; 
							maxlikeindex = iperm;
							bestfitvector = tempvector; 
						}
					
					// debug
					std::cout << " Calculate p-value" << std::endl; 
					std::cout << (fModelTop -> GetBestFitParameters()).size() << std::endl; 
					fModelTop -> CalculatePValue(fModelTop -> GetBestFitParameters(), true) -> Print(Form("output/pvalue_%i_%i.ps", ievent, iperm), 1, loglike); 
					//					fModelTop -> CalculatePValue(fModelTop -> GetBestFitParameters(), false);
					double pvalue = fModelTop -> GetPValue(); 
					std::cout << " pvalue = " << iperm << " " << pvalue << " (" << loglike << ")" << std::endl; 

					// print summary 
					fModelTop -> PrintResults(Form("output/summary_%i_%i.txt", ievent, iperm)); 

				} // end of loop over permutations 

			// increase counter of events with correct jet assignment 
			if (maxlikeindex == 0)
				nevents_correctcombination++; 

			// fill histograms for the best combintation, if there is one
			if (maxlikeindex >= 0)
				{
					// increase counter of accepted events 
					nevents_accepted++; 

					fModelTop -> SetPermutation(maxlikeindex); 
					fModelTop -> GetMinuit() -> SetPrintLevel(-1);
					//					fModelTop -> FindMode();
					fModelTop -> CalculateLorentzVectors(bestfitvector); 

					hist_combinations -> Fill(double(maxlikeindex)); 
					hist_best_comb -> Fill(double(maxlikeindex), maxloglike); 
					hist_sign_comb -> Fill(double(maxlikeindex), maxlike/sumlike); 
					hist_JES_all -> Fill(bestfitvector.at(6));
					hist_mass_Whad -> Fill(fModelTop -> fLV_Whad.M()); 
					hist_mass_Wlep -> Fill(fModelTop -> fLV_Wlep.M()); 
					hist_mass_Tophad -> Fill(fModelTop -> fLV_Tophad.M()); 
					hist_mass_Toplep -> Fill(fModelTop -> fLV_Toplep.M()); 

					hist_JES_all_mass_Tophad -> Fill(bestfitvector.at(6),
																					 fModelTop -> fLV_Tophad.M());
				}

			// delete model 
			delete fModelTop; 

		} // end of loop over events 

	// ---------------------------------------------------------
	// print to screen  
	// ---------------------------------------------------------

	std::cout << " ===================================== " << std::endl; 
	std::cout << " N (processed)  : " << nevents_processed << std::endl; 
	std::cout << " N (accepted)   : " << nevents_accepted << std::endl; 
	std::cout << " N (correct)    : " << nevents_correctcombination << std::endl; 
	std::cout << std::endl; 
	std::cout << " Eff (accepted) : " << double(nevents_accepted) / double(nevents_processed) * 100.0 << "%" << std::endl; 
	std::cout << " Eff (correct)  : " << double(nevents_correctcombination) / double(nevents_accepted)  * 100.0 << "%" << std::endl; 

	std::cout << std::endl;
	std::cout << " N (Minuit errors) : " << nevents_minuiterror << std::endl; 

	std::cout << " ===================================== " << std::endl; 



	// ---------------------------------------------------------
	// print histograms 
	// ---------------------------------------------------------

	gStyle -> SetPalette(1); 

	TCanvas * c1 = new TCanvas("c1"); 
	c1 -> Divide(2, 2); 
	c1 -> cd(1); 
	hist_mass_Whad -> Draw(); 
	c1 -> cd(2); 
	hist_mass_Wlep -> Draw(); 
	c1 -> cd(3); 
	hist_mass_Tophad -> Draw(); 
	c1 -> cd(4); 
	hist_mass_Toplep -> Draw(); 

	c1 -> Print("c1.eps"); 

	TCanvas * c2 = new TCanvas("c2"); 
	c2 -> cd(); 
	hist_combinations -> Draw(); 

	c2 -> Print("c2.eps"); 

	TCanvas * c3 = new TCanvas("c3"); 
	c3 -> cd(); 
	hist_best_comb -> Draw("COLZ"); 

	c3 -> Print("c3.eps"); 

	TCanvas * c4 = new TCanvas("c4"); 
	c4 -> cd(); 
	hist_sign_comb -> Draw("COLZ"); 

	c4 -> Print("c4.eps"); 

	TCanvas * c5 = new TCanvas("c5"); 
	c5 -> cd(); 
	hist_JES_all -> Draw(); 
	
	c5 -> Print("c5.eps"); 

	TCanvas * c6 = new TCanvas("c6"); 
	c6 -> cd(); 
	hist_JES_all_mass_Tophad -> Draw("COLZ"); 
	
	c6 -> Print("c6.eps"); 

	// ---------------------------------------------------------
	// close log file 
	// ---------------------------------------------------------

	// closes the log file 
	BCLog::CloseLog(); 

	return 0; 

}

// ---------------------------------------------------------
  
