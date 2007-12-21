#include "BCEngineMCMC.h" 

#include <TH1D.h> 
#include <TH2D.h> 
#include <TCanvas.h> 

int main()
{

	int nchains = 10; 

	BCEngineMCMC * enginemcmc = new BCEngineMCMC(nchains); 

	enginemcmc -> MCMCAddParameter( -20.0, 200.0); 
	enginemcmc -> MCMCAddParameter( -20.0, 200.0); 
	enginemcmc -> MCMCAddParameter( -20.0, 200.0); 
	enginemcmc -> MCMCAddParameter( -20.0, 200.0); 
	enginemcmc -> MCMCAddParameter( -20.0, 200.0); 

	enginemcmc -> MCMCSetFlagInitialPosition(1); 
	//	enginemcmc -> MCMCSetTrialFunctionScale(1.0); 
	enginemcmc -> MCMCSetTrialFunctionScale(0.05); 

	enginemcmc -> MCMCInitialize(); 

	enginemcmc -> MCMCMetropolis(); 
	
	//	cout << enginemcmc -> MCMCGetNIterationsConvergenceGlobal() << endl; 
	//	return 0; 

	std::vector <TH1D*> histcontainer_x1; 
	std::vector <TH1D*> histcontainer_y1; 
	std::vector <TH1D*> histcontainer_z; 
	std::vector <TH1D*> histcontainer_x2; 
	std::vector <TH1D*> histcontainer_y2; 


	for (int i = 0; i < nchains; ++i)
		{
			TH1D * hist_x1 = new TH1D(Form("hist_x1_%i", i), "", 120, -20.0, 100.0); 
			hist_x1 -> SetLineColor(2 + i); 

			histcontainer_x1.push_back(hist_x1); 

			TH1D * hist_y1 = new TH1D(Form("hist_y1_%i", i), "", 120,   0.0, 200.0); 
			hist_y1 -> SetLineColor(2 + i); 

			histcontainer_y1.push_back(hist_y1); 

			TH1D * hist_z = new TH1D(Form("hist_z_%i", i), "", 120, -20.0, 100.0); 
			hist_z -> SetLineColor(2 + i); 

			histcontainer_z.push_back(hist_z); 

			TH1D * hist_x2 = new TH1D(Form("hist_x2_%i", i), "", 120, -20.0, 100.0); 
			hist_x2 -> SetLineColor(2 + i); 

			histcontainer_x2.push_back(hist_x2); 

			TH1D * hist_y2 = new TH1D(Form("hist_y2_%i", i), "", 120,   0.0, 200.0); 
			hist_y2 -> SetLineColor(2 + i); 

			histcontainer_y2.push_back(hist_y2); 
		}

	int dn = int((enginemcmc -> MCMCGetNIterationsMax()) / 100); 
	
	TH1D * hist_r = new TH1D("hist_r", "", 100, 0.0, double(enginemcmc -> MCMCGetNIterationsMax())); 

	// perform burn-in run 

	for (int i = 0; i < enginemcmc -> MCMCGetNIterationsBurnIn(); ++i)
		for (int j = 0; j < enginemcmc -> MCMCGetNChains(); ++j)
			enginemcmc -> MCMCGetNewPointMetropolis(j, false);

	// reset run statistics 

	enginemcmc -> MCMCResetRunStatistics(); 
	
	// perform PCA run 

	//	enginemcmc -> MCMCPCARun(); 
	
	// perform run 

	for (int i = 0; i < enginemcmc -> MCMCGetNIterationsMax(); ++i)
		{
			for (int j = 0; j < enginemcmc -> MCMCGetNChains(); ++j)
				enginemcmc -> MCMCGetNewPointMetropolis(j, false);

			enginemcmc -> MCMCUpdateConvergence(); 

			if ((i % dn) == 0)
				hist_r -> SetBinContent(i / dn + 1, log(enginemcmc -> MCMCGetRValue())); 

			for (int j = 0; j < nchains; ++j)
				{
					histcontainer_x1[j] -> Fill((enginemcmc -> MCMCGetx()).at(5 * j + 0)); 
					histcontainer_y1[j] -> Fill((enginemcmc -> MCMCGetx()).at(5 * j + 1)); 
					histcontainer_z[j]  -> Fill((enginemcmc -> MCMCGetx()).at(5 * j + 2)); 
					histcontainer_x2[j] -> Fill((enginemcmc -> MCMCGetx()).at(5 * j + 3)); 
					histcontainer_y2[j] -> Fill((enginemcmc -> MCMCGetx()).at(5 * j + 4)); 
				}
		}

	cout << enginemcmc -> MCMCGetNIterationsConvergenceGlobal() << endl; 
	
	for (int i = 0; i < nchains; ++i)
		cout << " efficiency : " << double((enginemcmc -> MCMCGetNTrialsTrue()).at(i)) / double((enginemcmc -> MCMCGetNIterations()).at(i)) << endl; 

	TCanvas * canvas = new TCanvas(); 
	canvas -> Divide(3,2); 

	canvas -> cd(1); 
	histcontainer_x1[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_x1[i] -> Draw("SAME"); 

	canvas -> cd(2); 
	histcontainer_y1[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_y1[i] -> Draw("SAME"); 

	canvas -> cd(3); 
	histcontainer_z[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_z[i] -> Draw("SAME"); 

	canvas -> cd(4); 
	histcontainer_x2[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_x2[i] -> Draw("SAME"); 

	canvas -> cd(5); 
	histcontainer_y2[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_y2[i] -> Draw("SAME"); 

	canvas -> cd(6); 
	hist_r -> Draw(""); 

	canvas -> Print("canvas.ps"); 

	delete enginemcmc; 

	return 0; 

}
