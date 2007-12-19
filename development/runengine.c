#include "BCEngineMCMC.h" 

#include <TH1D.h> 
#include <TH2D.h> 
#include <TCanvas.h> 

int main()
{

	int nchains = 2; 

	BCEngineMCMC * enginemcmc = new BCEngineMCMC(nchains); 

	enginemcmc -> MCMCAddParameter( -20.0, 100.0); 
	enginemcmc -> MCMCAddParameter(   0.0, 200.0); 

	enginemcmc -> MCMCSetFlagInitialPosition(1); 
	enginemcmc -> MCMCSetTrialFunctionScale(0.01); 

	enginemcmc -> MCMCInitialize(); 

	std::vector <TH1D*> histcontainer_x; 
	std::vector <TH1D*> histcontainer_y; 
	std::vector <TH2D*> histcontainer_xy; 

		enginemcmc -> MCMCMetropolis(); 
		return 0; 

	for (int i = 0; i < nchains; ++i)
		{
			TH1D * hist_x = new TH1D(Form("hist_x%i", i), "", 120, -20.0, 100.0); 
			hist_x -> SetLineColor(2 + i); 

			histcontainer_x.push_back(hist_x); 

			TH1D * hist_y = new TH1D(Form("hist_y%i", i), "", 120,   0.0, 200.0); 
			hist_y -> SetLineColor(2 + i); 

			histcontainer_y.push_back(hist_y); 

			TH2D * hist_xy = new TH2D(Form("hist_xy%d", i), "", 120, -20.0, 100.0, 200, 0.0, 200.0); 
			hist_xy -> SetStats(kFALSE); 
			hist_xy -> SetMarkerColor(2 + i); 

			histcontainer_xy.push_back(hist_xy); 
		}

	int dn = int((enginemcmc -> MCMCGetNIterationsMax()) / 100); 
	
	TH1D * hist_r = new TH1D("hist_r", "", 100, 0.0, double(enginemcmc -> MCMCGetNIterationsMax())); 

	// perform burn-in run 

	for (int i = 0; i < enginemcmc -> MCMCGetNIterationsBurnIn(); ++i)
		for (int j = 0; j < enginemcmc -> MCMCGetNChains(); ++j)
			enginemcmc -> MCMCGetNewPointMetropolis(j, false);

	// reset run statistics 

	enginemcmc -> MCMCResetRunStatistics(); 

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
					histcontainer_x[j] -> Fill((enginemcmc -> MCMCGetx()).at(2 * j + 0)); 
					histcontainer_y[j] -> Fill((enginemcmc -> MCMCGetx()).at(2 * j + 1)); 
					histcontainer_xy[j] -> Fill((enginemcmc -> MCMCGetx()).at(2 * j + 0), (enginemcmc -> MCMCGetx()).at(2 * j + 1)); 
				}
		}

	TCanvas * canvas = new TCanvas(); 
	canvas -> Divide(2,2); 
	canvas -> cd(1); 
	histcontainer_x[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_x[i] -> Draw("SAME"); 
	canvas -> cd(2); 
	histcontainer_y[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_y[i] -> Draw("SAME"); 
	canvas -> cd(3); 
	histcontainer_xy[0] -> Draw(); 
	for (int i = 1; i < nchains; ++i) 
		histcontainer_xy[i] -> Draw("SAME"); 
	canvas -> cd(4); 
	hist_r -> Draw(""); 
	canvas -> Print("canvas.ps"); 

	delete enginemcmc; 

	return 0; 

}
