#include <fstream.h> 

void CreateData()
{
  
  // parameters 

	int nevents = 1000; 

	float M     = 100.0; 
  float Gamma =   2.0; 
	float sigma =   0.5; 
	
	float E_min =   0.0; 
	float E_max = 200.0; 

  // initialize random number generator 

  TRandom3* fRandom = new TRandom3(0); 

  // open file 
  
  std::fstream file_data;
  
  char filename[200]; 
  
  file_data.open("./data/data.txt", std::fstream::out);

	// histograms 

	TH1D * hist_M = new TH1D("hist_M", ";M;N", 101, 49.5, 150.5); 
	hist_M -> SetLineColor(kBlue); 

	TH1D * hist_Esumtrue = new TH1D("hist_Esumtrue", ";E_{1} + E_{2};N", 101, 49.5, 150.5); 
	hist_Esumtrue -> SetLineColor(kRed); 

	TH1D * hist_Esum = new TH1D("hist_Esum", ";E_{1} + E_{2};N", 101, 49.5, 150.5); 

  // loop over points 
	  
  for (int i = 0; i < nevents; i++)
    {

			// add resolution 

			double E1 = E_min-1.0; 
			double E2 = E_min-1.0; 

			while (!(E1>E_min&&E1<E_max&&E2>E_min&&E2<E_max))
				{
					// get center-of-mass energy 

					double Ecom = gRandom -> BreitWigner(M, Gamma); 
					double E     = Ecom / 2.0; 

					// add resolution 

					E1 = gRandom -> Gaus(E, sigma * sqrt(E)); 
					E2 = gRandom -> Gaus(E, sigma * sqrt(E)); 
	      }
					// write to file 
      
					file_data << E1 << " " << E2 << endl; 

					hist_M -> Fill(M); 
					hist_Esumtrue -> Fill(Ecom); 
					hist_Esum -> Fill(E1 + E2); 
		}

	// plot histograms 

	TCanvas * canvas = new TCanvas(); 
	
	canvas -> cd(); 

	hist_M -> Draw(); 
	hist_Esumtrue -> Draw("SAME"); 
	hist_Esum -> Draw("SAME"); 

  // close file 
	  
  file_data.close(); 
	  
}
