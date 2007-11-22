#include <fstream.h> 

void CreateData()
{
  
  // parameters 

	int nevents = 1000; 

	float M     =   50.0; 
  float Gamma =   5.0; 
	float sigma =   5.0; 
	
	float E_min =   0.0; 
	float E_max = 100.0; 

  // initialize random number generator 

  TRandom3* fRandom = new TRandom3(0); 

  // open file 
  
  std::fstream file_data;
  
  char filename[200]; 
  
  file_data.open("./data/data.txt", std::fstream::out);

	// histograms 

	TH1D * hist_M = new TH1D("hist_M", ";M;N", 101, -0.5, 100.5); 
	hist_M -> SetLineColor(kBlue); 

	TH1D * hist_E = new TH1D("hist_E", ";E_{1} + E_{2};N", 101, -0.5, 100.5); 
	hist_E -> SetLineColor(kRed); 

  // loop over points 
	  
  for (int i = 0; i < nevents; i++)
    {

			// add resolution 

			double E = E_min-1.0; 
			double Mevent = gRandom -> Gaus(M, Gamma); 

			while (!(E>E_min&&E<E_max))
				{
					// get center-of-mass energy 

					Mevent = gRandom -> Gaus(M, Gamma); 

					// add resolution 

					E = gRandom -> Gaus(Mevent, sigma); 
	      }
					// write to file 
      
			file_data << E << " " << endl; 

			hist_M -> Fill(Mevent); 
			hist_E -> Fill(E); 
		}

	// plot histograms 

	TCanvas * canvas = new TCanvas(); 
	
	canvas -> cd(); 

	hist_M -> Draw(); 
	hist_E -> Draw("SAME"); 

  // close file 
	  
  file_data.close(); 
	  
}
