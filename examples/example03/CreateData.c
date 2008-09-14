#include <fstream.h> 

int CreateSpectrum(); 

void CreateData()
{

  int nspectra = 1; 
  int sum = 0; 

  for (int i = 0; i < nspectra; i++)
    sum += CreateSpectrum(); 

  cout << nspectra << " " << sum/double(nspectra) << " +- " << TMath::Sqrt(sum/double(nspectra)) / TMath::Sqrt(double(nspectra)) << endl; 

}

int CreateSpectrum()
{
  
  // parameters 

  double background =   10.0; 
  double signal     =   10.0; 
  double mean       = 2039.0; 
  double sigma      =    2.0; 

  int nbins = 21; 

  double Emin = mean - 10.5; 
  double Emax = mean + 10.5; 
  double dE   = (Emax - Emin) / double(nbins); 

  // initialize random number generator 
  TRandom3* fRandom = new TRandom3(1000); 

  // loop over models  
  char filename[200]; 

  sprintf(filename, "./data/data.txt"); 
  
  // open filelist file 
  std::fstream data_stream; 

  data_stream.open(filename, std::fstream::out);

  // write energy interval to file 
  data_stream << Emin << " " << Emax << endl; 
  
  // reset counter 
  int nevents = 0; 

  // loop over bins
  for (int i = 0; i < nbins; i++)
    {
      double E = Emin + i * dE; 

      double t1 = (E - mean) / (TMath::Sqrt(2.0) * sigma); 
      double t2 = (E + dE - mean) / (TMath::Sqrt(2.0) * sigma); 

      double expectation = (background / double(nbins) + 
			    signal * 0.5 * (TMath::Erf(t2) - TMath::Erf(t1))); 

      int n = fRandom -> Poisson(expectation); 

      // increase counter 
      nevents += n; 

      data_stream << E << " " << double(n); 

      if (i < nbins - 1) 
	data_stream << endl; 
    }

  // close file 
  
  data_stream.close(); 

  return nevents; 

}
