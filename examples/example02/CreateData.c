#include <fstream.h> 

void CreateData()
{
  
  // parameters 

  int npoints = 10; 
  float par0; 
  float par1; 
  float par2; 

  float sigmay = 0.2; 

  float xmin =   0.0; 
  float xmax = 100.0; 

  bool random = false; 

  int ndatasets = 100 + 1; 

  // initialize random number generator 

  TRandom* fRandom = new TRandom(0); 

  // loop over polynomials 

  for (int ipol = 0; ipol < 3; ipol++)
    {
      if (ipol == 0) 
	{
	  par0 = 1.0; 
	  par1 = 0.0; 
	  par2 = 0.0; 
	}

      else if (ipol == 1)
	{
	  par0 = 2.0; 
	  par1 = 0.01; 
	  par2 = 0.0; 
	}

      else
	{
	  par0 = 1.0; 
	  par1 = 0.0; 
	  par2 = 0.0003; 
	}

      // open filelist file 

      std::fstream filelist_stream; 

      filelist_stream.open(Form("./data/filelist_pol%i.txt", ipol), std::fstream::out);

      // histogram for conditional probability 

      TCanvas* canvas = new TCanvas(); 

      TH1F* hist_pcond = new TH1F(Form("hist_pcond_pol%i", ipol), "", 100, -30.0, -15.0); 

      // loop over data sets 
      
      for (int idataset = 0; idataset < ndatasets; idataset++)
	{
	  
	  // open file 
	  
	  std::fstream file_stream;

	  char filename[200]; 
	  
	  if (idataset == 0) 
	    sprintf(filename, "./data/data_ModelPol%i.txt", ipol); 
	  else
	    sprintf(filename, "./data/data_ModelPol%i_%i.txt", ipol, idataset); 
	  
	  file_stream.open(filename, std::fstream::out);

	  // write file to filelist 

	  if (idataset > 0)
	    {
	      filelist_stream << filename; 
	      
	      if (idataset < ndatasets - 1) 
		filelist_stream << std::endl; 
	    }

	  double pcond = 1.0; 

	  // loop over points 
	  
	  for (int i = 0; i < npoints; i++)
	    {
	      // get x value 
	      
	      float x; 
	      
	      if (random == true)
		x = fRandom -> Uniform(xmin, xmax); 
	      
	      else
		x = (xmax - xmin) / float(npoints) * (float(i) + 0.5); 
	      
	      // get y value 
	      
	      float y; 
	      
	      float mean = par0 + par1 * x + par2 * x * x; 
	      
	      y = fRandom -> Gaus(mean * 100, sigmay * 100) / 100.0; 
	      
	      // conditional probability 

	      pcond *= TMath::Gaus(y * 100.0, (par0 + x * par1 + x * x * par2) * 100.0, sigmay * 100.0, true);

	      // write to file 
	      
	      file_stream << x << " " << y << " " << sigmay; 
	      
	      if (i < npoints - 1)
		file_stream << std::endl; 
	    }

	  // fill histogram 

	  hist_pcond -> Fill(TMath::Log10(pcond)); 

	  if (idataset == 0)
	    cout << " conditional probability = " << pcond << endl; 

	  // close file 
	  
	  file_stream.close(); 

	}
      
      // print histogram 
      
      canvas -> cd(); 
      hist_pcond -> Draw(); 

      // close list file 

      filelist_stream.close(); 

    }
  
}
