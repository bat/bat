#include <fstream.h> 

void CreateData()
{
  
  // parameters 

  int npoints = 10; 

  float par0 = 2.0; 
  float par1 = 0.01; 
  
  float sigmay = 0.2; 
  
  float xmin =   0.0; 
  float xmax = 100.0; 

  // initialize random number generator 

  TRandom3* fRandom = new TRandom3(0); 

  // open file 
  
  std::fstream file_data;
  
  char filename[200]; 
  
  file_data.open("./data/data.txt", std::fstream::out);

  // loop over points 
	  
  for (int i = 0; i < npoints; i++)
    {
      // get x value 
      
      float x; 
      
      x = (xmax - xmin) / float(npoints) * (float(i) + 0.5); 
      
      // get y value 
      
      float y; 
      
      float mean = par0 + par1 * x; 
	      
      y = fRandom -> Gaus(mean, sigmay); 
	      
      // write to file 
      
      file_data << x << " " << y << " " << sigmay; 
	      
      if (i < npoints - 1)
	file_data << std::endl; 
    }

  // close file 
	  
  file_data.close(); 
	  
}
