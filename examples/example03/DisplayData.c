#include <fstream.h> 
  
void DisplayData()
{

  int nfiles = 100; 

  // histogram 

  TH1D* hist = new TH1D("hist", "", 101, -0.5, 100.5); 
  hist -> SetXTitle("x"); 
  hist -> SetYTitle("Entries"); 
  
  // loop over files 

  for (int i = 0; i < nfiles; i++)
    {
      char filename[200]; 
	
      //      sprintf(filename, "data/data_Background-only-model_%i.txt", i); 
      sprintf(filename, "data/data_Simple-signal+background-model_%i.txt", i); 

      std::fstream stream; 
      
      stream.open(filename, std::fstream::in);
      
      // check if file is open 

      if (!stream.is_open())
	{
	  cout << " Couldn't open file " << filename << endl; 
	  return; 
	}

      // read from file 
      
      while (!stream.eof())
	{
	  double data; 
	  
	  stream >> data; 
	  
	  hist -> Fill(data); 
	}
      
      // close file 
      
      stream.close(); 
    }

  // print histogram 

  TCanvas* canvas = new TCanvas(); 
  canvas -> cd(); 

  hist -> Draw(); 

}
