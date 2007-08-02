#include "BCH1D.h" 

#include <TCanvas.h> 
#include <TLine.h> 
#include <TMarker.h> 
#include <TLegend.h> 

// --------------------------------------------------------- 

BCH1D::BCH1D() 
{
	fHistogram = 0; 

	// suppress the ROOT Info printouts
	gErrorIgnoreLevel=1;
}

// --------------------------------------------------------- 

BCH1D::~BCH1D()
{
	if (fHistogram) 
		delete fHistogram; 
}

// --------------------------------------------------------- 

double BCH1D::GetMode()
{
	return fHistogram -> GetBinCenter(fHistogram -> GetMaximumBin()); 
}

// --------------------------------------------------------- 

double BCH1D::GetQuantile(double probabilitysum)
{

  int nquantiles = 1; 
  double quantiles[1]; 
  double probsum[1]; 

  probsum[0] = probabilitysum; 

  fHistogram -> GetQuantiles(nquantiles, quantiles, probsum); 

  return quantiles[0]; 

}

// --------------------------------------------------------- 

double BCH1D::GetIntegral(double valuemin, double valuemax) 
{

  double integral = 0; 

  int binmin = fHistogram -> FindBin(valuemin); 
  int binmax = fHistogram -> FindBin(valuemax); 

  integral = fHistogram -> Integral(binmin, binmax); 

  return integral; 

}

// --------------------------------------------------------- 

double BCH1D::GetPValue(double probability)
{
  
  double valuemax = fHistogram -> GetBinCenter(fHistogram -> GetNbinsX()); 
  
  double integral = this -> GetIntegral(probability, valuemax); 

  return integral; 

}

// --------------------------------------------------------- 
 
void BCH1D::Print(char* filename, int options, double value) 
{

  // create temporary canvas 

  TCanvas* canvas = new TCanvas(); 

  canvas -> cd(); 

  // draw histogram 
  
  fHistogram -> Draw(); 

  // options: draw lines 

  double maximum = fHistogram -> GetBinContent(fHistogram -> FindBin(this -> GetMode())); 

  if (options == 1) 
    {
      // draw line for mode 
      
      TLine* line = new TLine(); 
      
      line -> DrawLine(this -> GetMode(), 0.0, 
		       this -> GetMode(), 1.05 * maximum); 
      
      // draw line for 16% 
      
      line -> DrawLine(this -> GetQuantile(0.16), 0.0, 
		       this -> GetQuantile(0.16), 1.05 * maximum); 
      
      // draw line for 84% 
      
      line -> DrawLine(this -> GetQuantile(0.84), 0.0, 
		       this -> GetQuantile(0.84), 1.05 * maximum); 
    } 

  else if (options == 2) 
    {
      // draw line at value 
      
      TLine* line = new TLine();  
      line -> SetLineColor(kRed); 
      
      line -> DrawLine(value, 0.0, 
		       value, 1.05 * maximum); 
    }

  // print to file 

  canvas -> Print(filename); 
  
}

// --------------------------------------------------------- 
 
void BCH1D::Print(char* filename, int options) 
{
  
  double value = 0; 

  this -> Print(filename, options, value); 

}

// --------------------------------------------------------- 
 
void BCH1D::Print(char* filename) 
{
  
  int options = 0; 
  double value = 0; 

  this -> Print(filename, options, value); 

}

// --------------------------------------------------------- 
 





