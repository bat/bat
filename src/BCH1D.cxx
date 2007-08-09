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
  
  TCanvas * canvas = new TCanvas();
  canvas -> cd();
  
  // draw histogram

  fHistogram -> Draw();
  
  // options: draw lines

  double maximum = fHistogram -> GetBinContent(fHistogram -> FindBin(this -> GetMode()));
  
  TLine * line = new TLine();
  
  switch(options)
    {
      // draw lines
      
    case 1:

      // draw line for mode

      line -> DrawLine(this -> GetMode(), 0.0, this -> GetMode(), 1.05 * maximum);

      // draw line for 16%

      line -> DrawLine(this -> GetQuantile(0.16), 0.0, this -> GetQuantile(0.16), 1.05 * maximum);

      // draw line for 84%

      line -> DrawLine(this -> GetQuantile(0.84), 0.0, this -> GetQuantile(0.84), 1.05 * maximum);
      
      break;
      
      // draw line at value
      
    case 2:

      line -> SetLineColor(kRed);
      line -> DrawLine(value, 0.0, value, 1.05 * maximum);
      
      break;

    case 3: 
      {
	// draw yellow shaded region for the 16% and 84% quantiles with a red triangle on top of the mode 
	
	TH1D * hist_shaded =  (TH1D *) fHistogram -> Clone(); 
	hist_shaded -> SetFillStyle(1001); 
	hist_shaded -> SetFillColor(kYellow); 
	
	// remove entries on the sides 
	
	int nbins = hist_shaded -> GetNbinsX(); 
	int lowerborder = hist_shaded -> FindBin(this -> GetQuantile(0.16)); 
	int upperborder = TMath::Min(hist_shaded -> FindBin(this -> GetQuantile(0.84)) + 1, nbins); 
	
	for (int i = 1; i < lowerborder; i++)
	  hist_shaded -> SetBinContent(i, 0.0); 
	
	for (int i = upperborder; i <= nbins; i++)
	  hist_shaded -> SetBinContent(i, 0.0); 
	
	// draw shaded histogram 
	
	hist_shaded -> Draw("SAME"); 
	gPad->RedrawAxis(); 
	
	// draw triangle 
	
	TMarker * marker = new TMarker(); 
	marker -> SetMarkerStyle(23); 
	marker -> SetMarkerSize(3); 
	marker -> SetMarkerColor(kRed); 
	
	marker -> DrawMarker(this -> GetMode(), 
			     1.0 * maximum); 
	
	// draw legend 
	
	TLegend * legend = new TLegend(0.20, 0.75, 0.55, 0.95); 
	legend -> SetFillColor(kWhite); 
	legend -> SetFillStyle(0); 
	legend -> SetBorderSize(0); 
	legend -> AddEntry(fHistogram,  "probability density", "L"); 
	legend -> AddEntry(hist_shaded, "68% prob. region", "F"); 
	legend -> AddEntry(marker,      "mode", "P"); 
	
	legend -> Draw(); 
      }

      break; 

      // bad options
      
    default:

      BCLog::Out(BCLog::warning, BCLog::warning, Form("BCH1D::Print. Invalid option %d",options));

      break;
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
 





