#include "BCH2D.h" 
#include "BCMath.h" 

#include <TCanvas.h> 
#include <TLine.h> 
#include <TMarker.h> 
#include <TLegend.h> 

// --------------------------------------------------------- 

BCH2D::BCH2D() 
{
	fHistogram = 0; 
	fIntegratedHistogram = 0; 
}

// --------------------------------------------------------- 

BCH2D::~BCH2D()
{
	if (fHistogram) 
		delete fHistogram; 

	if (fIntegratedHistogram) 
		delete fIntegratedHistogram; 
}

// --------------------------------------------------------- 

void BCH2D::GetMode(double& mode)
{
//	double mode[2]; 
//	int binx, biny, binz; 

//	fHistogram -> GetMaximumBin(binx, biny, binz); 
//	mode = fHistogram -> GetBinCenter(); 

//	return mode; 
}

// --------------------------------------------------------- 

void BCH2D::Print(char* filename, int options) 
{

  // create temporary canvas 
  TCanvas* canvas = new TCanvas(); 

  canvas -> cd(); 

  // draw histogram 
  fHistogram -> SetLineColor(kBlack); 
  fHistogram -> SetLineWidth(4); 

  if (options == 0) 
    fHistogram -> Draw("CONT0"); 

  if (options == 1)
    {
      fHistogram -> Draw("CONT3"); 

      // set contours 
      this -> CalculateIntegratedHistogram(); 
      
      double levels[4]; 
            
      levels[0] = 0.0; 
      levels[1] = this -> GetLevel(1.0 - 0.6827); 
      levels[2] = this -> GetLevel(1.0 - 0.9545); 
      levels[3] = this -> GetLevel(1.0 - 0.9973); 
  
      fHistogram -> SetContour(4, levels); 

      // best fit value 
      int maximumbin = fHistogram -> GetMaximumBin(); 

      int binx = maximumbin % (fHistogram -> GetNbinsX() + 2); 
      int biny = maximumbin / (fHistogram -> GetNbinsX() + 2); 

      double x = fHistogram -> GetXaxis() -> GetBinCenter(binx); 
      double y = fHistogram -> GetYaxis() -> GetBinCenter(biny); 

      TMarker* marker = new TMarker(x, y, 24); 
      marker -> Draw(); 

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95); 
      legend -> SetBorderSize(0); 
      legend -> SetFillColor(kWhite); 
      legend -> AddEntry(fHistogram, "68% prob. region", "L"); 
      legend -> AddEntry(marker, "Best fit", "P"); 
      legend -> Draw(); 
    }

  if (options == 2)
    {
      fHistogram -> Draw("CONT3"); 

      // set contours 
      this -> CalculateIntegratedHistogram(); 
      
      double levels[2]; 
      double level32 = this -> GetLevel(0.32); 
      
      levels[0] = 0.0; 
      levels[1] = level32; 
  
      fHistogram -> SetContour(2, levels); 

      // best fit value 
      int maximumbin = fHistogram -> GetMaximumBin(); 

      int binx = maximumbin % (fHistogram -> GetNbinsX() + 2); 
      int biny = maximumbin / (fHistogram -> GetNbinsX() + 2); 

      double x = fHistogram -> GetXaxis() -> GetBinCenter(binx); 
      double y = fHistogram -> GetYaxis() -> GetBinCenter(biny); 

      TMarker* marker = new TMarker(x, y, 24); 
      marker -> Draw(); 

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95); 
      legend -> SetBorderSize(0); 
      legend -> SetFillColor(kWhite); 
      legend -> AddEntry(fHistogram, "68% prob. region", "L"); 
      legend -> AddEntry(marker, "Best fit", "P"); 
      legend -> Draw(); 
    }

  if (options == 3) 
    {
      fHistogram -> Draw("CONT3"); 

      // set contours 
      this -> CalculateIntegratedHistogram(); 
      
      double levels[2]; 
      double level10 = this -> GetLevel(0.10); 
      
      levels[0] = 0.0; 
      levels[1] = level10; 
  
      fHistogram -> SetContour(2, levels); 

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95); 
      legend -> SetBorderSize(0); 
      legend -> SetFillColor(kWhite); 
      legend -> AddEntry(fHistogram, "90% prob. region", "L"); 
      legend -> Draw(); 
    }

  if (options == 4) 
    {
      fHistogram -> Draw("CONT3"); 

      // set contours 
      this -> CalculateIntegratedHistogram(); 
      
      double levels[2]; 
      double level5 = this -> GetLevel(0.05); 
      
      levels[0] = 0.0; 
      levels[1] = level5; 
  
      fHistogram -> SetContour(2, levels); 

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95); 
      legend -> SetBorderSize(0); 
      legend -> SetFillColor(kWhite); 
      legend -> AddEntry(fHistogram, "95% prob. region", "L"); 
      legend -> Draw(); 
    }


  // print to file 
  canvas -> Print(filename); 

}

// --------------------------------------------------------- 

void BCH2D::CalculateIntegratedHistogram()
{
	int nz = 100;

	double zmax = fHistogram -> GetMaximum();
	double dz   = zmax / double(nz);

	double nx = fHistogram -> GetNbinsX();
	double ny = fHistogram -> GetNbinsY();

	// create histogram
	if (fIntegratedHistogram)
		delete fIntegratedHistogram;

	fIntegratedHistogram = new TH1D(Form("%s_int_prob",fHistogram->GetName()), "", nz, 0.0, zmax);
	fIntegratedHistogram -> SetXTitle("z");
	fIntegratedHistogram -> SetYTitle("Integrated probability");
	fIntegratedHistogram -> SetStats(kFALSE);

	// loop over histogram
	for (int ix = 1; ix <= nx; ix++)
		for (int iy = 1; iy <= ny; iy++)
		{
			int binmin = BCMath::Nint(fHistogram -> GetBinContent(ix, iy) / dz);
			for (int i = binmin; i <= nz; i++)
				fIntegratedHistogram -> SetBinContent(i,
					fIntegratedHistogram -> GetBinContent(i) +
					fHistogram -> GetBinContent(ix, iy));
		}
}

// --------------------------------------------------------- 

double BCH2D::GetLevel(double p) 
{
	// calculate level 
	double level = 0.0; 

	int nquantiles = 1; 
	double quantiles[1]; 
	double probsum[1]; 

	probsum[0] = p; 

	fIntegratedHistogram -> GetQuantiles(nquantiles, quantiles, probsum); 

	level = quantiles[0]; 

	return level; 
}

// --------------------------------------------------------- 
