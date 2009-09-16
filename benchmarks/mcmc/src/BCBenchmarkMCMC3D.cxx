/*
 * Copyright (C) 2009, 
 * Daniel Kollar, Kevin Kroeninger and Jing Liu.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

//=============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "BAT/BCLog.h"

#include "BCBenchmarkMCMC3D.h"

//=============================================================================

BCBenchmarkMCMC3D::BCBenchmarkMCMC3D(
		TF3* testFunction,
		const char* outputFile,
		const char* modelName)
:BCModel(modelName), BCModelOutput(), fNbinx(100), fNbiny(100), fNbinz(100)
{
	BCLog::Out(BCLog::summary,BCLog::summary," setup test function ...");
	if (testFunction==NULL) {
		BCLog::Out(BCLog::error,BCLog::error,
				" The test function doesn't exist!");
		abort();
	}

	fTestFunction = testFunction;
	fXmin=fTestFunction->GetXmin(); 
	fXmax=fTestFunction->GetXmax();
	fYmin=fTestFunction->GetYmin(); 
	fYmax=fTestFunction->GetYmax();
	fZmin=fTestFunction->GetZmin();
	fZmax=fTestFunction->GetZmax();
	
	// set no. of points for histogram representation of test function
	// (used for Kolmogorov-Smirnov test)
	fTestFunction->SetNpx(fNbinx);
	fTestFunction->SetNpy(fNbiny);
	fTestFunction->SetNpz(fNbinz);


	BCLog::Out(BCLog::summary,BCLog::summary," setup fitting function ...");
	fFitFunction = new TF3("fFitFunction",
			fTestFunction->GetExpFormula(), fXmin, fXmax, fYmin, fYmax,
			fZmin, fZmax);
	for (int i=1; i<fTestFunction->GetNpar(); i++)
		fFitFunction->FixParameter(i,fTestFunction->GetParameter(i));


	BCLog::Out(BCLog::summary,BCLog::summary," add parameter x ...");
	this->AddParameter("x",fXmin,fXmax);
	BCLog::Out(BCLog::summary,BCLog::summary," add parameter y ...");
	this->AddParameter("y",fYmin,fYmax);
	BCLog::Out(BCLog::summary,BCLog::summary," add parameter z ...");
	this->AddParameter("z",fZmin,fZmax);


	BCLog::Out(BCLog::summary,BCLog::summary," setup trees for chains ...");
	SetModel(this);
	BCModelOutput::WriteMarkovChain(true);
	SetFile(outputFile);
	

	BCLog::Out(BCLog::summary,BCLog::summary," book histograms for analysis ...");
	for (int i=0; i<fMaxChains; i++) {
		for (int j=1; j<fMaxLags; j++)
			fHistXYZLags[i][j] = new TH3F(Form("fHistXYZChain%dLag%d",i,j),
					Form("lag of %d",j),
					fNbinx, fXmin, fXmax,
					fNbiny, fYmin, fYmax,
					fNbinz, fZmin, fZmax);

		for (int j=1; j<fMax10thOfIters; j++)
			fHistXYZIter[i][j] = new TH3F(Form("fHistXYZChain%dIter%d0",i,j),
					Form("after %d%% of total iterations",j*10),
					fNbinx, fXmin, fXmax,
					fNbiny, fYmin, fYmax,
					fNbinz, fZmin, fZmax);

		fHChi2vsLags[i] = new TH1F(Form("HChi2vsLagsChain%d",i),"",
				fMaxLags-1,1,fMaxLags);
		fHChi2vsIter[i] = new TH1F(Form("HChi2vsIterChain%d",i),"",
				fMax10thOfIters-1,1,fMax10thOfIters);
		
		fHKolmogorovProbVsLags[i] = new TH1F(
				Form("HKolmogorovProbVsLagsChain%d", i), "",
				fMaxLags-1, 1, fMaxLags);
		fHKolmogorovProbVsIter[i] = new TH1F(
				Form("HKolmogorovProbVsIterChain%d", i), "",
				fMax10thOfIters-1, 1, fMax10thOfIters);
	}
}

//=============================================================================

BCBenchmarkMCMC3D::~BCBenchmarkMCMC3D()
{
	if (fFitFunction) delete fFitFunction;

	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		for (int j=1; j<fMaxLags; j++)
			if (fHistXYZLags[i][j]) delete fHistXYZLags[i][j];

		for (int j=1; j<fMax10thOfIters; j++)
			if (fHistXYZIter[i][j]) delete fHistXYZIter[i][j];

		if (fHChi2vsLags[i]) delete fHChi2vsLags[i];
		if (fHChi2vsIter[i]) delete fHChi2vsIter[i];
		
		if (fHKolmogorovProbVsLags[i]) delete fHKolmogorovProbVsLags[i];
		if (fHKolmogorovProbVsIter[i]) delete fHKolmogorovProbVsIter[i];
	}
}

//=============================================================================

void BCBenchmarkMCMC3D::ProcessMCTrees()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) ProcessMCTree(i);
}

//=============================================================================

void BCBenchmarkMCMC3D::ProcessMCTree(int chainID)
{
	TTree *chain = this->MCMCGetMarkovChainTree(chainID);

	Double_t fMeanX, fMeanY, fMeanZ, fVarianceX, fVarianceY, fVarianceZ,
			fSkewnessX, fSkewnessY, fSkewnessZ;
	chain->Branch("fMeanX",&fMeanX,"fMeanX/D");
	chain->Branch("fMeanY",&fMeanY,"fMeanY/D");
	chain->Branch("fMeanZ",&fMeanZ,"fMeanZ/D");
	chain->Branch("fVarianceX",&fVarianceX,"fVarianceX/D");
	chain->Branch("fVarianceY",&fVarianceY,"fVarianceY/D");
	chain->Branch("fVarianceZ",&fVarianceZ,"fVarianceZ/D");
	chain->Branch("fSkewnessX",&fSkewnessX,"fSkewnessX/D");
	chain->Branch("fSkewnessY",&fSkewnessY,"fSkewnessY/D");
	chain->Branch("fSkewnessZ",&fSkewnessZ,"fSkewnessZ/D");

	Double_t par0, par1, par2;
	chain->SetBranchAddress("fParameter0",&par0);
	chain->SetBranchAddress("fParameter1",&par1);
	chain->SetBranchAddress("fParameter2",&par2);

	Double_t sum1X=0, sum1Y=0, sum1Z=0, sum2X=0, sum2Y=0, sum2Z=0,
			sum3X=0, sum3Y=0, sum3Z=0;
	int Niters = chain->GetEntries();
	for (int i=0; i<Niters; i++) {
		chain->GetEntry(i);
		
		sum1X += par0;
		fMeanX = sum1X/Double_t(i+1);
		
		sum1Y += par1;
		fMeanY = sum1Y/Double_t(i+1);
		
		sum1Z += par2;
		fMeanZ = sum1Z/Double_t(i+1);
		
		sum2X += (par0-fMeanX)*(par0-fMeanX);
		fVarianceX = sqrt(sum2X/Double_t(i+1));
		
		sum2Y += (par1-fMeanY)*(par1-fMeanY);
		fVarianceY = sqrt(sum2Y/Double_t(i+1));

		sum2Z += (par2-fMeanZ)*(par2-fMeanZ);
		fVarianceZ = sqrt(sum2Z/Double_t(i+1));

		sum3X += (par0-fMeanX)*(par0-fMeanX)*(par0-fMeanX);
		if (fVarianceX!=0) 
			fSkewnessX = (sum3X/Double_t(i+1))/(fVarianceX*fVarianceX*fVarianceX);
		
		sum3Y += (par1-fMeanY)*(par1-fMeanY)*(par1-fMeanY);
		if (fVarianceY!=0) 
			fSkewnessY = (sum3Y/Double_t(i+1))/(fVarianceY*fVarianceY*fVarianceY);
		
		sum3Z += (par2-fMeanZ)*(par2-fMeanZ)*(par2-fMeanZ);
		if (fVarianceZ!=0) 
			fSkewnessZ = (sum3Z/Double_t(i+1))/(fVarianceZ*fVarianceZ*fVarianceZ);
		
		chain->Fill();

		for (int j=1; j<fMaxLags; j++) {
			if (i%j==0) fHistXYZLags[chainID][j]->Fill(par0, par1, par2);
		}

		for (int j=1; j<fMax10thOfIters; j++) {
			if (i<j*Niters/10) 
				fHistXYZIter[chainID][j]->Fill(par0, par1, par2);
		}
	}
}

//=============================================================================

void BCBenchmarkMCMC3D::PerformLagsTest()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		Chi2vsLagsOfChain(i);
		KolmogorovVsLagsOfChain(i);
	}
}

//=============================================================================

void BCBenchmarkMCMC3D::Chi2vsLagsOfChain(int chainID)
{
	for (int i=1; i<fMaxLags; i++) {
		fHistXYZLags[chainID][i]->Fit(fFitFunction,"bq");
		fHChi2vsLags[chainID]->SetBinContent(i,
				fFitFunction->GetChisquare()/fFitFunction->GetNDF());
	}
}

//=============================================================================

void BCBenchmarkMCMC3D::PerformIterationsTest()
{
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		Chi2vsIterOfChain(i);
		KolmogorovVsIterOfChain(i);
	}
}

//=============================================================================

void BCBenchmarkMCMC3D::Chi2vsIterOfChain(int chainID)
{
	for (int i=1; i<fMax10thOfIters; i++) {
		fHistXYZIter[chainID][i]->Fit(fFitFunction,"bq");
		fHChi2vsIter[chainID]->SetBinContent(i,
				fFitFunction->GetChisquare()/fFitFunction->GetNDF());
	}
}

//=============================================================================

void BCBenchmarkMCMC3D::WriteResults()
{
	TFile *file = this->GetFile();

	file->cd();
	for (int i=0; i<BCEngineMCMC::fMCMCNChains; i++) {
		file->mkdir(Form("HistsForMCTree%d",i),
				Form("Histograms for Markov chain tree %d",i));
		file->cd(Form("HistsForMCTree%d",i));

		for (int j=1; j<fMaxLags; j++) {
			fHistXYZLags[i][j]->SetXTitle("x");
			fHistXYZLags[i][j]->SetYTitle("y");
			fHistXYZLags[i][j]->SetZTitle("z");
			fHistXYZLags[i][j]->Write();
		}

		for (int j=1; j<fMax10thOfIters; j++) {
			fHistXYZIter[i][j]->SetXTitle("x");
			fHistXYZIter[i][j]->SetYTitle("y");
			fHistXYZIter[i][j]->SetYTitle("z");
			fHistXYZIter[i][j]->Write();
		}

		fHChi2vsLags[i]->SetXTitle("Lags");
		fHChi2vsLags[i]->SetYTitle("#chi^{2}/NDF");
		fHChi2vsLags[i]->Write();

		fHChi2vsIter[i]->SetXTitle("10% of total iteration");
		fHChi2vsIter[i]->SetYTitle("#chi^{2}/NDF");
		fHChi2vsIter[i]->Write();
		
		fHKolmogorovProbVsLags[i]->SetXTitle("Lags");
		fHKolmogorovProbVsLags[i]->SetYTitle("Kolmogorov-Smirnov Probability");
		fHKolmogorovProbVsLags[i]->Write();

		fHKolmogorovProbVsIter[i]->SetXTitle("10% of total iteration");
		fHKolmogorovProbVsIter[i]->SetYTitle("Kolmogorov-Smirnov Probability");
		fHKolmogorovProbVsIter[i]->Write();
	}
	file->cd();

	this->Close();
}

//=============================================================================

/**
 * Perform a Kolmogorov-Smirnov test for all available lagging histograms
 * against a histogram of corresponding bin size containing the
 *"real" data obtained from the test function.
 */
void BCBenchmarkMCMC3D::KolmogorovVsLagsOfChain(int chainID)
{
	TH3F *testHisto = this->GetTestFunctionHistogram();
	
	for (int i=1; i<fMaxLags; i++) {
		double kolmogorovProb = 0.;
		kolmogorovProb = this->KolmogorovTest(fHistXYZLags[chainID][i],
				testHisto, "");

		fHKolmogorovProbVsLags[chainID]->SetBinContent(i,
				kolmogorovProb);
	}
	delete testHisto;
}

//=============================================================================

/**
 * Perform a Kolmogorov-Smirnov tests to see if a histogram of the "real"
 * data matches the distribution after 10%, 20%, 30%, ... of the total
 * number of iterations.
 */
void BCBenchmarkMCMC3D::KolmogorovVsIterOfChain(int chainID)
{
	TH3F *testHisto = this->GetTestFunctionHistogram();
	
	for (int i=1; i<fMax10thOfIters; i++) {
		double kolmogorovProb = 0.;
		kolmogorovProb = this->KolmogorovTest(fHistXYZIter[chainID][i],
				testHisto, "");

		fHKolmogorovProbVsIter[chainID]->SetBinContent(i,
				kolmogorovProb);
	}
	delete testHisto;
}

//=============================================================================

/**
 * Implementation of the Kolmogorov-Smirnov test for 3D histograms while
 * the ROOT implementation is still faulty.
 * This implementation is just taken from ROOT's TH3 class and
 * modified so it works correctly.
 */
Double_t BCBenchmarkMCMC3D::KolmogorovTest(TH1 *h1, TH1 *h2,
		Option_t *option)
{
   //  Statistical test of compatibility in shape between
   //  THIS histogram and h2, using Kolmogorov test.
   //     Default: Ignore under- and overflow bins in comparison
   //
   //     option is a character string to specify options
   //         "U" include Underflows in test
   //         "O" include Overflows
   //         "N" include comparison of normalizations
   //         "D" Put out a line of "Debug" printout
   //         "M" Return the Maximum Kolmogorov distance instead of prob
   //
   //   The returned function value is the probability of test
   //       (much less than one means NOT compatible)
   //
   //   The KS test uses the distance between the pseudo-CDF's obtained 
   //   from the histogram. Since in more than 1D the order for
   //   generating the pseudo-CDF is 
   //   arbitrary, we use the pseudo-CDF's obtained from all the possible
   //   6 combinatons of the 3 axis. 
   //   The average of all the maximum  distances obtained is used in the tests.  

   TString opt = option;
   opt.ToUpper();

   Double_t prb = 0;
//   TH1 *h1 = (TH1*)this;
   if (h1 == 0) return 0;
   if (h2 == 0) return 0;
   TAxis *xaxis1 = h1->GetXaxis();
   TAxis *xaxis2 = h2->GetXaxis();
   TAxis *yaxis1 = h1->GetYaxis();
   TAxis *yaxis2 = h2->GetYaxis();
   TAxis *zaxis1 = h1->GetZaxis();
   TAxis *zaxis2 = h2->GetZaxis();
   Int_t ncx1   = xaxis1->GetNbins();
   Int_t ncx2   = xaxis2->GetNbins();
   Int_t ncy1   = yaxis1->GetNbins();
   Int_t ncy2   = yaxis2->GetNbins();
   Int_t ncz1   = zaxis1->GetNbins();
   Int_t ncz2   = zaxis2->GetNbins();

   // Check consistency of dimensions
   if (h1->GetDimension() != 3 || h2->GetDimension() != 3) {
      printf("KolmogorovTest: Histograms must be 3-D\n");
      return 0;
   }

   // Check consistency in number of channels
   if (ncx1 != ncx2) {
      printf("KolmogorovTest: Number of channels in X is different, %d and %d\n",ncx1,ncx2);
      return 0;
   }
   if (ncy1 != ncy2) {
      printf("KolmogorovTest: Number of channels in Y is different, %d and %d\n",ncy1,ncy2);
      return 0;
   }
   if (ncz1 != ncz2) {
      printf("KolmogorovTest: Number of channels in Z is different, %d and %d\n",ncz1,ncz2);
      return 0;
   }

   // Check consistency in channel edges
   Bool_t afunc1 = kFALSE;
   Bool_t afunc2 = kFALSE;
   Double_t difprec = 1e-5;
   Double_t diff1 = TMath::Abs(xaxis1->GetXmin() - xaxis2->GetXmin());
   Double_t diff2 = TMath::Abs(xaxis1->GetXmax() - xaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      printf("KolmogorovTest: histograms with different binning along X");
      return 0;
   }
   diff1 = TMath::Abs(yaxis1->GetXmin() - yaxis2->GetXmin());
   diff2 = TMath::Abs(yaxis1->GetXmax() - yaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      printf("KolmogorovTest: histograms with different binning along Y");
      return 0;
   }
   diff1 = TMath::Abs(zaxis1->GetXmin() - zaxis2->GetXmin());
   diff2 = TMath::Abs(zaxis1->GetXmax() - zaxis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      printf("KolmogorovTest: histograms with different binning along Z");
      return 0;
   }

   //   Should we include Uflows, Oflows?
   Int_t ibeg = 1, jbeg = 1, kbeg = 1;
   Int_t iend = ncx1, jend = ncy1, kend = ncz1;
   if (opt.Contains("U")) {ibeg = 0; jbeg = 0; kbeg = 0;}
   if (opt.Contains("O")) {iend = ncx1+1; jend = ncy1+1; kend = ncz1+1;}

   Int_t i,j,k,bin;
   Double_t sum1  = 0;
   Double_t sum2  = 0;
   Double_t w1    = 0;
   Double_t w2    = 0;
   for (i = ibeg; i <= iend; i++) {
      for (j = jbeg; j <= jend; j++) {
         for (k = kbeg; k <= kend; k++) {
            bin = h1->GetBin(i,j,k);
            sum1 += h1->GetBinContent(bin);
            sum2 += h2->GetBinContent(bin);
            Double_t ew1   = h1->GetBinError(bin);
            Double_t ew2   = h2->GetBinError(bin);
            w1   += ew1*ew1;
            w2   += ew2*ew2;
         }
      }
   }


   //    Check that both scatterplots contain events
   if (sum1 == 0) {
      printf("KolmogorovTest: Integral is zero for h1=%s\n",h1->GetName());
      return 0;
   }
   if (sum2 == 0) {
      printf("KolmogorovTest: Integral is zero for h2=%s\n",h2->GetName());
      return 0;
   }
   // calculate the effective entries.  
   // the case when errors are zero (w1 == 0 or w2 ==0) are equivalent to 
   // compare to a function. In that case the rescaling is done only on sqrt(esum2) or sqrt(esum1) 
   Double_t esum1 = 0, esum2 = 0; 
   if (w1 > 0) 
      esum1 = sum1 * sum1 / w1; 
   else 
      afunc1 = kTRUE;    // use later for calculating z
   
   if (w2 > 0) 
      esum2 = sum2 * sum2 / w2; 
   else 
      afunc2 = kTRUE;    // use later for calculating z
   
   if (afunc2 && afunc1) { 
      printf("KolmogorovTest: Errors are zero for both histograms\n");
      return 0;
   }

   //   Find Kolmogorov distance
   //   order is arbitrary take average of all possible 6 starting orders x,y,z 
   int order[3] = {0,1,2};
   int binbeg[3]; 
   int binend[3]; 
   int ibin[3];
   binbeg[0] = ibeg; binbeg[1] = jbeg; binbeg[2] = kbeg; 
   binend[0] = iend; binend[1] = jend; binend[2] = kend; 
   Double_t vdfmax[6]; // there are in total 6 combinations 
   int icomb = 0; 
   Double_t s1 = 1/sum1;
   Double_t s2 = 1/sum2;
   Double_t rsum1=0, rsum2=0;
   do { 
      // loop on bins
      Double_t dmax = 0;
      for (i = binbeg[order[0] ]; i <= binend[order[0] ]; i++) {
         for ( j = binbeg[order[1] ]; j <= binend[order[1] ]; j++) {
            for ( k = binbeg[order[2] ]; k <= binend[order[2] ]; k++) {
                  ibin[ order[0] ] = i;
                  ibin[ order[1] ] = j;
                  ibin[ order[2] ] = k;
                  bin = h1->GetBin(ibin[0],ibin[1],ibin[2]);
                  rsum1 += s1*h1->GetBinContent(bin);
                  rsum2 += s2*h2->GetBinContent(bin);
                  dmax   = TMath::Max(dmax, TMath::Abs(rsum1-rsum2));
            }
         }
      }
      vdfmax[icomb] = dmax; 
      icomb++;
   } while (TMath::Permute(3,order)  );


   // get average of distances 
   Double_t dfmax = TMath::Mean(6,vdfmax);
   
   //    Get Kolmogorov probability
   Double_t factnm;
   if (afunc1)      factnm = TMath::Sqrt(sum2);
   else if (afunc2) factnm = TMath::Sqrt(sum1);
   else             factnm = TMath::Sqrt(sum1*sum2/(sum1+sum2));
   Double_t z  = dfmax*factnm;

   prb = TMath::KolmogorovProb(z); 

   Double_t prb1 = 0, prb2 = 0; 
   // option N to combine normalization makes sense if both afunc1 and afunc2 are false
   if (opt.Contains("N")  && !(afunc1 || afunc2 ) ) { 
      // Combine probabilities for shape and normalization
      prb1   = prb;
      Double_t d12    = esum1-esum2;
      Double_t chi2   = d12*d12/(esum1+esum2);
      prb2   = TMath::Prob(chi2,1);
      //     see Eadie et al., section 11.6.2
      if (prb > 0 && prb2 > 0) prb = prb*prb2*(1-TMath::Log(prb*prb2));
      else                     prb = 0;
   }

   //    debug printout
   if (opt.Contains("D")) {
      printf(" Kolmo Prob  h1 = %s, sum1=%g\n",h1->GetName(),sum1);
      printf(" Kolmo Prob  h2 = %s, sum2=%g\n",h2->GetName(),sum2);
      printf(" Kolmo Probabil = %f, Max Dist = %g\n",prb,dfmax);
      if (opt.Contains("N"))
         printf(" Kolmo Probabil = %f for shape alone, =%f for normalisation alone\n",prb1,prb2);
   }
   // This numerical error condition should never occur:
   if (TMath::Abs(rsum1-6) > 0.004) printf("KolmogorovTest: Numerical problems with h1=%s\n",h1->GetName());
   if (TMath::Abs(rsum2-6) > 0.004) printf("KolmogorovTest: Numerical problems with h2=%s\n",h2->GetName());

   if(opt.Contains("M"))      return dfmax;  // return avergae of max distance

   return prb;
}

//=============================================================================

/**
 * ROOT's built-in method to create a histogram (TH3) from a 3-dimensional
 * function does not work correctly. This method does the same job
 * and creates a histogram from the test function fTestFunction.
 */
TH3F* BCBenchmarkMCMC3D::GetTestFunctionHistogram()
{
	TH3F* testHisto = (TH3F*)fTestFunction->CreateHistogram();
	int nbinsX = testHisto->GetXaxis()->GetNbins();
	int nbinsY = testHisto->GetYaxis()->GetNbins();
	int nbinsZ = testHisto->GetZaxis()->GetNbins();
	
	double Xmin = testHisto->GetXaxis()->GetXmin();
	double Ymin = testHisto->GetYaxis()->GetXmin();
	double Zmin = testHisto->GetZaxis()->GetXmin();
	double Xmax = testHisto->GetXaxis()->GetXmax();
	double Ymax = testHisto->GetYaxis()->GetXmax();
	double Zmax = testHisto->GetZaxis()->GetXmax();
	
	for (int i = 0; i < nbinsX; i++) {
		double coord_x = Xmin + (Xmax - Xmin) / (double)nbinsX * (i + 0.5);
		for (int j = 0; j < nbinsY; j++) {
			double coord_y = Ymin + (Ymax - Ymin) / (double)nbinsY * (j + 0.5);
			for (int k = 0; k < nbinsZ; k++) {
				double coord_z = Zmin + (Zmax - Zmin) / (double)nbinsZ * (k + 0.5);
				testHisto->Fill(coord_x, coord_y, coord_z,
					fTestFunction->Eval(coord_x, coord_y, coord_z));
			}
		}
	}
	return testHisto;
}

