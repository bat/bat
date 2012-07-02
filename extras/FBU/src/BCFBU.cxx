// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "BCFBU.h"

#include <BAT/BCMath.h>
//#include <BAT/BCDataSet.h>
//#include <BAT/BCDataPoint.h>
#include <BAT/BCLog.h>

#include <TFile.h>
#include <TVector.h>
#include <TMatrix.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <map>

// ---------------------------------------------------------
BCFBU::BCFBU() 
	: BCModel("BFU")
	, fHistTruth(0)
	, fNBinsTruth(0)
	, fNBinsTruth1(0)
	, fNBinsTruth2(0)
	, fTruthMin(0)
	, fTruthMax(0)
	, fNBinsReco(0)
	, fNBinsReco1(0)
	, fNBinsReco2(0)
	, fResponseMatrix(0)
	, fMigrationMatrix(0)
	, fHistEfficiency(0)
	, fHistRatio(0)
	, fAcTotal(0)
	, fNDim(-1)
	, fNSyst(0)
	, fNSamples(0)
	, fSystNames(0)
	, fSampleNames(0)
	, fNominalHisto(0)
	, fNNormSyst(0)
	, fNormSystNames(0)
	, fNormSystSizes(0)
	, fInterpolationType(1)
{
};

// ---------------------------------------------------------
BCFBU::BCFBU(const char * name) 
	: BCModel(name)
	, fHistTruth(0)
	, fNBinsTruth(0)
	, fNBinsTruth1(0)
	, fNBinsTruth2(0)
	, fTruthMin(0)
	, fTruthMax(0)
	, fNBinsReco(0)
	, fNBinsReco1(0)
	, fNBinsReco2(0)
	, fResponseMatrix(0)
	, fMigrationMatrix(0)
	, fHistEfficiency(0)
	, fHistRatio(0)
	, fAcTotal(0)
	, fNDim(-1)
	, fNSyst(0)
	, fNSamples(0)
	, fSystNames(0)
	, fSampleNames(0)
	, fNominalHisto(0)
	, fNNormSyst(0)
	, fNormSystNames(0)
	, fNormSystSizes(0)
	, fInterpolationType(1)
{
};

// ---------------------------------------------------------
BCFBU::~BCFBU()
{
	if (fHistTruth)
		delete fHistTruth;

	if (fHistEfficiency)
		delete fHistEfficiency;
};

// ---------------------------------------------------------
int BCFBU::RebinHistograms(int rebin, TH2* h_migration, TH1 *h_truth, TH1 *h_background) 
{
	// check if migration matrix histogram exists 
  if (!h_migration) {
		BCLog::OutWarning("BCFBU::RebinHistograms: migration matrix histogram not found. Exit.");
		exit(1);
	}
  
	// check if truth histogram exists 
  if (!h_truth){
		BCLog::OutWarning("BCFBU::RebinHistograms: truth histogram not found. Exit.");
		exit(1);
	}

	// do not check if background histogram exists since it might not be
	// necessary to have
  if (rebin!=1) {
		h_migration->Rebin2D(rebin,rebin);
		h_truth->Rebin(rebin);
		if (h_background)
			h_background->Rebin(rebin);
	}

	// no error
	return 1;
}

// ---------------------------------------------------------
int BCFBU::PrepareResponseMatrix(TH2* h_migration, TH1* h_truth, TH1* h_background, std::vector<double> parmin, std::vector<double> parmax)
{
	// check if migration matrix histogram exists 
  if (!h_migration) {
		BCLog::OutWarning("BCFBU::PrepareResponseMatrix: migration matrix histogram not found.");
		return 0;
	}
  
	// check if truth histogram exists 
  if (!h_truth){
		BCLog::OutWarning("BCFBU::PrepareResponseMatrix: truth histogram not found.");
		return 0;
	}

	// check dimension
	if( (fNDim > 0) && (fNDim != h_truth->GetDimension()) ) {
		BCLog::OutWarning("BCBFU::PrepareResponseMatrix. Dimensions of truth histogram does not match previously defined dimension. Truth histogram rejected.");
		return 0; 
	}
	else if (fNDim <= 0) {
		// set dimension
		fNDim = h_truth->GetDimension();
	}

	// check number of reco bins
	if ( (fNBinsReco > 0) && (fNBinsReco != h_migration->GetYaxis()->GetNbins()) ) {
		BCLog::OutWarning("BCBFU::SetDataHistogram. Number of reco bins doens't match number of bins in the migration matrix. Migration matrix rejected.");
		return 0; 
		}

	// clone truth distribution
	if (fNDim==1) {
		fHistTruth = (TH1D*) h_truth->Clone();
		int nbins = fHistTruth->GetNbinsX();
		fRecoMin = fHistTruth->GetXaxis()->GetBinLowEdge(1);
		fRecoMax = fHistTruth->GetXaxis()->GetBinUpEdge(nbins);
	}
	else if (fNDim==2)
		fHistTruth = (TH2D*) h_truth->Clone();

	// get number of bins
  fNBinsTruth = h_migration->GetXaxis()->GetNbins();
  fNBinsReco  = h_migration->GetYaxis()->GetNbins();

	// prepare helper varible for log likelihood calculation
	fVectorTruth.clear();
	fVectorTruth.assign(fNBinsTruth, 0);
	
	fVectorReco.clear();
	fVectorReco.assign(fNBinsReco, 0);
  
	// print statements
  BCLog::OutDetail(Form("Dimension of problem : %i", fNDim));
  BCLog::OutDetail(Form("Number of truth bins : %i", fNBinsTruth));
  BCLog::OutDetail(Form("Number of reco bins  : %i", fNBinsReco));

	// bring response matrix into matrix form (TMatrixD)
  fResponseMatrix = new TMatrixD(fNBinsTruth, fNBinsReco);

	// calculate migration matrix 
  fMigrationMatrix = new TMatrixD(fNBinsTruth, fNBinsReco);
  
	// bring unnormalized migration matrix into matrix form (TMatrixD)
	// and switch axes
  TMatrix m = TMatrixD(fNBinsTruth, fNBinsReco);

  for (int i_tru = 0; i_tru < fNBinsTruth; ++i_tru) 
    for (int i_rec = 0; i_rec < fNBinsReco; ++i_rec) 
      m(i_tru,i_rec) = h_migration->GetBinContent(i_tru+1, i_rec+1);
	
	// normalization factor of h_migration
	double norm = 0;

	// normalize migration matrix
  for (int i_tru=0;i_tru<fNBinsTruth;++i_tru) {
		// sum over all reco bins
		double sum = 0;

		// increase sum over all reco bins
		for (int i_rec=0;i_rec<fNBinsReco;++i_rec)  
			sum += m(i_tru,i_rec);

		// increase norm		
		norm+= sum; 

		// normalize response matrix
		for (int i_rec=0;i_rec<fNBinsReco;++i_rec)
			(*fResponseMatrix)(i_tru,i_rec) = m(i_tru,i_rec)/sum;
	}
	
	// normalize migration matrix
  for (int i_tru=0;i_tru<fNBinsTruth;++i_tru) 
		for (int i_rec=0;i_rec<fNBinsReco;++i_rec)
			(*fMigrationMatrix)(i_tru,i_rec) = m(i_tru,i_rec)/norm;
	
	// debugKK: really print them? 
	//  m.Print();
	//  fMigrationMatrix->Print();

	// calculate efficiency for 1D case
  if (fNDim==1) {
		// create new efficiency histogram based on the same binning as
		// the truth distribution
		fHistEfficiency = (TH1D*) fHistTruth->Clone();
		fHistEfficiency->GetYaxis()->SetTitle("efficiency");
		
		// loop over all truth bins
		for (int i_tru=0;i_tru<fNBinsTruth;++i_tru) {
			double sum = 0;
			
			// sum all reco bins
			for (int i_rec=0;i_rec<fNBinsReco;++i_rec)
				sum += m(i_tru,i_rec);
			
			// calculate efficiency
			double eff = sum/fHistTruth->GetBinContent(i_tru+1); 

			// set efficiency
			fHistEfficiency->SetBinContent(i_tru+1, eff);

			// check parameter ranges
			double mini = 0.;
			double maxi = 20e4;

			if ( (parmin.size()) > 0 && (parmin.size() == parmax.size()) ) {
				mini = parmin[i_tru];
				maxi = parmax[i_tru];
			}

			// add parameter
			AddParameter(("T"+IntToString(i_tru+1)).c_str(), mini, maxi);
		} 
	}
  
	// debugKK: need to work on 2D part
  if (fNDim==2)
    {

      fHistEfficiency = (TH2D*) fHistTruth->Clone();
  
      std::cout << "dimensions of h_migration x = " << h_migration->GetXaxis()->GetNbins() << " y = " << h_migration->GetYaxis()->GetNbins() << " response(3,3) " << h_migration->GetBinContent(3,3) << std::endl;
      
      std::cout << "dimensions of h_migration, after rebinning x = " << h_migration->GetXaxis()->GetNbins() << " y = " << h_migration->GetYaxis()->GetNbins() << " reponse(3,3) " << h_migration->GetBinContent(3,3) << std::endl;
      
      fNBinsTruth1 = fHistTruth->GetXaxis()->GetNbins();
      fNBinsTruth2 = fHistTruth->GetYaxis()->GetNbins();
      
      fNBinsReco1 = h_background->GetXaxis()->GetNbins();
      fNBinsReco2 = h_background->GetYaxis()->GetNbins();
      
      std::cout << " truth bins X " << fNBinsTruth1 << " Y " << fNBinsTruth2 << std::endl;
      std::cout << " reco bins X " << fNBinsReco1 << " Y " << fNBinsReco2 << std::endl;
      
      for (int i_tru1=0;i_tru1<fNBinsTruth1;i_tru1++)
				for (int i_tru2=0;i_tru2<fNBinsTruth2;i_tru2++)
					{
						double sum = 0;
	    
						for (int i_rec1=0;i_rec1<fNBinsReco1;i_rec1++)
							for (int i_rec2=0;i_rec2<fNBinsReco2;i_rec2++)  
								{
									std::cout << "index x y " << i_tru1 << " " << i_tru2 << " " << i_rec1 << " " << i_rec2 << " " << Get2DIndex(i_tru1,i_tru2,i_rec1, i_rec2).first << " " << 	      Get2DIndex(i_tru1,i_tru2,i_rec1, i_rec2).second << std::endl;
		  
									sum += m(Get2DIndex(i_tru1,i_tru2,i_rec1, i_rec2).first,
													 Get2DIndex(i_tru1,i_tru2,i_rec1, i_rec2).second);
		  
								}
	    
						fHistEfficiency->SetBinContent(i_tru1+1,i_tru2+1, sum/fHistTruth->GetBinContent(i_tru1+1,i_tru2+1));
	    
						std::cout << "total truth " << fHistTruth->GetBinContent(i_tru1+1,i_tru2+1) << " efficiency bin " << i_tru1 << " " << i_tru2 << " " << fHistEfficiency->GetBinContent(i_tru1+1,i_tru2+1) << std::endl;
	    
					}  
            
      
      for (int i_tru1=0;i_tru1<fNBinsTruth1;i_tru1++)
				for (int i_tru2=0;i_tru2<fNBinsTruth2;i_tru2++)
					{
						std::cout << "truth bin " << i_tru1 << " " << i_tru2 << " " << fHistTruth->GetBinContent(i_tru1+1,i_tru2+1) << std::endl;
					}
      
      for (int j=0;j<fNBinsTruth2;++j)
				{
					TH1D* hist = new TH1D("", "A_C", 500, -0.35, 0.35);
	  
					BCH1D *histptr = new BCH1D();
	  
					histptr->SetHistogram(hist);
	  
					fHistRatio.push_back(histptr);
	  
				}
    }

	// debugKK: can the lines below be removed?
	//----
  TH1D* totalac = new TH1D("", "A_C", 500, -0.35, 0.35);
  
  fAcTotal = new BCH1D();
  
  fAcTotal->SetHistogram(totalac);
	//----

	// no error 
	return 1;
}

// ---------------------------------------------------------
void BCFBU::DefineParameters()
{
	// add parameters for 1D case
  if (fNDim==1) {
		// add one parameter per truth bin
		for (int i=1; i<=fNBinsTruth; ++i) {
			// debugKK: how to set the upper range?
			AddParameter(("T"+IntToString(i)).c_str(), 0, 20e4);
		}

		// add one parameter per source of systematic shape uncertainty
		for (int i=0; i<fNSyst; ++i) {
			AddParameter(fSystNames[i].c_str(), -5, 5);
			// debugKK: what is this?
			fSystParamMap[fSystNames[i]] = GetNBinsTruth() + i; 
		}
      
		// add one parameter per source of systematic normalization
		// uncertainty
		for (int i=0;i<fNNormSyst;++i) {
			AddParameter(fNormSystNames[i].c_str(), -5, 5);
			// debugKK: what is this?
			fSystParamMap[fNormSystNames[i]] = GetNBinsTruth() + fNSyst + i;
		}
	}

	// debugKK: did not check the 2D case yet
  if (fNDim==2)
    {
      for (int i=1;i<=GetNBinsTruthX();++i) // parameters corresponding to bins of sought after truth histogram
				for (int j=1;j<=GetNBinsTruthY();++j)
					{      
						AddParameter(("X"+IntToString(i)+"Y"+IntToString(j)).c_str(), 0, 20e4);
					}
  
      
      for (int i=0;i<fNSyst;++i)  // shape systematics
				{
					AddParameter(fSystNames[i].c_str(), -5, 5);
					fSystParamMap[fSystNames[i]] = GetNBinsTruthX()*GetNBinsTruthY() + i;
				}
      
      
      for (int i=0;i<fNNormSyst;++i) // normalisation systematics
				{
					AddParameter(fNormSystNames[i].c_str(), -5, 5);
					fSystParamMap[fNormSystNames[i]] = GetNBinsTruthX()*GetNBinsTruthY() + fNSyst + i;
	  
				}
    }
}

// ---------------------------------------------------------
double BCFBU::LinearInterpolate(double alpha, double nominal, double up, double down)
{
  double res = 0;

  if (alpha>=0)
    res = alpha*(up-nominal);
  else
    res = alpha*(nominal-down);

  return res;
}

// ---------------------------------------------------------
double BCFBU::ExponentialInterpolate(double alpha, double nominal, double up, double down)
{
  double res = 0;

  if (alpha>=0)
    res = exp(alpha*log(up/nominal));
  else
    res = exp(-alpha*log(down/nominal));

  return res;
}

// ---------------------------------------------------------
double BCFBU::NoInterpolation(double alpha, double nominal, double up, double down)
{
  double res = 0;

  if (alpha>=1)
    res = up;
  else if (alpha<=-1)
    res = down;
  else
    res = nominal;

  return res;
}

// ---------------------------------------------------------
double BCFBU::LogLikelihood(const std::vector<double> & parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

  double logprob = 0.;

  if (fNDim==2)
    {
      TMatrixD t = TMatrixD(fNBinsTruth1,fNBinsTruth2);
      
      for (int i=0;i<fNBinsTruth1;++i)
				for (int j=0;j<fNBinsTruth2;++j)
					t(i,j) = parameters.at(Get1DIndex(i,j));
      
      
      /*
				vector<double> alpha;
	
				for (int i=0;i<fNSyst;++i)
				{
				alpha.push_back(parameters.at(GetNBinsTruthX()*GetNBinsTruthY()+i));
				}
      */
      
  
      TMatrixD d = TMatrixD(fNBinsReco1,fNBinsReco2);
      
      TMatrixD r = TMatrixD(fNBinsReco1, fNBinsReco2);
      
      for (int i=0;i<fNBinsReco1;++i)
				for (int j=0;j<fNBinsReco2;++j)
					{
						d(i,j) = GetDataPoint(0)->GetValue(Get1DIndex(i,j));
	
					}
      
      
      
      for (int i_rec1=0;i_rec1<fNBinsReco1;i_rec1++) // loop over all reco bins
				for (int i_rec2=0;i_rec2<fNBinsReco2;i_rec2++)
					{
						double sum = 0;
	    
						for (int i_tru1=0;i_tru1<fNBinsTruth1;i_tru1++) // loop over truth bins
							for (int i_tru2=0;i_tru2<fNBinsTruth2;i_tru2++)
								{
		  
									double nominalresp = (*fMigrationMatrix)(Get2DIndex(i_tru1,i_tru2,i_rec1,i_rec2).first, Get2DIndex
																								(i_tru1,i_tru2,i_rec1,i_rec2).second ); // find appropriate entry of response matrix
		  
									double varresp = nominalresp;
		  
		  
									// loop through the systematics, compute shifts in response matrix, add them to nominal response
		  
									for( std::map<std::string, TH2*>::iterator it_map2=fSystResponseUpMap.begin(); it_map2!=fSystResponseUpMap.end(); it_map2++)
										{
											//	std::cout << "Now doing systematic with name " << (*it_map2).first << std::endl;
		      
											std::string curr_syst = (*it_map2).first;
		      
											int syst_param_num = GetSystParamNumber(curr_syst);
		      
											double responseup = fSystResponseUpMap[curr_syst]->GetBinContent(Get2DIndex(i_tru1,i_tru2,i_rec1,i_rec2).first+1, Get2DIndex
																																											 (i_tru1,i_tru2,i_rec1,i_rec2).second+1); // find shifted up response for current syst
		      
		      
											double responsedown = fSystResponseDownMap[curr_syst]->GetBinContent(Get2DIndex(i_tru1,i_tru2,i_rec1,i_rec2).first+1, Get2DIndex
																																													 (i_tru1,i_tru2,i_rec1,i_rec2).second+1); // find shifted down response for current syst
		      
											//		  std::cout << "nominal response " << nominalresp << " up " << responseup << " down " << responsedown << std::endl;
		      
											// find interpolated response matrix
		      
											if (fInterpolationType==0) 
												varresp += LinearInterpolate(parameters.at(syst_param_num),nominalresp,responseup,responsedown);
											else if (fInterpolationType==1)
												varresp *= ExponentialInterpolate(parameters.at(syst_param_num),nominalresp,responseup,responsedown);
											else if (fInterpolationType==2)
												varresp = NoInterpolation(parameters.at(syst_param_num),nominalresp,responseup,responsedown);
											else
												{
													std::cout << "Unknown extrapolation type, exiting" << std::endl;
													exit(1);
												}
		      
										}
	      
	      
									sum += t(i_tru1,i_tru2)*fHistEfficiency->GetBinContent(i_tru1+1,i_tru2+1)*varresp; // this is the contribution to reco bin (i_rec1,i_rec2) from truth bin (i_tru1,i_tru2)
	      
								}

						for (int i_sam=0;i_sam<fNSamples; i_sam++) // loop over all background samples
							{	    
	    
								std::string curr_sample = fSampleNames[i_sam];


								// there must be a nominal histo for each sample
	    
								double nominal = fNominalHisto[i_sam]->GetBinContent(i_rec1+1,i_rec2+1);

								//	    std::cout << curr_sample << " " << nominal << std::endl;
	    
								// check if affected by normalisation systematic
	    
								double norm = 1.0;
			
								if (fNNormSyst>0) {
									for( std::map<std::string, std::string >::iterator it_map=fNormUncer.begin(); it_map!=fNormUncer.end(); it_map++)
										{
											if ((*it_map).first==curr_sample)
												{
								
													int paramnumber = GetSystParamNumber(fNormUncer[curr_sample]);

													// std::cout << " normsyst for sample " << curr_sample << " parameter number " << paramnumber << " norm value " << 1 + parameters.at(paramnumber) << std::endl;
			
													norm *= 1 + parameters.at(paramnumber);
			
												}		     		    
										}
		  
								}
	    
								//	    std::cout << curr_sample << " nominal " << nominal << " nominal*norm " << nominal*norm << std::endl;

								r(i_rec1,i_rec2) += norm * nominal;
	    
	   
 
								for( std::map<std::string, std::map<std::string, TH1*> >::iterator it_map=fSystUpHisto.begin(); it_map!=fSystUpHisto.end(); it_map++)
									{
										//		std::cout << (*it_map).first;

										std::map<std::string,TH1*> map_sys = (*it_map).second;
		
										if ( ((*it_map).first)==curr_sample)  // for each sample, loop over the systematics affecting it
											{
		    
												for( std::map<std::string, TH1*>::iterator it_map2=map_sys.begin(); it_map2!=map_sys.end(); it_map2++)
													{
														//			std::cout << "Now doing systematic with name " << (*it_map2).first << " for sample " << curr_sample << " param number " << GetSystParamNumber((*it_map2).first) << 			  " param value " << parameters.at(GetSystParamNumber((*it_map2).first)) << std::endl;
			
														std::string curr_syst = (*it_map2).first;
			
														double up = fSystUpHisto[curr_sample][curr_syst]->GetBinContent(i_rec1+1,i_rec2+1);
														double down = fSystDownHisto[curr_sample][curr_syst]->GetBinContent(i_rec1+1,i_rec2+1);
			
														// find parameter number this systematic corresponds to
			
														int syst_param_num = GetSystParamNumber(curr_syst);


														if (fInterpolationType==0)
															r(i_rec1,i_rec2) += norm * LinearInterpolate(parameters.at(syst_param_num),nominal,up,down);
														else if (fInterpolationType==1)
															{
																r(i_rec1,i_rec2) *= ExponentialInterpolate(parameters.at(syst_param_num),nominal,up,down);
															}
														else if (fInterpolationType==2)
															{
																r(i_rec1,i_rec2) += NoInterpolation(parameters.at(syst_param_num),nominal,up,down);
															}
													}		  		  
											}
									}
	    
								//	    std::cout << " contrib from sample " << curr_sample <<  " " <<  r(i_rec1,i_rec2)-contrib << std::endl;
	   
							}

						r(i_rec1,i_rec2) += sum;

					}
      
      
      for (int i=0;i<fNBinsReco1;++i)
				for (int j=0;j<fNBinsReco2;++j)
					{
						if (r(i,j)>0)
							logprob += -r(i,j)+d(i,j)*log(r(i,j)) - BCMath::LogFact(d(i,j));
						else
							{
								logprob = -9999.0;
								//	    std::cout << " likelihood contains log(-1) " << std::endl;
								//	    exit(1);
							}
					}
    }
  

  if (fNDim==1)
    {
      for (int i=0;i<fNBinsTruth;++i) 
 				fVectorTruth[i] = parameters[i];

      for (int i_rec=0;i_rec<fNBinsReco;++i_rec)
				{
					double sum = 0;
	  
					fVectorReco[i_rec] = 0;

					for (int i_tru=0;i_tru<fNBinsTruth;++i_tru)
						{
	  
							double nominalresp = (*fMigrationMatrix)(i_tru, i_rec); // find appropriate entry of response matrix
	      
							double varresp = nominalresp;
	      
	      
							// loop through the systematics, compute shifts in response matrix, add them to nominal response
	      
							for( std::map<std::string, TH2*>::iterator it_map2=fSystResponseUpMap.begin(); it_map2!=fSystResponseUpMap.end(); it_map2++)
								{
									std::string curr_syst = (*it_map2).first;
		  
									int syst_param_num = GetSystParamNumber(curr_syst);
		  
									double responseup = fSystResponseUpMap[curr_syst]->GetBinContent(i_tru+1, i_rec+1); // find shifted up response for current syst
		  
		  
									double responsedown = fSystResponseDownMap[curr_syst]->GetBinContent(i_tru+1, i_rec+1); // find shifted down response for current syst
		  
									// find interpolated response matrix
		  
									if (fInterpolationType==0) 
										varresp += LinearInterpolate(parameters.at(syst_param_num),nominalresp,responseup,responsedown);
									else if (fInterpolationType==1)
										varresp *= ExponentialInterpolate(parameters.at(syst_param_num),nominalresp,responseup,responsedown);
									else if (fInterpolationType==2)
										varresp = NoInterpolation(parameters.at(syst_param_num),nominalresp,responseup,responsedown);
									else
										{
											std::cout << "Unknown extrapolation type, exiting" << std::endl;
											exit(1);
										}
		
								}

							sum += fVectorTruth[i_tru]*fHistEfficiency->GetBinContent(i_tru+1)*varresp;

							//	      std::cout << " bin " << i_rec << " truth contrib " << sum << " truth " << t(i_tru) << " eff " << fHistEfficiency->GetBinContent(i_tru+1) << " response " << varresp << std::endl;
						}
	      
					for (int i_sam=0;i_sam<fNSamples; i_sam++) // loop over all background samples
						{	    
		  
							std::string curr_sample = fSampleNames[i_sam];
		  
							// there must be a nominal histo for each sample
		  
							double nominal = fNominalHisto[i_sam]->GetBinContent(i_rec+1);
		  
							// check if affected by normalisation systematic
		  
							double norm = 1.0;
		  
							if (fNNormSyst>0)
								{
									for( std::map<std::string, std::string >::iterator it_map=fNormUncer.begin(); it_map!=fNormUncer.end(); it_map++)
										{
											if ((*it_map).first==curr_sample)
												{
			      
													int paramnumber = GetSystParamNumber(fNormUncer[curr_sample]);
			      
													norm *= 1 + parameters.at(paramnumber);
			      
												}		     		    
										}
		      
								}
		  
							fVectorReco[i_rec] += norm * nominal;

							for( std::map<std::string, std::map<std::string, TH1*> >::iterator it_map=fSystUpHisto.begin(); it_map!=fSystUpHisto.end(); it_map++)
								{
		      
									std::map<std::string,TH1*> map_sys = (*it_map).second;
		      
									if ( ((*it_map).first)==curr_sample)  // for each sample, loop over the systematics affecting it
										{
			  
											for( std::map<std::string, TH1*>::iterator it_map2=map_sys.begin(); it_map2!=map_sys.end(); it_map2++)
												{
			      
													std::string curr_syst = (*it_map2).first;
			      
													double up = fSystUpHisto[curr_sample][curr_syst]->GetBinContent(i_rec+1);
													double down = fSystDownHisto[curr_sample][curr_syst]->GetBinContent(i_rec+1);
			      
													// find parameter number this systematic corresponds to
			      
													int syst_param_num = GetSystParamNumber(curr_syst);
			      
													if (fInterpolationType==0)
														fVectorReco[i_rec] += norm * LinearInterpolate(parameters.at(syst_param_num),nominal,up,down);
													else if (fInterpolationType==1)
														{
															fVectorReco[i_rec] *= ExponentialInterpolate(parameters.at(syst_param_num),nominal,up,down);
														}
													else if (fInterpolationType==2)
														{
															fVectorReco[i_rec] += NoInterpolation(parameters.at(syst_param_num),nominal,up,down);
														}
												}		  		  
										}
								}
						}
	      
					fVectorReco[i_rec] += sum;
				}
     

  
      for (int i=0;i<fNBinsReco;++i)
				{
					//	    std::cout << i << " " << r(i) << " " << d(i) << std::endl;
					double r = fVectorReco[i];
					if (r>0) {
						double d = fHistData->GetBinContent(i+1);
						logprob += BCMath::LogPoisson(d, r);
					}
					else
						logprob = -9999.0;
				}

    }
      
  return logprob;
}

// ---------------------------------------------------------
double BCFBU::LogAPrioriProbability(const std::vector <double> & parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	if (fNDim==2)
		{
			for (int i = 0; i < GetNBinsTruthX()*GetNBinsTruthY(); ++i)
				logprob -= log(GetParameter(i)->GetRangeWidth());

			// gaussian profile for systematics
    
			for (int i = GetNBinsTruthX()*GetNBinsTruthY(); i < GetNBinsTruthX()*GetNBinsTruthY() + fNSyst; ++i)
				logprob += BCMath::LogGaus(parameters.at(i), 0., 1.0);


			// std::cout <<  GetNBinsTruthX() << " " << GetNBinsTruthY() << " " << fNSyst << " " << GetNParameters() << std::endl;
   
			for (unsigned int i = GetNBinsTruthX()*GetNBinsTruthY() + fNSyst; i < GetNParameters(); ++i)
				{   
					// std::cout << " in loop " << i << " " << fNormSystSizes[i-(GetNBinsTruthX()*GetNBinsTruthY() + fNSyst)] << std::endl;
       
					if ((1+parameters.at(i))>0)
						logprob += BCMath::LogGaus(parameters.at(i), 0., fNormSystSizes[i - (GetNBinsTruthX()*GetNBinsTruthY() + fNSyst)]);
					else
						logprob = -9999.0;
				}

		}


	if (fNDim==1)
		{
			for (int i = 0; i < GetNBinsTruth(); ++i)
				logprob -= log(GetParameter(i)->GetRangeWidth());
	
			// gaussian profile for systematics
	

			for (int i = GetNBinsTruth(); i < GetNBinsTruth() + fNSyst; ++i)
				logprob += BCMath::LogGaus(parameters.at(i), 0., 1.0);
	
       
			for (unsigned int i = GetNBinsTruth() + fNSyst; i < GetNParameters(); ++i)
				{   	    
					if ((1+parameters.at(i))>0)
						logprob += BCMath::LogGaus(parameters.at(i), 0., fNormSystSizes[i - (GetNBinsTruth() + fNSyst)]);
					else
						logprob = -9999.0;
				}
     
		}


	return logprob;
}
// ---------------------------------------------------------
void BCFBU::MCMCIterationInterface()
{
  // get number of chains
  int nchains = MCMCGetNChains();
  
  // get number of parameters
  int npar = GetNParameters();
     
  // loop over all chains and fill histogram
  for (int i = 0; i < nchains; ++i) {
    // get the current values of the parameters x and y. These are
    // stored in fMCMCx.
    double sum1, sum2;
    
    double sum1tot, sum2tot;

    sum1tot = 0; sum2tot = 0;

    for (int j=0;j<fNBinsTruth2;++j)
      {
				sum1 = 0;
				sum2 = 0;

				/*
					for (int i_tru=0;i_tru<(fNBinsTruth1/2);++i_tru)
					{
					std::cout << i_tru + j << std::endl;
					sum1 += fMCMCx.at(i * npar + i_tru + j );
					}
	
					for (int i_tru=(fNBinsTruth1/2);i_tru<fNBinsTruth1;++i_tru)
					{
					std::cout << i_tru + fNBinsTruth2 << std::endl;
					sum2 += fMCMCx.at(i * npar + i_tru + fNBinsTruth2);
					}
				*/
	
				sum1 = fMCMCx.at(i * npar + j );
				sum2 = fMCMCx.at(i * npar + j + fNBinsTruth2);
	
				sum1tot+=sum1;
				sum2tot+=sum2;

				// fill the ratio histogram
	
				double ac = (sum2-sum1)/(sum1+sum2);
	
				//    std::cout << " ac = " << ac << std::endl;
	
				(fHistRatio.at(j))->GetHistogram()->Fill(ac);
      }
    
    fAcTotal->GetHistogram()->Fill((sum2tot-sum1tot)/(sum2tot+sum1tot));
  }
}

// ---------------------------------------------------------
void BCFBU::PrintHistogram()
{  
  // print the BAT histogram to an eps file
  for (int j=0;j<fNBinsTruth2;++j)
    {

      (fHistRatio.at(j))->Print(("ac"+IntToString(j)+".eps").c_str());

      (fHistRatio.at(j))->Print(("ac"+IntToString(j)+".png").c_str());
    }

  fAcTotal->Print("acTotal.eps");
  fAcTotal->Print("acTotal.png");
}

// ---------------------------------------------------------
Double_t BCFBU::GetCurvature(const TVectorD& vec, const TMatrixD& curv) 
{      
	// Compute curvature of vector
	return vec*(curv*vec);
}

// ---------------------------------------------------------
void BCFBU::FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC, int fDdim )
{
	Double_t eps = 0.00001;

	Int_t ndim = tCurv.GetNrows();

	// Init
	tCurv *= 0;
	tC    *= 0;

	if (fDdim == 0) 
		for (Int_t i=0; i<ndim; ++i) tC(i,i) = 1;
	else if (fDdim == 1) {
		for (Int_t i=0; i<ndim; ++i) {
			if (i < ndim-1) tC(i,i+1) = 1.0;
			tC(i,i) = 1.0;
		}
	}
	else if (fDdim == 2) {
		for (Int_t i=0; i<ndim; ++i) {
			if (i > 0)      tC(i,i-1) = 1.0;
			if (i < ndim-1) tC(i,i+1) = 1.0;
			tC(i,i) = -2.0;
		}
		tC(0,0) = -1.0;
		tC(ndim-1,ndim-1) = -1.0;
	}
   
	// Add epsilon to avoid singularities
	for (Int_t i=0; i<ndim; ++i) tC(i,i) = tC(i,i) + eps;

	//Get curvature matrix
	for (Int_t i=0; i<ndim; ++i) {
		for (Int_t j=0; j<ndim; ++j) {
			for (Int_t k=0; k<ndim; k++) {
				tCurv(i,j) = tCurv(i,j) + tC(k,i)*tC(k,j);
			}
		}
	}


}

// ---------------------------------------------------------
void BCFBU::DefineSample(std::string samplename, TH1 *nominal)
{
  
  fNSamples++;
  fSampleNames.push_back(samplename);
  fNominalHisto.push_back(nominal);
  
}

// ---------------------------------------------------------
void BCFBU::DefineSystematic(std::string samplename, std::string parname, TH1 *up, TH1 *down)
{
   
  // has this systematic appeared before?

  bool found = false;
  
  for (unsigned int i=0;i<fSystNames.size();++i)
    if (fSystNames[i]==parname)
      found = true;
  
  // if not, add it to the list of systematics

  if (!found)
    {
      fNSyst++;
      fSystNames.push_back(parname);
    }
  
  fSystUpHisto[samplename][parname] = up;
  fSystDownHisto[samplename][parname] = down;

  std::cout << " systematic " << parname << " UP: " << std::endl;

  for (int i_rec1=0;i_rec1<fNBinsReco1;i_rec1++)
    for (int i_rec2=0;i_rec2<fNBinsReco2;i_rec2++)
      {
				std::cout << " bin " << i_rec1 << " " << i_rec2 << " " << up->GetBinContent(i_rec1+1,i_rec2+1) << std::endl;
			}

  std::cout << " systematic " << parname << " DOWN: " << std::endl;

  for (int i_rec1=0;i_rec1<fNBinsReco1;i_rec1++)
    for (int i_rec2=0;i_rec2<fNBinsReco2;i_rec2++)
      {
				std::cout << " bin " << i_rec1 << " " << i_rec2 << " " << down->GetBinContent(i_rec1+1,i_rec2+1) << std::endl;
			}
   
}

// ---------------------------------------------------------
int BCFBU::GetSystParamNumber(std::string systname)
{
  int result = -1;
  
  for( std::map<std::string, int>::iterator it_map2=fSystParamMap.begin(); it_map2!=fSystParamMap.end(); it_map2++)
    {
      if ( (*it_map2).first == systname)
				result = (*it_map2).second;
    }

  if (result==-1)
    {
      std::cout << " Systematic with name " << systname << " not found - exiting. " << std::endl;
      exit(1);
    }

  return result;
}

// ---------------------------------------------------------
void BCFBU::DefineSystematicResponse(std::string parname, TH2 *responseup, TH2 *responsedown)
{
  
  int fNBinsTruth = responseup->GetYaxis()->GetNbins();
  int fNBinsReco = responseup->GetXaxis()->GetNbins();
  
  TH2D *cl_responseup = (TH2D*) responseup->Clone();
  TH2D *cl_responsedown = (TH2D*) responsedown->Clone();

  TH2D *norm_responseup = (TH2D*) responseup->Clone();
  TH2D *norm_responsedown = (TH2D*) responsedown->Clone();

  for (int i_tru=0;i_tru<fNBinsTruth;++i_tru)
    for (int i_rec=0;i_rec<fNBinsReco;++i_rec)
      {
				cl_responseup->SetBinContent(i_tru+1,i_rec+1, responseup->GetBinContent(i_rec+1,i_tru+1));
				cl_responsedown->SetBinContent(i_tru+1,i_rec+1, responsedown->GetBinContent(i_rec+1,i_tru+1)) ;
      }
  
  for (int i_tru=0;i_tru<fNBinsTruth;++i_tru)
    {
      double sum = 0;

      for (int i_rec=0;i_rec<fNBinsReco;++i_rec)  
				sum += cl_responseup->GetBinContent(i_tru+1,i_rec+1);
      
      for (int i_rec=0;i_rec<fNBinsReco;++i_rec)
				norm_responseup->SetBinContent(i_tru+1,i_rec+1, cl_responseup->GetBinContent(i_tru+1,i_rec+1)/sum);
    }

  for (int i_tru=0;i_tru<fNBinsTruth;++i_tru)
    {
      double sum = 0;

      for (int i_rec=0;i_rec<fNBinsReco;++i_rec)  
				sum += cl_responsedown->GetBinContent(i_tru+1,i_rec+1);
      
      for (int i_rec=0;i_rec<fNBinsReco;++i_rec)
				norm_responsedown->SetBinContent(i_tru+1,i_rec+1, cl_responsedown->GetBinContent(i_tru+1,i_rec+1)/sum);
    }

  std::cout << "Response up " << std::endl;

  Dump(norm_responseup);

  std::cout << "Response down " << std::endl;

  Dump(norm_responsedown);

  //  exit(1);

  fSystResponseUpMap[parname]=norm_responseup;
  fSystResponseDownMap[parname]=norm_responsedown;
}

// ---------------------------------------------------------
void BCFBU::Dump(TH2D *bla)
{
  for (int i=1; i<bla->GetYaxis()->GetNbins();++i)
    {
      for (int j=1; j<bla->GetXaxis()->GetNbins();++j)
				{
					std::cout << bla->GetBinContent(i,j) << " ";
				}
      std::cout << std::endl;
    }

}


// ---------------------------------------------------------
void BCFBU::DefineNormalisationSystematic(std::string samplename, std::string parname, double uncertainty)
{
  fNNormSyst++;
  fNormSystNames.push_back(parname);  
  fNormSystSizes.push_back(uncertainty);

  fNormUncer[samplename]=parname;
  
}

// ---------------------------------------------------------
std::string BCFBU::IntToString(int x)
{
  std::stringstream ss;// create stream
  ss << x;//add number to the stream

  return ss.str();

}

// ---------------------------------------------------------
TH1* BCFBU::GetUnfoldedResult()
{
  TH1* unfolded = (TH1D*) fHistTruth->Clone();

  if (fNDim==1)
    {
      for (int i=1;i<=GetNBinsReco();++i)
				{
					std::stringstream ss;//create a std::stringstream
					ss << i;//add number to the stream
	  
					BCH1D *temph = GetMarginalized(("T"+ss.str()).c_str());
	  
					temph->Print(("marginalised"+ss.str()+".eps").c_str());
	  
					double lowbound = 0;
					double upbound = 0;
	  
	  
					temph->GetSmallestInterval (lowbound, upbound);
					std::cout << " bin " << i << " mode " << temph->GetMode() << " mean " << temph->GetMean() << " median " << temph->GetMedian() << " low " << lowbound << " up " << upbound << std::endl;
	  
					unfolded->SetBinContent(i,(upbound+lowbound)/2);
					unfolded->SetBinError(i,(upbound-lowbound)/2);
				}
    }
  
  
  if (fNDim==2)
    {           
      for (int i=1;i<=GetNBinsTruthX();++i)
				for (int j=1;j<=GetNBinsTruthY();++j)
					{       
	    
						BCH1D *temph = GetMarginalized(("X"+IntToString(i)+"Y"+IntToString(j)).c_str());
	    
						double lowbound = 0;
						double upbound = 0;
	    
						temph->Print(("marginalised"+IntToString(i)+"_"+IntToString(j)+".eps").c_str());
	    
						temph->GetSmallestInterval (lowbound, upbound);
	    
						std::cout << " unfolded result: truth bin " << i << " " << j << " mode " << temph->GetMode() << " mean " << temph->GetMean() << " median " << temph->GetMedian() << " low " << lowbound << " up " << upbound << std::endl;
	    
						unfolded->SetBinContent(i,j,(upbound+lowbound)/2);
						unfolded->SetBinError(i,j,(upbound-lowbound)/2);
	    
					}

      for (int j=1;j<=GetNBinsTruthY();++j)
				{ 
					TCanvas *bla = new TCanvas("bla","bla",600,600);
					bla->SetLogy();
	  
					TH1D *truthX = new TH1D("","",2,-3,3);
	  
					TH1D *unfoldedX = new TH1D("","",2,-3,3);
	  
					for (int i=1;i<=GetNBinsTruthX();++i)
						{
	      
							truthX->SetBinContent(i,fHistTruth->GetBinContent(i,j));
	      
							unfoldedX->SetBinContent(i,unfolded->GetBinContent(i,j));
	      
							unfoldedX->SetBinError(i,unfolded->GetBinError(i,j));
	      
							std::cout << "Y bin " << j << " X bin " << i << " truth value = " << fHistTruth->GetBinContent(i,j) << " unfolded value " << unfolded->GetBinContent(i,j) << std::endl;
	      
						}
	  
					truthX->Draw();
	  
					unfoldedX->SetMarkerColor(kBlue);
					unfoldedX->SetLineColor(kBlue);
					unfoldedX->Draw("E same");
	  
	  
					bla->SaveAs(("unfolded_Ybin"+IntToString(j)+".eps").c_str());
					bla->SaveAs(("unfolded_Ybin"+IntToString(j)+".png").c_str());
   
				}
    }

  return unfolded;

}

// ---------------------------------------------------------
void BCFBU::PrintSystematicsPosteriors()
{
  
  // save plots of posteriors for shape systematics
   
  for (int i=0;i<GetNSyst();++i)
    {
      
      BCH1D *temph = GetMarginalized((GetSystName(i)).c_str());
      temph->Print(("marginalised_"+GetSystName(i)+".eps").c_str());
    }
  
  // save plots of posteriors for normalisation systematics
  
  for (int i=0;i<GetNNormSyst();++i)
    {
      
      BCH1D *temph = GetMarginalized((GetNormSystName(i)).c_str());
      temph->Print(("marginalised_"+GetNormSystName(i)+".eps").c_str());
    }
}

// ---------------------------------------------------------
int BCFBU::SetDataHistogram(TH1 *h_data)
{
	// check dimension
	if( (fNDim > 0) && fNDim != h_data->GetDimension()) {
		BCLog::OutWarning("BCBFU::SetDataHistogram. Dimensions of data histogram does not match previously defined dimension. Data histogram rejected.");
		return 0; 
	}
	else if (fNDim <= 0) {
		// set dimension
		fNDim = h_data->GetDimension();
	}

	// handle 1-D case
  if (fNDim==1) {
		// check number of bins
		if ( (fNBinsReco > 0) && (fNBinsReco != h_data->GetNbinsX()) ) {
			BCLog::OutWarning("BCBFU::SetDataHistogram. Number of reco bins doens't match number of bins in the data histogram. Data histogram rejected.");
			return 0; 
		}
		
		// clone histogram
		fHistData = (TH1D*) h_data->Clone(); 
		
		// (re)set number of reco bins
		fNBinsReco = fHistData->GetNbinsX();
	}

	// handle 2-D case
  if (fNDim==2) {
		// check number of bins
		if (fNBinsReco != (h_data->GetNbinsX() * h_data->GetNbinsY()) ) {
			BCLog::OutWarning("BCBFU::SetDataHistogram. Number of reco bins doens't match number of bins in the data histogram. Data histogram rejected.");
			return 0; 
		}

		fHistData = (TH2D*) h_data->Clone();

		// (re)set number of reco bins
		fNBinsReco = fHistData->GetNbinsX() * h_data->GetNbinsY();
	}

	// no error
	return 1;
}

// ---------------------------------------------------------
void BCFBU::PrintResponseMatrix(const char * filename)
{
	// create new canvas
	TCanvas * c1 = new TCanvas();
	c1->cd();
	
	// create histogram
	TH2D* hist = new TH2D("hist", ";truth;reco", 
												fNBinsTruth, fHistTruth->GetXaxis()->GetBinLowEdge(1), fHistTruth->GetXaxis()->GetBinUpEdge(fNBinsTruth),
												fNBinsReco, fRecoMin, fRecoMax);
	hist->SetStats(kFALSE);

	// copy bins
  for (int i_tru = 0; i_tru < fNBinsTruth; ++i_tru) 
    for (int i_rec = 0; i_rec < fNBinsReco; ++i_rec)  
      hist->SetBinContent(i_tru+1,i_rec+1, (*fResponseMatrix)(i_tru,i_rec));
	
	// print
	hist->Draw("COLZ");
	c1->Print(filename);

	// free memory
	delete hist;
	delete c1;
}

// ---------------------------------------------------------
void BCFBU::PrintMigrationMatrix(const char * filename)
{
	// create new canvas
	TCanvas * c1 = new TCanvas();
	c1->cd();
	
	// create histogram
	TH2D* hist = new TH2D("hist", ";truth;reco", 
												fNBinsTruth, fHistTruth->GetXaxis()->GetBinLowEdge(1), fHistTruth->GetXaxis()->GetBinUpEdge(fNBinsTruth),
												fNBinsReco, fRecoMin, fRecoMax);
	hist->SetStats(kFALSE);

	// copy bins
  for (int i_tru = 0; i_tru < fNBinsTruth; ++i_tru) 
    for (int i_rec = 0; i_rec < fNBinsReco; ++i_rec)  
      hist->SetBinContent(i_tru+1,i_rec+1, (*fMigrationMatrix)(i_tru,i_rec));
	
	// print
	hist->Draw("COLZ");
	c1->Print(filename);

	// free memory
	delete hist;
	delete c1;
}

// ---------------------------------------------------------
void BCFBU::PrintEfficiencyHistogram(const char * filename)
{
	// check if efficiency histogram exists
  if (!fHistEfficiency) {
		BCLog::OutWarning("BCFBU::PrintEfficiency: efficiency histogram not found. Exit.");
		exit(1);
	}
 
	// create new canvas
	TCanvas * c1 = new TCanvas();
	c1->cd();
	fHistEfficiency->Draw();
	c1->Print(filename);

	// free memory
	delete c1;
}

// ---------------------------------------------------------
void BCFBU::PrintTruthHistogram(const char * filename)
{
	// check if truth histogram exists
  if (!fHistTruth) {
		BCLog::OutWarning("BCFBU::PrintTruth: truth histogram not found. Exit.");
		exit(1);
	}
 
	// create new canvas
	TCanvas * c1 = new TCanvas();
	c1->cd();
	fHistTruth->Draw();
	c1->Print(filename);

	// free memory
	delete c1;
}

// ---------------------------------------------------------
void BCFBU::PrintDataHistogram(const char * filename)
{
	// check if data histogram exists
  if (!fHistData) {
		BCLog::OutWarning("BCFBU::PrintData: data histogram not found. Exit.");
		exit(1);
	}
 
	// create new canvas
	TCanvas * c1 = new TCanvas();
	c1->cd();
	fHistData->Draw();
	c1->Print(filename);

	// free memory
	delete c1;
}

// ---------------------------------------------------------
