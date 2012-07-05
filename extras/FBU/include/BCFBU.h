#ifndef __BAT__BCFBU__H
#define __BAT__BCFBU__H

#include <BAT/BCModel.h>
#include <BAT/BCH1D.h>

#include <BCFBUBackground.h>

#include <TMatrixD.h>

#include <map>
#include <utility>

class TH1;
class TH2;
class TH1D;
class TH2D;

// This is a BCFBU header file.
// Model source code is located in file src/BCFBU.cxx

// ---------------------------------------------------------
class BCFBU : public BCModel
{
 public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      BCFBU();

      /**
       * A constructor.
       * @param name The name of the model */
      BCFBU(const char * name);

      /**
       * The default destructor. */
      ~BCFBU();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

			/** 
			 * @return the number of bins of the truth histogram in the 1D
			 * case */
      int GetNBinsTruth() 
			{return fNBinsTruth;};

			/** 
			 * @return the number of bins of the reco histogram in the 1D
			 * case */
      int GetNBinsReco() 
			{return fNBinsReco;};

			/** 
			 * @return the number of bins of the truth histogram in the 2D
			 * case (axis 1) */
      int GetNBinsTruthX() 
			{return fNBinsTruth1;};

			/** 
			 * @return the number of bins of the truth histogram in the 2D
			 * case (axis 2) */
      int GetNBinsTruthY() 
			{return fNBinsTruth2;};

			/** 
			 * @return the number of bins of the reco histogram in the 2D
			 * case (axis 1) */
      int GetNBinsRecoX() 
			{return fNBinsReco1;};

			/** 
			 * @return the number of bins of the reco histogram in the 2D
			 * case (axis 2) */
      int GetNBinsRecoY() 
			{return fNBinsReco2;};

      int Get1DIndex(int i,int j) 
			{return (fNBinsTruth2*i+j);};

      int GetNSyst() 
			{return this->fNSyst;};

      int GetNNormSyst() 
			{return this->fNNormSyst;};

      std::pair<int,int> Get2DIndex(int i,int j, int p, int q) 
			{ std::pair<int,int> res(fNBinsTruth1*j+i,fNBinsReco1*q+p); 
				return res;};

      std::string GetSystName(int i) 
				{return fSystNames[i];};

      std::string GetNormSystName(int i) 
				{return fNormSystNames[i];};

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      int SetDataHistogram(TH1 *h_data);

       /** @} */
      /** \name Member functions (miscellaneous methods) */
      /** @{ */

     // Methods to overload, see file src/BCFBU.cxx
      void DefineParameters(int mode=0, double min=0, double max=1.0e6);
      double LogAPrioriProbability(const std::vector <double> & params);
      double LogLikelihood(const std::vector<double> & params);
			void MCMCIterationInterface();

      /**
       * Prepares the response matrix and calculates the efficiency
       * @param h_migration the migration matrix
			 * @param h_truth the truth distribution
			 * @param h_background the background distribution
			 * @return an error code
			 */
      int PrepareResponseMatrix(TH2* h_migration, TH1* h_truth, TH1 *h_background, std::vector<double> paramin = std::vector<double>(0), std::vector<double> parmax = std::vector<double>(0));

      /**
       * Rebins all histograms
       * @param h_migration the migration matrix
			 * @param h_truth the truth distribution
			 * @param h_background the background distribution
			 * @param rebin a rebin factor (default 1, no rebinning) 
			 * @return an error code
			 */
			int RebinHistograms(int rebin, TH2* h_migration, TH1 *h_truth, TH1 *h_background=0); 
      
      Double_t GetCurvature(const TVectorD& vec, const TMatrixD& curv);
      void FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC, int fDdim );

      void AddBackgroundProcess(std::string backgroundname, TH1 *h_background);

      void DefineSystematic(std::string samplename, std::string parname, TH1 *up, TH1 *down);
      
      void DefineSystematicResponse(std::string parname, TH2 *responseup, TH2 *responsedown);
      
      void DefineNormalisationSystematic(std::string samplename, std::string parname, double uncertainty);

      TH1* GetUnfoldedResult();

       /** @} */
      /** \name Member functions (printing) */
      /** @{ */

			/**
			 * Prints the response matrix. 
			 * @param filename a filename
			 */
			void PrintResponseMatrix(const char * filename = "response.eps");

			/**
			 * Prints the migration matrix. 
			 * @param filename a filename
			 */
			void PrintMigrationMatrix(const char * filename = "migration.eps");

			/**
			 * Prints the efficiency histogram
			 * @param filename a filename
			 */
			void PrintEfficiencyHistogram(const char * filename = "efficiency.eps");

			/**
			 * Prints the truth histogram
			 * @param filename a filename
			 */
			void PrintTruthHistogram(const char * filename = "truth.eps");

			/**
			 * Prints the data histogram
			 * @param filename a filename
			 */
			void PrintDataHistogram(const char * filename = "truth.eps");

      void PrintSystematicsPosteriors();
            
      void PrintHistogram();
      
   /** @} */

 private:
			
			/**
       * Histogram of the truth distribution. Entries are number of
       * events. */
      TH1* fHistTruth;
			
     /**
       * Number of bins of the truth histogram in the 1D case. */
      int fNBinsTruth; 

     /**
       * Number of bins of the truth histogram in the 2D case (axis 1). */
      int fNBinsTruth1; 
			
			/**
       * Number of bins of the truth histogram in the 2D case (axis 2). */
			int fNBinsTruth2;
			
			/**
       * Lower range of truth histogram in the 1D case. */			
			double fTruthMin;

			/**
       * Upper range of truth histogram in the 1D case. */			
			double fTruthMax;

			/**
       * Histogram of the observed distribution. Entries are number of
       * events. */
			TH1* fHistData;

			int fNBinsReco;
      int fNBinsReco1;
			int fNBinsReco2;

			double fRecoMin;
			double fRecoMax;
      
			/**
			 * Response matrix: describes conditional probability of observing r if t was produced, P(r;t).
			 * x: truth
			 * y: reco
			 * normalization: sum over all reco <= 1 for a given truth value */
      TMatrixD *fResponseMatrix;     
      
			/**
			 * Migration matrix: describes joint probability of observing r and having produced t, P(r,t).
			 * x: truth
			 * y: reco
			 * normalization: <= 1 */
      TMatrixD *fMigrationMatrix;     
      
			/**
			 * Efficiency for a certain truth value t to be reconstructed at
			 * all. */
			TH1 *fHistEfficiency;

      std::vector<BCH1D*> fHistRatio;
      
      BCH1D *fAcTotal;
      
			/**
			 * The dimension of the data spectrum (1D or 2D). */
      int fNDim;

			// a container of background processes
			std::vector<BCFBUBackground*> fBackgroundProcesses;
      
      int fNSyst, fNSamples;
      std::vector<std::string> fSystNames;
      std::vector<std::string> fSampleNames;
      std::vector<TH1*> fNominalHisto;

      int fNNormSyst;
      std::vector<std::string> fNormSystNames;
      std::vector<double> fNormSystSizes;
      std::map<std::string, std::string> fNormUncer;

      std::map<std::string, std::map<std::string, TH1*> > fSystUpHisto;
      std::map<std::string, std::map<std::string, TH1*> > fSystDownHisto;
      std::map<std::string, std::map<std::string, TH2*> > fSystResponseUp;
      std::map<std::string, std::map<std::string, TH2*> > fSystResponseDown;
      
      std::map<std::string, int> fSystParamMap;
      std::map<std::string, TH2*> fSystResponseUpMap;
      std::map<std::string, TH2*> fSystResponseDownMap;
      
      BCDataSet * fDataSet;
      
      double NoInterpolation(double alpha, double nominal, double up, double down);
      double LinearInterpolate(double alpha, double nominal, double up, double down);
      double ExponentialInterpolate(double alpha, double nominal, double up, double down);
      int GetSystParamNumber(std::string systname);

      int fInterpolationType;

      std::string IntToString(int x);
      
      void Dump(TH2D *bla);

			// helper variables
			std::vector<double> fVectorTruth;
			std::vector<double> fVectorReco;

};
// ---------------------------------------------------------

#endif

